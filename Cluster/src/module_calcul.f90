MODULE module_calcul
  USE module_structure
  USE module_solve
  USE module_embed
CONTAINS

  !*****************************************
  !calcul du sigma
  SUBROUTINE calculsigma(dataw,sigma)
    IMPLICIT NONE
    TYPE(type_data) :: dataw
    DOUBLE PRECISION :: sigma,norme,norme1,sigma1
    INTEGER :: i1,j1,k1

    sigma=0.0
    sigma1=0.0
    DO i1=1,dataw%nb
       DO j1=i1+1,dataw%nb
          norme=0.0
          DO k1=1,dataw%dim
             norme=norme+&
                  (dataw%point(i1)%coord(k1)-dataw%point(j1)%coord(k1))**2
          ENDDO
          sigma=max(sigma,sqrt(norme))
       ENDDO
    ENDDO
    sigma=sigma/(2*exp(log(float(dataw%nb))*(1.0/float(dataw%dim))))
    !securite
    IF (sigma==0.0) sigma=1.0
    RETURN
  END SUBROUTINE calculsigma

  !*****************************************
  !calcul du sigma pour l'interface
  SUBROUTINE calculsigmainterface(numproc,dataw,sigma,bornes,decoupe,epsilon)
    IMPLICIT NONE
    TYPE(type_data) :: dataw
    DOUBLE PRECISION :: sigma
    INTEGER :: i,j,k,numproc,nb,ok
    DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: bornes
    INTEGER,DIMENSION(:),POINTER :: decoupe,decoupe0
    INTEGER,DIMENSION(:,:),POINTER :: tableau
    DOUBLE PRECISION :: epsilon
    DOUBLE PRECISION :: long,sigma0
    DOUBLE PRECISION :: volext,volint
    !nb de decoupes
    nb=1
    DO i=1,dataw%dim
       nb=nb*decoupe(i)
    ENDDO
    !creation des decoupes
    ALLOCATE(tableau(nb,0:dataw%dim))
    ALLOCATE(decoupe0(dataw%dim)); decoupe0(:)=1
    DO i=1,nb
       DO j=1,dataw%dim
          tableau(i,j)=decoupe0(j)
       ENDDO
       decoupe0(1)=decoupe0(1)+1
       ok=0; k=1
       DO WHILE(ok==0)
          IF (decoupe0(k)>decoupe(k)) THEN
             decoupe0(k)=1
             IF (k<dataw%dim) decoupe0(k+1)=decoupe0(k+1)+1
          ELSE
             ok=1
          ENDIF
       ENDDO
    ENDDO
    DEALLOCATE(decoupe0)
    !valeur de sigma
    sigma0=0.0
    DO i=1,nb
       volext=1.0; volint=1.0
       DO j=1,dataw%dim
          k=tableau(i,j)
          long=bornes(j,k,2)-bornes(j,k,1)
          volext=volext*long
          volint=volint*max(0.0,long-2.0*epsilon)
       ENDDO
       sigma0=sigma0+volext-volint
    ENDDO
    DEALLOCATE(tableau)
    !calcul du grandeur equivalente
    sigma0=exp(1.0/float(dataw%dim)*log(sigma0))
    !calcul du sigma
    sigma0=sigma0/(2.0*exp(log(float(dataw%nb))*(1.0/float(dataw%dim))))
    !calcul du sigma formule globale
    CALL calculsigma(dataw,sigma)
#if aff
    PRINT *,numproc,'valeur de sigma calculee pour interface:',sigma0
#endif
#if aff
    PRINT *,numproc,'valeur sigma interface',sigma
#endif
    RETURN
  END SUBROUTINE calculsigmainterface


  !*****************************************
  !calcul des clusters
  SUBROUTINE calculclusters(numproc,nblimit,nbideal,dataw,sigma)
    IMPLICIT INTEGER(i,j,q)

    INCLUDE 'mpif.h'
    TYPE(type_data) :: dataw
    INTEGER :: numproc,nbproc
    DOUBLE PRECISION :: sigma
    DOUBLE PRECISION,DIMENSION(:,:),POINTER :: A,Z,A2,cluster_center
    DOUBLE PRECISION,DIMENSION(:),POINTER :: W
    INTEGER :: n,k,nbcluster
    DOUBLE PRECISION, DIMENSION(:),POINTER :: D
    DOUBLE PRECISION,DIMENSION(:),POINTER :: ratiomax,cluster_energy,&
         ratiomin,ratiomoy,ratiorii,ratiorij
    INTEGER,DIMENSION(:),POINTER ::cluster,cluster_population,nbinfo
    INTEGER :: nblimit, nbideal, nb, nbvp
    DOUBLE PRECISION :: norme,ratio,ratio1,ratio2,seuilrij
    CHARACTER*30 :: num,files
    DOUBLE PRECISION :: t1, t2, t_cons_vp

    ! deux valeurs qui quand elles ne sont pas declarees
    ! et donc implicitement des REAL font que ca marche mieux
    DOUBLE PRECISION :: val, value

    ! solveur au valeur propre => parametre de controle
    INTEGER :: solver

    !creation de la matrice
    PRINT *,numproc,'valeur du sigma',sigma
#if aff
    PRINT *,numproc,'valeur du sigma',sigma
#endif
    n=dataw%nb
    ALLOCATE(A(n,n));A(:,:)=0.0

    ! A(i,i) = 0
    DO i=1,n-1
       DO j=i+1,n
          norme=0.0
          DO k=1,dataw%dim
             norme=norme+(dataw%point(i)%coord(k)-dataw%point(j)%coord(k))**2
          ENDDO
          value=exp(-norme/sigma)
          ! partie triangulaire sup
          A(i,j) = value
          ! partie triangulaire inf
          A(j,i)=A(i,j)
       ENDDO
    ENDDO

    !Normalisation de la matrice affinite
    ALLOCATE(D(n)); D(:)=0.0
    DO i=1,n
       DO j=1,i-1
          D(i)=D(i)+A(i,j)
       ENDDO
       DO j=i+1, n
          D(i)=D(i)+A(j,i)
       ENDDO
    ENDDO

    DO i=1,n
       DO j=1,n
          ! la matrice A n'est plus symetrique
          A(i,j)=A(i,j)/D(i)
       ENDDO
    ENDDO
    DEALLOCATE(D)

    solver = 0

    IF(solver == 0) THEN
      PRINT *, numproc, 'solveur lapack'

      nbvp = n

      ALLOCATE(A2(n,n)); A2(:,:)=0.0
      DO i=1,n
        DO j=1,n
          A2(i,j)=A(i,j)
        ENDDO
      ENDDO

      ! lapack
      ! n DIMENSION in
      ! A2 matrice in
      ! Z(N,N) vecteurs propres out
      ! W(N) valeurs propres out

      t1 = MPI_WTIME()
      CALL solve_dgeev(n,A2,Z,W)
    ELSE
      PRINT *, numproc, 'solveur arpack'

      ! arpack
      nb = 2*nblimit
      nbvp = nb
      CALL solve_arpack_full(A, n, nb, W, Z)

    ENDIF

    t2 = MPI_WTIME()

    t_cons_vp = t2 - t1
    PRINT *, numproc, 'cout construction vp', t_cons_vp

    DO i=1,nbvp-1
       DO j=i+1,nbvp
          IF (W(i)<W(j)) THEN
             val=W(i); W(i)=W(j); W(j)=val
             DO k=1,n
                val=Z(k,i); Z(k,i)=Z(k,j); Z(k,j)=val
             ENDDO
          ENDIF
       ENDDO
    ENDDO

    !Test spectral embedding avec different nbcluster   
    !***********************
    ! Spectral embedding

    IF ((nbideal==0).AND.(n>2)) THEN
       !** recherche du meilleur decoupage
       ALLOCATE(ratiomax(nblimit)); ratiomax(:)=0
       ALLOCATE(ratiomin(nblimit)); ratiomin(:)=0
       ALLOCATE(ratiomoy(nblimit)); ratiomoy(:)=0
       ALLOCATE(ratiorii(nblimit)); ratiorii(:)=0
       ALLOCATE(ratiorij(nblimit)); ratiorij(:)=0

       ALLOCATE(nbinfo(nblimit)); nbinfo(:)=0
       DO nbcluster=2,min(n,nblimit)

          ALLOCATE(cluster(n));cluster(:)=0.0
          ALLOCATE(cluster_center(nbcluster,nbcluster)); cluster_center(:,:)=0.0
          ALLOCATE(cluster_population(nbcluster));cluster_population(:)=0.0
          ALLOCATE(cluster_energy(nbcluster));cluster_energy(:)=0.0

          CALL spectral_embedding(nbcluster,n,Z,A,&
               ratiomax(nbcluster),cluster,cluster_center,cluster_population,&
               cluster_energy,nbinfo(nbcluster),numproc,ratiomoy(nbcluster), &
               ratiorij(nbcluster),ratiorii(nbcluster))


          DEALLOCATE(cluster);DEALLOCATE(cluster_center);
          DEALLOCATE(cluster_energy)
          DEALLOCATE(cluster_population);
       ENDDO

#if aff
PRINT *, 'ratio de frobenius'
#endif
       !*******************************
       ! Ratio de norme de frobenius
       ratio=ratiomax(nblimit)
       dataw%nbclusters=nblimit
       ratio1=0.0;ratio2=1e+10
       DO i=2,nblimit
          IF ((numproc==0).AND.(nbproc>1)) THEN 
             seuilrij=1e-1
          ELSE
             seuilrij=1e-4
          ENDIF
          IF ((ratiorii(i)>=0.95*ratio1).AND.(ratiorij(i)-ratio2<=seuilrij)) THEN  
             dataw%nbclusters=i
             ratio1=ratiorii(i)
             ratio2=ratiorij(i)
          ENDIF
       ENDDO

    ELSEIF ((nbideal==1).AND.(n>nbideal)) THEN
       !** test avec un cluster impose
       ALLOCATE(nbinfo(nbideal)); nbinfo(:)=0
       ALLOCATE(ratiomin(1)); ratiomin(:)=0.0
       dataw%nbclusters=nbideal
    ELSE
       !** cas d'un domaine avec moins de points que nbideal ou 1 seul point
       ALLOCATE(nbinfo(n)); nbinfo(:)=0
       ALLOCATE(ratiomin(1)); ratiomin(:)=0.0
       dataw%nbclusters=n
       ALLOCATE(ratiomax(n)); ratiomax(:)=0
       ALLOCATE(ratiomoy(n)); ratiomoy(:)=0
       ALLOCATE(ratiomin(n)); ratiomin(:)=0
       ALLOCATE(ratiorii(n)); ratiorii(:)=0
       ALLOCATE(ratiorij(n)); ratiorij(:)=0
    ENDIF
    ! cas ou nbcluster==1
    IF (dataw%nbclusters==2) THEN
       PRINT *, 'difference ratio',ratiorij(2)/ratiorii(2)
       IF (ratiomax(2)>=0.6) THEN 
          dataw%nbclusters=1
       ELSE 
          dataw%nbclusters=2
       ENDIF
    ENDIF
#if aff
    PRINT *,numproc,'cluster final obtenu : ',dataw%nbclusters
#endif

    !** calcul du clustering final
    IF (dataw%nbclusters>1) THEN
       CALL spectral_embedding(dataw%nbclusters,n,Z,A,ratio,cluster,&
            cluster_center,cluster_population,cluster_energy,&
            nbinfo(dataw%nbclusters),numproc,ratiomin(1),ratiorij(1),ratiorii(1))
       DO i=1,dataw%nb
          dataw%point(i)%cluster=cluster(i)
       ENDDO
       DEALLOCATE(cluster)
       DEALLOCATE(cluster_population)
       DEALLOCATE(ratiomax)
       DEALLOCATE(cluster_energy)
       DEALLOCATE(ratiomin)
       DEALLOCATE(ratiomoy)
       DEALLOCATE(ratiorii)
       DEALLOCATE(ratiorij)
       DEALLOCATE(A)
       DEALLOCATE(Z)
       IF(solver == 0) DEALLOCATE(A2)
       DEALLOCATE(cluster_center)
       DEALLOCATE(W)
    ELSE 
#if aff
       PRINT *, numproc, 'ok'
#endif
       DO i=1,dataw%nb
          dataw%point(i)%cluster=1
       ENDDO
#if aff
       PRINT *,numproc,'cluster'
#endif
    ENDIF
	
    RETURN
  END SUBROUTINE calculclusters

END MODULE module_calcul
