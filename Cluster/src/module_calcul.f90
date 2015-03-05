MODULE module_calcul
  USE module_structure
  USE module_solve
  USE module_embed
CONTAINS


  SUBROUTINE get_sigma(dataw, sigma)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: dataw

    !====  OUT ====
    DOUBLE PRECISION :: sigma

    !#### Variables  ####
    DOUBLE PRECISION :: norme
    DOUBLE PRECISION :: sigma1
    INTEGER :: i1
    INTEGER :: j1
    INTEGER :: k1

    !###########################################
    ! INSTRUCTIONS
    !###########################################
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
    ! Safety
    IF (sigma==0.0) sigma=1.0
    RETURN
  END SUBROUTINE get_sigma


  SUBROUTINE get_sigma_interface(numproc, dataw, sigma, bounds, decoupe, epsilon)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: dataw
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bounds
    DOUBLE PRECISION :: epsilon
    INTEGER, DIMENSION(:), POINTER :: decoupe
    INTEGER :: numproc

    !====  OUT ====
    DOUBLE PRECISION :: sigma

    !#### Variables  ####
    INTEGER, DIMENSION(:,:), POINTER :: tableau
    INTEGER, DIMENSION(:), POINTER :: decoupe0
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: nb
    DOUBLE PRECISION :: long
    DOUBLE PRECISION :: sigma0
    DOUBLE PRECISION :: volext
    DOUBLE PRECISION :: volint

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    !nb de decoupes
    nb=1
    DO i=1,dataw%dim
       nb=nb*decoupe(i)
    ENDDO
    ! Creation of partitionning
    ALLOCATE(tableau(nb,0:dataw%dim))
    ALLOCATE(decoupe0(dataw%dim))
    decoupe0(:)=1
    DO i=1,nb
       DO j=1,dataw%dim
          tableau(i,j)=decoupe0(j)
       ENDDO
       decoupe0(1)=decoupe0(1)+1
       k=1
       DO WHILE(decoupe0(k)>decoupe(k))
          decoupe0(k)=1
          IF (k<dataw%dim) decoupe0(k+1)=decoupe0(k+1)+1
       ENDDO
    ENDDO
    DEALLOCATE(decoupe0)
!!======================TODO : debut de #if aff ??????
    ! Value of sigma
    sigma0=0.0
    DO i=1,nb
       volext=1.0
       volint=1.0
       DO j=1,dataw%dim
          k=tableau(i,j)
          long=bounds(j,k,2)-bounds(j,k,1)
          volext=volext*long
          volint=volint*max(0.0D1,long-2.0*epsilon)
       ENDDO
       sigma0=sigma0+volext-volint
    ENDDO
    DEALLOCATE(tableau)
    ! Computing of scale length
    sigma0=exp(1.0/float(dataw%dim)*log(sigma0))
    ! Sigma computing
    sigma0=sigma0/(2.0*exp(log(float(dataw%nb))*(1.0/float(dataw%dim))))
!!======================TODO : fin de #if aff ??????
#if aff
    PRINT *,numproc,'valeur de sigma calculee pour interface:',sigma0
#endif
    ! Sigma computing, global formula
    CALL get_sigma(dataw,sigma)
#if aff
    PRINT *,numproc,'valeur sigma interface',sigma
#endif
    RETURN
  END SUBROUTINE get_sigma_interface

  SUBROUTINE apply_spectral_clustering(numproc, nblimit, nbideal, dataw, sigma)
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    DOUBLE PRECISION :: sigma
    INTEGER :: nbideal
    INTEGER :: nblimit
    INTEGER :: numproc

    !=== IN/OUT ===
    TYPE(type_data) :: dataw

    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A2
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: cluster_center
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z
    DOUBLE PRECISION, DIMENSION(:), POINTER :: cluster_energy
    DOUBLE PRECISION, DIMENSION(:), POINTER :: D
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratiomax
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratiomin
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratiomoy
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratiorii
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratiorij
    DOUBLE PRECISION, DIMENSION(:), POINTER :: W
    DOUBLE PRECISION :: norme
    DOUBLE PRECISION :: ratio
    DOUBLE PRECISION :: ratio1
    DOUBLE PRECISION :: ratio2
    DOUBLE PRECISION :: seuilrij
    DOUBLE PRECISION :: t_cons_vp
    DOUBLE PRECISION :: t1
    DOUBLE PRECISION :: t2
    DOUBLE PRECISION :: val
    DOUBLE PRECISION :: value
    INTEGER, DIMENSION(:), POINTER :: cluster
    INTEGER, DIMENSION(:), POINTER :: cluster_population
    INTEGER, DIMENSION(:), POINTER :: nbinfo
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: n
    INTEGER :: nb
    INTEGER :: nbcluster
    INTEGER :: nbproc !TODO : mettre en parametre et WTF faut-il faire car on lit une variable vide
    INTEGER :: nbvp
    INTEGER :: solver ! solveur au valeur propre => parametre de controle

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ! Matrix creation
    PRINT *,numproc,'valeur du sigma',sigma
    n=dataw%nb
    ! Forall i, A(i,i) = 0
    ALLOCATE(A(n,n))
    A(:,:)=0.0

    DO i=1,n-1
       DO j=i+1,n
          norme=0.0
          DO k=1,dataw%dim
             norme=norme+(dataw%point(i)%coord(k)-dataw%point(j)%coord(k))**2
          ENDDO
          value=exp(-norme/sigma)
          ! Upper triangular part
          A(i,j) = value
          ! Lower triangular part
          A(j,i)=A(i,j)
       ENDDO
    ENDDO

    ! Normalizing of affinity matrix
    ALLOCATE(D(n))
    D(:)=0.0
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
          ! Matrix A is not symmetric anymore
          A(i,j)=A(i,j)/D(i)
       ENDDO
    ENDDO
    DEALLOCATE(D)

    solver = 0

    IF(solver == 0) THEN
      PRINT *, numproc, 'solveur lapack'

      nbvp = n

      ALLOCATE(A2(n,n))
      A2(:,:)=0.0
      DO i=1,n
        DO j=1,n
          A2(i,j)=A(i,j)
        ENDDO
      ENDDO

      t1 = MPI_WTIME()
      CALL solve_dgeev(n,A2,Z,W)
    ELSE
      PRINT *, numproc, 'solveur arpack'

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
             val=W(i)
             W(i)=W(j)
             W(j)=val
             DO k=1,n
                val=Z(k,i)
                Z(k,i)=Z(k,j)
                Z(k,j)=val
             ENDDO
          ENDIF
       ENDDO
    ENDDO

    ! Spectral embedding
    IF ((nbideal==0).AND.(n>2)) THEN
       ! Search of the best partitionning
       ALLOCATE(ratiomax(nblimit))
       ratiomax(:)=0
       ALLOCATE(ratiomin(nblimit))
       ratiomin(:)=0
       ALLOCATE(ratiomoy(nblimit))
       ratiomoy(:)=0
       ALLOCATE(ratiorii(nblimit))
       ratiorii(:)=0
       ALLOCATE(ratiorij(nblimit))
       ratiorij(:)=0

       ALLOCATE(nbinfo(nblimit))
       nbinfo(:)=0
       DO nbcluster=2,min(n,nblimit)

          ALLOCATE(cluster(n))
          cluster(:)=0
          ALLOCATE(cluster_center(nbcluster,nbcluster))
          cluster_center(:,:)=0.0
          ALLOCATE(cluster_population(nbcluster))
          cluster_population(:)=0
          ALLOCATE(cluster_energy(nbcluster))
          cluster_energy(:)=0.0

          CALL spectral_embedding(nbcluster,n,Z,A,&
               ratiomax(nbcluster),cluster,cluster_center,cluster_population,&
               cluster_energy,nbinfo(nbcluster),numproc,ratiomoy(nbcluster), &
               ratiorij(nbcluster),ratiorii(nbcluster))


          DEALLOCATE(cluster)
          DEALLOCATE(cluster_center)
          DEALLOCATE(cluster_energy)
          DEALLOCATE(cluster_population)
       ENDDO

#if aff
PRINT *, 'ratio de frobenius'
#endif
       ! Ratio of frobenius norm
       ratio=ratiomax(nblimit)
       dataw%nbclusters=nblimit
       ratio1=0.0
       ratio2=1e+10
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
       ! Test with an imposed cluster
       ALLOCATE(nbinfo(nbideal))
       nbinfo(:)=0
       ALLOCATE(ratiomin(1))
       ratiomin(:)=0.0
       dataw%nbclusters=nbideal
    ELSE
       ! Case of a domain with less points than nbideal or only one point
       ALLOCATE(nbinfo(n))
       nbinfo(:)=0
       ALLOCATE(ratiomin(1))
       ratiomin(:)=0.0
       dataw%nbclusters=n
       ALLOCATE(ratiomax(n))
       ratiomax(:)=0
       ALLOCATE(ratiomoy(n))
       ratiomoy(:)=0
       ALLOCATE(ratiomin(n))
       ratiomin(:)=0
       ALLOCATE(ratiorii(n))
       ratiorii(:)=0
       ALLOCATE(ratiorij(n))
       ratiorij(:)=0
    ENDIF
    ! Case of nbcluster==1
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

    ! Final clustering computing
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
  END SUBROUTINE apply_spectral_clustering

END MODULE module_calcul
