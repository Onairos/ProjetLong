MODULE module_embed
  USE module_structure
CONTAINS



  SUBROUTINE spectral_embedding(nbcluster,n,Z,A,ratio,cluster,&
       cluster_center,cluster_population,cluster_energy,nbinfo,numproc,&
       ratiomoy,ratiorij,ratiorii)

    !*****************************************
    ! spectral embedding
    !
    ! nbcluster = nbre de cluster
    ! dataw : points
    ! Z : matrice des vecteurs propres
    ! M : nbre de vp trouvÃÂÃÂ©esx
    ! ratio : max des ration de frob sur matrice aff rÃÂÃÂ©ordonnancÃÂÃÂ©e suivant
    ! les clusters
    ! cluster : appartenance des clusters
    ! cluster_center : centre des nbclusters clusters
    ! cluster_population : nbre de points par cluster
    ! cluster_energy : somme des ÃÂÃÂ©nergies par cluster
    !

    IMPLICIT NONE
    !TYPE(type_data) :: dataw
    DOUBLE PRECISION,DIMENSION(:,:),POINTER:: Z,A,cluster_center
    INTEGER ::nbcluster,n,nbinfo,numproc
    DOUBLE PRECISION ::ratio,test,ratiomin,ratiorii,ratiorij, ratiomoy
    DOUBLE PRECISION,DIMENSION(:),POINTER :: cluster_energy,Z3
    INTEGER,DIMENSION(:),POINTER ::cluster,cluster_population
    !INTEGER,DIMENSION(:),POINTER::ordaffperclus
    DOUBLE PRECISION, DIMENSION(:,:),POINTER :: Frob
    DOUBLE PRECISION,DIMENSION(:,:),POINTER::Z1,Z2
    INTEGER :: it_max,it_num,i,j,k
    INTEGER,DIMENSION(:,:),POINTER :: clustercorresp
    INTEGER :: ki,kj,ni,nj,ok,nbmax
    ALLOCATE(cluster(n));
    ALLOCATE(cluster_center(nbcluster,nbcluster));
    ALLOCATE(cluster_population(nbcluster));
    ALLOCATE(cluster_energy(nbcluster));
    ALLOCATE(Z1(n,nbcluster));  ALLOCATE(Z2(nbcluster,n));
    ALLOCATE(Z3(n));Z3(:)=0.0

    DO i=1,n
       DO j=1,nbcluster
          Z1(i,j)=Z(i,j)
          !if (i==1) PRINT *,numproc,'matrice vp',j,W(j)
          Z3(i)=Z3(i)+Z1(i,j)**2
       ENDDO
    ENDDO

    DO i=1,n
       test=0.0
       DO j=1,nbcluster
          Z2(j,i)=Z1(i,j)/(sqrt(Z3(i)))
          test=test+Z2(j,i)**2
       ENDDO
    ENDDO

!    PRINT *, numproc,'methode kmeans'

    it_max=n*n !1000.

    CALL kmeans_01 ( nbcluster, n, nbcluster, it_max, it_num,Z2,&
         cluster, cluster_center, cluster_population, cluster_energy, &
         numproc)
    !PRINT *,numproc,'fin de kmeans. nb d iterations effectuees : ',it_num

    !PRINT *,numproc,'Nombre points par cluster', cluster_population
    ! PRINT *,'vecteur cluster',cluster(1:5) 

    !*****************************
    ! Mesure de qualitÃÂÃÂ©
    !PRINT *,'Indexation'
    nbmax=0
    DO i=1,nbcluster
       nbmax=max(nbmax,cluster_population(i))
    ENDDO
    ALLOCATE(clustercorresp(nbcluster,nbmax)); clustercorresp(:,:)=0
    DO i=1,n
       j=cluster(i)
       ok=0;k=1
       DO WHILE(ok==0)
          IF (clustercorresp(j,k)==0) THEN
             ok=1
          ELSE
             k=k+1
          ENDIF
       ENDDO
       clustercorresp(j,k)=i
    ENDDO
    ALLOCATE(Frob(nbcluster,nbcluster)); Frob(:,:)=0.0 
    DO i=1,nbcluster
       DO j=1,nbcluster
          DO ki=1,cluster_population(i)
             ni=clustercorresp(i,ki)
             DO kj=1,cluster_population(j)
                nj=clustercorresp(j,kj)
                Frob(i,j)=Frob(i,j)+A(ni,nj)**2
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(clustercorresp)
    ratio=0.0; ratiomin=1.D+16;ratiorii=0.0;ratiorij=0.0
    nbinfo=nbcluster
    DO i=1,nbcluster
       IF ((cluster_population(i)/=0).AND.(Frob(i,i)/=0)) THEN
          DO j=1,nbcluster
             IF (i/=j) THEN
                ratio=ratio+Frob(i,j)/Frob(i,i)
                !ratio=max(ratio,Frob(i,j)/Frob(i,i))
                ratiomoy=ratiomoy+Frob(i,j)/Frob(i,i)
                ratiorij=ratiorij+Frob(i,j)
                ratiorii=ratiorii+Frob(i,i)
                ratiomin=min(ratiomin,Frob(i,j)/Frob(i,i))
             ENDIF
          ENDDO
       ELSE
          nbinfo=nbinfo-1
       ENDIF
       ratiorij=ratiorij*2/(nbcluster*(nbcluster-1))
       ratiorii=ratiorii!/nbcluster
    ENDDO
    DEALLOCATE(Frob)

!!$    p2=1
!!$    dataw%point(:)%cluster=cluster(:)
!!$    ALLOCATE(P(n,n)); P(:,:)=0.0
!!$    !PRINT *,'petit test', dataw%point(1:5)%cluster
!!$    DO i=1,nbcluster
!!$       DO j=1,n
!!$          ! P(j,:)=0.0
!!$          IF (dataw%point(j)%cluster==i) THEN 
!!$             P(p2,j)=1.0
!!$             p2=p2+1
!!$          ENDIF
!!$       ENDDO
!!$    ENDDO
!!$    !PRINT *,
!!$    !PRINT *, 'matrice permut',P
!!$    ALLOCATE(P1(n,n)); P1(:,:)=0.0
!!$    P1=MATMUL(P,A)
!!$    PRINT *,numproc,'Mesure de qualite : Frobenius norm'
!!$    ALLOCATE(ordaffperclus(nbcluster+1));ordaffperclus(:)=0.0
!!$    ordaffperclus(1)=0.0
!!$    DO i=2,nbcluster+1
!!$       ordaffperclus(i)=ordaffperclus(i-1)+cluster_population(i-1)
!!$    ENDDO
!!$    ! PRINT *, 'ordreaffperclus',  ordaffperclus
!!$    ALLOCATE(Frob(nbcluster,nbcluster)); Frob(:,:)=0.0 
!!$    nbinfo=nbcluster
!!$    DO i=1,nbcluster
!!$       IF (cluster_population(i)/=0) THEN
!!$          DO p3=i,nbcluster
!!$             normfro=0.0
!!$             DO j=ordaffperclus(i)+1,ordaffperclus(i+1)
!!$                DO q=ordaffperclus(p3)+1,ordaffperclus(p3+1)
!!$                   normfro=normfro+P1(j,q)**2
!!$                ENDDO
!!$             ENDDO
!!$             IF (i/=p3) THEN                
!!$                Frob(i,p3)=normfro/Frob(i,i)
!!$                ratio=ratio+Frob(i,p3)/nbcluster
!!$             ELSE
!!$                Frob(i,i)=normfro
!!$             ENDIF
!!$          ENDDO
!!$       ELSE
!!$          nbinfo=nbinfo-1
!!$       ENDIF
!!$    ENDDO
!!$    ! PRINT *, 
!!$    ! PRINT *,'norme de frobenius ',Frob
#if aff
    PRINT *,numproc,'nbinfo=', nbinfo,' nbcluster=',nbcluster
#endif
    !DEALLOCATE(ordaffperclus);DEALLOCATE(Frob)
    !DEALLOCATE(Z1), DEALLOCATE(Z2),DEALLOCATE(P);DEALLOCATE(P1)

    RETURN 
  END SUBROUTINE spectral_embedding


  !********************************************************
  !K-means
  SUBROUTINE kmeans_01 (dim_num,point_num,cluster_num,it_max, it_num, point, &
       cluster, cluster_center, cluster_population, cluster_energy, numproc)

    !*****************************************************************************80
    !
    !! KMEANS_01 applies the K-Means algorithm.
    !
    !  Discussion:
    !
    !    Given a matrix of POINT_NUM observations on DIM_NUM variables, the
    !    observations are to be ALLOCATEd to CLUSTER_NUM clusters in such 
    !    a way that the within-cluster sum of squares is minimized.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    18 November 2004
    !
    !  Author:
    !
    !    FORTRAN77 original version by David Sparks
    !    FORTRAN90 version by John Burkardt
    !
    !  Reference:
    !
    !    David Sparks,
    !    Algorithm AS 58: 
    !    Euclidean Cluster Analysis,
    !    Applied Statistics,
    !    Volume 22, Number 1, 1973, pages 126-130.
    !
    !  Parameters:
    !
    !    Input, INTEGER ( KIND = 4 ) DIM_NUM, the number of spatial dimensions.
    !
    !    Input, INTEGER ( KIND = 4 ) POINT_NUM, the number of points.
    !
    !    Input, INTEGER ( KIND = 4 ) CLUSTER_NUM, the number of clusters.
    !
    !    Input, INTEGER ( KIND = 4 ) IT_MAX, the maximum number of iterations.
    !
    !    Output, INTEGER ( KIND = 4 ) IT_NUM, the number of iterations taken.
    !
    !    Input, DOUBLE PRECISION POINT(DIM_NUM,POINT_NUM), the points.
    !
    !    Output, INTEGER ( KIND = 4 ) CLUSTER(POINT_NUM), indicates which cluster
    !    each point belongs to.
    !
    !    Input/output, DOUBLE PRECISION CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
    !    the cluster centers.
    !
    !    Output, INTEGER ( KIND = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number 
    !    of points in each cluster.
    !
    !    Output, DOUBLE PRECISION CLUSTER_ENERGY(CLUSTER_NUM), the 
    !    cluster energies.
    !
    IMPLICIT NONE
    INTEGER ( KIND = 4 ) cluster_num
    INTEGER ( KIND = 4 ) dim_num
    INTEGER ( KIND = 4 ) point_num
    INTEGER ( KIND = 4 ) cluster(point_num)
    DOUBLE PRECISION cluster_center(dim_num,cluster_num)
    DOUBLE PRECISION stockcenter(dim_num,cluster_num)
    !DOUBLE PRECISION diffcenter(dim_num)
    DOUBLE PRECISION listnorm(point_num,cluster_num)
    DOUBLE PRECISION cluster_energy(cluster_num),stockenergy(cluster_num)
    INTEGER ( KIND = 4 ) cluster_population(cluster_num)
    INTEGER ( KIND = 4 ) stockpopulation(cluster_num)
    !DOUBLE PRECISION dc,de,epsilon
    !DOUBLE PRECISION f(point_num)
    INTEGER ( KIND = 4 ) i
    !INTEGER ( KIND = 4 ) il
    !INTEGER ( KIND = 4 ) diffpopulation
    INTEGER ( KIND = 4 ) it_max
    INTEGER ( KIND = 4 ) it_num
    INTEGER ( KIND = 4 ) j
    INTEGER ( KIND = 4 ) k
    !INTEGER ( KIND = 4 ) list(1)
    DOUBLE PRECISION point(dim_num,point_num)
    INTEGER ( KIND = 4 ) swap
    INTEGER  :: ok,p,ok2
    DOUBLE PRECISION :: val,valmax,seuil,norme !,diffenergy
    INTEGER :: cluster_id(cluster_num)

    INTEGER :: numproc
    CHARACTER*30 :: files, num

    it_num = 0
    !
    !  Idiot checks.
    !
    IF ( cluster_num < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  CLUSTER_NUM < 1.0'
       STOP
    ENDIF

    IF ( dim_num < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  DIM_NUM < 1.0'
       STOP
    ENDIF

    IF ( point_num < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  POINT_NUM < 1.0'
       STOP
    ENDIF
    !
    !  For each observation, calculate the distance from each cluster
    !  center, and assign to the nearest.
    !

    !
    !  Assign one point to each cluster center.
    !
    cluster_center(:,1) = point(:,1)
    cluster_id(:)=0; cluster_id(1)=1
    p=2
    seuil=0.4
    !do i=2,cluster_num
    !   ok=0
    !   DO WHILE (ok==0)
    !      valmax=2.0*seuil
    !      DO j=1,i-1
    !         val=0.0
    !         DO k=1,dim_num
    !            val=max(val,abs(cluster_center(k,j)-point(k,p)))
    !         ENDDO
    !         valmax=min(valmax,val)
    !      ENDDO
    !      IF (valmax>seuil) ok=1
    !      p=p+1
    !   ENDDO
    !   cluster_center(:,i)=point(:,p-1)
    !ENDDO
#if aff
PRINT *, 'recherche des centres'
#endif
    DO i = 2, cluster_num
       ok=0
       DO WHILE(ok==0)
          valmax=2.0*seuil
          !recherche si le point est deja utilise dans comme centre
          ok2=0
          DO j=1,i-1
             IF (cluster_id(j)==p) ok2=1
          ENDDO
          !si point pas centre, teste par rapport au seuil
          IF (ok2==0) THEN
             DO j=1,i-1
                val=0.0; norme=0.0
                DO k=1,dim_num
                   val=max(val,abs(cluster_center(k,j)-point(k,p)))
                  !norme=norme+abs(cluster_center(k,j)-point(k,p))**2
                ENDDO
               !norme=sqrt(norme)
                !valmax=min(valmax,val/norme)
                valmax=min(val,valmax)
               !if (valmax>seuil) ok=1
             ENDDO
             IF (valmax>=seuil) ok=1
          ENDIF
         p=p+1
!if (seuil<=1e-4) THEN
!do i1=p+1,cluster_num

         !abaisse le seuil si pas assez de centre sont trouves
         IF ((p>point_num).AND.(ok==0)) THEN 
            seuil=0.9*seuil
#if aff
            PRINT *,'abaisse seuil :',seuil
#endif
            p=1
          ENDIF
       ENDDO
       p=p-1
       cluster_center(:,i)=point(:,p)
       cluster_id(i)=p
    ENDDO
#if aff
   PRINT *,'centres initiaux',p
#endif
!cluster_center

!!! boucle            
    it_num = 0
    swap=1
    cluster(:)=1
    DO WHILE ((it_num<it_max).AND.(swap/=0))
       it_num = it_num + 1
       swap=0
       DO i=1,cluster_num
          stockenergy(i)=cluster_energy(i);
          stockpopulation(i)=cluster_population(i);
          DO j=1,dim_num
             stockcenter(j,i)=cluster_center(j,i)
          ENDDO
       ENDDO

       !! Calcul de toutes les distances
       cluster_population(1:cluster_num) = 1
       ! ALLOCATE(listnorm(point_num,cluster_num));
       listnorm(:,:)=0.0
       DO i=1,point_num
          DO j=1,cluster_num
             DO k=1,dim_num
                listnorm(i,j)=listnorm(i,j)+(point(k,i)-cluster_center(k,j))**2
             ENDDO
          ENDDO
       ENDDO

       !!assignation par rapport au min des distances
       cluster_population(:)=0.0
       DO i=1,point_num
          DO j=1,cluster_num
             IF (listnorm(i,j)<listnorm(i,cluster(i))) THEN
                cluster(i)=j
                swap=swap+1
             ENDIF
          ENDDO
          cluster_energy(cluster(i))=cluster_energy(cluster(i))&
               +listnorm(i,cluster(i))
          cluster_population(cluster(i))=cluster_population(cluster(i))+1
       ENDDO

       !WRITE(num,*) numproc
       !files = 'kmeans.'//trim(adjustl(num))
       !OPEN(FILE=files, UNIT=15, ACCESS="APPEND")
       !WRITE (15,*) numproc,'iteration',it_num,'test clustering',cluster(1:3)
       !CALL FLUSH(15)
       !CLOSE(15)

       !PRINT *,numproc,'nbre de permut',swap

       !! mise ÃÂÃÂ  jour des centres
       cluster_center(:,:)=0.0
       DO j=1,point_num
          i=cluster(j) 
          DO k=1,dim_num
             cluster_center(k,i)=cluster_center(k,i)+point(k,j)
          ENDDO
       ENDDO
       DO i=1,cluster_num
          cluster_center(:,i)=cluster_center(:,i)/cluster_population(i)
       ENDDO



    ENDDO










!!$
!!$    DO i = 1, point_num
!!$       DO j = 1, cluster_num
!!$          cluster_energy(j) = sum ( &
!!$               ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
!!$       ENDDO
!!$       list = minloc ( cluster_energy(1:cluster_num) )
!!$       cluster(i) = list(1)
!!$    ENDDO
!!$    !
!!$    !  Determine the cluster population counts.
!!$    !
!!$    cluster_population(1:cluster_num) = 0
!!$
!!$    DO i = 1, point_num
!!$       j = cluster(i)
!!$       cluster_population(j) = cluster_population(j) + 1
!!$    ENDDO
!!$    !
!!$    !  Calculate the mean and sum of squares for each cluster.
!!$    !
!!$
!!$    DO i = 1, point_num
!!$       j = cluster(i)
!!$       cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
!!$            + point(1:dim_num,i)
!!$    ENDDO
!!$
!!$    DO i = 1, cluster_num
!!$       IF ( 0 < cluster_population(i) ) THEN
!!$          cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
!!$               / REAL ( cluster_population(i), KIND = 8 )
!!$       ENDIF
!!$    ENDDO
!!$    !
!!$    !  Set the point energies.
!!$    !
!!$    f(1:point_num) = 0.0
!!$
!!$    DO i = 1, point_num
!!$       j = cluster(i)
!!$       f(i) = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
!!$    ENDDO
!!$    !
!!$    !  Set the cluster energies.
!!$    !
!!$    cluster_energy(1:cluster_num) = 0.0
!!$
!!$    DO i = 1, point_num
!!$       j = cluster(i)
!!$       cluster_energy(j) = cluster_energy(j) + f(i)
!!$    ENDDO
!!$    !
!!$    !  Adjust the point energies by a weight factor.
!!$    !
!!$    DO i = 1, point_num
!!$       j = cluster(i)
!!$       IF ( 1 < cluster_population(j) ) THEN
!!$          f(i) = f(i) * REAL ( cluster_population(j), KIND = 8 ) &
!!$               / REAL ( cluster_population(j) - 1, KIND = 8 )
!!$       ENDIF
!!$    ENDDO
!!$    !
!!$    !  Examine each observation in turn to see IF it should be
!!$    !  reassigned to a different cluster.
!!$    !
!!$    it_num = 0
!!$    DO WHILE ( it_num < it_max )
!!$       it_num = it_num + 1
!!$       swap = 0
!!$       DO i = 1, point_num
!!$          il = cluster(i)
!!$          ir = il
!!$          IF ( cluster_population(il) <= 1 ) THEN
!!$             cycle
!!$          ENDIF
!!$          dc = f(i)
!!$          DO j = 1, cluster_num
!!$             IF ( j /= il ) THEN
!!$                de = sum ( &
!!$                     ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 ) &
!!$                     * REAL ( cluster_population(j), KIND = 8 ) &
!!$                     / REAL ( cluster_population(j) + 1, KIND = 8 )
!!$                IF ( de < dc ) THEN
!!$                   dc = de
!!$                   ir = j
!!$                ENDIF
!!$             ENDIF
!!$          ENDDO
!!$          !
!!$          !  If the lowest value was obtained by staying in the current cluster,
!!$          !  THEN cycle.
!!$          !
!!$          IF ( ir == il ) THEN
!!$             cycle
!!$          ENDIF
!!$          !
!!$          !  Reassign the point from cluster IL to cluster IR.
!!$          !
!!$          cluster_center(1:dim_num,il) = &
!!$               ( cluster_center(1:dim_num,il) &
!!$               * REAL ( cluster_population(il), KIND = 8 ) &
!!$               - point(1:dim_num,i) ) / REAL ( cluster_population(il) - 1, KIND = 8 )
!!$
!!$          cluster_center(1:dim_num,ir) = &
!!$               ( cluster_center(1:dim_num,ir) &
!!$               * REAL ( cluster_population(ir), KIND = 8 ) &
!!$               + point(1:dim_num,i) ) / REAL ( cluster_population(ir) + 1, KIND = 8 )
!!$
!!$          cluster_energy(il) = cluster_energy(il) - f(i)
!!$          cluster_energy(ir) = cluster_energy(ir) + dc
!!$          cluster_population(ir) = cluster_population(ir) + 1
!!$          cluster_population(il) = cluster_population(il) - 1
!!$
!!$          cluster(i) = ir
!!$          !
!!$          !  Adjust the value of F for points in clusters IL and IR.
!!$          !
!!$          DO j = 1, point_num
!!$             k = cluster(j)
!!$             IF ( k == il .OR. k == ir ) THEN
!!$                f(j) = sum ( &
!!$                     ( point(1:dim_num,j) - cluster_center(1:dim_num,k) )**2 )
!!$                IF ( 1 < cluster_population(k) ) THEN
!!$                   f(j) = f(j) * REAL ( cluster_population(k), KIND = 8 ) &
!!$                        / ( REAL ( cluster_population(k) - 1, KIND = 8 ) )
!!$                ENDIF
!!$             ENDIF
!!$          ENDDO
!!$          swap = swap + 1
!!$       ENDDO
!!$       !
!!$       !  Exit IF no reassignments were made during this iteration.
!!$       !
!!$       IF ( swap == 0 ) THEN
!!$          exit
!!$       ENDIF
!!$    ENDDO
!!$    !
!!$    !  Compute the cluster energies.
!!$    !
!!$    cluster_energy(1:cluster_num) = 0.0
!!$    DO i = 1, point_num
!!$       j = cluster(i)
!!$       cluster_energy(j) = cluster_energy(j) + sum ( &
!!$            ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
!!$    ENDDO
!!$

    RETURN
  END SUBROUTINE kmeans_01



END MODULE module_embed
