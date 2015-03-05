MODULE module_embed
  USE module_structure
CONTAINS



  SUBROUTINE spectral_embedding(nbcluster, n, Z, A, ratio,cluster, &
       cluster_center, cluster_population, cluster_energy, nbinfo, numproc, &
       ratiomoy, ratiorij, ratiorii)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z ! matrice des vecteurs propres
    INTEGER :: n
    INTEGER :: nbcluster ! nbre de cluster
    INTEGER :: numproc
    
    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: cluster_center ! centre des nbclusters clusters
    DOUBLE PRECISION, DIMENSION(:), POINTER :: cluster_energy ! somme des energies par cluster
    DOUBLE PRECISION :: ratio ! max des ration de frob sur matrice aff reordonnancee suivant
    DOUBLE PRECISION :: ratiomoy
    DOUBLE PRECISION :: ratiorii
    DOUBLE PRECISION :: ratiorij
    INTEGER, DIMENSION(:), POINTER :: cluster ! appartenance des clusters
    INTEGER, DIMENSION(:), POINTER :: cluster_population ! nbre de points par cluster
    INTEGER :: nbinfo
    
    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Frob
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z1
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z2
    DOUBLE PRECISION, DIMENSION(:), POINTER :: Z3
    DOUBLE PRECISION :: ratiomin
    DOUBLE PRECISION :: test
    INTEGER, DIMENSION(:,:), POINTER :: clustercorresp
    INTEGER :: i
    INTEGER :: it_max
    INTEGER :: it_num
    INTEGER :: j
    INTEGER :: k
    INTEGER :: ki
    INTEGER :: kj
    INTEGER :: nbmax
    INTEGER :: ni
    INTEGER :: nj
    LOGICAL :: ok
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
    ALLOCATE(cluster(n))
    ALLOCATE(cluster_center(nbcluster,nbcluster))
    ALLOCATE(cluster_population(nbcluster))
    ALLOCATE(cluster_energy(nbcluster))
    ALLOCATE(Z1(n,nbcluster))
    ALLOCATE(Z2(nbcluster,n))
    ALLOCATE(Z3(n))
    Z3(:)=0.0

    DO i=1,n
       DO j=1,nbcluster
          Z1(i,j)=Z(i,j)
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

    it_max=n*n

    CALL apply_kmeans( nbcluster, n, nbcluster, it_max, it_num,Z2,&
         cluster, cluster_center, cluster_population, cluster_energy, &
         numproc)

    !*****************************
    ! Mesure de qualite
    nbmax=0
    DO i=1,nbcluster
       nbmax=max(nbmax,cluster_population(i))
    ENDDO
    ALLOCATE(clustercorresp(nbcluster,nbmax))
    clustercorresp(:,:)=0
    DO i=1,n
       j=cluster(i)
       ok=.FALSE.
       k=1
       DO WHILE(.NOT. ok)
          IF (clustercorresp(j,k)==0) THEN
             ok=.TRUE.
          ELSE
             k=k+1
          ENDIF
       ENDDO
       clustercorresp(j,k)=i
    ENDDO
    ALLOCATE(Frob(nbcluster,nbcluster))
    Frob(:,:)=0.0 
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
    ratio=0.0
    ratiomin=1.D+16
    ratiorii=0.0
    ratiorij=0.0
    nbinfo=nbcluster
    DO i=1,nbcluster
       IF ((cluster_population(i)/=0).AND.(Frob(i,i)/=0)) THEN
          DO j=1,nbcluster
             IF (i/=j) THEN
                ratio=ratio+Frob(i,j)/Frob(i,i)
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
    ENDDO
    DEALLOCATE(Frob)

#if aff
    PRINT *,numproc,'nbinfo=', nbinfo,' nbcluster=',nbcluster
#endif

    RETURN 
  END SUBROUTINE spectral_embedding



  SUBROUTINE apply_kmeans(dim_num, point_num, cluster_num, it_max, it_num, point, &
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

    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: point_num ! the number of points
                !TODO : a reflechir sur l'ordre de declaration
    DOUBLE PRECISION :: point (dim_num, point_num) ! the points
    INTEGER :: cluster_num ! the number of clusters
    INTEGER :: dim_num ! the number of spatial dimensions
    INTEGER :: it_max ! the maximum number of iterations

    !=== IN/OUT ===
    DOUBLE PRECISION :: cluster_center (dim_num, cluster_num) ! the cluster centers

    !====  OUT ====
    DOUBLE PRECISION :: cluster_energy (cluster_num) ! the cluster energies
    INTEGER :: it_num ! the number of iterations taken
    INTEGER :: cluster (point_num) ! indicates which cluster each point belongs to
    INTEGER :: cluster_population (cluster_num) ! the number of points in each cluster

    !== USELESS ===
    INTEGER :: numproc ! TODO : etudier si garder ou pas
    
    !#### Variables  ####
    DOUBLE PRECISION :: listnorm (point_num, cluster_num)
    DOUBLE PRECISION :: stockcenter (dim_num, cluster_num)
    DOUBLE PRECISION :: stockenergy (cluster_num)
    DOUBLE PRECISION :: norme
    DOUBLE PRECISION :: seuil
    DOUBLE PRECISION :: val
    DOUBLE PRECISION :: valmax
    INTEGER :: cluster_id (cluster_num)
    INTEGER :: stockpopulation (cluster_num)
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: swap
    INTEGER :: p
    LOGICAL :: ok
    LOGICAL :: ok2
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################   
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
    cluster_id(:)=0
    cluster_id(1)=1
    p=2
    seuil=0.4
#if aff
PRINT *, 'recherche des centres'
#endif
    DO i = 2, cluster_num
       ok=.FALSE.
       DO WHILE(.NOT. ok)
          valmax=2.0*seuil
          ! Test if the point is already used as center
          ok2=.FALSE.
          DO j=1,i-1
             IF (cluster_id(j)==p) ok2=.TRUE.
          ENDDO
          ! If the point is not a center, test against the threshold
          IF (.NOT. ok2) THEN
             DO j=1,i-1
                val=0.0
                norme=0.0
                DO k=1,dim_num
                   val=max(val,abs(cluster_center(k,j)-point(k,p)))
                ENDDO
                valmax=min(val,valmax)
             ENDDO
             IF (valmax>=seuil) ok=.TRUE.
          ENDIF
         p=p+1

         ! Lower the threshold if not enough centers found
         IF ((p>point_num).AND.(.NOT. ok)) THEN 
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
    it_num = 0
    swap=1
    cluster(:)=1
    DO WHILE ((it_num<it_max).AND.(swap/=0))
       it_num = it_num + 1
       swap=0
       DO i=1,cluster_num
          stockenergy(i)=cluster_energy(i)
          stockpopulation(i)=cluster_population(i)
          DO j=1,dim_num
             stockcenter(j,i)=cluster_center(j,i)
          ENDDO
       ENDDO

       ! Computing of the distances
       cluster_population(1:cluster_num) = 1
       listnorm(:,:)=0.0
       DO i=1,point_num
          DO j=1,cluster_num
             DO k=1,dim_num
                listnorm(i,j)=listnorm(i,j)+(point(k,i)-cluster_center(k,j))**2
             ENDDO
          ENDDO
       ENDDO

       ! Allocation related to the minimum of the distances
       cluster_population(:)=0
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

       ! Update of centers
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

    RETURN
  END SUBROUTINE apply_kmeans



END MODULE module_embed
