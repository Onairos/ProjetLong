MODULE module_embed
  USE module_structure
CONTAINS



  SUBROUTINE apply_spectral_embedding(nb_clusters, n, Z, A, ratio,clusters, &
       clusters_centers, points_by_clusters, clusters_energies, nb_info, proc_id, &
       ratio_moy, ratio_rij, ratio_rii)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z ! matrice des vecteurs propres
    INTEGER :: n
    INTEGER :: nb_clusters ! nbre de clusters
    INTEGER :: proc_id
    
    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: clusters_centers ! centre des nbclusters clusters
    DOUBLE PRECISION, DIMENSION(:), POINTER :: clusters_energies ! somme des energies par cluster
    DOUBLE PRECISION :: ratio ! max des ration de frob sur matrice aff reordonnancee suivant
    DOUBLE PRECISION :: ratio_moy
    DOUBLE PRECISION :: ratio_rii
    DOUBLE PRECISION :: ratio_rij
    INTEGER, DIMENSION(:), POINTER :: clusters ! appartenance des clusters
    INTEGER, DIMENSION(:), POINTER :: points_by_clusters ! nbre de points par cluster
    INTEGER :: nb_info
    
    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Frob
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z1
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z2
    DOUBLE PRECISION, DIMENSION(:), POINTER :: Z3
    DOUBLE PRECISION :: ratiomin
    DOUBLE PRECISION :: test
    INTEGER, DIMENSION(:,:), POINTER :: clustercorresp
    INTEGER :: i
    INTEGER :: nb_iter_max
    INTEGER :: nb_iter
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
    ALLOCATE(clusters(n))
    ALLOCATE(clusters_centers(nb_clusters,nb_clusters))
    ALLOCATE(points_by_clusters(nb_clusters))
    ALLOCATE(clusters_energies(nb_clusters))
    ALLOCATE(Z1(n,nb_clusters))
    ALLOCATE(Z2(nb_clusters,n))
    ALLOCATE(Z3(n))
    Z3(:)=0.0

    DO i=1,n
       DO j=1,nb_clusters
          Z1(i,j)=Z(i,j)
          Z3(i)=Z3(i)+Z1(i,j)**2
       ENDDO
    ENDDO

    DO i=1,n
       test=0.0
       DO j=1,nb_clusters
          Z2(j,i)=Z1(i,j)/(sqrt(Z3(i)))
          test=test+Z2(j,i)**2
       ENDDO
    ENDDO

    nb_iter_max=n*n

    CALL apply_kmeans( nb_clusters, n, nb_clusters, nb_iter_max, nb_iter,Z2,&
         clusters, clusters_centers, points_by_clusters, clusters_energies, &
         proc_id)

    !*****************************
    ! Mesure de qualite
    nbmax=0
    DO i=1,nb_clusters
       nbmax=max(nbmax,points_by_clusters(i))
    ENDDO
    ALLOCATE(clustercorresp(nb_clusters,nbmax))
    clustercorresp(:,:)=0
    DO i=1,n
       j=clusters(i)
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
    ALLOCATE(Frob(nb_clusters,nb_clusters))
    Frob(:,:)=0.0 
    DO i=1,nb_clusters
       DO j=1,nb_clusters
          DO ki=1,points_by_clusters(i)
             ni=clustercorresp(i,ki)
             DO kj=1,points_by_clusters(j)
                nj=clustercorresp(j,kj)
                Frob(i,j)=Frob(i,j)+A(ni,nj)**2
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(clustercorresp)
    ratio=0.0
    ratiomin=1.D+16
    ratio_rii=0.0
    ratio_rij=0.0
    nb_info=nb_clusters
    DO i=1,nb_clusters
       IF ((points_by_clusters(i)/=0).AND.(Frob(i,i)/=0)) THEN
          DO j=1,nb_clusters
             IF (i/=j) THEN
                ratio=ratio+Frob(i,j)/Frob(i,i)
                ratio_moy=ratio_moy+Frob(i,j)/Frob(i,i)
                ratio_rij=ratio_rij+Frob(i,j)
                ratio_rii=ratio_rii+Frob(i,i)
                ratiomin=min(ratiomin,Frob(i,j)/Frob(i,i))
             ENDIF
          ENDDO
       ELSE
          nb_info=nb_info-1
       ENDIF
       ratio_rij=ratio_rij*2/(nb_clusters*(nb_clusters-1))
    ENDDO
    DEALLOCATE(Frob)

#if aff
    PRINT *, 'DEBUG : process n', proc_id,' : nb_info=', nb_info, ' nb_clusters=', nb_clusters
#endif

    RETURN 
  END SUBROUTINE apply_spectral_embedding



  SUBROUTINE apply_kmeans(dim, nb_points, nb_clusters, nb_iter_max, nb_iter, points, &
       clusters, clusters_centers, points_by_clusters, clusters_energies, proc_id)

    !*****************************************************************************80
    !
    !! KMEANS_01 applies the K-Means algorithm.
    !
    !  Discussion:
    !
    !    Given a matrix of NB_POINTS observations on DIMENSION variables, the
    !    observations are to be ALLOCATEd to NB_CLUSTERS clusters in such 
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
    INTEGER :: dim ! the number of spatial dimensions
    INTEGER :: nb_points ! the number of points + TODO reorganisation
    DOUBLE PRECISION :: points (dim, nb_points) ! the points
    INTEGER :: nb_clusters ! the number of clusters
    INTEGER :: nb_iter_max ! the maximum number of iterations

    !=== IN/OUT ===
    DOUBLE PRECISION :: clusters_centers (dim, nb_clusters) ! the cluster centers

    !====  OUT ====
    DOUBLE PRECISION :: clusters_energies (nb_clusters) ! the cluster energies
    INTEGER :: nb_iter ! the number of iterations taken
    INTEGER :: clusters (nb_points) ! indicates which cluster each point belongs to
    INTEGER :: points_by_clusters (nb_clusters) ! the number of points in each cluster

    !== USELESS ===
    INTEGER :: proc_id ! TODO : etudier si garder ou pas : le garder finalement
    
    !#### Variables  ####
    DOUBLE PRECISION :: listnorm (nb_points, nb_clusters)
    DOUBLE PRECISION :: stockcenter (dim, nb_clusters)
    DOUBLE PRECISION :: stockenergy (nb_clusters)
    DOUBLE PRECISION :: norme
    DOUBLE PRECISION :: seuil
    DOUBLE PRECISION :: val
    DOUBLE PRECISION :: valmax
    INTEGER :: cluster_id (nb_clusters)
    INTEGER :: stockpopulation (nb_clusters)
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
    nb_iter = 0
    !
    !  Idiot checks.
    !
    IF ( nb_clusters < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  NB_CLUSTERS < 1.0'
       STOP
    ENDIF

    IF ( dim < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  DIMENSION < 1.0'
       STOP
    ENDIF

    IF ( nb_points < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  NB_POINTS < 1.0'
       STOP
    ENDIF
    !
    !  For each observation, calculate the distance from each cluster
    !  center, and assign to the nearest.
    !

    !
    !  Assign one point to each cluster center.
    !
    clusters_centers(:,1) = points(:,1)
    cluster_id(:)=0
    cluster_id(1)=1
    p=2
    seuil=0.4
#if aff
PRINT *, 'DEBUG : searching centers'
#endif
    DO i = 2, nb_clusters
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
                DO k=1,dim
                   val=max(val,abs(clusters_centers(k,j)-points(k,p)))
                ENDDO
                valmax=min(val,valmax)
             ENDDO
             IF (valmax>=seuil) ok=.TRUE.
          ENDIF
         p=p+1

         ! Lower the threshold if not enough centers found
         IF ((p>nb_points).AND.(.NOT. ok)) THEN 
            seuil=0.9*seuil
#if aff
            PRINT *, 'DEBUG : Lower threshold : ', seuil
#endif
            p=1
          ENDIF
       ENDDO
       p=p-1
       clusters_centers(:,i)=points(:,p)
       cluster_id(i)=p
    ENDDO
#if aff
   PRINT *, 'DEBUG : initial centers : ', p
#endif
    nb_iter = 0
    swap=1
    clusters(:)=1
    DO WHILE ((nb_iter<nb_iter_max).AND.(swap/=0))
       nb_iter = nb_iter + 1
       swap=0
       DO i=1,nb_clusters
          stockenergy(i)=clusters_energies(i)
          stockpopulation(i)=points_by_clusters(i)
          DO j=1,dim
             stockcenter(j,i)=clusters_centers(j,i)
          ENDDO
       ENDDO

       ! Computing of the distances
       points_by_clusters(1:nb_clusters) = 1
       listnorm(:,:)=0.0
       DO i=1,nb_points
          DO j=1,nb_clusters
             DO k=1,dim
                listnorm(i,j)=listnorm(i,j)+(points(k,i)-clusters_centers(k,j))**2
             ENDDO
          ENDDO
       ENDDO

       ! Allocation related to the minimum of the distances
       points_by_clusters(:)=0
       DO i=1,nb_points
          DO j=1,nb_clusters
             IF (listnorm(i,j)<listnorm(i,clusters(i))) THEN
                clusters(i)=j
                swap=swap+1
             ENDIF
          ENDDO
          clusters_energies(clusters(i))=clusters_energies(clusters(i))&
               +listnorm(i,clusters(i))
          points_by_clusters(clusters(i))=points_by_clusters(clusters(i))+1
       ENDDO

       ! Update of centers
       clusters_centers(:,:)=0.0
       DO j=1,nb_points
          i=clusters(j) 
          DO k=1,dim
             clusters_centers(k,i)=clusters_centers(k,i)+points(k,j)
          ENDDO
       ENDDO
       DO i=1,nb_clusters
          clusters_centers(:,i)=clusters_centers(:,i)/points_by_clusters(i)
       ENDDO



    ENDDO

    RETURN
  END SUBROUTINE apply_kmeans



END MODULE module_embed
