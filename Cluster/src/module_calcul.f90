MODULE module_calcul
  USE module_structure
  USE module_solve
  USE module_embed
CONTAINS


  SUBROUTINE get_sigma(partitioned_data, sigma)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: partitioned_data

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
    DO i1=1,partitioned_data%nb
       DO j1=i1+1,partitioned_data%nb
          norme=0.0
          DO k1=1,partitioned_data%dim
             norme=norme+&
                  (partitioned_data%point(i1)%coord(k1)-partitioned_data%point(j1)%coord(k1))**2
          ENDDO
          sigma=max(sigma,sqrt(norme))
       ENDDO
    ENDDO
    sigma=sigma/(2*exp(log(float(partitioned_data%nb))*(1.0/float(partitioned_data%dim))))
    ! Safety
    IF (sigma==0.0) sigma=1.0
    RETURN
  END SUBROUTINE get_sigma


  SUBROUTINE get_sigma_interface(proc_id, partitioned_data, sigma, bounds, partitioning, epsilon)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: partitioned_data
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bounds
    DOUBLE PRECISION :: epsilon
    INTEGER, DIMENSION(:), POINTER :: partitioning
    INTEGER :: proc_id

    !====  OUT ====
    DOUBLE PRECISION :: sigma

    !#### Variables  ####
    INTEGER, DIMENSION(:,:), POINTER :: tableau
    INTEGER, DIMENSION(:), POINTER :: partitioning_tmp
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
#if aff
    ! Number of partitionings
    nb=1
    DO i=1,partitioned_data%dim
       nb=nb*partitioning(i)
    ENDDO
    ! Creation of partitioning
    ALLOCATE(tableau(nb,0:partitioned_data%dim))
    ALLOCATE(partitioning_tmp(partitioned_data%dim))
    partitioning_tmp(:)=1
    DO i=1,nb
       DO j=1,partitioned_data%dim
          tableau(i,j)=partitioning_tmp(j)
       ENDDO
       partitioning_tmp(1)=partitioning_tmp(1)+1
       k=1
       DO WHILE(partitioning_tmp(k)>partitioning(k))
          partitioning_tmp(k)=1
          IF (k<partitioned_data%dim) partitioning_tmp(k+1)=partitioning_tmp(k+1)+1
       ENDDO
    ENDDO
    DEALLOCATE(partitioning_tmp)
    ! Value of sigma
    sigma0=0.0
    DO i=1,nb
       volext=1.0
       volint=1.0
       DO j=1,partitioned_data%dim
          k=tableau(i,j)
          long=bounds(j,k,2)-bounds(j,k,1)
          volext=volext*long
          volint=volint*max(0.0D1,long-2.0*epsilon)
       ENDDO
       sigma0=sigma0+volext-volint
    ENDDO
    DEALLOCATE(tableau)
    ! Computing of scale length
    sigma0=exp(1.0/float(partitioned_data%dim)*log(sigma0))
    ! Sigma computing
    sigma0=sigma0/(2.0*exp(log(float(partitioned_data%nb))*(1.0/float(partitioned_data%dim))))
    PRINT *, 'DEBUG : process n', proc_id, ' : value of computed sigma for interfacing : ', sigma0
#endif
    ! Sigma computing, global formula
    CALL get_sigma(partitioned_data,sigma)
#if aff
    PRINT *, 'DEBUG : process n', proc_id,' : value of sigma for interfacing', sigma
#endif
    RETURN
  END SUBROUTINE get_sigma_interface



FUNCTION poly_kernel( dataw, gam, delta )
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: dataw
    DOUBLE PRECISION ::gam 
    DOUBLE PRECISION :: delta

    !====  OUT  ====
    DOUBLE PRECISION, DIMENSION(dataw%nb,dataw%nb) :: poly_kernel
    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: K
    INTEGER :: n
    INTEGER :: i
    INTEGER :: j
    INTEGER :: d

    n=dataw%nb
    ALLOCATE(K(n,n))
    !ALLOCATE(polyKernel)
    K(:,:)=0.0

    DO i=1,n-1
      DO j=1,n-1
        DO d=1,dataw%dim
        K(i,j)=K(i,j)+dataw%point(i)%coord(d)*dataw%point(j)%coord(d)
        ENDDO 
        K(i,j)=(K(i,j)+gam)**delta
      ENDDO
    ENDDO
    poly_kernel=K
    RETURN
  END

    FUNCTION gaussian_kernel( dataw, sigma )
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: dataw
    DOUBLE PRECISION sigma


    !====  OUT  ====
    DOUBLE PRECISION, DIMENSION(dataw%nb,dataw%nb) :: gaussian_kernel

    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: K
    INTEGER n
    INTEGER i
    INTEGER j
    INTEGER d

    n=dataw%nb
    ALLOCATE(K(n,n))
    !ALLOCATE(gaussian_kernel)
    K(:,:)=0.0

    DO i=1,n-1
      DO j=i+1,n
        DO d=1,dataw%dim
        K(i,j)=K(i,j)+(dataw%point(i)%coord(d)-dataw%point(j)%coord(d))**2
        ENDDO
        K(i,j)=exp(- K(i,j)/(2*sigma**2))
        ! Symetry
        K(j,i)=K(i,j)
      ENDDO
    ENDDO
    gaussian_kernel=K
    RETURN
  END


!Stop when converged compute E = sum_N(sum_M( Indicatrice (xi E Ck)*||phi(xi)-mk||²))



SUBROUTINE apply_kernel_k_means(proc_id,nblimit,nb_clusters_opt,dataw,clust_param)
    IMPLICIT NONE

   INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
  
    INTEGER :: nb_clusters_opt
    INTEGER :: nblimit
    INTEGER :: proc_id
    TYPE(type_clustering_param) :: clust_param

    !=== IN/OUT ===
    TYPE(type_data) :: dataw


    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Variables  ####

    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Ker
    INTEGER :: it_max ! the maximum number of iterations
    DOUBLE PRECISION :: cluster_center (dataw%dim, dataw%nbclusters) ! the cluster centers
    DOUBLE PRECISION :: cluster_energy (dataw%nbclusters) ! the cluster energies
    INTEGER :: it_num ! the number of iterations taken
    INTEGER :: cluster (dataw%nb) ! indicates which cluster each point belongs to
    INTEGER :: cluster_population (dataw%nbclusters) ! the number of points in each cluster
    DOUBLE PRECISION :: listnorm (dataw%nb, dataw%nbclusters)
    DOUBLE PRECISION :: stockcenter (dataw%dim, dataw%nbclusters)
    DOUBLE PRECISION :: stockenergy (dataw%nbclusters)
    DOUBLE PRECISION :: norme
    DOUBLE PRECISION :: seuil
    DOUBLE PRECISION :: val
    DOUBLE PRECISION :: valmax
    INTEGER :: cluster_id (dataw%nbclusters)
    INTEGER :: stockpopulation (dataw%nbclusters)
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: l

   ! LOGICAL :: found
    DOUBLE PRECISION :: num1
    DOUBLE PRECISION :: den1
    DOUBLE PRECISION :: num2
    DOUBLE PRECISION :: den2

    INTEGER :: ok ! TODO: booleen ?
    INTEGER :: ok2 ! TODO: booleen ?
    INTEGER :: swap
    INTEGER :: p
    

    ALLOCATE(Ker(dataw%nb,dataw%nb))
    Ker(:,:)=0.0


    !###########################################      
    ! INSTRUCTIONS
    !###########################################   
    it_num = 0
    !
    !  Idiot checks.
    !
    IF ( dataw%nbclusters < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  CLUSTER_NUM < 1.0'
       STOP
    ENDIF

    IF ( dataw%dim < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  DIM_NUM < 1.0'
       STOP
    ENDIF

    IF ( dataw%nb < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  POINT_NUM < 1.0'
       STOP
    ENDIF

 
    IF (clust_param%kernelfunindex==0 .AND. (clust_param%gam <0.0 .OR. clust_param%delta<0.0) )  THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  GAMMA AND DELTA NOT INITIALIZED IN POLYNOMIAL KERNEL'
       STOP  
    

    ELSEIF (clust_param%kernelfunindex==1 .AND. clust_param%sigma<0.0) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  SIGMA NOT INITIALIZED IN GAUSSIAN KERNEL'
       STOP 
    ENDIF


    IF (clust_param%kernelfunindex==0) THEN
        Ker=poly_kernel( dataw, clust_param%gam, clust_param%delta)
    ELSEIF (clust_param%kernelfunindex==1) THEN
        Ker=gaussian_kernel(dataw, clust_param%sigma)
    ENDIF



    !
    !  For each observation, calculate the distance from each cluster
    !  center, and assign to the nearest.
    !

    !
    !  Assign one point to each cluster center.
    !
    cluster_center(:,1) = dataw%point(1)%coord(:) !point(:,1) %%
    cluster_id(:)=0; cluster_id(1)=1
    p=2
    seuil=0.4
!#if aff
PRINT *, 'recherche des centres'
!#endif
    DO i = 2, dataw%nbclusters
       ok=0
       DO WHILE(ok==0)
          valmax=2.0*seuil
          !recherche si le point est deja utilise dans comme centre
          ok2=0
          DO j=1,i-1
             IF (dataw%point(j)%cluster==p) ok2=1
          ENDDO
          !si point pas centre, teste par rapport au seuil
          IF (ok2==0) THEN
             DO j=1,i-1
                val=0.0; norme=0.0
                DO k=1,dataw%dim
                   val=max(val,abs(cluster_center(k,j)-dataw%point(p)%coord(k))) 
!VOIR SI CELA DOIT ÊTRE MODIFIE EN FONCTION DES KERNEL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ENDDO
                valmax=min(val,valmax)
             ENDDO
             IF (valmax>=seuil) ok=1
          ENDIF
         p=p+1

         !abaisse le seuil si pas assez de centre sont trouves
         IF ((p>dataw%nb).AND.(ok==0)) THEN 
            seuil=0.9*seuil
!#if aff
            PRINT *,'abaisse seuil :',seuil
!#endif
            p=1
          ENDIF
       ENDDO
       p=p-1
       cluster_center(:,i)= dataw%point(P)%coord(:) !point(:,p) 
       cluster_id(i)=p
    ENDDO
!#if aff
   PRINT *,'centres initiaux',p
!#endif
!cluster_center

!!! boucle            
    it_num = 0
    swap=1
    dataw%point(:)%cluster=1 !  cluster(:)=1
    DO WHILE ((it_num<it_max).AND.(swap/=0))
       it_num = it_num + 1
       swap=0
       DO i=1,dataw%nbclusters
          stockenergy(i)=cluster_energy(i);
          stockpopulation(i)=cluster_population(i);
          DO j=1,dataw%dim
             stockcenter(j,i)=cluster_center(j,i)
          ENDDO
       ENDDO

       !! Calcul de toutes les distances
       cluster_population(1:dataw%nbclusters) = 1
       listnorm(:,:)=0.0
       num1=0.0
       den1=0.0
       num2=0.0
       den2=0.0
       DO k=1,dataw%nbclusters
           DO i=1,dataw%nb
               DO j=1,dataw%nb
                   IF ( dataw%point(j)%cluster.EQ.k) THEN
                   num1=num1 + 2*(Ker(i,j))
                   den1=den1+1 
                   ENDIF
                   DO l=1,dataw%nb
                       IF ( dataw%point(j)%cluster.EQ.k .AND. dataw%point(l)%cluster.EQ.k) THEN
                       num2=num2 + Ker(j,l)
                       den2=den2+1
                       ENDIF
                   ENDDO
               ENDDO
           listnorm(i,k)= Ker(i,i)
               IF (den1.NE.0.0) THEN
               listnorm(i,k)=  listnorm(i,k) - (num1/den1)
               ENDIF
               IF (den2.NE.0.0) THEN
               listnorm(i,k)=  listnorm(i,k)  + (num2/den2)
               ENDIF
           ENDDO
       ENDDO
     

       !!assignation par rapport au min des distances
       cluster_population(:)=0
       DO i=1,dataw%nb
          DO j=1,dataw%nbclusters
             IF (listnorm(i,j)<listnorm(i,dataw%point(i)%cluster)) THEN
                dataw%point(i)%cluster=j
                swap=swap+1
             ENDIF
          ENDDO
          cluster_energy(dataw%point(i)%cluster)=cluster_energy(dataw%point(i)%cluster)&
               +listnorm(i,dataw%point(i)%cluster)
          cluster_population(dataw%point(i)%cluster)=cluster_population(dataw%point(i)%cluster)+1
       ENDDO

       !! mise a jour des centres
       cluster_center(:,:)=0.0
       DO j=1,dataw%nb
          i=dataw%point(j)%cluster 
          DO k=1,dataw%dim
             cluster_center(k,i)=cluster_center(k,i)+dataw%point(j)%coord(k)
          ENDDO
       ENDDO
       DO i=1,dataw%nbclusters
          cluster_center(:,i)=cluster_center(:,i)/cluster_population(i)
       ENDDO



    ENDDO

    RETURN
  END SUBROUTINE apply_kernel_k_means











  SUBROUTINE apply_spectral_clustering(proc_id, nb_clusters_max, nb_clusters_opt, partitioned_data, sigma,clust_param)
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    DOUBLE PRECISION :: sigma
    INTEGER :: nb_clusters_opt
    INTEGER :: nb_clusters_max
    INTEGER :: proc_id
    TYPE(type_clustering_param) :: clust_param

    !=== IN/OUT ===
    TYPE(type_data) :: partitioned_data

    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A2
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: clusters_centers
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z
    DOUBLE PRECISION, DIMENSION(:), POINTER :: clusters_energies
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
    INTEGER, DIMENSION(:), POINTER :: points_by_clusters
    INTEGER, DIMENSION(:), POINTER :: nb_info
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: n
    INTEGER :: nb
    INTEGER :: nb_clusters
    INTEGER :: nbproc !TODO : mettre en parametre et WTF faut-il faire car on lit une variable vide
    INTEGER :: nbvp
    INTEGER :: solver ! solveur au valeur propre => parametre de controle

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ! Matrix creation
    PRINT *, 'Process n', proc_id, ' : value of sigma : ', sigma
    n=partitioned_data%nb
    ! Forall i, A(i,i) = 0
    ALLOCATE(A(n,n))
    A(:,:)=0.0

    DO i=1,n-1
       DO j=i+1,n
          norme=0.0
          DO k=1,partitioned_data%dim
             norme=norme+(partitioned_data%point(i)%coord(k)-partitioned_data%point(j)%coord(k))**2
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
      PRINT *, 'Process n', proc_id, ' : Lapack solver'

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
      PRINT *, 'Process n', proc_id, ' : Arpack solver'

      nb = 2*nb_clusters_max
      nbvp = nb
      CALL solve_arpack_full(A, n, nb, W, Z)

    ENDIF

    t2 = MPI_WTIME()

    t_cons_vp = t2 - t1
    PRINT *, 'Process n', proc_id, ' : Time for eigen values construction : ', t_cons_vp

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
    IF ((nb_clusters_opt==0).AND.(n>2)) THEN
       ! Search of the best partitioning
       ALLOCATE(ratiomax(nb_clusters_max))
       ratiomax(:)=0
       ALLOCATE(ratiomin(nb_clusters_max))
       ratiomin(:)=0
       ALLOCATE(ratiomoy(nb_clusters_max))
       ratiomoy(:)=0
       ALLOCATE(ratiorii(nb_clusters_max))
       ratiorii(:)=0
       ALLOCATE(ratiorij(nb_clusters_max))
       ratiorij(:)=0

       ALLOCATE(nb_info(nb_clusters_max))
       nb_info(:)=0
       DO nb_clusters=2,min(n,nb_clusters_max)

          ALLOCATE(cluster(n))
          cluster(:)=0
          ALLOCATE(clusters_centers(nb_clusters,nb_clusters))
          clusters_centers(:,:)=0.0
          ALLOCATE(points_by_clusters(nb_clusters))
          points_by_clusters(:)=0
          ALLOCATE(clusters_energies(nb_clusters))
          clusters_energies(:)=0.0

          CALL apply_spectral_embedding(nb_clusters,n,Z,A,&
               ratiomax(nb_clusters),cluster,clusters_centers,points_by_clusters,&
               clusters_energies,nb_info(nb_clusters),proc_id,ratiomoy(nb_clusters), &
               ratiorij(nb_clusters),ratiorii(nb_clusters))


          DEALLOCATE(cluster)
          DEALLOCATE(clusters_centers)
          DEALLOCATE(clusters_energies)
          DEALLOCATE(points_by_clusters)
       ENDDO

#if aff
PRINT *, 'DEBUG : Frobenius ratio'
#endif
       ! Norm of frobenius ratio
       ratio=ratiomax(nb_clusters_max)
       partitioned_data%nbclusters=nb_clusters_max
       ratio1=0.0
       ratio2=1e+10
       DO i=2,nb_clusters_max
          IF ((proc_id==0).AND.(nbproc>1)) THEN 
             seuilrij=1e-1
          ELSE
             seuilrij=1e-4 !TODO : garder que celui-ci
          ENDIF
          IF ((ratiorii(i)>=0.95*ratio1).AND.(ratiorij(i)-ratio2<=seuilrij)) THEN  
             partitioned_data%nbclusters=i
             ratio1=ratiorii(i)
             ratio2=ratiorij(i)
          ENDIF
       ENDDO

    ELSEIF ((nb_clusters_opt==1).AND.(n>nb_clusters_opt)) THEN
       ! Test with an imposed cluster
       ALLOCATE(nb_info(nb_clusters_opt))
       nb_info(:)=0
       ALLOCATE(ratiomin(1))
       ratiomin(:)=0.0
       partitioned_data%nbclusters=nb_clusters_opt
    ELSE
       ! Case of a domain with less points than nb_clusters_opt or only one point
       ALLOCATE(nb_info(n))
       nb_info(:)=0
       ALLOCATE(ratiomin(1))
       ratiomin(:)=0.0
       partitioned_data%nbclusters=n
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
    ! Case of nb_clusters==1
    IF (partitioned_data%nbclusters==2) THEN
       PRINT *, 'Ratio difference : ', ratiorij(2)/ratiorii(2)
       IF (ratiomax(2)>=0.6) THEN 
          partitioned_data%nbclusters=1
       ELSE 
          partitioned_data%nbclusters=2
       ENDIF
    ENDIF
#if aff
    PRINT *, 'DEBUG : Process n', proc_id,' : final cluster got : ', partitioned_data%nbclusters
#endif

    ! Final clustering computing
    IF (partitioned_data%nbclusters>1) THEN
       CALL apply_spectral_embedding(partitioned_data%nbclusters,n,Z,A,ratio,cluster,&
            clusters_centers,points_by_clusters,clusters_energies,&
            nb_info(partitioned_data%nbclusters),proc_id,ratiomin(1),ratiorij(1),ratiorii(1))
       DO i=1,partitioned_data%nb
          partitioned_data%point(i)%cluster=cluster(i)
       ENDDO
       DEALLOCATE(cluster)
       DEALLOCATE(points_by_clusters)
       DEALLOCATE(ratiomax)
       DEALLOCATE(clusters_energies)
       DEALLOCATE(ratiomin)
       DEALLOCATE(ratiomoy)
       DEALLOCATE(ratiorii)
       DEALLOCATE(ratiorij)
       DEALLOCATE(A)
       DEALLOCATE(Z)
       IF(solver == 0) DEALLOCATE(A2)
       DEALLOCATE(clusters_centers)
       DEALLOCATE(W)
    ELSE 
#if aff
       PRINT *, 'Process n', proc_id, ' : OK'
#endif
       DO i=1,partitioned_data%nb
          partitioned_data%point(i)%cluster=1
       ENDDO
#if aff
       PRINT *, 'Process n', proc_id,' : Cluster'
#endif
    ENDIF

    RETURN
  END SUBROUTINE apply_spectral_clustering



SUBROUTINE mean_shift(proc_id,nblimit,nb_clusters_opt,dataw,bandWidth)

   INCLUDE 'mpif.h'
    !IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
  
    INTEGER :: nb_clusters_opt
    INTEGER :: nblimit
    INTEGER :: proc_id
    INTEGER :: bandWidth !bandwidth parameter

    !=== IN/OUT ===
    TYPE(type_data) :: dataw
    
    !###########################################
    ! DECLARATIONS
    !###########################################      
    INTEGER ::point_num												!number of points
    INTEGER ::dim_num													!number of dimensions
    INTEGER ::cluster_num												!number of clusters
    
    !#### Variables  ####
    INTEGER :: numClust 												!the cluster number  
    DOUBLE PRECISION :: bandSq										!square of bandWidth
    DOUBLE PRECISION :: stopThresh									!when mean has converged
    INTEGER :: beenVisitedFlag(dataw%nb)								!track if a point has been seen already
    INTEGER :: numInitPts												!number of points to possibly use as initialization points
    INTEGER :: thisClusterVotes(dataw%nb)								!used to resolve conflicts on cluster membership
    INTEGER :: stInd													!start point of mean
    DOUBLE PRECISION :: myMean(dataw%dim)								!mean of this cluster
    DOUBLE PRECISION :: myOldMean(dataw%dim)							!old mean computed for this cluster
    INTEGER :: myMembers(dataw%nb)										!1 if the point belongs to the cluster, else 0
    INTEGER :: mergeWith												!used to merge clusters
    DOUBLE PRECISION :: clustCent(dataw%dim,dataw%nbclusters)			!centers of each cluster
	INTEGER :: clusterVotes(dataw%nbclusters,dataw%nb)					!number of votes for each point for each cluster
    INTEGER :: i
    INTEGER :: j
    INTEGER :: num
    DOUBLE PRECISION :: sqDist
    INTEGER :: cN
    
    
    
    !INITIALIZE STUFF    
    point_num=dataw%nb
    dim_num=dataw%dim
    numClust = 1
    bandSq = bandWidth**2
    stopThresh = 1e-3*bandWidth
    beenVisitedFlag(:) = 0
    numInitPts = point_num
    clusterVotes(:,:) = 0		
    
    DO WHILE (numInitPts>0)
    
		!take the first point as start of mean
		DO i=1, point_num
			IF (beenVisitedFlag(i)==0) THEN
				stInd = i
				EXIT
			ENDIF
		ENDDO
		myMean = dataw%point(stInd)%coord	!initialize mean to this points location
		DO j=1, dim_num
			myMean(j) = dataw%point(i)%coord(j)
		ENDDO
		myMembers(:) = 0
		thisClusterVotes(:) = 0	!used to resolve conflicts on cluster membership
		
		DO
		
			DO i=1, point_num
				!dist squared from mean to all points still active
				sqDist = 0
				DO j=1, dim_num
					sqDist = sqDist + (dataw%point(i)%coord(j) - myMean(j))**2
				ENDDO
				IF (sqDist < bandSq) THEN
					thisClusterVotes(i) = thisClusterVotes(i) + 1	!add a vote for all the in points belonging to this cluster
					myMembers(i) = 1								!add any point within bandWidth to the cluster
					beenVisitedFlag(i) = 1							!mark that these points have been visited
				ENDIF
			ENDDO
			
			myOldMean = myMean
			
			!compute the new mean
			DO i=1, point_num
				num = 0
				IF (myMembers(i)==1) THEN
					DO j=1, dim_num
						myMean(j) = myMean(j) + dataw%point(i)%coord(j)
					ENDDO
					num = num + 1
				ENDIF
			ENDDO
			myMean = myMean/num
			
			!compute the distance from myMean to myOldMean
			sqDist = 0
			DO j=1, dim_num
				sqDist = sqDist + (myOldMean(j) - myMean(j))**2
			ENDDO
			
			!if mean doesn't move much stop this cluster
			IF (sqDist < stopThresh**2) THEN
			
				!check for merge posibilities
				mergeWith = 0
				DO cN=1, numclust-1
					!compute the distance from possible new clust max to old clust max
					sqDist = 0
					DO j=1, dim_num
						sqDist = sqDist + (clustCent(j,cN) - myMean(j))**2
					ENDDO
					IF (sqDist < (bandWidth/2)**2) THEN
						mergeWith = cN
						EXIT
					ENDIF
				ENDDO
				
				IF (mergeWith > 0) THEN		!something to merge
				
					clustCent(:,mergeWith) = (myMean+clustCent(:,mergeWith))/2					!mean of centers
					clusterVotes(mergeWith,:) = clusterVotes(mergeWith,:) + thisClusterVotes;	!add these votes to the merged cluster
					
				ELSE
				
					numClust = numClust + 1
					clustCent(:,numClust) = myMean
					clusterVotes(numClust,:) = thisClusterVotes
				ENDIF
				EXIT
				
			ENDIF
		
		ENDDO
			
		DO i=1, point_num
			numInitPts = 0
			IF (beenVisitedFlag(i)==0) THEN
				numInitPts = numInitPts + 1
			ENDIF
		ENDDO
			
	ENDDO
	
	DO i=1, point_num
		dataw%point(i)%cluster = MAXLOC(clusterVotes(:,i), DIM=1)
	ENDDO
	
END SUBROUTINE mean_shift

END MODULE module_calcul
