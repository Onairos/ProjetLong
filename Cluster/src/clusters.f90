PROGRAM clusters
  USE module_structure
  USE module_entree
  USE module_decoupe
  USE module_MPI
  USE module_calcul
  USE module_sortie

  IMPLICIT NONE
 
  ! MPI library
  INCLUDE 'mpif.h'
  !###########################################
  ! DECLARATIONS
  !###########################################
  !#### Variables  ####
  TYPE(type_clustering_param) :: clust_param
  TYPE(type_clusters), DIMENSION(:), POINTER :: array_clusters
  TYPE(type_data) :: data
  TYPE(type_data) :: partitioned_data
  CHARACTER (LEN=80) :: proc_name ! MPI variable
  CHARACTER (LEN=30) :: input
  CHARACTER (LEN=30) :: input_file
  INTEGER :: i
  INTEGER :: ierr ! MPI variable
  INTEGER :: j
  INTEGER :: length ! MPI variable
  INTEGER :: n_max
  INTEGER :: nb_clusters
  INTEGER :: nb_clusters_max
  INTEGER :: nb_clusters_opt
  INTEGER :: nb_proc ! MPI variable
  INTEGER :: proc_id ! MPI variable
  INTEGER :: status(MPI_STATUS_SIZE) ! MPI variable
  INTEGER :: tag ! MPI variable
  INTEGER,DIMENSION(:), POINTER :: list_nb_clusters
  INTEGER,DIMENSION(:), POINTER :: partitioning
  INTEGER,DIMENSION(:), POINTER :: points_by_cluster
  INTEGER,DIMENSION(:), POINTER :: points_by_domain
  INTEGER,DIMENSION(:,:), POINTER :: assignments
  INTEGER,DIMENSION(:,:), POINTER :: cluster_map
  DOUBLE PRECISION :: end_time
  DOUBLE PRECISION :: epsilon
  DOUBLE PRECISION :: sigma
  DOUBLE PRECISION :: start_time
  DOUBLE PRECISION :: t1
  DOUBLE PRECISION :: t2
  DOUBLE PRECISION :: t_parall
  DOUBLE PRECISION :: t_parallg
  DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
  DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min
  DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bounds
  LOGICAL :: exist_bool

PRINT*, 'L55 clusters :DEBUG BANDWIDTH BEFORE READ :', clust_param%bandwidth


  !###########################################
  ! INSTRUCTIONS
  !###########################################
  ! Timers init
  t1 = 0.0
  t2 = 0.0
  start_time = 0.0
  end_time = 0.0
  ! MPI init
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,proc_id,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nb_proc,ierr)
  CALL MPI_GET_PROCESSOR_NAME(proc_name,length,ierr)
  PRINT *, 'Process ID : ', proc_id, '. Process name : ', proc_name

  IF(proc_id==0) THEN
    start_time = MPI_WTIME()
  ENDIF




#if aff
  PRINT *,'DEBUG : Launching process ',proc_id,' of ',nb_proc
#endif
  PRINT *,'----------------------------------------'

  ! Data reading
  IF (proc_id==0) THEN
     t1 = MPI_WTIME()
#if aff
     PRINT *
     PRINT *,'------------------------------'
     PRINT *,'DEBUG : Spectral clustering of data'
     PRINT *,'------------------------------'
     PRINT *
#endif
     IF (iargc()>0) THEN
        ! Gets the input file name
        CALL getarg(1,input)
        INQUIRE(FILE=input,EXIST=exist_bool)
        IF (.NOT.exist_bool) CALL help
     ELSE
       ! If no input file, default=fort.1
        CALL help
     ENDIF
     ! File reading
#if aff
     PRINT *
     PRINT *,'DEBUG : Reading data : ',input
#endif
     OPEN(FILE=input,UNIT=1)

     CALL read_params(nb_proc, data, clust_param, input_file, nb_clusters_max, &
                   list_nb_clusters, partitioning, epsilon, sigma, coord_max, &
                   coord_min)

     t2 = MPI_WTIME()
     PRINT *, 'Time for reading data : ', t2-t1
     CLOSE(1)
  ENDIF

  ! Partitioning and sigma computing
  IF (proc_id==0) THEN
     ! Partitioning
#if aff
     PRINT *
     PRINT *,'DEBUG : Partitioning data...'
#endif
    t1 = MPI_WTIME()
    CALL partition_data(data, nb_proc, partitioning, epsilon, coord_max, coord_min, points_by_domain, assignments, bounds)
    t2 = MPI_WTIME()
    PRINT *,'Time for partitioning data : ', t2-t1
  ENDIF
     
  IF(proc_id==0) THEN
    t1 = MPI_WTIME()
  ENDIF

  ! Exchanges
  IF (nb_proc>1) THEN
     ! Case of several proc 

     ! Sigma computing if auto global mode
     CALL MPI_BCAST(sigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     IF ((sigma==0.0).AND.(proc_id==0)) THEN
     PRINT *,'DEBUG : clustmethid : ', clust_param%clustering_method_id
        CALL get_sigma_interface(data, proc_id, partitioning, epsilon, bounds, sigma)
     ENDIF

     ! Sigma sending
     IF (sigma>=0.0) THEN
          CALL MPI_BCAST(sigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     ENDIF

     ! Nblimit sending
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     CALL MPI_BCAST(nb_clusters_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

     ! Nbideal sending
     IF (proc_id==0) THEN
        DO i=1,nb_proc-1
           tag=i
           CALL MPI_SEND(list_nb_clusters(i),1,MPI_INTEGER,i,tag,MPI_COMM_WORLD,ierr)
        ENDDO
        nb_clusters_opt=list_nb_clusters(0)
     ELSE
        tag=proc_id
        CALL MPI_RECV(nb_clusters_opt,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
     ENDIF
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

     ! Data transferring
     IF (proc_id==0) THEN
        ! Data sending
#if aff


        PRINT *
        PRINT *,'DEBUG : Transferring partitioned data...'
#endif
        CALL send_partitioning(nb_proc, points_by_domain, assignments, data, partitioned_data)
#if aff
        PRINT *
        PRINT *,'DEBUG : Computing clusters...'
#endif
     ELSE
        ! Data receiving
        CALL receive_partitioning(proc_id, partitioned_data)
     ENDIF

  ELSE
     ! Case of 1 proc alone
     partitioned_data%nb_points=data%nb_points
     partitioned_data%dim=data%dim
     partitioned_data%nb_clusters=0
     ALLOCATE(partitioned_data%points(data%nb_points))
     DO i=1,data%nb_points
        ALLOCATE(partitioned_data%points(i)%coords(data%dim))
        partitioned_data%points(i)%coords=data%points(i)%coords
        partitioned_data%points(i)%cluster=0
     ENDDO
     nb_clusters_opt=list_nb_clusters(1)
  ENDIF

  ! Sigma computing if auto individual mode
  IF (sigma<0.0) THEN
#if aff
     PRINT *
#endif
     IF (proc_id==0) CALL get_sigma(data, sigma)
     CALL MPI_BCAST(sigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     IF (proc_id==0) THEN
        PRINT *, proc_id,' : computing global sigma :',sigma
        IF (data%is_interfacing==1) THEN 
           CALL get_sigma_interface(partitioned_data, proc_id, partitioning, epsilon, bounds, sigma) 
           PRINT *, proc_id,' : computing interface sigma interface :',sigma
        ENDIF
     ENDIF
  ENDIF

  IF(proc_id==0) THEN
    t2 = MPI_WTIME()
    PRINT *,'Time for sending data and computing sigma : ', t2-t1
  ENDIF

  ! Clusters computing
  ! Parallel part
  t1 = MPI_WTIME()
  IF (partitioned_data%nb_points>0) THEN
#if aff
     PRINT *,'DEBUG : ', proc_id, ' : computing clusters...'
#endif



        CALL send_clustering_param( clust_param )


    !SELECT CASE (clust_param%clustering_method_id)
    !CASE (1)
     ! CALL apply_spectral_clustering(clust_param, nb_clusters_max, nb_clusters_opt, proc_id, sigma, partitioned_data)
    !CASE (2)
     ! CALL mean_shift(proc_id,nb_clusters_max,nb_clusters_opt,partitioned_data,clust_param)
    !CASE (3)
      CALL apply_kernel_k_means(proc_id,nb_clusters_max,nb_clusters_opt,partitioned_data,clust_param)
    !END SELECT

  ENDIF
  t2 = MPI_WTIME()
  t_parall = t2 - t1
  PRINT *, proc_id, ' : computing parallel cluster : ', t_parall

  CALL MPI_REDUCE(t_parall, t_parallg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  IF(proc_id==0) PRINT *, 'Time for computing global parallel cluster :', t_parallg

  IF(proc_id==0) THEN
    t1 = MPI_WTIME()
  ENDIF

  ! Saves the partial clusters
  CALL write_partial_clusters(proc_id, partitioned_data)

  ! Exchanges part
  IF (nb_proc>1) THEN
     ! Clusters grouping
     IF (proc_id==0) THEN
#if aff
        PRINT *
        PRINT *,'DEBUG : Grouping clusters...'
#endif
        ! Receiving of the number of clusters with duplications
        CALL receive_number_clusters(partitioned_data, nb_proc, points_by_domain, nb_clusters, array_clusters)
#if aff
        PRINT *,'DEBUG : number of clusters with duplications found : ', nb_clusters
#endif
        ! Receiving of clusters info
        ALLOCATE(cluster_map(nb_clusters,data%nb_points))
        CALL receive_clusters(array_clusters, partitioned_data, nb_clusters, nb_proc, points_by_domain,  &
                  assignments, points_by_cluster, cluster_map)
     ELSE
        ! Sends the number of clusters
        CALL send_number_clusters(partitioned_data, proc_id)
        ! Sends the clusters
        CALL send_clusters(partitioned_data, proc_id)
     ENDIF

     ! End of post-process
     IF (proc_id==0) THEN
        ! Groups the clusters and removes duplicates from the set of found clusters
        CALL group_clusters(nb_clusters, points_by_cluster, cluster_map, data)
     ENDIF

  ELSE
     ! Case of 1 proc alone
     nb_clusters=partitioned_data%nb_clusters
     ALLOCATE(points_by_cluster(nb_clusters))
     points_by_cluster(:)=0
     n_max=0
     DO i=1,partitioned_data%nb_points
        j=partitioned_data%points(i)%cluster
        points_by_cluster(j)=points_by_cluster(j)+1
        n_max=max(n_max,points_by_cluster(j))
     ENDDO
     ALLOCATE(cluster_map(nb_clusters,n_max))
     cluster_map(:,:)=0
     points_by_cluster(:)=0
     DO i=1,partitioned_data%nb_points
        j=partitioned_data%points(i)%cluster
        points_by_cluster(j)=points_by_cluster(j)+1
        cluster_map(j,points_by_cluster(j))=i
     ENDDO
  ENDIF

  IF(proc_id==0) THEN
    t2 = MPI_WTIME()
    PRINT *,'Time for grouping clusters : ', t2-t1
  ENDIF

  IF(proc_id==0) THEN
    t1 = MPI_WTIME()
  ENDIF
  ! Outputs
  IF (proc_id==0) THEN
     ! Writing of the cluster.final files
PRINT*, 'Calling write_final_clusters'
     CALL write_final_clusters(points_by_cluster, cluster_map, nb_clusters)

     ! Information writing
     CALL write_metadata(data, input_file, nb_clusters, nb_proc)
  ENDIF
  IF(proc_id==0) THEN
    t2 = MPI_WTIME()
    PRINT *,'Time for writing clusters : ', t2-t1
  ENDIF
  
  ! Ending of MPI
  IF (proc_id==0) THEN
    end_time = MPI_WTIME()
    PRINT *,'End of computing'
    PRINT *,'Overall time : ', end_time-start_time
  ENDIF
  
  CALL MPI_FINALIZE(ierr)
  STOP

END PROGRAM clusters
