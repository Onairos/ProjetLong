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
  TYPE(type_clusters), DIMENSION(:), POINTER :: nclust
  TYPE(type_data) :: data
  TYPE(type_data) :: partitioned_data
  CHARACTER (LEN=80) :: procname ! MPI variable
  CHARACTER (LEN=30) :: entree
  CHARACTER (LEN=30) :: input_file
  DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bounds
  DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
  DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min
  DOUBLE PRECISION :: endtime
  DOUBLE PRECISION :: epsilon
  DOUBLE PRECISION :: sigma
  DOUBLE PRECISION :: starttime
  DOUBLE PRECISION :: t_parall
  DOUBLE PRECISION :: t_parallg
  DOUBLE PRECISION :: t1
  DOUBLE PRECISION :: t2
  INTEGER,DIMENSION(:,:), POINTER :: cluster_map
  INTEGER,DIMENSION(:,:) ,POINTER :: ddat
  INTEGER,DIMENSION(:), POINTER :: partitionning
  INTEGER,DIMENSION(:), POINTER :: iclust
  INTEGER,DIMENSION(:), POINTER :: ldat
  INTEGER,DIMENSION(:), POINTER :: list_nb_clusters
  INTEGER :: i
  INTEGER :: ierr ! MPI variable
  INTEGER :: j
  INTEGER :: len ! MPI variable
  INTEGER :: nbclust
  INTEGER :: nbideal
  INTEGER :: nblimit
  INTEGER :: nbproc ! MPI variable
  INTEGER :: nmax
  INTEGER :: numproc ! MPI variable
  INTEGER :: status(MPI_STATUS_SIZE) ! MPI variable
  INTEGER :: tag ! MPI variable
  LOGICAL :: existe

  !###########################################
  ! INSTRUCTIONS
  !###########################################
  ! Timers init
  t1 = 0.0
  t2 = 0.0
  starttime = 0.0
  endtime = 0.0
  ! MPI init
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,numproc,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nbproc,ierr)
  CALL MPI_GET_PROCESSOR_NAME(procname,len,ierr)
  PRINT *, 'rank', numproc, ' procname',procname

  IF(numproc==0) THEN
    starttime = MPI_WTIME()
  ENDIF
#if aff
  PRINT *,'lancement du proc',numproc,' de ',nbproc
#endif
  PRINT *,'----------------------------------------'

  ! Data reading
  IF (numproc==0) THEN
     t1 = MPI_WTIME()
#if aff
     PRINT *
     PRINT *,'------------------------------'
     PRINT *,'spectral clustering de datas'
     PRINT *,'------------------------------'
     PRINT *
#endif
     IF (iargc()>0) THEN
        ! Gets the input file name
        CALL getarg(1,entree)
        INQUIRE(FILE=entree,EXIST=existe)
        IF (.NOT.existe) CALL help
     ELSE
       ! If no input file, default=fort.1
        CALL help
     ENDIF
     ! File reading
#if aff
     PRINT *
     PRINT *,'lecture des data... ',entree
#endif
     OPEN(FILE=entree,UNIT=1)
     CALL read_file(data,epsilon,coord_min,coord_max,nbproc,partitionning,input_file,&
          sigma,nblimit,list_nb_clusters)
     t2 = MPI_WTIME()
     PRINT *, 'temps lecture data ', t2-t1
     CLOSE(1)
  ENDIF

  ! Partitionning and sigma computing
  IF (numproc==0) THEN
     ! Partitionning
#if aff
     PRINT *
     PRINT *,'decoupage des datas...'
#endif
    t1 = MPI_WTIME()
    CALL partition_data(data,epsilon,nbproc,coord_min,coord_max,partitionning,&
         ldat,ddat,bounds)
    t2 = MPI_WTIME()
    PRINT *,'temps decoupage des datas...', t2-t1
  ENDIF
     
  IF(numproc==0) THEN
    t1 = MPI_WTIME()
  ENDIF

  ! Exchanges
  IF (nbproc>1) THEN
     ! Case of several proc 

     ! Sigma computing if auto global mode
     CALL MPI_BCAST(sigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     IF ((sigma==0.0).AND.(numproc==0)) THEN
        CALL get_sigma_interface(numproc,data,sigma,bounds,partitionning,epsilon)
     ENDIF

     ! Sigma sending
     IF (sigma>=0.0) THEN
          CALL MPI_BCAST(sigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     ENDIF

     ! Nblimit sending
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     CALL MPI_BCAST(nblimit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

     ! Nbideal sending
     IF (numproc==0) THEN
        DO i=1,nbproc-1
           tag=i
           CALL MPI_SEND(list_nb_clusters(i),1,MPI_INTEGER,i,tag,MPI_COMM_WORLD,ierr)
        ENDDO
        nbideal=list_nb_clusters(0)
     ELSE
        tag=numproc
        CALL MPI_RECV(nbideal,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
     ENDIF
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

     ! Data transferring
     IF (numproc==0) THEN
        ! Data sending
#if aff
        PRINT *
        PRINT *,'transfert des datas decoupees...'
#endif
        CALL send_partitionning(nbproc,data,ldat,ddat,partitioned_data)
#if aff
        PRINT *
        PRINT *,'calcul des clusters...'
#endif
     ELSE
        ! Data receiving
        CALL receive_partitionning(numproc,partitioned_data)
     ENDIF

  ELSE
     ! Case of 1 proc alone
     partitioned_data%nb=data%nb
     partitioned_data%dim=data%dim
     partitioned_data%nbclusters=0
     ALLOCATE(partitioned_data%point(data%nb))
     DO i=1,data%nb
        ALLOCATE(partitioned_data%point(i)%coord(data%dim))
        partitioned_data%point(i)%coord=data%point(i)%coord
        partitioned_data%point(i)%cluster=0
     ENDDO
     nbideal=list_nb_clusters(1)
  ENDIF

  ! Sigma computing if auto individual mode
  IF (sigma<0.0) THEN
#if aff
     PRINT *
#endif
     IF (numproc==0) CALL get_sigma(data,sigma)
     CALL MPI_BCAST(sigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     IF (numproc==0) THEN
        PRINT *,numproc,'calcule le sigma global :',sigma
#if aff
        PRINT *,numproc,'calcule le sigma global :',sigma
#endif
        IF (data%interface==1) THEN 
           CALL get_sigma_interface(numproc,partitioned_data,sigma,bounds,partitionning,epsilon) 
           PRINT *,numproc,'calcule le sigma interface :',sigma
#if aff
           PRINT *,numproc,'calcule le sigma interface :',sigma
#endif
        ENDIF
     ENDIF
  ENDIF

  IF(numproc==0) THEN
    t2 = MPI_WTIME()
    PRINT *,'temps envoi donnees et calcul sigma', t2-t1
  ENDIF

  ! Clusters computing
  ! Parallel part
  t1 = MPI_WTIME()
  IF (partitioned_data%nb>0) THEN
#if aff
     PRINT *,numproc,'calcul des clusters...'
#endif
     CALL apply_spectral_clustering(numproc,nblimit,nbideal,partitioned_data,sigma)
  ENDIF
  t2 = MPI_WTIME()
  t_parall = t2 - t1
  PRINT *, numproc, 'calcul cluster //', t_parall

  CALL MPI_REDUCE(t_parall, t_parallg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  IF(numproc==0) PRINT *, 'temps calcul cluster global //', t_parallg

  IF(numproc==0) THEN
    t1 = MPI_WTIME()
  ENDIF

  ! Saves the partial clusters
  CALL write_partial_clusters(numproc,partitioned_data)

  ! Exchanges part
  IF (nbproc>1) THEN
     ! Clusters grouping
     IF (numproc==0) THEN
#if aff
        PRINT *
        PRINT *,'regroupement des clusters...'
#endif
        ! Receiving of the number of clusters with duplications
        CALL receive_number_clusters(nbproc,nbclust,ldat,partitioned_data,nclust)
#if aff
        PRINT *,'  > nb de clusters avec doublons obtenus :',nbclust
#endif
        ! Receiving of clusters info
        ALLOCATE(cluster_map(nbclust,data%nb))
        CALL receive_clusters(nbproc,nbclust,ldat,ddat,partitioned_data,&
             cluster_map,nclust,iclust)
     ELSE
        ! Sends the number of clusters
        CALL send_number_clusters(numproc,partitioned_data)
        ! Sends the clusters
        CALL send_clusters(numproc,partitioned_data)
     ENDIF

     ! End of post-process
     IF (numproc==0) THEN
        ! Groups the clusters and removes duplicates from the set of found clusters
        CALL group_clusters(nbclust,iclust,cluster_map,data)
     ENDIF

  ELSE
     ! Case of 1 proc alone
     nbclust=partitioned_data%nbclusters
     ALLOCATE(iclust(nbclust))
     iclust(:)=0
     nmax=0
     DO i=1,partitioned_data%nb
        j=partitioned_data%point(i)%cluster
        iclust(j)=iclust(j)+1
        nmax=max(nmax,iclust(j))
     ENDDO
     ALLOCATE(cluster_map(nbclust,nmax))
     cluster_map(:,:)=0
     iclust(:)=0
     DO i=1,partitioned_data%nb
        j=partitioned_data%point(i)%cluster
        iclust(j)=iclust(j)+1
        cluster_map(j,iclust(j))=i
     ENDDO
  ENDIF

  IF(numproc==0) THEN
    t2 = MPI_WTIME()
    PRINT *,'temps regroupement des clusters', t2-t1
  ENDIF

  IF(numproc==0) THEN
    t1 = MPI_WTIME()
  ENDIF
  ! Outputs
  IF (numproc==0) THEN
     ! Writing of the cluster.final files
     CALL write_final_clusters(nbclust,iclust,cluster_map)

     ! Information writing
     CALL write_metadata(input_file,data,nbproc,nbclust)
  ENDIF
  IF(numproc==0) THEN
    t2 = MPI_WTIME()
    PRINT *,'temps ecriture des clusters', t2-t1
  ENDIF
  
  ! Ending of MPI
  IF (numproc==0) THEN
    endtime = MPI_WTIME()
    PRINT *,' Fin du calcul'
    PRINT *,'  > temps total=', endtime-starttime
  ENDIF
  
  CALL MPI_FINALIZE(ierr)
  STOP

END PROGRAM clusters
