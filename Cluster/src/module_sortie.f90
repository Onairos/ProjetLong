!>Contains methods enabling writing results in specific formatted files
MODULE module_sortie
  USE module_structure
CONTAINS


!>Writes the file containing the domain definitions
!!@details This methods writes the domain definitions using the
!!following formatting: each line contains the two point coordinates
!!of a bound separated by the character "|".
!!@note The written file is <em>fort.2</em>.S
!! @param[in] data the entire data for computing
!! @param[in] domains the domains constructed from the bounds
!! @param[in] nb_proc the number of processors used
  SUBROUTINE write_domains(data, nb_proc, domains)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: domains  
    INTEGER :: nb_proc
    
    !#### Variables  ####
    INTEGER :: i
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    DO i=1,nb_proc-data%interface
       WRITE(2,*) domains(i,:,1),'|', domains(i,:,2)
    ENDDO
    CALL flush(2)
    CLOSE(2)
    RETURN
  END SUBROUTINE write_domains


!>Writes the files containing the partitionning
!!@details Each file correspond to a domain. The first number
!!is the number of points in the domain. Then the following
!!lines are simply the coordinates of the points (in case of
!!coordinates format) or the color of the pixels (in case of
!!picture format) separated by blank spaces.
!!@note The written files are <em>decoupe.x</em> with x the ids of the processes.
!! @param[in] data the entire data for computing
!! @param[in] assignements the assignement of each point in a partition
!! @param[in] nb_proc the number of processors used
!! @param[in] points_by_domain the number of points in each partition
  SUBROUTINE write_partitioning(nb_proc, data, points_by_domain, assignements)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    INTEGER,DIMENSION(:,:),POINTER :: assignements
    INTEGER,DIMENSION(:),POINTER :: points_by_domain
    INTEGER :: nb_proc
    
    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb_domains
    INTEGER :: offset
    
    !###########################################      
    ! INSTRUCTIONS
    !########################################### 
    PRINT *, '> Partitioning review :'
    offset=1
    nb_domains=nb_proc
    IF ((data%interface==1).AND.(nb_proc>1)) THEN
       offset=0
       nb_domains=nb_proc-1
    ENDIF
    IF (data%recouvrement==1) THEN
       offset=0
       nb_domains=nb_proc-1
    ENDIF
    DO i=offset,nb_domains
       PRINT *, '> Zone ', i, ' : ', points_by_domain(i)
       ! File name
       WRITE(num,*) i
       files='decoupe.'//trim(adjustl(num))
       OPEN(FILE=files,UNIT=10)
       WRITE(10,*) points_by_domain(i)
       DO j=1,points_by_domain(i)
          IF (data%coord==1) THEN
             ! Writing in coordinates
             WRITE(10,*) data%point(assignements(i,j))%coord(:)
          ELSEIF ((data%image==1).OR.(data%seuil==1).OR.(data%geom==1)) THEN
             ! Writing in picture format
             WRITE(10,*) assignements(i,j)
          ENDIF
       ENDDO
       CALL flush(10)
       CLOSE(10)
    ENDDO
    CALL flush(6)
    RETURN
  END SUBROUTINE write_partitioning


!>Writes the file containing the clusters computed by the process
!!@details The file starts with the number of points in the domain
!!and the number of dimensions separated by a blank space. The
!!following lines are :
!!<ul>
!!<li> <b> For coordinates format </b> : the coordinates of a point
!!and the id of the cluster which it belongs to separated by a comma </li>
!!<li> <b> For picture formats </b> : the id of the point and the 
!!id of the cluster which it belongs to separated by a comma</li>
!!</ul>
!!@note The written file is <em>cluster.partiel.x</em> with x the id of the process
!! @param[in] partitioned_data the partitioned data for computing
!! @param[in] proc_id the processus identifier
  SUBROUTINE write_partial_clusters(proc_id, partitioned_data)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====  
    TYPE(type_data) :: partitioned_data
    INTEGER :: proc_id
    
    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    INTEGER :: i
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################
    ! File name
    WRITE(num,*),proc_id
    num=adjustl(num)
    files='cluster.partiel.'//trim(num)
    PRINT *, 'Process n', proc_id, ' : clusters writing : ', files
    OPEN(FILE=files,UNIT=10)
    WRITE(10,*) partitioned_data%nb,partitioned_data%dim
    DO i=1,partitioned_data%nb
       IF (partitioned_data%coord==1) THEN
          WRITE(10,*) partitioned_data%point(i)%coord(:), partitioned_data%point(i)%cluster
       ELSE
          WRITE(10,*) i, partitioned_data%point(i)%cluster
       ENDIF
    ENDDO
    CALL flush(10)
    CLOSE(10)
    RETURN
  END SUBROUTINE write_partial_clusters


!>Writes the files containing the final clusters after grouping
!!@details Each file correspond to a cluster. The first number
!!is the number of points in the cluster. Then all the following
!!numbers are the indices of the points that belongs to the cluster.
!!@note The written files are <em>cluster.final.x</em> with x the ids of the clusters.
!! @param[in] cluster_map the cluster indices and the number of points in each cluster
!! @param[in] points_by_cluster the number of points in each cluster
!! @param[in,out] nb_clusters the number of clusters
!! @param[in,out] nb_clusters the number of clusters
!! @param[in,out] nb_clusters the number of clusters
  SUBROUTINE write_final_clusters(nb_clusters, points_by_cluster, cluster_map)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER,DIMENSION(:,:),POINTER :: cluster_map
    INTEGER,DIMENSION(:),POINTER :: points_by_cluster

    !=== IN/OUT === 
    INTEGER :: nb_clusters

    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    
    !###########################################      
    ! INSTRUCTIONS
    !########################################### 
    PRINT *, '> Result writing...'
    k=0
    DO i=1,nb_clusters
       IF (points_by_cluster(i)>0) THEN
          k=k+1
          WRITE(num,*) k
          files='cluster.final.'//trim(adjustl(num))
          OPEN(FILE=files,UNIT=20)     
          PRINT *, '> Cluster n', k, ' : ', points_by_cluster(i), ' -> ', files
          WRITE(20,*) points_by_cluster(i)
          DO j=1,points_by_cluster(i)
             WRITE(20,*) cluster_map(i,j)
          ENDDO
          CALL flush(20)
          CLOSE(20)
       ENDIF
    ENDDO
    nb_clusters=k
    RETURN
  END SUBROUTINE write_final_clusters


!>Writes a file containing various information
!!@details The following information is written in 
!!<em>fort.3</em> file :
!!<ol>
!!<li> The data file name </li>
!!<li> The number of the points in the entire data set </li>
!!<li> The number of processes used </li>
!!<li> The partitioning mode (by interface or overlapping) </li>
!!<li> The number of clusters found </li>
!!<li> The data file format </li>
!!</ol>
!!In the case of a picture format: this extra information is
!!written :
!!<ol>
!!<li> The image dimension </li>
!!<li> The image partitioning </li>
!!<li> The number of attributes </li>
!!<li> The number of steps (only in geometric format) </li>
!!</ol> 
!!@see read_metadata()
!! @param[in] data the entire data for computing
!! @param[in] input_file the name of the text file where input data is written
!! @param[in] nb_clusters the number of clusters
!! @param[in] nb_clusters the number of clusters
!! @param[in] nb_clusters the number of clusters
!! @param[in] nb_proc the number of processors used
  SUBROUTINE write_metadata(mesh, data, nb_proc, nb_clusters)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    CHARACTER (LEN=30) :: mesh
    INTEGER :: nb_clusters
    INTEGER :: nb_proc
    
    !###########################################      
    ! INSTRUCTIONS
    !########################################### 
    WRITE(3,*) '# Mesh file : '
    WRITE(3,*) mesh
    WRITE(3,*) '# Number of points : '
    WRITE(3,*) data%nb
    WRITE(3,*) '# DIMENSION : '
    WRITE(3,*) data%dim
    WRITE(3,*) '# Number of process : '
    WRITE(3,*) nb_proc
    WRITE(3,*) '# Partitioning by interfacing : '
    WRITE(3,*) data%interface
    WRITE(3,*) '# Partitioning by overlapping : '
    WRITE(3,*) data%recouvrement
    WRITE(3,*) '# Number of clusters : '
    WRITE(3,*) nb_clusters
    WRITE(3,*) '# Coord format : '
    WRITE(3,*) data%coord
    WRITE(3,*) '# Image format : '
    WRITE(3,*) data%image
    WRITE(3,*) '# Geom format : '
    WRITE(3,*) data%geom
    WRITE(3,*) '# Threshold format : '
    WRITE(3,*) data%seuil
    IF ((data%image==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
       WRITE(3,*) '# DIMENSION : '
       WRITE(3,*) data%imgdim
       WRITE(3,*) '# Partitioning : '
       WRITE(3,*) data%imgmap(:)
       WRITE(3,*) '# Number of time : '
       WRITE(3,*) data%imgt
       IF (data%geom==1) THEN
          WRITE(3,*) '## No mesh : '
          WRITE(3,*) data%pas(:)
       ENDIF
    ENDIF
    CALL flush(3)
    CLOSE(3)
    RETURN
  END SUBROUTINE write_metadata

END MODULE module_sortie
