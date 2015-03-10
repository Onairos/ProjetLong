MODULE module_sortie
  USE module_structure
CONTAINS


  SUBROUTINE write_domains(data, nbproc, domains)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: domains  
    INTEGER :: nbproc
    
    !#### Variables  ####
    INTEGER :: i
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    DO i=1,nbproc-data%interface
       WRITE(2,*) domains(i,:,1),'|', domains(i,:,2)
    ENDDO
    CALL flush(2)
    CLOSE(2)
    RETURN
  END SUBROUTINE write_domains


  SUBROUTINE write_partitioning(nbproc, data, points_by_domain, assignements)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    INTEGER,DIMENSION(:,:),POINTER :: assignements
    INTEGER,DIMENSION(:),POINTER :: points_by_domain
    INTEGER :: nbproc
    
    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nbdom
    INTEGER :: offset
    
    !###########################################      
    ! INSTRUCTIONS
    !########################################### 
    PRINT *, '> Partitioning review :'
    offset=1
    nbdom=nbproc
    IF ((data%interface==1).AND.(nbproc>1)) THEN
       offset=0
       nbdom=nbproc-1
    ENDIF
    IF (data%recouvrement==1) THEN
       offset=0
       nbdom=nbproc-1
    ENDIF
    DO i=offset,nbdom
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


  SUBROUTINE write_final_clusters(nbclust, points_by_cluster, cluster_map)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER,DIMENSION(:,:),POINTER :: cluster_map
    INTEGER,DIMENSION(:),POINTER :: points_by_cluster

    !=== IN/OUT === 
    INTEGER :: nbclust

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
    DO i=1,nbclust
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
    nbclust=k
    RETURN
  END SUBROUTINE write_final_clusters


  SUBROUTINE write_metadata(mesh, data, nbproc, nbclust)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    CHARACTER (LEN=30) :: mesh
    INTEGER :: nbclust
    INTEGER :: nbproc
    
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
    WRITE(3,*) nbproc
    WRITE(3,*) '# Partitioning by interfacing : '
    WRITE(3,*) data%interface
    WRITE(3,*) '# Partitioning by overlapping : '
    WRITE(3,*) data%recouvrement
    WRITE(3,*) '# Number of clusters : '
    WRITE(3,*) nbclust
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
