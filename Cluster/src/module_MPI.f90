MODULE module_MPI
  USE module_structure
CONTAINS



  SUBROUTINE send_partitioning(nb_proc, data, points_by_domain, ddat, partitioned_data)
    IMPLICIT NONE    
    ! MPI library
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER, DIMENSION(:,:), POINTER :: ddat
    INTEGER, DIMENSION(:), POINTER :: points_by_domain
    INTEGER :: nb_proc

    !=== IN/OUT ===
    TYPE(type_data) :: data

    !====  OUT ====
    TYPE(type_data) :: partitioned_data
    
    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coord
    INTEGER :: i
    INTEGER :: j
    INTEGER :: ierr
    INTEGER :: m
    INTEGER :: n
    INTEGER :: tag
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
    DO i=1,nb_proc-1
       m=points_by_domain(i)
       n=data%dim
       tag=i
       CALL MPI_SEND(m,1,MPI_INTEGER,i,tag,MPI_COMM_WORLD,ierr)
       CALL MPI_SEND(n,1,MPI_INTEGER,i,tag,MPI_COMM_WORLD,ierr)   
       IF (m>0) THEN
          ! Creation of coordinates arrays
          ALLOCATE(coord(m,n))
          coord=0.0
          DO j=1,m
             coord(j,1:n)=data%point(ddat(i,j))%coord(1:n)
          ENDDO
          ! Sending arrays
          tag=i*10
          CALL MPI_SEND(coord,m*n,MPI_DOUBLE_PRECISION,i,tag,MPI_COMM_WORLD,ierr)
          DEALLOCATE(coord)
       ENDIF
    ENDDO
    ! Creation of TYPE partitioned_data of interface
    m=points_by_domain(0)
    n=data%dim
    partitioned_data%nb=m
    partitioned_data%dim=n
    partitioned_data%nbclusters=0
    IF (m>0) THEN
       ALLOCATE(partitioned_data%point(m))
       DO i=1,m
          ALLOCATE(partitioned_data%point(i)%coord(n))
          partitioned_data%point(i)%coord(:)=data%point(ddat(0,i))%coord(:)
          partitioned_data%point(i)%cluster=0
       ENDDO
    ENDIF
    ! Sending flags picture, threshold, geometric...
    n=data%coord
    partitioned_data%coord=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%image
    partitioned_data%image=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%geom
    partitioned_data%geom=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%seuil
    partitioned_data%seuil=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%recouvrement
    partitioned_data%recouvrement=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%interface
    partitioned_data%interface=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%dim
    partitioned_data%dim=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    RETURN
  END SUBROUTINE send_partitioning


  SUBROUTINE receive_partitioning(proc_id, partitioned_data)
    IMPLICIT NONE
    ! MPI library
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: proc_id

    !====  OUT ====
    TYPE(type_data) :: partitioned_data
    
    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coord
    INTEGER status(MPI_STATUS_SIZE)
    INTEGER :: i
    INTEGER :: ierr
    INTEGER :: m
    INTEGER :: n
    INTEGER :: tag
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################   
    ! Receiving dimensions
    tag=proc_id
    CALL MPI_RECV(m,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
    CALL MPI_RECV(n,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
    partitioned_data%nb=m
    partitioned_data%dim=n
    partitioned_data%nbclusters=0
    IF (m>0) THEN
       ALLOCATE(coord(m,n))
       coord=0.0
       ! Receiving arrays
       tag=proc_id*10
       CALL MPI_RECV(coord,m*n,MPI_DOUBLE_PRECISION,0,tag,&
            MPI_COMM_WORLD,status,ierr)
       ! Creation of TYPE partitioned_data of subdomain
       ALLOCATE(partitioned_data%point(m))
       DO i=1,m
          ALLOCATE(partitioned_data%point(i)%coord(n))
          partitioned_data%point(i)%coord=coord(i,:)
          partitioned_data%point(i)%cluster=0
       ENDDO
       DEALLOCATE(coord)
    ENDIF
    ! Sending flags picture, threshold, geometric...
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    partitioned_data%coord=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    partitioned_data%image=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    partitioned_data%geom=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    partitioned_data%seuil=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    partitioned_data%recouvrement=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    partitioned_data%interface=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    partitioned_data%dim=n
    RETURN
  END SUBROUTINE receive_partitioning



  SUBROUTINE receive_number_clusters(nb_proc, nbclust, points_by_domain, partitioned_data, array_clust)
    IMPLICIT NONE
    ! MPI library
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ==== 
    TYPE(type_data) ::partitioned_data
    INTEGER, DIMENSION(:), POINTER :: points_by_domain
    INTEGER :: nb_proc

    !====  OUT ====
    TYPE(type_clusters), DIMENSION(:), POINTER :: array_clust
    INTEGER :: nbclust
    
    !#### Variables  ####    
    INTEGER status(MPI_STATUS_SIZE)
    INTEGER :: i
    INTEGER :: ierr
    INTEGER :: nb
    INTEGER :: tag
    
    !###########################################      
    ! INSTRUCTIONS
    !########################################### 
    IF (partitioned_data%nb>0) THEN
       nbclust=partitioned_data%nbclusters
    ELSE
       nbclust=0
    ENDIF
    ! Numebr of clusters
    ALLOCATE(array_clust(nb_proc))
    array_clust(:)%nb=0
    DO i=1,nb_proc-1
       IF (points_by_domain(i)>0) THEN
          tag=i*11
          CALL MPI_RECV(nb,1,MPI_INTEGER,i,tag,MPI_COMM_WORLD,status,ierr)
          nbclust=nbclust+nb
          array_clust(i)%nb=nb
       ENDIF
    ENDDO
    ! Number of points by cluster
    DO i=1,nb_proc-1
       IF (points_by_domain(i)>0) THEN
          tag=i*11+1
          ALLOCATE(array_clust(i)%nbelt(array_clust(i)%nb))
          CALL MPI_RECV(array_clust(i)%nbelt,array_clust(i)%nb,MPI_INTEGER,i,tag,MPI_COMM_WORLD,status,ierr)
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE receive_number_clusters


  SUBROUTINE send_number_clusters(proc_id, partitioned_data)
    IMPLICIT NONE
    ! MPI library
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) ::partitioned_data
    INTEGER :: proc_id

    !#### Variables  ####
    INTEGER, DIMENSION(:), POINTER :: list
    INTEGER :: tag
    INTEGER :: i
    INTEGER :: ierr
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################
    IF (partitioned_data%nb>0) THEN
       ! Number of clusters
       tag=proc_id*11
       CALL MPI_SEND(partitioned_data%nbclusters,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,ierr)
       ! Number of points by cluster
       ALLOCATE(list(partitioned_data%nbclusters))
       list(:) = 0
       DO i=1,partitioned_data%nb
          list(partitioned_data%point(i)%cluster)=list(partitioned_data%point(i)%cluster)+1
       ENDDO
       tag=tag+1
       CALL MPI_SEND(list,partitioned_data%nbclusters,MPI_INTEGER,0,tag,MPI_COMM_WORLD,ierr)
    ENDIF
    RETURN
  END SUBROUTINE send_number_clusters



  SUBROUTINE send_clusters(proc_id, partitioned_data)
    IMPLICIT NONE
    ! MPI library
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: partitioned_data
    INTEGER :: proc_id
    
    !#### Variables  ####  
    INTEGER, DIMENSION(:), POINTER :: lclust
    INTEGER :: i
    INTEGER :: ierr
    INTEGER :: tag
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    IF (partitioned_data%nb>0) THEN
       ALLOCATE(lclust(partitioned_data%nb))
       DO i=1,partitioned_data%nb
          lclust(i)=partitioned_data%point(i)%cluster
       ENDDO
       tag=proc_id*12
       CALL MPI_SEND(lclust,partitioned_data%nb,MPI_INTEGER,0,tag,MPI_COMM_WORLD,ierr)
       DEALLOCATE(lclust)
    ENDIF
    RETURN
  END SUBROUTINE send_clusters



  SUBROUTINE receive_clusters(nb_proc, nbclust, points_by_domain, ddat, partitioned_data, cluster_map, &
       array_clust, points_by_cluster)
    IMPLICIT NONE
    ! MPI library
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_clusters), DIMENSION(:), POINTER :: array_clust
    TYPE(type_data) ::partitioned_data
    INTEGER, DIMENSION(:,:), POINTER :: ddat
    INTEGER, DIMENSION(:), POINTER :: points_by_domain 
    INTEGER :: nb_proc
    INTEGER :: nbclust

    !====  OUT ====
    INTEGER, DIMENSION(:,:), POINTER :: cluster_map
    INTEGER, DIMENSION(:), POINTER :: points_by_cluster
    
    !#### Variables  ####
    INTEGER, DIMENSION(:), POINTER :: lclust
    INTEGER status(MPI_STATUS_SIZE)
    INTEGER :: i
    INTEGER :: i0
    INTEGER :: ierr
    INTEGER :: j
    INTEGER :: k
    INTEGER :: maxldat
    INTEGER :: p
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    i0=0
    ALLOCATE(points_by_cluster(nbclust))
    points_by_cluster(:)=0
    IF (partitioned_data%nb>0) THEN
       ! Storage of local clusters in the global array
       DO i=1,partitioned_data%nb
          j=partitioned_data%point(i)%cluster
          points_by_cluster(j)=points_by_cluster(j)+1
          cluster_map(j,points_by_cluster(j))=ddat(0,i)
       ENDDO
       i0=i0+partitioned_data%nbclusters
    ENDIF
    maxldat = maxval(points_by_domain)
    ALLOCATE(lclust(maxldat))
    DO i=1,nb_proc-1
       IF (points_by_domain(i)>0) THEN
          ! Receiving local allocations of subdomain points
          CALL MPI_RECV(lclust,maxldat,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
#if aff
          PRINT *, 'DEBUG : MPI_RECV ', i, status(1), status(2), status(3), status(4) 
#endif
          p = status(MPI_SOURCE)
          ! Storage of local clusters in the global array
          DO j=1,points_by_domain(p)
             k=lclust(j)+i0
             points_by_cluster(k)=points_by_cluster(k)+1
             cluster_map(k,points_by_cluster(k))=ddat(p,j)
          ENDDO
          i0=i0+array_clust(p)%nb
       ENDIF
    ENDDO
    DEALLOCATE(lclust)
    RETURN
  END SUBROUTINE receive_clusters


END MODULE module_MPI
