MODULE module_MPI
  USE module_structure
CONTAINS


  !****************************************
  !envoi des decoupages
  SUBROUTINE send_partitionning(nbproc,data,ldat,ddat,dataw)
    IMPLICIT NONE    
    ! librairie MPI
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER, DIMENSION(:,:), POINTER :: ddat
    INTEGER, DIMENSION(:), POINTER :: ldat
    INTEGER :: nbproc

    !=== IN/OUT ===
    TYPE(type_data) :: data

    !====  OUT ====
    TYPE(type_data) :: dataw
    
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
    DO i=1,nbproc-1
       m=ldat(i); n=data%dim
       tag=i
       CALL MPI_SEND(m,1,MPI_INTEGER,i,tag,MPI_COMM_WORLD,ierr)
       CALL MPI_SEND(n,1,MPI_INTEGER,i,tag,MPI_COMM_WORLD,ierr)   
       IF (m>0) THEN
          !creation des tableaux de coordonnees
          ALLOCATE(coord(m,n)); coord=0.0
          DO j=1,m
             coord(j,1:n)=data%point(ddat(i,j))%coord(1:n)
          ENDDO
          !envoi des tableaux
          tag=i*10
          CALL MPI_SEND(coord,m*n,MPI_DOUBLE_PRECISION,i,tag,MPI_COMM_WORLD,ierr)
          DEALLOCATE(coord)
       ENDIF
    ENDDO
    !creation du TYPE dataw de l'interface
    m=ldat(0); n=data%dim
    dataw%nb=m; dataw%dim=n; dataw%nbclusters=0
    IF (m>0) THEN
       ALLOCATE(dataw%point(m))
       DO i=1,m
          ALLOCATE(dataw%point(i)%coord(n))
          dataw%point(i)%coord(:)=data%point(ddat(0,i))%coord(:)
          dataw%point(i)%cluster=0
       ENDDO
    ENDIF
    !envoi des flags image, seuil, geom,...
    n=data%coord; dataw%coord=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%image; dataw%image=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%geom; dataw%geom=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%seuil; dataw%seuil=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%recouvrement; dataw%recouvrement=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%interface; dataw%interface=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%dim; dataw%dim=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    RETURN
  END SUBROUTINE send_partitionning

  !***************************************
  !reception des decoupages
  SUBROUTINE receive_partitionning(numproc,dataw)
    IMPLICIT NONE
    ! librairie MPI
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: numproc

    !====  OUT ====
    TYPE(type_data) :: dataw
    
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
    !reception des dimensions
    tag=numproc
    CALL MPI_RECV(m,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
    CALL MPI_RECV(n,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
    dataw%nb=m; dataw%dim=n; dataw%nbclusters=0
    IF (m>0) THEN
       ALLOCATE(coord(m,n)); coord=0.0
       !reception des tableaux
       tag=numproc*10
       CALL MPI_RECV(coord,m*n,MPI_DOUBLE_PRECISION,0,tag,&
            MPI_COMM_WORLD,status,ierr)
       !creation du TYPE dataw du sous-domaine
       ALLOCATE(dataw%point(m))
       DO i=1,m
          ALLOCATE(dataw%point(i)%coord(n))
          dataw%point(i)%coord=coord(i,:)
          dataw%point(i)%cluster=0
       ENDDO
       DEALLOCATE(coord)
    ENDIF
    !envoi des flags image, seuil, geom,...
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    dataw%coord=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    dataw%image=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    dataw%geom=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    dataw%seuil=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    dataw%recouvrement=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    dataw%interface=n
    CALL MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    dataw%dim=n
    RETURN
  END SUBROUTINE receive_partitionning


  !***************************************
  !compte des clusters avec doublons
  SUBROUTINE receive_number_clusters(nbproc,nbclust,ldat,dataw,nclust)
    IMPLICIT NONE
    ! librairie MPI
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ==== 
    TYPE(type_data) ::dataw
    INTEGER, DIMENSION(:), POINTER :: ldat
    INTEGER :: nbproc

    !====  OUT ====
    TYPE(type_clusters), DIMENSION(:), POINTER :: nclust
    INTEGER :: nbclust
    
    !#### Variables  ####    
    INTEGER status(MPI_STATUS_SIZE)
    INTEGER :: i
    INTEGER :: j
    INTEGER :: ierr
    INTEGER :: nb
    INTEGER :: tag
    
    !###########################################      
    ! INSTRUCTIONS
    !########################################### 
    IF (dataw%nb>0) THEN
       nbclust=dataw%nbclusters
    ELSE
       nbclust=0
    ENDIF
    !nb de clusters
    ALLOCATE(nclust(nbproc)); nclust(:)%nb=0
    DO i=1,nbproc-1
       IF (ldat(i)>0) THEN
          tag=i*11
          CALL MPI_RECV(nb,1,MPI_INTEGER,i,tag,MPI_COMM_WORLD,status,ierr)
          nbclust=nbclust+nb
          nclust(i)%nb=nb;
       ENDIF
    ENDDO
    !nb de points par cluster
    DO i=1,nbproc-1
       IF (ldat(i)>0) THEN
          tag=i*11+1
          ALLOCATE(nclust(i)%nbelt(nclust(i)%nb))
          CALL MPI_RECV(nclust(i)%nbelt,nclust(i)%nb,MPI_INTEGER,i,tag,MPI_COMM_WORLD,status,ierr)
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE receive_number_clusters

  !***************************************
  !envoi des nb de clusters
  SUBROUTINE send_number_clusters(numproc,dataw)
    IMPLICIT NONE
    ! librairie MPI
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) ::dataw
    INTEGER :: numproc

    !#### Variables  ####
    INTEGER, DIMENSION(:), POINTER :: list
    INTEGER :: tag
    INTEGER :: i
    INTEGER :: ierr
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################
    IF (dataw%nb>0) THEN
       !nb de clusters
       tag=numproc*11
       CALL MPI_SEND(dataw%nbclusters,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,ierr)
       !nb de points par cluster
       ALLOCATE(list(dataw%nbclusters))
       list(:) = 0
       DO i=1,dataw%nb
          list(dataw%point(i)%cluster)=list(dataw%point(i)%cluster)+1
       ENDDO
       tag=tag+1
       CALL MPI_SEND(list,dataw%nbclusters,MPI_INTEGER,0,tag,MPI_COMM_WORLD,ierr)
    ENDIF
    RETURN
  END SUBROUTINE send_number_clusters

  !***************************************
  !envoi des clusters
  SUBROUTINE send_clusters(numproc,dataw)
    IMPLICIT NONE
    ! librairie MPI
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: dataw
    INTEGER :: numproc
    
    !#### Variables  ####  
    INTEGER, DIMENSION(:), POINTER :: lclust
    INTEGER :: i
    INTEGER :: ierr
    INTEGER :: tag
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    IF (dataw%nb>0) THEN
       ALLOCATE(lclust(dataw%nb))
       DO i=1,dataw%nb
          lclust(i)=dataw%point(i)%cluster
       ENDDO
       tag=numproc*12
       CALL MPI_SEND(lclust,dataw%nb,MPI_INTEGER,0,tag,MPI_COMM_WORLD,ierr)
       DEALLOCATE(lclust)
    ENDIF
    RETURN
  END SUBROUTINE send_clusters

  !***************************************
  !reception des clusters
  SUBROUTINE receive_clusters(nbproc,nbclust,ldat,ddat,dataw,clustermap,&
       nclust,iclust)
    IMPLICIT NONE
    ! librairie MPI
    INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_clusters), DIMENSION(:), POINTER :: nclust
    TYPE(type_data) ::dataw
    INTEGER, DIMENSION(:,:), POINTER :: ddat
    INTEGER, DIMENSION(:), POINTER :: ldat 
    INTEGER :: nbproc
    INTEGER :: nbclust

    !====  OUT ====
    INTEGER, DIMENSION(:,:), POINTER :: clustermap
    INTEGER, DIMENSION(:), POINTER :: iclust
    
    !#### Variables  ####
    INTEGER, DIMENSION(:), POINTER :: lclust
    INTEGER, DIMENSION(:), POINTER :: listclust
    INTEGER status(MPI_STATUS_SIZE)
    INTEGER :: i
    INTEGER :: i0
    INTEGER :: ierr
    INTEGER :: j
    INTEGER :: k
    INTEGER :: longueur
    INTEGER :: m
    INTEGER :: maxldat
    INTEGER :: p
    INTEGER :: tag
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    i0=0; ALLOCATE(iclust(nbclust)); iclust(:)=0
    IF (dataw%nb>0) THEN
       !stockage des clusters locaux dans le tableau global
       DO i=1,dataw%nb
          j=dataw%point(i)%cluster
          iclust(j)=iclust(j)+1
          clustermap(j,iclust(j))=ddat(0,i)
       ENDDO
       i0=i0+dataw%nbclusters
    ENDIF
    maxldat = maxval(ldat)
    ALLOCATE(lclust(maxldat))
    DO i=1,nbproc-1
       IF (ldat(i)>0) THEN
          !reception des affectations locales des points du sous-domaine
          CALL MPI_RECV(lclust,maxldat,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
#if aff
          PRINT *, 'mpi_recv', i, status(1), status(2), status(3), status(4) 
#endif
          p = status(MPI_SOURCE)
          !stockage des clusters locaux dans le tableau global
          DO j=1,ldat(p)
             k=lclust(j)+i0
             iclust(k)=iclust(k)+1
             clustermap(k,iclust(k))=ddat(p,j)
          ENDDO
          i0=i0+nclust(p)%nb
       ENDIF
    ENDDO
    DEALLOCATE(lclust)
    RETURN
  END SUBROUTINE receive_clusters


END MODULE module_MPI
