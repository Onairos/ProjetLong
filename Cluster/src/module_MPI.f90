MODULE module_MPI
  USE module_structure
CONTAINS

  !****************************************
  !envoi des decoupages
  SUBROUTINE envoidecoupes(nbproc,data,ldat,ddat,dataw)
    IMPLICIT NONE
    ! librairie MPI
    INCLUDE 'mpif.h'
    !variables
    TYPE(type_data) :: data,dataw
    INTEGER :: nbproc
    INTEGER,DIMENSION(:),POINTER :: ldat
    INTEGER,DIMENSION(:,:),POINTER :: ddat
    INTEGER :: i,j,m,n,tag,ierr
    REAL*8,DIMENSION(:,:),POINTER :: coord
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
  END SUBROUTINE envoidecoupes

  !***************************************
  !reception des decoupages
  SUBROUTINE recoitdecoupes(numproc,dataw)
    IMPLICIT NONE
    ! librairie MPI
    INCLUDE 'mpif.h'
    !variables
    INTEGER :: numproc
    TYPE(type_data) :: dataw
    INTEGER :: m,n,tag,ierr,i
    REAL*8,DIMENSION(:,:),POINTER :: coord
    INTEGER status(MPI_STATUS_SIZE)
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
  END SUBROUTINE recoitdecoupes


  !***************************************
  !compte des clusters avec doublons
  SUBROUTINE preparecpclusters(nbproc,nbclust,ldat,dataw,nclust)
    IMPLICIT NONE
    ! librairie MPI
    INCLUDE 'mpif.h'
    !variables
    INTEGER :: nbproc,nbclust
    TYPE(type_data) ::dataw
    TYPE(type_clusters),DIMENSION(:),POINTER :: nclust
    INTEGER status(MPI_STATUS_SIZE)
    INTEGER,DIMENSION(:),POINTER :: ldat
    INTEGER :: i,j,nb,tag,ierr
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
  END SUBROUTINE preparecpclusters

  !***************************************
  !envoi des nb de clusters
  SUBROUTINE prepaenvclusters(numproc,dataw)
    IMPLICIT NONE
    ! librairie MPI
    INCLUDE 'mpif.h'
    !variables
    INTEGER :: numproc
    TYPE(type_data) ::dataw
    INTEGER :: tag,ierr,i
    INTEGER,DIMENSION(:),POINTER :: list
    IF (dataw%nb>0) THEN
       !nb de clusters
       tag=numproc*11
       CALL MPI_SEND(dataw%nbclusters,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,ierr)
       !nb de points par cluster
       ALLOCATE(list(dataw%nbclusters))
       list(:) = 0
       DO i=1,dataw%nb
         ! PRINT *, dataw%point(i)%cluster
          list(dataw%point(i)%cluster)=list(dataw%point(i)%cluster)+1
       ENDDO
       tag=tag+1
       CALL MPI_SEND(list,dataw%nbclusters,MPI_INTEGER,0,tag,MPI_COMM_WORLD,ierr)
    ENDIF
    RETURN
  END SUBROUTINE prepaenvclusters

  !***************************************
  !envoi des clusters
  SUBROUTINE envoiclusters(numproc,dataw)
    IMPLICIT NONE
    ! librairie MPI
    INCLUDE 'mpif.h'
    !variables
    INTEGER :: numproc
    TYPE(type_data) ::dataw
    INTEGER,DIMENSION(:),POINTER :: lclust
    INTEGER :: i,tag,ierr
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
  END SUBROUTINE envoiclusters

  !***************************************
  !reception des clusters
  SUBROUTINE recepclusters(nbproc,nbclust,ldat,ddat,dataw,clustermap,&
       nclust,iclust)
    IMPLICIT NONE
    ! librairie MPI
    INCLUDE 'mpif.h'
    !variables
    INTEGER :: nbproc,nbclust
    TYPE(type_data) ::dataw
    INTEGER,DIMENSION(:),POINTER :: ldat
    INTEGER,DIMENSION(:,:),POINTER :: ddat
    INTEGER,DIMENSION(:),POINTER :: lclust,iclust,listclust
    TYPE(type_clusters),DIMENSION(:),POINTER :: nclust
    INTEGER :: i,j,k,tag,ierr,i0
    INTEGER,DIMENSION(:,:),POINTER :: clustermap
    INTEGER status(MPI_STATUS_SIZE)
    INTEGER :: p, longueur, maxldat, m
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
          !tag=i*12
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
  END SUBROUTINE recepclusters


END MODULE module_MPI
