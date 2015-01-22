module module_MPI
  use module_structure
contains

  !****************************************
  !envoi des decoupages
  subroutine envoidecoupes(nbproc,data,ldat,ddat,dataw)
    implicit none
    ! librairie MPI
    include 'mpif.h'
    !variables
    type(type_data) :: data,dataw
    integer :: nbproc
    integer,dimension(:),pointer :: ldat
    integer,dimension(:,:),pointer :: ddat
    integer :: i,j,m,n,tag,ierr
    real*8,dimension(:,:),pointer :: coord
    do i=1,nbproc-1
       m=ldat(i); n=data%dim
       tag=i
       call MPI_SEND(m,1,MPI_INTEGER,i,tag,MPI_COMM_WORLD,ierr)
       call MPI_SEND(n,1,MPI_INTEGER,i,tag,MPI_COMM_WORLD,ierr)   
       if (m>0) then
          !creation des tableaux de coordonnees
          allocate(coord(m,n)); coord=0.0
          do j=1,m
             coord(j,1:n)=data%point(ddat(i,j))%coord(1:n)
          enddo
          !envoi des tableaux
          tag=i*10
          call MPI_SEND(coord,m*n,MPI_DOUBLE_PRECISION,i,tag,MPI_COMM_WORLD,ierr)
          deallocate(coord)
       endif
    end do
    !creation du type dataw de l'interface
    m=ldat(0); n=data%dim
    dataw%nb=m; dataw%dim=n; dataw%nbclusters=0
    if (m>0) then
       allocate(dataw%point(m))
       do i=1,m
          allocate(dataw%point(i)%coord(n))
          dataw%point(i)%coord(:)=data%point(ddat(0,i))%coord(:)
          dataw%point(i)%cluster=0
       enddo
    endif
    !envoi des flags image, seuil, geom,...
    n=data%coord; dataw%coord=n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%image; dataw%image=n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%geom; dataw%geom=n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%seuil; dataw%seuil=n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%recouvrement; dataw%recouvrement=n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%interface; dataw%interface=n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    n=data%dim; dataw%dim=n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    return
  end subroutine envoidecoupes

  !***************************************
  !reception des decoupages
  subroutine recoitdecoupes(numproc,dataw)
    implicit none
    ! librairie MPI
    include 'mpif.h'
    !variables
    integer :: numproc
    type(type_data) :: dataw
    integer :: m,n,tag,ierr,i
    real*8,dimension(:,:),pointer :: coord
    integer status(MPI_STATUS_SIZE)
    !reception des dimensions
    tag=numproc
    call MPI_RECV(m,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
    call MPI_RECV(n,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
    dataw%nb=m; dataw%dim=n; dataw%nbclusters=0
    if (m>0) then
       allocate(coord(m,n)); coord=0.0
       !reception des tableaux
       tag=numproc*10
       call MPI_RECV(coord,m*n,MPI_DOUBLE_PRECISION,0,tag,&
            MPI_COMM_WORLD,status,ierr)
       !creation du type dataw du sous-domaine
       allocate(dataw%point(m))
       do i=1,m
          allocate(dataw%point(i)%coord(n))
          dataw%point(i)%coord=coord(i,:)
          dataw%point(i)%cluster=0
       enddo
       deallocate(coord)
    endif
    !envoi des flags image, seuil, geom,...
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    dataw%coord=n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    dataw%image=n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    dataw%geom=n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    dataw%seuil=n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    dataw%recouvrement=n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    dataw%interface=n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    dataw%dim=n
    return
  end subroutine recoitdecoupes


  !***************************************
  !compte des clusters avec doublons
  subroutine preparecpclusters(nbproc,nbclust,ldat,dataw,nclust)
    implicit none
    ! librairie MPI
    include 'mpif.h'
    !variables
    integer :: nbproc,nbclust
    type(type_data) ::dataw
    type(type_clusters),dimension(:),pointer :: nclust
    integer status(MPI_STATUS_SIZE)
    integer,dimension(:),pointer :: ldat
    integer :: i,j,nb,tag,ierr
    if (dataw%nb>0) then
       nbclust=dataw%nbclusters
    else
       nbclust=0
    endif
    !nb de clusters
    allocate(nclust(nbproc)); nclust(:)%nb=0
    do i=1,nbproc-1
       if (ldat(i)>0) then
          tag=i*11
          call MPI_RECV(nb,1,MPI_INTEGER,i,tag,MPI_COMM_WORLD,status,ierr)
          nbclust=nbclust+nb
          nclust(i)%nb=nb;
       endif
    enddo
    !nb de points par cluster
    do i=1,nbproc-1
       if (ldat(i)>0) then
          tag=i*11+1
          allocate(nclust(i)%nbelt(nclust(i)%nb))
          call MPI_RECV(nclust(i)%nbelt,nclust(i)%nb,MPI_INTEGER,i,tag,MPI_COMM_WORLD,status,ierr)
       endif
    enddo
    return
  end subroutine preparecpclusters

  !***************************************
  !envoi des nb de clusters
  subroutine prepaenvclusters(numproc,dataw)
    implicit none
    ! librairie MPI
    include 'mpif.h'
    !variables
    integer :: numproc
    type(type_data) ::dataw
    integer :: tag,ierr,i
    integer,dimension(:),pointer :: list
    if (dataw%nb>0) then
       !nb de clusters
       tag=numproc*11
       call MPI_SEND(dataw%nbclusters,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,ierr)
       !nb de points par cluster
       allocate(list(dataw%nbclusters))
       list(:) = 0
       do i=1,dataw%nb
         ! print *, dataw%point(i)%cluster
          list(dataw%point(i)%cluster)=list(dataw%point(i)%cluster)+1
       enddo
       tag=tag+1
       call MPI_SEND(list,dataw%nbclusters,MPI_INTEGER,0,tag,MPI_COMM_WORLD,ierr)
    endif
    return
  end subroutine prepaenvclusters

  !***************************************
  !envoi des clusters
  subroutine envoiclusters(numproc,dataw)
    implicit none
    ! librairie MPI
    include 'mpif.h'
    !variables
    integer :: numproc
    type(type_data) ::dataw
    integer,dimension(:),pointer :: lclust
    integer :: i,tag,ierr
    if (dataw%nb>0) then
       allocate(lclust(dataw%nb))
       do i=1,dataw%nb
          lclust(i)=dataw%point(i)%cluster
       enddo
       tag=numproc*12
       call MPI_SEND(lclust,dataw%nb,MPI_INTEGER,0,tag,MPI_COMM_WORLD,ierr)
       deallocate(lclust)
    endif
    return
  end subroutine envoiclusters

  !***************************************
  !reception des clusters
  subroutine recepclusters(nbproc,nbclust,ldat,ddat,dataw,clustermap,&
       nclust,iclust)
    implicit none
    ! librairie MPI
    include 'mpif.h'
    !variables
    integer :: nbproc,nbclust
    type(type_data) ::dataw
    integer,dimension(:),pointer :: ldat
    integer,dimension(:,:),pointer :: ddat
    integer,dimension(:),pointer :: lclust,iclust,listclust
    type(type_clusters),dimension(:),pointer :: nclust
    integer :: i,j,k,tag,ierr,i0
    integer,dimension(:,:),pointer :: clustermap
    integer status(MPI_STATUS_SIZE)
    integer :: p, longueur, maxldat, m
    i0=0; allocate(iclust(nbclust)); iclust(:)=0
    if (dataw%nb>0) then
       !stockage des clusters locaux dans le tableau global
       do i=1,dataw%nb
          j=dataw%point(i)%cluster
          iclust(j)=iclust(j)+1
          clustermap(j,iclust(j))=ddat(0,i)
       enddo
       i0=i0+dataw%nbclusters
    endif
    maxldat = maxval(ldat)
    allocate(lclust(maxldat))
    do i=1,nbproc-1
       if (ldat(i)>0) then
          !reception des affectations locales des points du sous-domaine
          !tag=i*12
          call MPI_RECV(lclust,maxldat,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
#if aff
          print *, 'mpi_recv', i, status(1), status(2), status(3), status(4) 
#endif
          p = status(MPI_SOURCE)
          !stockage des clusters locaux dans le tableau global
          do j=1,ldat(p)
             k=lclust(j)+i0
             iclust(k)=iclust(k)+1
             clustermap(k,iclust(k))=ddat(p,j)
          enddo
          i0=i0+nclust(p)%nb
       endif
    enddo
    deallocate(lclust)
    return
  end subroutine recepclusters


end module module_MPI
