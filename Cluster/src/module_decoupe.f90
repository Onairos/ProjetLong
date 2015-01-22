module module_decoupe
  use module_structure
  use module_sortie
contains

  !****************************************
  !decoupage pour transfert MPI
  subroutine decoupedata(data,epsilon,nbproc,coordmin,coordmax,decoupe,&
       ldat,ddat,bornes)
    implicit none
    type(type_data) :: data
    real*8 :: epsilon
    integer :: nbproc, ierr
    integer,dimension(:),pointer :: decoupe
    real*8,dimension(:),pointer :: coordmax,coordmin
    integer,dimension(:),pointer :: ldat
    integer,dimension(:,:),pointer :: ddat
    real*8,dimension(:,:,:),pointer :: bornes,domaines
 

    !definition des bornes
    call definit_bornes(data,coordmin,coordmax,bornes,decoupe,epsilon,nbproc)

    !definition des domaines
    call definit_domaines(nbproc,data,domaines,bornes,decoupe)

    !ecriture des domaines decoupes
    call ecrit_domaines(data,nbproc,domaines)

    !definition des decoupages
    if ((data%interface==1).or.(nbproc==1)) then
       !decoupage avec interface
       call decoupe_interface(nbproc,data,ldat,ddat,domaines,epsilon)
    else
       !decoupage avec recouvrement
       call decoupe_recouvrement(nbproc,data,ldat,ddat,domaines)
    end if
    deallocate(domaines)

    !sauvegarde des decoupages
    call ecrit_decoupages(nbproc,data,ldat,ddat)

    !CALL MPI_ABORT(ierr)

    return
  end subroutine decoupedata

  !****************************************
  !definition des bornes avec interface
  subroutine definit_bornes(data,coordmin,coordmax,bornes,decoupe,epsilon,nbproc)
    implicit none
    integer :: nbproc
    real*8 :: epsilon
    type(type_data) :: data
    integer,dimension(:),pointer :: decoupe
    real*8,dimension(:),pointer :: coordmax,coordmin
    real*8,dimension(:,:,:),pointer :: bornes
    integer :: i,j
    double precision :: prod1,prod2,som1,prod
    character*30 :: num,files

    prod1=1.0
    prod2=0.0
    som1=1.0
    prod=1.0
    !volume max
    do i=1,data%dim
       !print *,'coordmax',coordmax(i)
       !print *,'coordmin',coordmin(i)
       prod=prod*(coordmax(i)-coordmin(i))
       !print *,'difference',(coordmax(i)-coordmin(i))
    end do
    !print *,'surface globale',data%dim
    files='diminterface'
    write(num,*),0
    num=adjustl(num)
    files=trim(files)//'.'//trim(num)
    ! len=len(trim(num))
    ! print *,numproc,'ecriture dim interface : '
    open(file=files,unit=20)
    do i=1,data%dim
       som1=som1*(decoupe(i)-1)
       prod2=prod2+(decoupe(i)-1)*prod/(coordmax(i)-coordmin(i))    
    end do
    !print *,'surface bande',prod2     
    write(20,*)  prod,epsilon*prod2-som1*(epsilon)**data%dim
    close(20)
 



 if ((data%coord==1).or.(data%geom==1).or.(data%seuil==1)) then
       !traitement en coordonnees ou image en coodonnees ou image seuillee
       allocate(bornes(data%dim,max(nbproc-data%interface,1),2));bornes(:,:,:)=0.0
       do i=1,data%dim
          coordmin(i)=coordmin(i)-epsilon*1.1
          coordmax(i)=coordmax(i)+epsilon*1.1
          do j=1,decoupe(i)
             bornes(i,j,1)=coordmin(i)+(j-1)*(coordmax(i)-coordmin(i))/decoupe(i)
             bornes(i,j,2)=coordmin(i)+j*(coordmax(i)-coordmin(i))/decoupe(i)
          enddo
          if (data%recouvrement==1) then
             !decoupage en mode interface
             do j=1,decoupe(i)
                bornes(i,j,1)=bornes(i,j,1)-epsilon
                bornes(i,j,2)=bornes(i,j,2)+epsilon
             enddo
          end if
          bornes(i,1,1)=coordmin(i)-0.01*abs(coordmin(i))
          bornes(i,decoupe(i),2)=coordmax(i)+0.01*abs(coordmax(i))
       enddo
    elseif (data%image==1) then
       !traitement pour decoupage en pixel d'images
       allocate(bornes(data%imgdim,max(nbproc-1,1),2));bornes(:,:,:)=0.0
       if ((data%imgdim/=2).and.(data%imgdim/=3)) then
#if aff
          print *
          print *,'format d images /= 2d, 3d pas supporte !!!!'
#endif
          stop
       end if
       do i=1,data%imgdim
          coordmin(i)=1.0-epsilon*1.1
          coordmax(i)=data%imgmap(i)+epsilon*1.1
          do j=1,decoupe(i)
             bornes(i,j,1)=coordmin(i)+(j-1)*(coordmax(i)-coordmin(i))/decoupe(i)
             bornes(i,j,2)=coordmin(i)+j*(coordmax(i)-coordmin(i))/decoupe(i)
          enddo
          if (data%recouvrement==1) then
             !decoupage en mode interface
             do j=1,decoupe(i)
                bornes(i,j,1)=bornes(i,j,1)-epsilon
                bornes(i,j,2)=bornes(i,j,2)+epsilon
             enddo
          end if
          bornes(i,1,1)=coordmin(i)-0.01*abs(coordmin(i))
          bornes(i,decoupe(i),2)=coordmax(i)+0.01*abs(coordmax(i))
       enddo
    end if
    return
  end subroutine definit_bornes

  !****************************************
  !definition des domaines de decoupages
  subroutine definit_domaines(nbproc,data,domaines,bornes,decoupe)
    implicit none
    integer :: nbproc
    type(type_data) :: data
    integer,dimension(:),pointer :: decoupe
    integer,dimension(:),pointer :: list
    real*8,dimension(:,:,:),pointer :: bornes,domaines
    integer :: k,n,ok
    if ((data%coord==1).or.(data%geom==1).or.(data%seuil==1)) then
       !traitement en coordonnees ou image en coodonnees ou image seuillee
       allocate(domaines(max(1,nbproc-data%interface),data%dim,2))
       domaines(:,:,:)=0.0
       allocate(list(data%dim)); list(:)=1
       if (nbproc>1) then
          !** mode >1 proc
          do n=1,nbproc-data%interface
             do k=1,data%dim
                domaines(n,k,:)=bornes(k,list(k),:)
             enddo
             ok=1
             do k=data%dim,1,-1
                if (ok==1) then
                   list(k)=list(k)+1
                   if (list(k)>decoupe(k)) then
                      list(k)=1
                   else
                      ok=0
                   endif
                end if
             end do
          enddo
       else
          !** mode 1 proc
          do k=1,data%dim
             domaines(1,k,:)=bornes(k,1,:)
          enddo
       endif
       deallocate(list)
    elseif (data%image==1) then
       !traitement pour decoupage en pixel d'images
       allocate(domaines(max(1,nbproc-data%interface),data%imgdim,2))
       domaines(:,:,:)=0.0
       allocate(list(data%imgdim)); list(:)=1
       if (nbproc>1) then
          !** mode >1 proc
          do n=1,nbproc-data%interface
             do k=1,data%imgdim
                domaines(n,k,:)=bornes(k,list(k),:)
             enddo
             ok=1
             do k=data%imgdim,1,-1
                if (ok==1) then
                   list(k)=list(k)+1
                   if (list(k)>decoupe(k)) then
                      list(k)=1
                   else
                      ok=0
                   endif
                end if
             end do
          enddo
       else
          !** mode 1 proc
          do k=1,data%imgdim
             domaines(1,k,:)=bornes(k,1,:)
          enddo
       endif
       deallocate(list)
    end if
    return
  end subroutine definit_domaines

  !****************************************
  !decoupage avec interface
  subroutine decoupe_interface(nbproc,data,ldat,ddat,domaines,epsilon)
    implicit none
    integer :: nbproc
    type(type_data) :: data
    real*8 :: epsilon
    integer,dimension(:),pointer :: ldat
    integer,dimension(:,:),pointer :: ddat
    real*8,dimension(:,:,:),pointer :: domaines
    integer :: i,j,n,ok,ierr
    allocate(ldat(0:max(1,nbproc-1))); ldat(:)=0
    allocate(ddat(0:max(1,nbproc-1),data%nb)); ddat(:,:)=0
    do i=1,data%nb
       !recherche des paquets
       n=0; ok=0
       do while(ok==0)
          n=n+1; ok=1
          if ((data%coord==1).or.(data%geom==1).or.(data%seuil==1)) then
             !traitement en coordonnees ou image en coodonnees ou image seuillee
             do j=1,data%dim
                if ((data%point(i)%coord(j)>domaines(n,j,2)).or.&
                     (data%point(i)%coord(j)<domaines(n,j,1))) ok=0
             end do
          elseif (data%image==1) then
             !traitement pour decoupage en pixel d'images
             do j=1,data%imgdim
                if ((data%refimg(i,j)>domaines(n,j,2)).or.&
                     (data%refimg(i,j)<domaines(n,j,1))) ok=0
             end do
          end if
          if ((n>nbproc-1).and.(nbproc>1)) then
#if aff
             print *,'bug dans le decoupage !',n,nbproc-1
#endif
             if (data%geom==0) then
#if aff
                print *,data%point(i)%coord(:)
#endif
             else
#if aff
                print *,i,data%refimg(i,:)
#endif
             end if
             call MPI_ABORT(ierr)
             stop
          end if
       end do
       ldat(n)=ldat(n)+1
       ddat(n,ldat(n))=i
       if (nbproc>1) then
          !** recherche de l'interface si plus de 1 proc
          ok=0
          if ((data%coord==1).or.(data%geom==1).or.(data%seuil==1)) then
             !traitement en coordonnees ou image en coodonnees ou image seuillee
             do j=1,data%dim
                if ((abs(data%point(i)%coord(j)-domaines(n,j,1))<epsilon).or.&
                     (abs(data%point(i)%coord(j)-domaines(n,j,2))<epsilon)) ok=1
             enddo
          elseif (data%image==1) then
             !traitement pour decoupage en pixel d'images
             do j=1,data%imgdim
                if ((abs(data%refimg(i,j)-domaines(n,j,1))<epsilon).or.&
                     (abs(data%refimg(i,j)-domaines(n,j,2))<epsilon)) ok=1
             enddo
          end if
          if (ok==1) then
             ldat(0)=ldat(0)+1
             ddat(0,ldat(0))=i
             write(7,*) ddat(0,ldat(0))
          endif
       end if
    end do
    write(7,*) ldat(0)
    return
  end subroutine decoupe_interface

  !****************************************
  !decoupage avec recouvrement
  subroutine decoupe_recouvrement(nbproc,data,ldat,ddat,domaines)
    implicit none
    integer :: nbproc
    type(type_data) :: data
    integer,dimension(:),pointer :: ldat
    integer,dimension(:,:),pointer :: ddat
    real*8,dimension(:,:,:),pointer :: domaines
    integer :: i,j,n,ok
    allocate(ldat(0:max(1,nbproc-1))); ldat(:)=0
    allocate(ddat(0:max(1,nbproc-1),data%nb)); ddat(:,:)=0
    do i=1,data%nb
       !recherche des paquets
       do n=1,nbproc
          ok=1
          if ((data%coord==1).or.(data%geom==1).or.(data%seuil==1)) then
             !traitement en coordonnees ou image en coodonnees ou image seuillee
             do j=1,data%dim
                if ((data%point(i)%coord(j)>domaines(n,j,2)).or.&
                     (data%point(i)%coord(j)<domaines(n,j,1))) ok=0
             end do
          elseif (data%image==1) then
             !traitement pour decoupage en pixel d'images
             do j=1,data%imgdim
                if ((data%refimg(i,j)>domaines(n,j,2)).or.&
                     (data%refimg(i,j)<domaines(n,j,1))) ok=0
             end do
          end if
          if (ok==1) then
             ldat(n-1)=ldat(n-1)+1
             ddat(n-1,ldat(n-1))=i
          end if
       end do
    end do
    return
  end subroutine decoupe_recouvrement

  !****************************************
  !elimine les doublons dans le clustering
  subroutine regroupe(nbclust,iclust,clustermap,data)
    implicit none
    integer :: nbclust
    type(type_data) :: data
    integer,dimension(:),pointer :: iclust
    integer,dimension(:,:),pointer :: clustermap
    integer :: ok,i,j,k,ok2,i2,j2,n,j3,ok3
    ok=0; i=1; j=0
#if aff
    print *,'  > elimination des doublons...'
    print *,'    > regroupement du sous-cluster ',1
#endif
    do while(ok==0)
       j=j+1 
       if (j>iclust(i)) then
          !la ligne i est completement testee
#if aff
          print *,'      > nb d elements apres regroupement :',iclust(i)
#endif
          i=i+1; j=1
#if aff
          print *,'    > regroupement du cluster ',i
#endif
       endif
       if (i>nbclust-1) then
          !plus de points a tester
          ok=1
       elseif (iclust(i)>0) then
          !stockage de l'indice
          data%point(clustermap(i,j))%cluster=i
          !test des recouvrement
          ok2=0; i2=i+1; j2=1
          do while(ok2==0)
             if (j2>iclust(i2)) then
                !ligne i2 completement teste pour le point (i,j)
                i2=i2+1; j2=1
             endif
             if (i2>nbclust) then
                !fin des tests pour le point (i,j)
                ok2=1
             else
                !test d'intersections :
                if (clustermap(i,j)==clustermap(i2,j2)) then
                   !intersection trouvee :
                   !ligne i2 ajoutee a la ligne i
                   n=0
                   do k=1,iclust(i2)
                      !test d'elimination de doublon
                      ok3=1
                      do j3=1,iclust(i)
                         if (clustermap(i2,k)==clustermap(i,j3)) ok3=0
                      enddo
                      if (ok3==1) then
                         n=n+1
                         clustermap(i,iclust(i)+n)=clustermap(i2,k)
                         clustermap(i2,k)=0                         
                      endif
                   enddo
                   iclust(i)=iclust(i)+n
                   iclust(i2)=0
                else
                   !test d'un nouveau point
                   j2=j2+1
                endif
             endif
          end do
       end if
    end do
#if aff
    print *,'      > nb d elements apres regroupement :',iclust(i)
#endif
    return
  end subroutine regroupe

end module module_decoupe
