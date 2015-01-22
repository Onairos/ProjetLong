module module_sortie
  use module_structure
contains

  !********************************
  !ecriture des domaines decoupes
  subroutine ecrit_domaines(data,nbproc,domaines)
    implicit none
    type(type_data) :: data
    integer :: nbproc
    real*8,dimension(:,:,:),pointer :: domaines
    integer :: i
    do i=1,nbproc-data%interface
       write(2,*) domaines(i,:,1),'|', domaines(i,:,2)
    end do
    call flush(2)
    close(2)
    return
  end subroutine ecrit_domaines

  !********************************
  !ecriture des decoupages
  subroutine ecrit_decoupages(nbproc,data,ldat,ddat)
    implicit none
    integer :: nbproc
    type(type_data) :: data
    integer,dimension(:),pointer :: ldat
    integer,dimension(:,:),pointer :: ddat
    integer :: i,j,offset,nbdom
    character*30 :: files,num
    print *,'  > bilan decoupage :'
    offset=1; nbdom=nbproc
    if ((data%interface==1).and.(nbproc>1)) then
       offset=0; nbdom=nbproc-1
    end if
    if (data%recouvrement==1) then
       offset=0; nbdom=nbproc-1
    end if
    do i=offset,nbdom
       print *,'    > zone ',i,':',ldat(i)
       !nom du fichier
       write(num,*) i
       files='decoupe.'//trim(adjustl(num))
       open(file=files,unit=10)
       write(10,*) ldat(i)
       do j=1,ldat(i)
          if (data%coord==1) then
             !ecriture en coordonnees
             write(10,*) data%point(ddat(i,j))%coord(:)
          elseif ((data%image==1).or.(data%seuil==1).or.(data%geom==1)) then
             !ecriture en image
             write(10,*) ddat(i,j)
          end if
       enddo
       call flush(10)
       close(10)
    enddo
    call flush(6)
    return
  end subroutine ecrit_decoupages

  !********************************
  !ecriture des clusters regroupes
  subroutine ecritcluster(numproc,dataw)
    implicit none
    integer :: numproc
    type(type_data) :: dataw
    integer :: i
    character*30 :: files,num
    !nom du fichier
    write(num,*),numproc
    num=adjustl(num)
    files='cluster.partiel.'//trim(num)
    !len=len(trim(num))
    print *,numproc,'ecriture des clusters : ',files
    open(file=files,unit=10)
    write(10,*) dataw%nb,dataw%dim
    do i=1,dataw%nb
       if (dataw%coord==1) then
          write(10,*) dataw%point(i)%coord(:),dataw%point(i)%cluster
       else
          write(10,*) i,dataw%point(i)%cluster
       end if
    enddo
    call flush(10)
    close(10)
    return
  end subroutine ecritcluster

  !****************************
  !ecriture de cluster.final.
  subroutine ecritclusterfinal(nbclust,iclust,clustermap)
    implicit none
    integer :: nbclust
    integer,dimension(:),pointer :: iclust
    integer,dimension(:,:),pointer :: clustermap
    integer :: i,j,k
    character*30 :: files,num
    print *,'  > Ecriture du resultat...'
    k=0
    do i=1,nbclust
       if (iclust(i)>0) then
          k=k+1
          write(num,*) k
          files='cluster.final.'//trim(adjustl(num))
          open(file=files,unit=20)     
          print *,'    > cluster ',k,' :',iclust(i),' -> ',files
          write(20,*) iclust(i)
          do j=1,iclust(i)
             write(20,*) clustermap(i,j)
          enddo
          call flush(20)
          close(20)
       end if
    end do
    nbclust=k
    return
  end subroutine ecritclusterfinal

  !***************************
  !ecriture des informations
  subroutine ecrit_info(mesh,data,nbproc,nbclust)
    implicit none
    character*30 :: mesh
    integer :: nbproc,nbclust
    type(type_data) :: data
    !ecriture du fichier fort.3
    write(3,*) '# fichier de maillage :'
    write(3,*) mesh
    write(3,*) '#nb de points :'
    write(3,*) data%nb
    write(3,*) '# dimension :'
    write(3,*) data%dim
    write(3,*) '# nb de proc :'
    write(3,*) nbproc
    write(3,*) '# decoupage par interface :'
    write(3,*) data%interface
    write(3,*) '# decoupage par recouvrement :'
    write(3,*) data%recouvrement
    write(3,*) '# nb de clusters :'
    write(3,*) nbclust
    write(3,*) '# format coord :'
    write(3,*) data%coord
    write(3,*) '# format image :'
    write(3,*) data%image
    write(3,*) '# format geom :'
    write(3,*) data%geom
    write(3,*) '# format seuil :'
    write(3,*) data%seuil
    if ((data%image==1).or.(data%geom==1).or.(data%seuil==1)) then
       write(3,*) '# dimension :'
       write(3,*) data%imgdim
       write(3,*) '# decoupage :'
       write(3,*) data%imgmap(:)
       write(3,*) '# nb de temps :'
       write(3,*) data%imgt
       if (data%geom==1) then
          write(3,*) '## pas de maillage :'
          write(3,*) data%pas(:)
       end if
    end if
    call flush(3)
    close(3)
    return
  end subroutine ecrit_info

end module module_sortie
