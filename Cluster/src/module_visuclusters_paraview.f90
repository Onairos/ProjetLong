module module_visuclusters_paraview
  use module_visuclusters_structure
contains

  !************************
  !ecriture de la geometrie du decoupage
  subroutine ecrit_decoupage_paraview(params)
    implicit none
    type(type_params) :: params
    integer :: i
    character*30 :: num
    real,pointer,dimension(:) :: xmin,ymin,zmin,xmax,ymax,zmax
    real :: x0,y0,z0,x1,y1,z1
    !lecture
    open(file='fort.2',unit=2)
    allocate(xmin(params%nbproc))
    allocate(xmax(params%nbproc))
    allocate(ymin(params%nbproc))
    allocate(ymax(params%nbproc))
    allocate(zmin(params%nbproc))
    allocate(zmax(params%nbproc))
    do i=1,params%nbproc-params%interface
       if (((params%coord==1).and.(params%dim==2)) &
          .or.((params%image==1).and.(params%imgdim==2)) &
          .or.((params%seuil==1).and.(params%imgdim==2)) &
          .or.((params%geom==1).and.(params%imgdim+params%imgt==2))) then
          !2D
          if (params%coord==1) then
             read(2,*) x0,y0,num,x1,y1
             xmin(i)=x0
             ymin(i)=y0
             xmax(i)=x1
             ymax(i)=y1
          elseif (params%seuil==1) then
             read(2,*) x0,num,x1
             xmin(i)=1
             ymin(i)=-1
             xmax(i)=params%imgmap(2)
             ymax(i)=-params%imgmap(1)
             zmin(i)=x0
             zmax(i)=x1
          elseif ((params%image==1).or.(params%geom==1)) then
             read(2,*) x0,y0,num,x1,y1
             xmin(i)=y0
             ymin(i)=-x0
             xmax(i)=y1
             ymax(i)=-x1
             zmin(i)=x0
          end if
       elseif (((params%coord==1).and.(params%dim==3)) &
          .or.((params%image==1).and.(params%imgdim==3)) &
          .or.((params%geom==1).and.(params%imgdim+params%imgt==3))) then
          !3D
          if (params%coord==1) then
             read(2,*) x0,y0,z0,num,x1,y1,z1
             xmin(i)=x0
             ymin(i)=y0
             zmin(i)=z0
             xmax(i)=x1
             ymax(i)=y1
             zmax(i)=z1
          elseif ((params%image==1).or.(params%geom==1)) then
             read(2,*) x0,y0,z0,num,x1,y1,z1
             !xmin(i)=-x0
             !ymin(i)=y0
             !zmin(i)=z0
             !xmax(i)=-x1
             !ymax(i)=y1
             !zmax(i)=z1
             xmin(i)=y0
             ymin(i)=-x0
             zmin(i)=z0
             xmax(i)=y1
             ymax(i)=-x1
             zmax(i)=z1
          end if
       end if
    enddo
    close(2)
    !ecriture
    print *,'-> visu/decoupe.geo'
    print *,'-> visu/decoupe.indices'
    open(file='visu/decoupe.geo',unit=10)
    open(file='visu/decoupe.indices',unit=11)
    write(10,*) '** sortie de visuclusters **'
    write(10,*) '** decoupage des sous-clusters **'
    write(10,'(a)') 'node id assign'
    write(10,'(a)') 'element id assign'
    write(11,*) '** indices des process **'
    write(10,'(a)') 'part'
    write(10,*) 1
    write(10,*) '** decoupages **'
    write(10,'(a)') 'coordinates'
    if (((params%coord==1).and.(params%dim==2)) &
         .or.((params%image==1).and.(params%imgdim==2)) &
         .or.((params%seuil==1).and.(params%imgdim==2)) &
         .or.((params%geom==1).and.(params%imgdim+params%imgt==2))) then
       !2D
       if ((params%coord==1).or.(params%geom==1).or.(params%image==1)) then
          write(10,*) 4*(params%nbproc-params%interface)
          do i=1,params%nbproc-params%interface
             write(10,*) xmin(i)
             write(10,*) xmax(i)
             write(10,*) xmax(i)
             write(10,*) xmin(i)
          end do
          do i=1,params%nbproc-params%interface
             write(10,*) ymin(i)
             write(10,*) ymin(i)
             write(10,*) ymax(i)
             write(10,*) ymax(i)
          end do
          do i=1,params%nbproc-params%interface
             write(10,*) 0.
             write(10,*) 0.
             write(10,*) 0.
             write(10,*) 0.
          end do
          write(10,'(a)') 'quad4'
          write(10,*) params%nbproc-params%interface
          write(11,'(a)') 'part'
          write(11,*) 1
          write(11,'(a)') 'quad4'
          do i=1,params%nbproc-params%interface
             write(10,*) 4*(i-1)+1,4*(i-1)+2,4*(i-1)+3,4*(i-1)+4
             write(11,*) i
          end do
       elseif (params%seuil==1) then
          write(10,*) 8*(params%nbproc-params%interface)
          do i=1,params%nbproc-params%interface
             write(10,*) xmin(i)
             write(10,*) xmax(i)
             write(10,*) xmax(i)
             write(10,*) xmin(i)
             write(10,*) xmin(i)
             write(10,*) xmax(i)
             write(10,*) xmax(i)
             write(10,*) xmin(i)
          end do
          do i=1,params%nbproc-params%interface
             write(10,*) ymin(i)
             write(10,*) ymin(i)
             write(10,*) ymax(i)
             write(10,*) ymax(i)
             write(10,*) ymin(i)
             write(10,*) ymin(i)
             write(10,*) ymax(i)
             write(10,*) ymax(i)
          end do
          do i=1,params%nbproc-params%interface
             write(10,*) zmin(i)
             write(10,*) zmin(i)
             write(10,*) zmin(i)
             write(10,*) zmin(i)
             write(10,*) zmax(i)
             write(10,*) zmax(i)
             write(10,*) zmax(i)
             write(10,*) zmax(i)
          end do
          write(10,'(a)') 'hexa8'
          write(10,*) params%nbproc-params%interface
          write(11,'(a)') 'part'
          write(11,*) 1
          write(11,'(a)') 'hexa8'
          do i=1,params%nbproc-params%interface
             write(10,*) 8*(i-1)+1,8*(i-1)+2,8*(i-1)+3,8*(i-1)+4,&
                  8*(i-1)+5,8*(i-1)+6,8*(i-1)+7,8*(i-1)+8
             write(11,*) i
          end do
       end if
    elseif (((params%coord==1).and.(params%dim==3)) &
         .or.((params%image==1).and.(params%imgdim==3)) &
         .or.((params%geom==1).and.(params%imgdim+params%imgt==3))) then
       !3D
       if ((params%coord==1).or.(params%geom==1).or.(params%image==1)) then
          write(10,*) 8*(params%nbproc-params%interface)
          do i=1,params%nbproc-params%interface
             write(10,*) xmin(i)
             write(10,*) xmax(i)
             write(10,*) xmax(i)
             write(10,*) xmin(i)
             write(10,*) xmin(i)
             write(10,*) xmax(i)
             write(10,*) xmax(i)
             write(10,*) xmin(i)
          end do
          do i=1,params%nbproc-params%interface
             write(10,*) ymin(i)
             write(10,*) ymin(i)
             write(10,*) ymax(i)
             write(10,*) ymax(i)
             write(10,*) ymin(i)
             write(10,*) ymin(i)
             write(10,*) ymax(i)
             write(10,*) ymax(i)
          end do
          do i=1,params%nbproc-params%interface
             write(10,*) zmin(i)
             write(10,*) zmin(i)
             write(10,*) zmin(i)
             write(10,*) zmin(i)
             write(10,*) zmax(i)
             write(10,*) zmax(i)
             write(10,*) zmax(i)
             write(10,*) zmax(i)
          end do
          write(10,'(a)') 'hexa8'
          write(10,*) params%nbproc-params%interface
          write(11,'(a)') 'part'
          write(11,*) 1
          write(11,'(a)') 'hexa8'
          do i=1,params%nbproc-params%interface
             write(10,*) 8*(i-1)+1,8*(i-1)+2,8*(i-1)+3,8*(i-1)+4,&
                  8*(i-1)+5,8*(i-1)+6,8*(i-1)+7,8*(i-1)+8
             write(11,*) i
          end do
       end if
    end if
    close(10)
    close(11)
    deallocate(xmin)
    deallocate(xmax)
    deallocate(ymin)
    deallocate(ymax)
    deallocate(zmin)
    deallocate(zmax)
    !fichier maitre
    print *,'-> decoupe.case'
    open(file='decoupe.case',unit=12)
    write(12,'(a)') 'FORMAT'
    write(12,'(a)') 'type: ensight gold'
    write(12,*) 
    write(12,'(a)') 'GEOMETRY'
    write(12,'(a)') 'model:   visu/decoupe.geo'
    write(12,*) 
    write(12,'(a)') 'VARIABLE'
    write(12,'(a)') 'scalar per element:   process   visu/decoupe.indices'
    close(12)
    return
  end subroutine ecrit_decoupage_paraview

  !***********************
  !initialisation du fichier de decoupage
  subroutine affectation_paraview(params)
    implicit none
    type(type_params) :: params
    real,dimension(:,:),pointer :: coord
    character*30 :: files,num
    integer :: i,j,nb,offset,totnum
    integer,dimension(:),pointer :: ind,indp
    !lecture des fichiers
    if (params%nbproc==1) then
       offset=1
       totnum=1
    else
       offset=0
       totnum=params%nbproc-1
    end if
    do i=offset,totnum
       if ((i==0).and.(params%interface==1).and.(params%nbproc>1)) then
          !fichier a part pour l'interface
          print *,'-> visu/affectation-interface.geo'
          print *,'-> visu/affectation-interface.indices'
          open(file='visu/affectation-interface.geo',unit=10)
          open(file='visu/affectation-interface.indices',unit=11)
          write(10,*) '** sortie de visuclusters **'
          write(10,*) '** points sur l interface **'
          write(10,'(a)') 'node id assign'
          write(10,'(a)') 'element id assign'
          write(11,*) '** indices des process **'
       elseif (((i==1).and.(params%interface==1)).or. &
            ((i==0).and.(params%interface==0))) then
          !fichier general pour les autres sous-domaines
          print *,'-> visu/affectation.geo'
          print *,'-> visu/affectation.indices'
          open(file='visu/affectation.geo',unit=10)
          open(file='visu/affectation.indices',unit=11)
          write(10,*) '** sortie de visuclusters **'
          write(10,*) '** decoupage des sous-clusters **'
          write(10,'(a)') 'node id assign'
          write(10,'(a)') 'element id assign'
          write(11,*) '** indices des process **'
       end if
       !nom du fichier
       write(num,*) i
       files='decoupe.'//trim(adjustl(num))
       open(file=files,unit=20)
       read(20,*) nb
       print *,'  > ',i,' :',nb
        allocate(coord(nb,params%dim))
       coord(:,:)=0.
       allocate(ind(nb)); ind(:)=0
       allocate(indp(nb)); indp(:)=0
       if (nb>0) then
          do j=1,nb
             if (params%coord==1) then
                !decoupage par coordonnees
                read(20,*) coord(j,:)
                ind(j)=i
             else
                !decoupage d'image 1D
                read (20,*) indp(j)
                ind(j)=i
             end if
          enddo
          !ecriture
          if (params%coord==1) then
             !decoupage par coordonnees
             call ecritpoint_paraview(10,11,nb,params%dim,coord,ind,1)
          else
             !decoupage d'image 1D
             call ecritpointimage_paraview(10,11,nb,params,ind,indp)
          end if
       end if
       deallocate(coord)
       deallocate(ind)
       close(20)
       if ((i==0).and.(params%interface==1)) then
          !fermeture des fichiers de l'interface
          close(10)
          close(11)
       end if
    end do
    close(10)
    close(11)
    !main de l'interface
    if ((params%interface==1).and.(params%nbproc>1)) then
       print *,'-> affectation-interface.case'
       open(file='affectation-interface.case',unit=12)
       write(12,'(a)') 'FORMAT'
       write(12,'(a)') 'type: ensight gold'
       write(12,*) 
       write(12,'(a)') 'GEOMETRY'
       write(12,'(a)') 'model:   visu/affectation-interface.geo'
       write(12,*) 
       write(12,'(a)') 'VARIABLE'
       write(12,'(a)') 'scalar per node:   process   visu/affectation-interface.indices'
       close(12)
    end if
    !main des autres sous-domaines
    print *,'-> affectation.case'
    open(file='affectation.case',unit=12)
    write(12,'(a)') 'FORMAT'
    write(12,'(a)') 'type: ensight gold'
    write(12,*) 
    write(12,'(a)') 'GEOMETRY'
    write(12,'(a)') 'model:   visu/affectation.geo'
    write(12,*) 
    write(12,'(a)') 'VARIABLE'
    write(12,'(a)') 'scalar per node:   process   visu/affectation.indices'
    close(12)
    return
  end subroutine affectation_paraview

  !***********************
  !ecriture des clusters avant regroupement
  subroutine sous_clusters_paraview(params)
    implicit none
    type(type_params) :: params
    integer :: i,j,k,nb,lenn,nbstar
    character*30 :: num,files,star
    real,dimension(:,:),pointer :: coord
    integer,dimension(:),pointer :: corresp,ind,indp
    !nb de characteres generiques a utiliser
    nbstar=floor(log(real(params%nbproc-1))/log(real(10)))+1
    do i=1,nbstar
       star(i:i)='0'
    end do
    do i=0,params%nbproc-1
       !nom du fichier
       !datas
       write(num,*) i
       num=adjustl(num)
       lenn=len(trim(num))
       files='cluster.partiel.'//trim(num)
       open(file=files,unit=20)
       !sortie
       !fichier general pour les autres sous-domaines
       files='cluster.partiel.'//star(1:nbstar-lenn)//num(1:lenn)
       print *,'-> visu/'//trim(files)//'.geo'
       print *,'-> visu/'//trim(files)//'.indices'
       open(file='visu/'//trim(files)//'.geo',unit=10)
       open(file='visu/'//trim(files)//'.indices',unit=11)
       write(10,*) '** sortie de visuclusters **'
       write(10,*) '** decoupage du sous-clusters '//trim(num)//' **'
       write(10,'(a)') 'node id assign'
       write(10,'(a)') 'element id assign'
       write(11,*) '** indices des elements clusterises **'
       !lecture
       read(20,*) nb,k
       print *,'  > ',i,' :',nb,' -> ',files
       allocate(coord(nb,k))
       allocate(ind(max(1,nb))); ind(:)=0
       allocate(indp(max(1,nb))); indp(:)=0
       if ((params%image==1).or.(params%geom==1).or.(params%seuil==1)) then
          !lecture des correspondances
          files='decoupe.'//trim(num)
          open(file=files,unit=21)
          read(21,*)
          allocate(corresp(nb))
          do j=1,nb
             read(21,*) corresp(j)
          end do
          close(21)
       end if
       do j=1,nb
          if (params%coord==1) then
             read(20,*) coord(j,:),ind(j)
          else
             read(20,*) indp(j),ind(j)
             indp(j)=corresp(indp(j))
          end if
       enddo
       close(20); 
       !ecriture
       if (params%coord==1) then
          !decoupage par coordonnees
          call ecritpoint_paraview(10,11,nb,params%dim,coord,ind,1)
       else
          !decoupage d'image 1D
          call ecritpointimage_paraview(10,11,nb,params,ind,indp)
       end if
       deallocate(coord)
       deallocate(ind)
       deallocate(indp)
       close(10)
       close(11)
       if ((params%image==1).or.(params%geom==1).or.(params%seuil==1)) then
          deallocate(corresp)
       end if
    enddo
    !main
    print *,'-> cluster.partiel.case'
    do i=1,nbstar
       star(i:i)='*'
    end do
    open(file='cluster.partiel.case',unit=12)
    write(12,'(a)') 'FORMAT'
    write(12,'(a)') 'type: ensight gold'
    write(12,*) 
    write(12,'(a)') 'GEOMETRY'
    write(12,'(a)') 'model: 1 visu/cluster.partiel.'//star(1:nbstar)//'.geo'
    write(12,*) 
    write(12,'(a)') 'VARIABLE'
    write(12,'(a)') 'scalar per node: 1 cluster visu/cluster.partiel.'//&
         star(1:nbstar)//'.indices'
    write(12,*) 
    write(12,'(a)') 'TIME'
    write(12,'(a)') 'time set:         1'
    write(num,*) params%nbproc
    write(12,'(a)') 'number of steps: '//trim(adjustl(num))
    write(12,'(a)') 'filename start number: 0'
    write(12,'(a)') 'filename increment: 1'
    write(12,'(a)') 'time values:'
    do i=0,params%nbproc-1
       write(12,*) i
    end do
    close(12)
    return
  end subroutine sous_clusters_paraview

  !***********************
  !ecriture des clusters apres regroupement
  subroutine cluster_final_paraview(params)
    implicit none
    type(type_params) :: params
    integer :: i,j,k,nb,nb0
    character*30 :: num,files
    real,dimension(:,:),pointer :: coord
    integer,dimension(:),pointer :: ind,indp
    if (params%coord==1) then
       open(file=params%mesh,unit=1)
       read(1,*) j,k
       print *,'lecture du fichier de maillage...',j,k
       allocate(coord(j,k))
       nb0=0
       do i=1,j
          read(1,*,end=100) coord(i,:)
          nb0=nb0+1
       enddo
100    print *,'nb de points ',nb0
       j=nb0
       close(1)
    end if

    nb0=k !dim des points
    !sortie
    print *,'-> visu/cluster.final.geo'
    print *,'-> visu/cluster.final.indices'
    open(file='visu/cluster.final.geo',unit=10)
    open(file='visu/cluster.final.indices',unit=11)
    write(10,*) '** sortie de visuclusters **'
    write(10,*) '** decoupage cluster final **'
    write(10,'(a)') 'node id assign'
    write(10,'(a)') 'element id assign'
    write(11,*) '** clusters **'
    allocate(ind(max(1,params%nbp))); ind(:)=0
    allocate(indp(max(1,params%nbp))); indp(:)=0
    !lecture des fichiers
    do i=1,params%nbclusters
       !nom du fichier
       write(num,*) i
       files='cluster.final.'//trim(adjustl(num))
       open(file=files,unit=20)
       read(20,*) nb
       print *,'  > ',i,' :',nb
       do j=1,nb
          read(20,*) k
          ind(k)=i
          indp(k)=k
       end do
       close(20)
    end do
    if (params%coord==1) then
       !coordonnees classiques
       call ecritpoint_paraview(10,11,params%nbp,params%dim,coord,ind,1)
    else
       !reassemblage des images
       call ecritpointimage_paraview(10,11,params%nbp,params,ind,indp)
    end if
    close(10)
    close(11)
    deallocate(ind)
    deallocate(indp)
    if (params%coord==1) deallocate(coord)
    !fichier maitre
    print *,'-> cluster.final.case'
    open(file='cluster.final.case',unit=12)
    write(12,'(a)') 'FORMAT'
    write(12,'(a)') 'type: ensight gold'
    write(12,*) 
    write(12,'(a)') 'GEOMETRY'
    write(12,'(a)') 'model:   visu/cluster.final.geo'
    write(12,*) 
    write(12,'(a)') 'VARIABLE'
    write(12,'(a)') 'scalar per element:   cluster   visu/cluster.final.indices'
    close(12)
    return
  end subroutine cluster_final_paraview

  !*************************
  !subroutine ecriture de points
  subroutine ecritpoint_paraview(unitgeo,unitind,nb,dim,coord,ind,k)
    implicit none
    integer :: unitgeo,unitind,k,nb,i,dim
    real,dimension(:,:),pointer :: coord
    integer,dimension(:),pointer :: ind
    write(unitgeo,'(a)') 'part'
    write(unitgeo,*) ind(k)
    write(unitgeo,*) '** decoupages **'
    write(unitgeo,'(a)') 'coordinates'
    write(unitgeo,*) nb
    write(unitind,'(a)') 'part'
    write(unitind,*) ind(k)
    write(unitind,'(a)') 'point'
    do i=1,nb
       write(unitgeo,*) coord(i,1)
       write(unitind,*) ind(i)
    end do
    do i=1,nb
       write(unitgeo,*) coord(i,2)
    end do
    do i=1,nb
       if (dim==2) then
          write(unitgeo,*) 0.
       else
          write(unitgeo,*) coord(i,3)
       end if
    end do
    write(unitgeo,'(a)') 'point'
    write(unitgeo,*)nb
    do i=1,nb
       write(unitgeo,*) i
    end do
    return
  end subroutine ecritpoint_paraview

  !*************************
  !subroutine ecriture de points en format image
  subroutine ecritpointimage_paraview(unitgeo,unitind,nbp,params,ind,indp)
    implicit none
    integer :: nbp,unitgeo,unitind,i,k,ix,iy
    type(type_params) :: params
    integer,dimension(:),pointer :: indp,ind
    real,dimension(:),pointer :: kx,ky,kz,data
    allocate(kx(nbp)); kx(:)=0
    allocate(ky(nbp)); ky(:)=0
    allocate(kz(nbp)); kz(:)=0
    if (((params%image==1).or.(params%geom==1).or.(params%seuil==1)) &
         .and.(params%imgdim==2)) then
       allocate(data(params%nbp))
       open(file=params%mesh,unit=50)
       read(50,*)
       read(50,*)
       do i=1,params%nbp
          read(50,*) data(i)
       end do
       close(50)
    end if
    !recherche des points
    do i=1,nbp
       k=indp(i)
       ix=params%refimg(k,1)
       iy=params%refimg(k,2)
       if (params%imgdim==2) then
          !points en 2D
          if (params%geom==1) then
             !coordonnees redimensionnees
             kx(i)=iy*params%pas(2)
             ky(i)=-ix*params%pas(1)
          else
             kx(i)=float(iy)
             ky(i)=-float(ix)
          end if
          if (((params%image==1).or.(params%geom==1).or.(params%seuil==1)) &
               .and.(params%imgdim==2)) then
             kz(i)=data(k)
          else
             kz(i)=0.
          end if
       elseif (params%imgdim==3) then
          !points en 3D
          if (params%geom==1) then
             kx(i)=iy*params%pas(2)
             ky(i)=-ix*params%pas(1)
             kz(i)=float(params%refimg(k,3))*params%pas(3)
          else
             !kx(i)=-float(ix)
             !ky(i)=float(iy)
             !kz(i)=float(params%refimg(k,3))
             kx(i)=float(iy)
             ky(i)=-float(ix)
             kz(i)=float(params%refimg(k,3))
          end if
       end if
    end do
    !ecriture
    write(unitgeo,'(a)') 'part'
    write(unitgeo,*) ind(1)
    write(unitgeo,*) '** decoupages **'
    write(unitgeo,'(a)') 'coordinates'
    write(unitgeo,*) nbp
    write(unitind,'(a)') 'part'
    write(unitind,*) ind(1)
    write(unitind,'(a)') 'point'
    do i=1,nbp
       write(unitgeo,*) kx(i)
       write(unitind,*) ind(i)
    end do
    do i=1,nbp
       write(unitgeo,*) ky(i)
    end do
    do i=1,nbp
       write(unitgeo,*) kz(i)
    end do
    write(unitgeo,'(a)') 'point'
    write(unitgeo,*)nbp
    do i=1,nbp
       write(unitgeo,*) i
    end do
    deallocate(kx)
    deallocate(ky)
    deallocate(kz)
    if (((params%image==1).or.(params%geom==1).or.(params%seuil==1)) &
         .and.(params%imgdim==2)) deallocate(data)
    return
  end subroutine ecritpointimage_paraview

  !************************
  !liste des commandes
  subroutine commandes_paraview
    print *,'paraview --data=decoupe.case'
    print *,'paraview --data=affectation.case'
    print *,'paraview --data=affectation-interface.case'
    print *,'paraview --data=cluster.partiel.case'
    print *,'paraview --data=cluster.final.case'
    return
  end subroutine commandes_paraview
  
end module module_visuclusters_paraview
