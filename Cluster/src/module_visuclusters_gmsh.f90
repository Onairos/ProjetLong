module module_visuclusters_gmsh
  use module_visuclusters_structure
contains

  !************************
  !ecriture de la geometrie du decoupage
  subroutine ecrit_decoupage_gmsh(params)
    implicit none
    type(type_params) :: params
    integer :: i
    character*30 :: num
    real :: xmin,ymin,zmin,xmax,ymax,zmax
    real :: x0,y0,z0,x1,y1,z1
    print *,'-> decoupe.geo'
    open(file='fort.2',unit=2)
    open(file='decoupe.geo',unit=10)
    do i=1,params%nbproc-params%interface
       if (((params%coord==1).and.(params%dim==2)) &
          .or.((params%image==1).and.(params%imgdim==2)) &
          .or.((params%seuil==1).and.(params%imgdim==2)) &
          .or.((params%geom==1).and.(params%imgdim+params%imgt==2))) then
          !2D
          if (params%coord==1) then
             read(2,*) x0,y0,num,x1,y1
             xmin=x0
             ymin=y0
             xmax=x1
             ymax=y1
          elseif (params%seuil==1) then
             read(2,*) x0,num,x1
             xmin=1
             ymin=-1
             xmax=params%imgmap(2)
             ymax=-params%imgmap(1)
             zmin=x0
             zmax=x1
          elseif ((params%image==1).or.(params%geom==1)) then
             read(2,*) x0,y0,num,x1,y1
             xmin=y0
             ymin=-x0
             xmax=y1
             ymax=-x1
             zmin=x0
          end if
          if (params%seuil==1) then
             write(10,*) 'Point(',8*(i-1)+1,')={',xmin,',',ymin,',',zmin,'};'
             write(10,*) 'Point(',8*(i-1)+2,')={',xmax,',',ymin,',',zmin,'};'
             write(10,*) 'Point(',8*(i-1)+3,')={',xmax,',',ymax,',',zmin,'};'
             write(10,*) 'Point(',8*(i-1)+4,')={',xmin,',',ymax,',',zmin,'};'
             write(10,*) 'Point(',8*(i-1)+5,')={',xmin,',',ymin,',',zmax,'};'
             write(10,*) 'Point(',8*(i-1)+6,')={',xmax,',',ymin,',',zmax,'};'
             write(10,*) 'Point(',8*(i-1)+7,')={',xmax,',',ymax,',',zmax,'};'
             write(10,*) 'Point(',8*(i-1)+8,')={',xmin,',',ymax,',',zmax,'};'
             write(10,*) 'Line(',12*(i-1)+1,')={',8*(i-1)+1,',',8*(i-1)+2,'};'
             write(10,*) 'Line(',12*(i-1)+2,')={',8*(i-1)+2,',',8*(i-1)+3,'};'
             write(10,*) 'Line(',12*(i-1)+3,')={',8*(i-1)+3,',',8*(i-1)+4,'};'
             write(10,*) 'Line(',12*(i-1)+4,')={',8*(i-1)+4,',',8*(i-1)+1,'};'
             write(10,*) 'Line(',12*(i-1)+5,')={',8*(i-1)+5,',',8*(i-1)+6,'};'
             write(10,*) 'Line(',12*(i-1)+6,')={',8*(i-1)+6,',',8*(i-1)+7,'};'
             write(10,*) 'Line(',12*(i-1)+7,')={',8*(i-1)+7,',',8*(i-1)+8,'};'
             write(10,*) 'Line(',12*(i-1)+8,')={',8*(i-1)+8,',',8*(i-1)+5,'};'
             write(10,*) 'Line(',12*(i-1)+9,')={',8*(i-1)+1,',',8*(i-1)+5,'};'
             write(10,*) 'Line(',12*(i-1)+10,')={',8*(i-1)+2,',',8*(i-1)+6,'};'
             write(10,*) 'Line(',12*(i-1)+11,')={',8*(i-1)+3,',',8*(i-1)+7,'};'
             write(10,*) 'Line(',12*(i-1)+12,')={',8*(i-1)+4,',',8*(i-1)+8,'};'
          else
             write(10,*) 'Point(',4*(i-1)+1,')={',xmin,',',ymin,',0.};'
             write(10,*) 'Point(',4*(i-1)+2,')={',xmax,',',ymin,',0.};'
             write(10,*) 'Point(',4*(i-1)+3,')={',xmax,',',ymax,',0.};'
             write(10,*) 'Point(',4*(i-1)+4,')={',xmin,',',ymax,',0.};'
             write(10,*) 'Line(',4*(i-1)+1,')={',4*(i-1)+1,',',4*(i-1)+2,'};'
             write(10,*) 'Line(',4*(i-1)+2,')={',4*(i-1)+2,',',4*(i-1)+3,'};'
             write(10,*) 'Line(',4*(i-1)+3,')={',4*(i-1)+3,',',4*(i-1)+4,'};'
             write(10,*) 'Line(',4*(i-1)+4,')={',4*(i-1)+4,',',4*(i-1)+1,'};'
          end if
       elseif (((params%coord==1).and.(params%dim==3)) &
          .or.((params%image==1).and.(params%imgdim==3)) &
          .or.((params%geom==1).and.(params%imgdim+params%imgt==3))) then
          !3D
          if (params%coord==1) then
             read(2,*) x0,y0,z0,num,x1,y1,z1
             xmin=x0
             ymin=y0
             zmin=z0
             xmax=x1
             ymax=y1
             zmax=z1
          elseif ((params%image==1).or.(params%geom==1)) then
             read(2,*) x0,y0,z0,num,x1,y1,z1
             !xmin=-x0
             !ymin=y0
             !zmin=z0
             !xmax=-x1
             !ymax=y1
             !zmax=z1
             xmin=y0
             ymin=-x0
             zmin=z0
             xmax=y1
             ymax=-x1
             zmax=z1
          end if
          write(10,*) 'Point(',8*(i-1)+1,')={',xmin,',',ymin,',',zmin,'};'
          write(10,*) 'Point(',8*(i-1)+2,')={',xmax,',',ymin,',',zmin,'};'
          write(10,*) 'Point(',8*(i-1)+3,')={',xmax,',',ymax,',',zmin,'};'
          write(10,*) 'Point(',8*(i-1)+4,')={',xmin,',',ymax,',',zmin,'};'
          write(10,*) 'Point(',8*(i-1)+5,')={',xmin,',',ymin,',',zmax,'};'
          write(10,*) 'Point(',8*(i-1)+6,')={',xmax,',',ymin,',',zmax,'};'
          write(10,*) 'Point(',8*(i-1)+7,')={',xmax,',',ymax,',',zmax,'};'
          write(10,*) 'Point(',8*(i-1)+8,')={',xmin,',',ymax,',',zmax,'};'
          write(10,*) 'Line(',12*(i-1)+1,')={',8*(i-1)+1,',',8*(i-1)+2,'};'
          write(10,*) 'Line(',12*(i-1)+2,')={',8*(i-1)+2,',',8*(i-1)+3,'};'
          write(10,*) 'Line(',12*(i-1)+3,')={',8*(i-1)+3,',',8*(i-1)+4,'};'
          write(10,*) 'Line(',12*(i-1)+4,')={',8*(i-1)+4,',',8*(i-1)+1,'};'
          write(10,*) 'Line(',12*(i-1)+5,')={',8*(i-1)+5,',',8*(i-1)+6,'};'
          write(10,*) 'Line(',12*(i-1)+6,')={',8*(i-1)+6,',',8*(i-1)+7,'};'
          write(10,*) 'Line(',12*(i-1)+7,')={',8*(i-1)+7,',',8*(i-1)+8,'};'
          write(10,*) 'Line(',12*(i-1)+8,')={',8*(i-1)+8,',',8*(i-1)+5,'};'
          write(10,*) 'Line(',12*(i-1)+9,')={',8*(i-1)+1,',',8*(i-1)+5,'};'
          write(10,*) 'Line(',12*(i-1)+10,')={',8*(i-1)+2,',',8*(i-1)+6,'};'
          write(10,*) 'Line(',12*(i-1)+11,')={',8*(i-1)+3,',',8*(i-1)+7,'};'
          write(10,*) 'Line(',12*(i-1)+12,')={',8*(i-1)+4,',',8*(i-1)+8,'};'
       end if
    enddo
    close(10)
    close(2)
    return
  end subroutine ecrit_decoupage_gmsh

  !***********************
  !initialisation du fichier de decoupage
  subroutine affectation_gmsh(params)
    implicit none
    type(type_params) :: params
    real,dimension(:,:),pointer :: coord
    character*30 :: files,num
    integer :: i,j,k,nb,offset,totnum
    open(file='decoupe.visu',unit=1)
    write(1,*) 'View "MPI" {'
    !lecture des fichiers
    if (params%nbproc==1) then
       offset=1
       totnum=1
    else
       offset=0
       totnum=params%nbproc-1
    end if
    allocate(coord(1,params%dim))
    do i=offset,totnum
       !nom du fichier
       write(num,*) i
       files='decoupe.'//trim(adjustl(num))
       open(file=files,unit=10)
       read(10,*) nb
       print *,'  > ',i,' :',nb
       do j=1,nb
          if (params%coord==1) then
             !decoupage par coordonnees
             coord(:,:)=0.
             read(10,*) coord(1,:)
             call ecritpoint_gmsh(1,params%dim,coord,i,1)
          else
             !decoupage d'image 1D
             read (10,*) k
             call ecritpointimage_gmsh(1,params,i,k)
          end if
       enddo
       close(10)
    end do
    deallocate(coord)
    write(1,*) '};'
    close(1)
    print *,'-> decoupe.visu'
    return
  end subroutine affectation_gmsh

  !***********************
  !ecriture des clusters avant regroupement
  subroutine sous_clusters_gmsh(params)
    implicit none
    type(type_params) :: params
    integer :: i,j,k,nb,ind,lenn
    character*30 :: num,files
    real,dimension(:,:),pointer :: coord
    integer,dimension(:),pointer :: corresp
    do i=0,params%nbproc-1
       !nom du fichier
       !datas
       write(num,*) i
       num=adjustl(num)
       lenn=len(trim(num))
       files='cluster.partiel.'//trim(num)
       open(file=files,unit=10)
       !sortie
       files='cluster.partiel.'//num(1:lenn)//'.visu'
       open(file=files,unit=1)
       write(1,*) 'View "Clusters for part '//num(1:lenn)//'" {'
       read(10,*) nb,k
       print *,'  > ',i,' :',nb,' -> ',files
       allocate(coord(1,k))
       if ((params%image==1).or.(params%geom==1).or.(params%seuil==1)) then
          !lecture des correspondances
          files='decoupe.'//trim(num)
          open(file=files,unit=11)
          read(11,*)
          allocate(corresp(nb))
          do j=1,nb
             read(11,*) corresp(j)
          end do
          close(11)
       end if
       do j=1,nb
          if (params%coord==1) then
             read(10,*) coord(1,:),k
             call ecritpoint_gmsh(1,params%dim,coord,k,1)
          else
             read(10,*) k,ind
             call ecritpointimage_gmsh(1,params,ind,corresp(k))
          end if
       enddo
       close(10); 
       write(1,*) '};'
       close(1)
       deallocate(coord)
       if ((params%image==1).or.(params%geom==1).or.(params%seuil==1)) then
          deallocate(corresp)
       end if
    enddo
    return
  end subroutine sous_clusters_gmsh

  !***********************
  !ecriture des clusters apres regroupement
  subroutine cluster_final_gmsh(params)
    implicit none
    type(type_params) :: params
    integer :: i,j,k,nb,nb0
    character*30 :: num,files
    real,dimension(:,:),pointer :: coord
    open(file=params%mesh,unit=1)
    if (params%coord==1) then
       !lecture au format coordonnes
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
    open(file='cluster.final.visu',unit=1)
    write(1,*) 'View "Clusters" {'
    !lecture des fichiers
    do i=1,params%nbclusters
       !nom du fichier
       write(num,*) i
       files='cluster.final.'//trim(adjustl(num))
       open(file=files,unit=10)
       read(10,*) nb
       print *,'  > ',i,' :',nb
       do j=1,nb
          read(10,*) k
          if (params%coord==1) then
             !coordonnees classiques
             call ecritpoint_gmsh(1,nb0,coord,i,k)
          else
             !reassemblage des images
             call ecritpointimage_gmsh(1,params,i,k)
          end if
       end do
       close(10)
    end do
    write(1,*) '};'
    close(1)
    print *,'-> cluster.final.visu'
    if (params%coord==1) deallocate(coord)
    return
  end subroutine cluster_final_gmsh

  !*************************
  !subroutine ecriture de point
  subroutine ecritpoint_gmsh(unit,dim,coord,ind,k)
    implicit none
    integer :: unit,dim,ind,k
    real,dimension(:,:),pointer :: coord
    if (dim==2) then
       !2D
       write(unit,*) 'SP(',coord(k,1),',',coord(k,2),',',0.,'){',ind,'};' 
    elseif (dim==3) then
       !3D
       write(unit,*) 'SP(',coord(k,1),',',coord(k,2),',',coord(k,3),'){',ind,'};'
    end if
    return
  end subroutine ecritpoint_gmsh

  !*************************
  !subroutine ecriture de point en format image
  subroutine ecritpointimage_gmsh(unit,params,ind,k)
    implicit none
    type(type_params) :: params
    integer :: unit,ind,k,ix,iy,i
    real :: kx,ky,kz
    real,dimension(:),pointer :: data
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
    !coordonnees
    ix=params%refimg(k,1)
    iy=params%refimg(k,2)
    if (params%imgdim==2) then
       !points en 2D
       if (params%geom==1) then
          !coordonnees redimensionnees
          kx=iy*params%pas(2)
          ky=-ix*params%pas(1)
       else
          kx=float(iy)
          ky=-float(ix)
       end if
       if (((params%image==1).or.(params%geom==1).or.(params%seuil==1)) &
            .and.(params%imgdim==2)) then
          kz=data(k)
       else
          kz=0.
       end if
    elseif (params%imgdim==3) then
       !points en 3D
       if (params%geom==1) then
          kx=iy*params%pas(2)
          ky=-ix*params%pas(1)
          kz=float(params%refimg(k,3))*params%pas(3)
       else
          !kx=-float(ix)
          !ky=float(iy)
          !kz=float(params%refimg(k,3))
          kx=float(iy)
          ky=-float(ix)
          kz=float(params%refimg(k,3))
       end if
    end if
    !ecriture
    write(unit,*) 'SP(',kx,',',ky,',',kz,'){',ind,'};'
    return
  end subroutine ecritpointimage_gmsh

  !************************
  !liste des commandes
  subroutine commandes_gmsh
    print *,'gmsh decoupe.visu'
    print *,'gmsh decoupe.geo'
    print *,'gmsh decoupe.visu decoupe.geo'
    print *,'gmsh cluster.partiel.*.visu'
    print *,'gmsh cluster.final.visu'
    print *,'gmsh decoupe.geo cluster.partiel.*.visu'
    print *,'gmsh decoupe.geo cluster.final.visu'
    return
  end subroutine commandes_gmsh


end module module_visuclusters_gmsh
