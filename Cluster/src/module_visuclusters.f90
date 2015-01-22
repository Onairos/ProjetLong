module module_visuclusters
  use module_visuclusters_structure
  use module_visuclusters_gmsh
  use module_visuclusters_paraview
contains

  !************************
  !lecture des fichiers d'entree
  subroutine lit(params)
    implicit none
    type(type_params) :: params
    integer :: i,j,k,n
    params%geom=0
    params%seuil=0
    print *
    print *,'lecture des infos...'
    read (3,*)
    read(3,*) params%mesh
    print *,'  > nom du fichier de maillage : ',params%mesh
    read (3,*)
    read(3,*) params%nbp
    print *,'  > nb de points :',params%nbp
    read (3,*)
    read (3,*) params%dim
    print *,'  > dimension : ',params%dim
    read (3,*)
    read(3,*) params%nbproc
    print *,'  > nb de proc :',params%nbproc
    read (3,*)
    read(3,*) params%interface
    print *,'  > decoupage par interface ?',params%interface
    read (3,*)
    read(3,*) params%recouvrement
    print *,'  > decoupage par recouvrement ?',params%recouvrement
    read (3,*)
    read(3,*) params%nbclusters  
    print *,'  > nb de clusters obtenus :',params%nbclusters
    read (3,*)
    read(3,*) params%coord
    print *,'  > format coord ?',params%coord
    read (3,*)
    read(3,*) params%image
    print *,'  > format image ?',params%image
    read (3,*)
    read(3,*) params%geom
    print *,'  > format geom ?',params%geom
    read (3,*)
    read(3,*) params%seuil
    print *,'  > format seuil ?',params%seuil
    if ((params%image==1).or.(params%geom==1).or.(params%seuil==1)) then
       read (3,*)
       read(3,*) params%imgdim
       print *,'  > dimension image :',params%imgdim
       allocate(params%imgmap(params%imgdim))
       read (3,*)
       read(3,*) params%imgmap(:)
       print *,'  > decoupage image :',params%imgmap
       read (3,*)
       read(3,*) params%imgt
       print *,'  > nb de temps :',params%imgt
       !reference des points
       allocate(params%refimg(params%nbp,params%imgdim))
       params%refimg(:,:)=0
       n=0
       do i=1,params%imgmap(1)
          do j=1,params%imgmap(2)
             if (params%imgdim==2) then
                n=n+1
                params%refimg(n,1)=i
                params%refimg(n,2)=j
             else
                do k=1,params%imgmap(3)
                   n=n+1
                   params%refimg(n,1)=i
                   params%refimg(n,2)=j
                   params%refimg(n,3)=k
                end do
             end if
          end do
       end do
       if (params%geom==1) then
          allocate(params%pas(params%imgdim))
          read (3,*)
          read (3,*) params%pas(:)
          print *,'    > pas de maillage :',params%pas(:)
       end if
    end if
    return
  end subroutine lit

  !*************************
  !ecriture de la geometrie du decoupage
  subroutine ecrit_decoupage(formato,params)
    implicit none
    type(type_params) :: params
    character*30 :: formato
    print *
    print *,'ecriture de la geometrie du decoupage...'
    select case(formato)
    case('gmsh')
       call ecrit_decoupage_gmsh(params)
    case('paraview')
       call ecrit_decoupage_paraview(params)
    end select
    return
  end subroutine ecrit_decoupage

  !***********************
  !ecriture des affectations de decoupage
  subroutine affectation(formato,params)
    implicit none
    type(type_params) :: params
    character*30 :: formato
    print *
    print *,'ecriture des affectations du decoupage...'
    select case(formato)
    case('gmsh')
       call affectation_gmsh(params)
    case('paraview')
       call affectation_paraview(params)
    end select
    return
  end subroutine affectation

  !***********************
  !ecriture des clusters avant regroupement
  subroutine sous_clusters(formato,params)
    implicit none
    type(type_params) :: params
    character*30 :: formato
    print *
    print *,'lecture des clusters avant regroupement...'
    select case(formato)
    case('gmsh')
       call sous_clusters_gmsh(params)
    case('paraview')
       call sous_clusters_paraview(params)
    end select
    return
  end subroutine sous_clusters

  !***********************
  !ecriture des clusters apres regroupement
  subroutine cluster_final(formato,params)
    implicit none
    type(type_params) :: params
    character*30 :: formato
    print *
    print *,'lecture des clusters apres regroupement...'
    select case(formato)
    case('gmsh')
       call cluster_final_gmsh(params)
    case('paraview')
       call cluster_final_paraview(params)
    end select
    return
  end subroutine cluster_final

  !***********************
  !liste des commandes
  subroutine commandes(formato)
    implicit none
    character*30 :: formato
    print *
    print *,'-------------------------------------'
    print *,'liste de commandes de visualisation :'
    print *
    select case(formato)
    case('gmsh')
       call commandes_gmsh
    case('paraview')
       call commandes_paraview
    end select
    return
  end subroutine commandes


end module module_visuclusters
