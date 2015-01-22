program visuclusters
  use module_visuclusters_structure
  use module_visuclusters

  implicit none
  type(type_params) :: params
  character*30 :: formato
  real :: elapsed(2)     ! For receiving user and system time
  real :: temps

  print *
  print *,'-----------------------------------'
  print *,'visualisation des clusters en 2D/3D'
  print *,'-----------------------------------'

  !lecture des infos
  call lit(params)

  !choix du format de sortie
  if (iargc()>0) then
     call getarg(1,formato)
     print *,' > format de sortie ? [gmsh,paraview]'
     print *,formato
     goto 11
  endif

10  print *
  print *,' > format de sortie ? [gmsh,paraview]'
  read *,formato
  formato=trim(adjustl(formato))

11 if ((formato/='gmsh').and.(formato/='paraview')) goto 10

  !creation de la sub-directory visu/ pour stocker les fichiers paraview

  if (formato=='paraview') call system('mkdir visu')

  !geometrie du decoupage
  !if (params%image==0) 
  call ecrit_decoupage(formato,params)

  !fichier de sortie
  call affectation(formato,params)

  !fichier de sortie des clusters avant regroupement
  if (params%nbproc>1) call sous_clusters(formato,params)

  !fichier de sortie des clusters apres regroupement
  call cluster_final(formato,params)

  !liste des commandes
  call commandes(formato)

  !fin du fichier
!  temps = etime(elapsed)
  select case(formato)
  case('gmsh')
     open(file='visuclusters.gmsh',unit=100)
  case('paraview')
     open(file='visuclusters.paraview',unit=100)
  end select
  write(100,*) '# temps total :'
  write(100,*) temps
  write(100,*) '# temps user :'
  write(100,*) elapsed(1)
  write(100,*) '# temps systeme :'
  write(100,*) elapsed(2)
  close(100)

  print *
  print *,'-----------------------------------'
  print *
  stop

end program visuclusters
