program gmsh2cluster
  implicit none
  integer :: nbnoe
  integer :: i,j,dim
  real :: x,y,z
  character(30) :: nom
  print *
  print *,'conversion gmsh vers entree cluster'
  print *
  print *,'nom du fichier gmsh ?'
  read *,nom
  open(file=nom,unit=1)
  print *,'nom du fichier du sortie ?'
  read *,nom
  open(file=nom,unit=2)
  print *,'dimension ?'
  read *,dim
  read (1,*)
  read (1,*)  
  read (1,*)
  read (1,*)  
  read (1,*)  nbnoe
  print *,'nb de noeuds :',nbnoe
  write(2,*) nbnoe,dim
  do i=1,nbnoe
     read (1,*) j,x,y,z
     if (dim==2) then
        write(2,*) x,y
     else
        write (2,*) x,y,z
     end if
  end do
  close(1); close(2)
  stop
end program gmsh2cluster
