program visudecoup
  implicit none
  integer :: n
  character*30 :: files,num
  integer :: i,j,nb
  real :: coord(2),xmin,xmax,ymin,ymax

  print *
  print *,'visualisation du decoupage parallele en 2D'
  print *

  print *,'nb de decoupages +1 ? (=nbproc)'
  read *,n
  print *
  
  !geometrie du decoupage
  open(file='fort.2',unit=2)
  open(file='decoupe.geo',unit=10)
  do i=1,n-1
     read(2,*) xmin,ymin,num,xmax,ymax
     write(10,*) 'Point(',4*(i-1)+1,')={',xmin,',',ymin,',0.};'
     write(10,*) 'Point(',4*(i-1)+2,')={',xmax,',',ymin,',0.};'
     write(10,*) 'Point(',4*(i-1)+3,')={',xmax,',',ymax,',0.};'
     write(10,*) 'Point(',4*(i-1)+4,')={',xmin,',',ymax,',0.};'
     write(10,*) 'Line(',4*(i-1)+1,')={',4*(i-1)+1,',',4*(i-1)+2,'};'
     write(10,*) 'Line(',4*(i-1)+2,')={',4*(i-1)+2,',',4*(i-1)+3,'};'
     write(10,*) 'Line(',4*(i-1)+3,')={',4*(i-1)+3,',',4*(i-1)+4,'};'
     write(10,*) 'Line(',4*(i-1)+4,')={',4*(i-1)+4,',',4*(i-1)+1,'};'
  enddo
  close(10)

  !fichier de sortie
  open(file='decoupe.visu',unit=1)
  write(1,*) 'View "MPI" {'
  !lecture des fichiers
  do i=0,n-1
     !nom du fichier
     if (i<10) then
        write(num,'(i1)'),i
     elseif (i<100) then
        write(num,'(i2)'),i
     elseif (i<1000) then
        write(num,'(i3)'),i
     endif
     files='decoupe.'//trim(adjustl(num))
     open(file=files,unit=10)
     read(10,*) nb
     print *,'  > ',i,' :',nb
     do j=1,nb
        read(10,*) coord(:)
        write(1,*) 'SP(',coord(1),',',coord(2),',',0.,'){',i,'};'
     enddo
     close(10)

  end do
  write(1,*) '};'
  close(1)
  print *
  print *,'gmsh decoupe.visu'
  print *,'gmsh decoupe.geo'
  print *,'gmsh decoupe.visu decoupe.geo'
  stop
end program visudecoup
