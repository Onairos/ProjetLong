PROGRAM visudecoup
  IMPLICIT NONE
  INTEGER :: n
  CHARACTER*30 :: files,num
  INTEGER :: i,j,nb
  DOUBLE PRECISION :: coord(2),xmin,xmax,ymin,ymax

  PRINT *
  PRINT *,'visualisation du decoupage parallele en 2D'
  PRINT *

  PRINT *,'nb de decoupages +1 ? (=nbproc)'
  READ *,n
  PRINT *
  
  !geometrie du decoupage
  OPEN(FILE='fort.2',UNIT=2)
  OPEN(FILE='decoupe.geo',UNIT=10)
  DO i=1,n-1
     READ(2,*) xmin,ymin,num,xmax,ymax
     WRITE(10,*) 'Point(',4*(i-1)+1,')={',xmin,',',ymin,',0.};'
     WRITE(10,*) 'Point(',4*(i-1)+2,')={',xmax,',',ymin,',0.};'
     WRITE(10,*) 'Point(',4*(i-1)+3,')={',xmax,',',ymax,',0.};'
     WRITE(10,*) 'Point(',4*(i-1)+4,')={',xmin,',',ymax,',0.};'
     WRITE(10,*) 'Line(',4*(i-1)+1,')={',4*(i-1)+1,',',4*(i-1)+2,'};'
     WRITE(10,*) 'Line(',4*(i-1)+2,')={',4*(i-1)+2,',',4*(i-1)+3,'};'
     WRITE(10,*) 'Line(',4*(i-1)+3,')={',4*(i-1)+3,',',4*(i-1)+4,'};'
     WRITE(10,*) 'Line(',4*(i-1)+4,')={',4*(i-1)+4,',',4*(i-1)+1,'};'
  ENDDO
  CLOSE(10)

  !fichier de sortie
  OPEN(FILE='decoupe.visu',UNIT=1)
  WRITE(1,*) 'View "MPI" {'
  !lecture des fichiers
  DO i=0,n-1
     !nom du fichier
     IF (i<10) THEN
        WRITE(num,'(i1)'),i
     ELSEIF (i<100) THEN
        WRITE(num,'(i2)'),i
     ELSEIF (i<1000) THEN
        WRITE(num,'(i3)'),i
     ENDIF
     files='decoupe.'//trim(adjustl(num))
     OPEN(FILE=files,UNIT=10)
     READ(10,*) nb
     PRINT *,'  > ',i,' :',nb
     DO j=1,nb
        READ(10,*) coord(:)
        WRITE(1,*) 'SP(',coord(1),',',coord(2),',',0.,'){',i,'};'
     ENDDO
     CLOSE(10)

  ENDDO
  WRITE(1,*) '};'
  CLOSE(1)
  PRINT *
  PRINT *,'gmsh decoupe.visu'
  PRINT *,'gmsh decoupe.geo'
  PRINT *,'gmsh decoupe.visu decoupe.geo'
  STOP
END PROGRAM visudecoup
