PROGRAM gmsh2cluster
  IMPLICIT NONE
  INTEGER :: nbnoe
  INTEGER :: i,j,dim
  REAL :: x,y,z
  CHARACTER(30) :: nom
  PRINT *
  PRINT *,'conversion gmsh vers entree cluster'
  PRINT *
  PRINT *,'nom du fichier gmsh ?'
  READ *,nom
  OPEN(FILE=nom,UNIT=1)
  PRINT *,'nom du fichier du sortie ?'
  READ *,nom
  OPEN(FILE=nom,UNIT=2)
  PRINT *,'DIMENSION ?'
  READ *,dim
  READ (1,*)
  READ (1,*)  
  READ (1,*)
  READ (1,*)  
  READ (1,*)  nbnoe
  PRINT *,'nb de noeuds :',nbnoe
  WRITE(2,*) nbnoe,dim
  DO i=1,nbnoe
     READ (1,*) j,x,y,z
     IF (dim==2) THEN
        WRITE(2,*) x,y
     ELSE
        WRITE (2,*) x,y,z
     ENDIF
  ENDDO
  CLOSE(1); CLOSE(2)
  STOP
END PROGRAM gmsh2cluster
