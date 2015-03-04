PROGRAM visuclusters
  USE module_visuclusters_structure
  USE module_visuclusters

  IMPLICIT NONE
  !###########################################
  ! DECLARATIONS
  !###########################################
  !#### Variables  ####
  TYPE(type_params) :: params
  CHARACTER (LEN=30) :: formato
  REAL :: elapsed(2) ! For receiving user and system time
  REAL :: temps

  !###########################################
  ! INSTRUCTIONS
  !###########################################
  PRINT *
  PRINT *,'-----------------------------------'
  PRINT *,'visualisation des clusters en 2D/3D'
  PRINT *,'-----------------------------------'

  ! Reads infos
  CALL read_params(params)

  ! Choice of output format
  IF (iargc()>0) THEN
     CALL getarg(1,formato)
     PRINT *,' > format de sortie ? [gmsh,paraview]'
     PRINT *,formato
     GOTO 11
  ENDIF

10  PRINT *
  PRINT *,' > format de sortie ? [gmsh,paraview]'
  READ *,formato
  formato=trim(adjustl(formato))

11 IF ((formato/='gmsh').AND.(formato/='paraview')) GOTO 10

  ! Builds a sub-directory visu/ to store paraview files

  IF (formato=='paraview') CALL system('mkdir visu')

  ! Geometry of partitionning
  CALL write_partionning(formato,params)

  ! Output file
  CALL affectation(formato,params)

  ! Outout file of clusters before regrouping
  IF (params%nbproc>1) CALL write_partial_clusters(formato,params)

  ! Outout file of clusters after regrouping
  CALL write_final_clusters(formato,params)

  ! Commands list
  CALL list_commands(formato)

  ! End file
  SELECT CASE(formato)
  CASE('gmsh')
     OPEN(FILE='visuclusters.gmsh',UNIT=100)
  CASE('paraview')
     OPEN(FILE='visuclusters.paraview',UNIT=100)
  END SELECT
  WRITE(100,*) '# temps total :'
  WRITE(100,*) temps
  WRITE(100,*) '# temps user :'
  WRITE(100,*) elapsed(1)
  WRITE(100,*) '# temps systeme :'
  WRITE(100,*) elapsed(2)
  CLOSE(100)

  PRINT *
  PRINT *,'-----------------------------------'
  PRINT *
  STOP

END PROGRAM visuclusters
