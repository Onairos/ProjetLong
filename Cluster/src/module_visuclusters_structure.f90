MODULE module_visuclusters_structure

  TYPE type_params
     CHARACTER (LEN=30) :: mesh
     DOUBLE PRECISION, DIMENSION(:), POINTER :: pas
     INTEGER, DIMENSION(:,:), POINTER :: refimg
     INTEGER, DIMENSION(:), POINTER :: imgmap
     INTEGER :: nbp
     INTEGER :: dim
     INTEGER :: nbproc
     INTEGER :: nbclusters
     INTEGER :: seuil
     INTEGER :: image
     INTEGER :: imgdim
     INTEGER :: imgt
     INTEGER :: coord
     INTEGER :: geom
     INTEGER :: interface
     INTEGER :: recouvrement
  END TYPE type_params

END MODULE module_visuclusters_structure
