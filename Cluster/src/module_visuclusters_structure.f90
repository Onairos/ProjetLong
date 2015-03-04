MODULE module_visuclusters_structure

  TYPE type_params
     CHARACTER (LEN=30) :: mesh
     DOUBLE PRECISION, DIMENSION(:), POINTER :: pas
     INTEGER, DIMENSION(:,:), POINTER :: refimg
     INTEGER, DIMENSION(:), POINTER :: imgmap
     INTEGER :: coord
     INTEGER :: dim
     INTEGER :: geom
     INTEGER :: image
     INTEGER :: imgdim
     INTEGER :: imgt
     INTEGER :: interface
     INTEGER :: nbclusters
     INTEGER :: nbp
     INTEGER :: nbproc
     INTEGER :: recouvrement
     INTEGER :: seuil
  END TYPE type_params

END MODULE module_visuclusters_structure
