!>Contains useful data structures
MODULE module_visuclusters_structure

  TYPE type_params
     CHARACTER (LEN=30) :: input_file
     DOUBLE PRECISION, DIMENSION(:), POINTER :: step
     INTEGER, DIMENSION(:,:), POINTER :: image_ref
     INTEGER, DIMENSION(:), POINTER :: partitioning
     INTEGER :: coords
     INTEGER :: dim
     INTEGER :: is_geom
     INTEGER :: is_image
     INTEGER :: image_dim
     INTEGER :: image_times
     INTEGER :: is_interfacing
     INTEGER :: nb_clusters
     INTEGER :: nb_points
     INTEGER :: nb_proc
     INTEGER :: is_overlapping
     INTEGER :: is_threshold
  END TYPE type_params

END MODULE module_visuclusters_structure
