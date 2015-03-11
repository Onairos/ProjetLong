!>Contains structure types required by the different modules
MODULE module_structure

  !#### Datatype for computing ####
  TYPE type_data
     ! Data parameters
     TYPE(type_points), DIMENSION(:), POINTER :: points
     INTEGER :: nb_points ! Number of points
     INTEGER :: dim ! Dimension of points
     INTEGER :: nb_clusters ! Number of clusters

     ! Input parameters : format + processing
     INTEGER :: coords ! Data format classic points
     INTEGER :: is_image ! If image mode activated
     INTEGER :: is_geom ! Image in mode geom ?
     INTEGER :: is_threshold ! Image in threshold mode ?
     INTEGER :: is_interfacing ! Partitioning + interfacing ?
     INTEGER :: is_overlapping ! Partitioning + overlapping ?

     ! Image processing parameters
     DOUBLE PRECISION, DIMENSION(:), POINTER :: step ! Step for geom mode
     INTEGER, DIMENSION(:,:), POINTER :: image_ref ! Point reference in pixel coordinates
     INTEGER, DIMENSION(:), POINTER :: partitioning   ! Pixel partitionings of picture
     INTEGER :: image_dim ! Dimension of images
     INTEGER :: image_times ! Number of "times"

  END TYPE type_data


  !#### Points description ####
  TYPE type_points
     DOUBLE PRECISION, DIMENSION(:), POINTER :: coords
     INTEGER :: cluster
  END TYPE type_points


  !#### Sub-clusters description ####
  TYPE type_clusters
     INTEGER, DIMENSION(:), POINTER :: nb_elements
     INTEGER :: nb
  END TYPE type_clusters
  !### Kernel parameter
  TYPE type_clustering_param
     ! Clustering method id
      INTEGER :: clustering_method_id
      INTEGER :: kernelfunindex
      DOUBLE PRECISION :: sigma
      DOUBLE PRECISION :: gamma
      DOUBLE PRECISION :: delta
     !mean shift
     INTEGER :: bandwidth !bandwidth for mean shift
  END TYPE type_clustering_param

END MODULE module_structure
