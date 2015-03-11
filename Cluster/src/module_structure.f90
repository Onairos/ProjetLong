!>Contains structure types required by the different modules
MODULE module_structure

  !#### Datatype for computing ####
  TYPE type_data
     ! Data parameters
     TYPE(type_points), DIMENSION(:), POINTER :: point
     INTEGER :: nb_points ! Number of points
     INTEGER :: dim ! Dimension of points
     INTEGER :: nb_clusters ! Number of clusters

     ! Input parameters : format + processing
     INTEGER :: coord ! Data format classic points
     INTEGER :: image ! If image mode activated
     INTEGER :: geom ! Image in mode geom ?
     INTEGER :: seuil ! Image in threshold mode ?
     INTEGER :: interface ! Partitioning + interfacing ?
     INTEGER :: recouvrement ! Partitioning + overlapping ?

     ! Image processing parameters
     DOUBLE PRECISION, DIMENSION(:), POINTER :: pas ! Step for geom mode
     INTEGER, DIMENSION(:,:), POINTER :: refimg ! Point reference in pixel coordinates
     INTEGER, DIMENSION(:), POINTER :: imgmap   ! Pixel partitionings of picture
     INTEGER :: imgdim ! Dimension of images
     INTEGER :: imgt ! Number of "times"

  END TYPE type_data


  !#### Points description ####
  TYPE type_points
     DOUBLE PRECISION, DIMENSION(:), POINTER :: coord
     INTEGER :: cluster
  END TYPE type_points


  !#### Sub-clusters description ####
  TYPE type_clusters
     INTEGER, DIMENSION(:), POINTER :: nbelt
     INTEGER :: nb
  END TYPE type_clusters
  !### Kernel parameter
  TYPE type_clustering_param
     ! Clustering method id
      INTEGER :: clustering_method_id
      INTEGER :: kernelfunindex
      DOUBLE PRECISION :: sigma
      DOUBLE PRECISION :: gam
      DOUBLE PRECISION :: delta
     !mean shift
     INTEGER :: bandwidth !bandwidth for mean shift
  END TYPE type_clustering_param

END MODULE module_structure
