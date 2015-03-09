MODULE module_structure

  !#### Datatype for computing ####
  TYPE type_data
     ! Data parameters
     TYPE(type_points), DIMENSION(:), POINTER :: point
     INTEGER :: nb ! Number of points
     INTEGER :: dim ! Dimension of points
     INTEGER :: nbclusters ! Number of clusters

     ! Input parameters : format + processing
     INTEGER :: coord ! Data format classic points
     INTEGER :: image ! If image mode activated
     INTEGER :: geom ! Image in mode geom ?
     INTEGER :: seuil ! Image in threshold mode ?
     INTEGER :: interface ! Partitionning + interfacing ?
     INTEGER :: recouvrement ! Partitionning + overlapping ?

     ! Image processing parameters
     DOUBLE PRECISION, DIMENSION(:), POINTER :: pas ! Step for geom mode
     INTEGER, DIMENSION(:,:), POINTER :: refimg ! Point reference in pixel coordinates
     INTEGER, DIMENSION(:), POINTER :: imgmap   ! Pixel partitionnings of picture
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
  END TYPE type_kernel

END MODULE module_structure
