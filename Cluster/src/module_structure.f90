MODULE module_structure

  !***********************
  !type de donnees pour le calcul
  TYPE type_data
     !** parametres des donnees
     INTEGER :: nb, &  !nb de points
          dim, &       !dimension des points
          nbclusters   !nb de clusters
     TYPE(type_points), DIMENSION(:), POINTER :: point
     !** parametres d'entree : format + traitement
     INTEGER :: coord, &  !format datas points classique
          image, &        !si mode image active
          geom, &         !image en mode geom ?
          seuil, &        !image en mode seuil ?
          interface, &    !decoupage + interface ?
          recouvrement    !decoupage + recouvrement ?
     !** parametres de traitement d'image
     INTEGER :: imgdim, & !dimension des images
          imgt            !nb de "temps" 
     INTEGER, DIMENSION(:), POINTER :: imgmap   !decoupages pixel de l'image
     INTEGER, DIMENSION(:,:), POINTER :: refimg !reference du point dans les coordonnes pixel
     DOUBLE PRECISION, DIMENSION(:), POINTER :: pas !pas utilise pour le mode geom
  END TYPE type_data

  !***********************
  !description des points
  TYPE type_points
     INTEGER :: cluster
     DOUBLE PRECISION, DIMENSION(:), POINTER :: coord
  END TYPE type_points

  !***********************
  !description des sous-clusters
  TYPE type_clusters
     INTEGER :: nb
     INTEGER, DIMENSION(:), POINTER :: nbelt
  END TYPE type_clusters

END MODULE module_structure
