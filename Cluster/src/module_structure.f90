module module_structure

  !***********************
  !type de donnees pour le calcul
  type type_data
     !** parametres des donnees
     integer :: nb, &  !nb de points
          dim, &       !dimension des points
          nbclusters   !nb de clusters
     type(type_points),dimension(:),pointer :: point
     !** parametres d'entree : format + traitement
     integer :: coord, &  !format datas points classique
          image, &        !si mode image active
          geom, &         !image en mode geom ?
          seuil, &        !image en mode seuil ?
          interface, &    !decoupage + interface ?
          recouvrement    !decoupage + recouvrement ?
     !** parametres de traitement d'image
     integer :: imgdim, & !dimension des images
          imgt            !nb de "temps" 
     integer,dimension(:),pointer :: imgmap   !decoupages pixel de l'image
     integer,dimension(:,:),pointer :: refimg !reference du point dans les coordonnes pixel
     real,dimension(:),pointer :: pas !pas utilise pour le mode geom
  end type type_data

  !***********************
  !description des points
  type type_points
     integer :: cluster
     double precision,dimension(:),pointer :: coord
  end type type_points

  !***********************
  !description des sous-clusters
  type type_clusters
     integer :: nb
     integer,dimension(:),pointer :: nbelt
  end type type_clusters

end module module_structure
