MODULE module_visuclusters_structure

  TYPE type_params
     CHARACTER*30 :: mesh
     INTEGER :: nbp,dim,nbproc,nbclusters,seuil,&
          image,imgdim,imgt,coord,geom,interface,recouvrement
     INTEGER,DIMENSION(:),POINTER :: imgmap
     INTEGER,DIMENSION(:,:),POINTER :: refimg
     DOUBLE PRECISION,DIMENSION(:),POINTER :: pas
  END TYPE type_params

END MODULE module_visuclusters_structure
