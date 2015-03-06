MODULE module_entree
 USE module_structure
CONTAINS

  SUBROUTINE help
    IMPLICIT NONE
    !###########################################
    ! INSTRUCTIONS
    !###########################################
    PRINT *
    PRINT *,'syntaxe d appel : clusters fichier_d_entree'
    PRINT *
    PRINT *,'mots cles du fichier d entree :'
    PRINT *
    PRINT *,'DATA'
    PRINT *,'COORD (si donnees par coordonnees)'
    PRINT *,'IMAGE (si fichier de maillage sous forme image + decoupage par pixel)'
    PRINT *,'GEOM  (si fichier de maillage sous forme image + decoupage geom)'
    PRINT *,'SEUIL  (si fichier de maillage sous forme image + decoupage par seuil)'
    PRINT *,'  fichier_de_maillage'
    PRINT *
    PRINT *,'EPAISSEUR'
    PRINT *,'  epaisseur_de_la_tranche'
    PRINT *
    PRINT *,'NBLIMIT'
    PRINT *,'  nb_max_de_clusters'
    PRINT *
    PRINT *,'NBCLUST'
    PRINT *,'  (facultatif)'
    PRINT *,'  nb_de_clusters_par_sous-domaine'
    PRINT *
    PRINT *,'SIGMA'
    PRINT *,'  (facultatif)'
    PRINT *,'  valeur_de_sigma_imposee'
    PRINT *
    PRINT *,'DECOUPAGE'
    PRINT *,'INTERFACE (decoupage + interface) '
    PRINT *,'RECOUVREMENT (decoupage avec recouvrement)'
    PRINT *,'  nb_de_sous-domaines_par_DIMENSION'
    PRINT *
    PRINT *,'END'
    PRINT *,'  (fin du fichier d entree)'
    STOP
    RETURN
  END SUBROUTINE help


  SUBROUTINE read_file(data, epsilon, coord_min, coord_max, nbproc, partitionning, &
       input_file, sigma, nblimit, list_nb_clusters)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: nbproc

    !=== IN/OUT ===
    TYPE(type_data) :: data

    !====  OUT ====
    CHARACTER (LEN=30) :: input_file
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min
    DOUBLE PRECISION :: epsilon
    DOUBLE PRECISION :: sigma
    INTEGER, DIMENSION(:), POINTER :: partitionning
    INTEGER, DIMENSION(:), POINTER :: list_nb_clusters
    INTEGER :: nblimit

    !#### Variables  ####
    CHARACTER (LEN=30) :: mot
    INTEGER :: decoupage
    INTEGER :: i
    INTEGER :: ierr
    INTEGER :: tot
    LOGICAL :: ok

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    epsilon=0.0
    sigma=-1.0
    data%coord=0
    data%image=0
    data%geom=0
    data%seuil=0
    data%interface=0
    data%recouvrement=0
    nblimit=4
    decoupage=0
    IF (nbproc>1) THEN
       ALLOCATE(list_nb_clusters(0:nbproc-1))
    ELSE
       ALLOCATE(list_nb_clusters(1))
    ENDIF
    list_nb_clusters(:)=0
    ! Reading
    ok=.FALSE.
    DO WHILE (.NOT. ok)
       ok=.TRUE.
       READ(1,*) mot
       PRINT *,mot
       SELECT CASE(mot)
       CASE('DATA')
          ok=.FALSE.
          READ(1,*) input_file
          IF (input_file=='IMAGE') THEN
             data%image=1
             READ (1,*) input_file
             PRINT *,'  > format d entree image + decoupage par pixel'
             PRINT *,'  > lecture du fichier de data : ',input_file
             CALL read_picture_data(input_file,data,coord_min,coord_max)
          ELSEIF (input_file=='GEOM') THEN
             data%geom=1
             READ (1,*) input_file
             PRINT *,'  > format d entree image + decoupage geometrique'
             PRINT *,'  > lecture du fichier de data : ',input_file
             CALL read_geometric_data(input_file,data,coord_min,coord_max)
          ELSEIF (input_file=='SEUIL') THEN
             data%seuil=1
             READ (1,*) input_file
             PRINT *,'  > format d entree image + decoupage par seuil'
             PRINT *,'  > lecture du fichier de data : ',input_file
             CALL read_threshold_data(input_file,data,coord_min,coord_max)
          ELSEIF (input_file=='COORD') THEN
             data%coord=1
             READ (1,*) input_file
             PRINT *,'  > format d entree image + decoupage par seuil'
             PRINT *,'  > lecture du fichier de data : ',input_file
             CALL read_coordinates_data(input_file,data,coord_min,coord_max)
          ELSE
             PRINT *
             PRINT *,'format de donnees non reconnu !!!'
             CALL help
          ENDIF
          IF ((data%image==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
             ! Creation of array pixels/coordinates
             PRINT *,'  > decodage du format image...'
             CALL tableau_image(data)
          ENDIF
       CASE('EPAISSEUR')
          ok=.FALSE.
          READ(1,*) epsilon
          PRINT *,'  > epaisseur de la tranche :',epsilon
       CASE('NBLIMIT')
          ok=.FALSE.
          READ(1,*) nblimit
          PRINT *,'  > nb maximal de clusters recherches :',nblimit
       CASE('NBCLUST')
          ok=.FALSE.
          READ(1,*) list_nb_clusters(:)
          PRINT *,'  > test pour nb de clusters=',list_nb_clusters
       CASE('SIGMA')
          ok=.FALSE.
          READ(1,*) sigma
          PRINT *,'  > valeur de sigma imposee :',sigma
          IF (data%image==1) THEN
             IF (sigma<1.0) THEN
                PRINT *,'epaisseur trop petite pour le mode image !!!!'
                STOP
             ENDIF
          ENDIF
       CASE('DECOUPAGE')
          decoupage=1
          ok=.FALSE.
          READ (1,*) mot
          SELECT CASE(mot)
          CASE('INTERFACE')
             data%interface=1
             PRINT *,'  > decoupage par interface active.'
          CASE('RECOUVREMENT')
             data%recouvrement=1
             PRINT *,'  > decoupage par recouvrement active.'
          CASE DEFAULT
             PRINT *
             PRINT *,'mauvais format de decoupage !!!'
             PRINT *
             CALL help
          END SELECT
          PRINT *, 'dim', data%dim
          IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
             ALLOCATE(partitionning(data%dim))
          ELSEIF (data%image==1) THEN
             ! Partitionning per pixel
             ALLOCATE(partitionning(data%imgdim))
          ENDIF
          READ(1,*) partitionning(:)
          PRINT *, 'partitionning', partitionning
          IF (nbproc>1) THEN
             tot=1
             IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
                DO i=1,data%dim
                   tot=tot*partitionning(i)
                ENDDO
             ELSEIF (data%image==1) THEN
                ! Partitionning per pixel
                DO i=1,data%imgdim
                   tot=tot*partitionning(i)
                ENDDO
             ENDIF
             IF (tot/=nbproc-data%interface) THEN
                PRINT *,'decoupage non valide !'
                PRINT *,'le nombre de proc doit etre egal a',tot+data%interface
                CALL MPI_ABORT(ierr)
                STOP
             ENDIF
          ELSE
             ! 1 proc
             IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
                DO i=1,data%dim
                   partitionning(i)=1
                ENDDO
                tot=1
             ELSEIF (data%image==1) THEN
                ! Partitionning per pixel
                DO i=1,data%imgdim
                   partitionning(i)=1
                ENDDO
                tot=1
             ENDIF
          ENDIF
          PRINT *,'  > decoupage :',partitionning
       CASE('END')
          ok=.TRUE.
       CASE DEFAULT
          ok=.FALSE.
          PRINT *,'mot cle inconnu :',mot
       END SELECT
    ENDDO
    ! Partitionning parameter
    IF ((nbproc>1).AND.(decoupage==0)) THEN
       PRINT *
       PRINT *,'mot cle DECOUPAGE absent !'
       CALL help 
    ENDIF
    ! 1 proc
    IF (nbproc==1) THEN
       ! Initialization to 1 by default of all the partitionning parameters
       IF (decoupage==1) DEALLOCATE(partitionning)
       IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
          ALLOCATE(partitionning(data%dim) )
       ELSEIF (data%image==1) THEN
          ALLOCATE(partitionning(data%imgdim))
       ENDIF
       partitionning(:)=1
       epsilon=1.0
    ENDIF   
    ! Validation of the combinations of input parameters
    tot=data%geom+data%seuil+data%coord+data%image
    IF (tot/=1) THEN
       PRINT *
       PRINT *,'pb dans les formats de data entres !'
       CALL help
    ENDIF
    tot=data%interface+data%recouvrement
    IF (tot/=1) THEN
       PRINT *
       PRINT *,'pb dans les formats de decoupage entres !'
       CALL help
    ENDIF
    RETURN
  END SUBROUTINE read_file


  SUBROUTINE read_coordinates_data(input_file, data, coord_min, coord_max)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    CHARACTER (LEN=30) :: input_file

    !=== IN/OUT ===
    TYPE(type_data) :: data

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ! Reading classic data
    OPEN(FILE=input_file,UNIT=2)
    READ(2,*) data%nb,data%dim
    data%nbclusters=0
    PRINT *,'    > nb de points :',data%nb
    PRINT *,'    > dimension :',data%dim
    ALLOCATE(data%point(data%nb))
    ALLOCATE(coord_max(data%dim))
    ALLOCATE(coord_min(data%dim))
    nb=0
    DO i=1,data%nb
       ALLOCATE(data%point(i)%coord(data%dim))
       READ(2,*,END=100) data%point(i)%coord(:)
       nb=nb+1
       data%point(i)%cluster=-1
       IF (i==1) THEN
          coord_max(:)=data%point(1)%coord(:)
          coord_min(:)=data%point(1)%coord(:)
       ELSE
          DO j=1,data%dim
             coord_min(j)=min(coord_min(j),data%point(i)%coord(j))
             coord_max(j)=max(coord_max(j),data%point(i)%coord(j))
          ENDDO
       ENDIF
    ENDDO
100 PRINT *,'nb de points',nb
    data%nb=nb
    CLOSE(2)
    PRINT *,'  > coordonnees min/max :'
    DO j=1,data%dim
       PRINT *,'    > ',j,':',coord_min(j),coord_max(j)
    ENDDO
    RETURN
  END SUBROUTINE read_coordinates_data


  SUBROUTINE read_picture_data(input_file, data, coord_min, coord_max)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    CHARACTER (LEN=30) :: input_file
    !=== IN/OUT ===
    TYPE(type_data) :: data
    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    OPEN(FILE=input_file,UNIT=2)
    READ(2,*) data%imgdim,data%imgt
    PRINT *,'    >dimension de l image:',data%imgdim
    PRINT *,'    >nb de temps:',data%imgt
    ALLOCATE(data%imgmap(data%imgdim))
    READ(2,*) data%imgmap(:)
    PRINT *,'    >decoupage spatial:',data%imgmap
    data%nb=1
    DO i=1,data%imgdim
       data%nb=data%nb*data%imgmap(i)
    ENDDO
    data%dim=data%imgt
    data%nbclusters=0
    PRINT *,'    > nb de points a lire :',data%nb
    ALLOCATE(data%point(data%nb))
    ALLOCATE(coord_max(data%imgdim))
    ALLOCATE(coord_min(data%imgdim))
    coord_min(:)=0.9
    DO i=1,data%imgdim
       coord_max(i)=data%imgmap(i)+0.1
    ENDDO
    nb=0
    DO i=1,data%nb
       ALLOCATE(data%point(i)%coord(data%dim))
       READ(2,*,END=200) data%point(i)%coord(:)
       nb=nb+1
       data%point(i)%cluster=-1
    ENDDO
200 PRINT *,'    > nb de points lus',nb       
    data%nb=nb
    CLOSE(2)
    PRINT *,'  > coordonnees min/max :'
    DO j=1,data%dim
       PRINT *,'    > ',j,':',coord_min(j),coord_max(j)
    ENDDO
    RETURN
  END SUBROUTINE read_picture_data
  

  SUBROUTINE read_geometric_data(input_file, data, coord_min, coord_max)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    CHARACTER (LEN=30) :: input_file

    !=== IN/OUT ===
    TYPE(type_data) :: data

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min

    !#### Variables  ####
    DOUBLE PRECISION :: pasmax
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    OPEN(FILE=input_file,UNIT=2)
    READ(2,*) data%imgdim,data%imgt
    PRINT *,'    >dimension de l image:',data%imgdim
    PRINT *,'    >nb de temps:',data%imgt
    ALLOCATE(data%pas(data%imgdim))
    data%pas(:)=0.0
    ALLOCATE(data%imgmap(data%imgdim))
    READ(2,*) data%imgmap(:)
    PRINT *,'    >decoupage spatial:',data%imgmap
    data%nb=1
    DO i=1,data%imgdim
       data%nb=data%nb*data%imgmap(i)
    ENDDO
    data%dim=data%imgdim+data%imgt
    data%nbclusters=0
    PRINT *,'    > nb de points a lire :',data%nb
    ALLOCATE(data%point(data%nb))
    ALLOCATE(coord_max(data%dim))
    ALLOCATE(coord_min(data%dim))
    nb=0
    DO i=1,data%nb
       ALLOCATE(data%point(i)%coord(data%dim))
       data%point(i)%coord(:)=0.0
       READ(2,*,END=300) data%point(i)%coord(data%imgdim+1:data%imgdim+data%imgt)
       nb=nb+1
       data%point(i)%cluster=-1
       IF (i==1) THEN
          coord_max(:)=data%point(1)%coord(:)
          coord_min(:)=data%point(1)%coord(:)
       ELSE
          DO j=1,data%dim
             coord_min(j)=min(coord_min(j),data%point(i)%coord(j))
             coord_max(j)=max(coord_max(j),data%point(i)%coord(j))
          ENDDO
       ENDIF
    ENDDO
300 PRINT *,'    > nb de points lus',nb
    data%nb=nb
    CLOSE(2)
    PRINT *,'  > coordonnees min/max :'
    pasmax=1.e-13
    DO j=data%imgdim+1,data%imgdim+data%imgt
       pasmax=max(pasmax,coord_max(j)-coord_min(j))
       PRINT *,'    > ',j,':',coord_min(j),coord_max(j)
    ENDDO
    PRINT *,'  > pas max :',pasmax
    ! Searching steps by picture dimension
    DO j=1,data%imgdim
       data%pas(j)=pasmax/data%imgmap(j)
       PRINT *,'  > pas :',j,data%pas(j)
       coord_min(j)=0.9*data%pas(j)
       coord_max(j)=(data%imgmap(j)+1)*data%pas(j)
    ENDDO
    RETURN
  END SUBROUTINE read_geometric_data



  SUBROUTINE read_threshold_data(input_file, data, coord_min, coord_max)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    CHARACTER (LEN=30) :: input_file

    !=== IN/OUT ===
    TYPE(type_data) :: data

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ! Reading classic data
    OPEN(FILE=input_file,UNIT=2)
    READ(2,*) data%imgdim,data%imgt
    PRINT *,'    >dimension de l image:',data%imgdim
    PRINT *,'    >nb de temps:',data%imgt
    ALLOCATE(data%imgmap(data%imgdim))
    READ(2,*) data%imgmap(:)
    PRINT *,'    >decoupage spatial:',data%imgmap
    data%nb=1
    DO i=1,data%imgdim
       data%nb=data%nb*data%imgmap(i)
    ENDDO
    data%dim=data%imgt
    data%nbclusters=0
    PRINT *,'    > nb de points a lire :',data%nb
    ALLOCATE(data%point(data%nb))
    ALLOCATE(coord_max(data%dim))
    ALLOCATE(coord_min(data%dim))
    nb=0
    DO i=1,data%nb
       ALLOCATE(data%point(i)%coord(data%dim))
       READ(2,*,END=400) data%point(i)%coord(:)
       nb=nb+1
       data%point(i)%cluster=-1
       IF (i==1) THEN
          coord_max(:)=data%point(1)%coord(:)
          coord_min(:)=data%point(1)%coord(:)
       ELSE
          DO j=1,data%dim
             coord_min(j)=min(coord_min(j),data%point(i)%coord(j))
             coord_max(j)=max(coord_max(j),data%point(i)%coord(j))
          ENDDO
       ENDIF
    ENDDO
400 PRINT *,'nb de points',nb
    data%nb=nb
    CLOSE(2)
    PRINT *,'  > coordonnees min/max :'
    DO j=1,data%dim
       PRINT *,'    > ',j,':',coord_min(j),coord_max(j)
    ENDDO
    RETURN
  END SUBROUTINE read_threshold_data



  SUBROUTINE tableau_image(data)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !=== IN/OUT ===
    TYPE(type_data) :: data

    !#### Variables  ####
    INTEGER, DIMENSION(:), POINTER :: plan
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    LOGICAL :: ok

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ! Creation of array points/image_coordinates
    ALLOCATE(data%refimg(data%nb,data%imgdim))
    ALLOCATE(plan(data%imgdim))
    plan(:)=1
    DO i=1,data%nb
       DO j=1,data%imgdim
          ! Index in the array points/pixel
          data%refimg(i,j)=plan(j)
          IF (data%geom==1) THEN
             ! Input of coordinates 1:imgdim for the geometric cluster
             data%point(i)%coord(j)=plan(j)*data%pas(j)
          ENDIF
       ENDDO
       ok=.FALSE.
       k=data%imgdim
       DO WHILE(.NOT. ok)
          IF (plan(k)<data%imgmap(k)) THEN
             plan(k)=plan(k)+1
             ok=.TRUE.
          ELSE
             plan(k)=1
             k=k-1
          ENDIF
          IF (k==0) ok=.TRUE.
       ENDDO
    ENDDO
    DEALLOCATE(plan)
    RETURN
  END SUBROUTINE tableau_image


END MODULE module_entree
