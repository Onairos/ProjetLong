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


  SUBROUTINE read_file(data, epsilon, coordmin, coordmax, nbproc, decoupe, &
       mesh, sigma, nblimit, listenbideal)
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
    CHARACTER (LEN=30) :: mesh
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmax
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmin
    DOUBLE PRECISION :: epsilon
    DOUBLE PRECISION :: sigma
    INTEGER, DIMENSION(:), POINTER :: decoupe
    INTEGER, DIMENSION(:), POINTER :: listenbideal
    INTEGER :: nblimit

    !#### Variables  ####
    CHARACTER (LEN=30) :: mot
    INTEGER :: decoupage
    INTEGER :: i
    INTEGER :: ierr
    INTEGER :: ok !TODO utilise comme un booleen, modifier en LOGICAL ??
    INTEGER :: tot

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
       ALLOCATE(listenbideal(0:nbproc-1))
    ELSE
       ALLOCATE(listenbideal(1))
    ENDIF
    listenbideal(:)=0
    !lecture
    ok=0
    DO WHILE (ok/=1)
       ok=1
       READ(1,*) mot
       PRINT *,mot
       SELECT CASE(mot)
       CASE('DATA')
          ok=0
          READ(1,*) mesh
          IF (mesh=='IMAGE') THEN
             data%image=1
             READ (1,*) mesh
             PRINT *,'  > format d entree image + decoupage par pixel'
             PRINT *,'  > lecture du fichier de data : ',mesh
             CALL read_picture_data(mesh,data,coordmin,coordmax)
          ELSEIF (mesh=='GEOM') THEN
             data%geom=1
             READ (1,*) mesh
             PRINT *,'  > format d entree image + decoupage geometrique'
             PRINT *,'  > lecture du fichier de data : ',mesh
             CALL read_geometric_data(mesh,data,coordmin,coordmax)
          ELSEIF (mesh=='SEUIL') THEN
             data%seuil=1
             READ (1,*) mesh
             PRINT *,'  > format d entree image + decoupage par seuil'
             PRINT *,'  > lecture du fichier de data : ',mesh
             CALL read_threshold_data(mesh,data,coordmin,coordmax)
          ELSEIF (mesh=='COORD') THEN
             data%coord=1
             READ (1,*) mesh
             PRINT *,'  > format d entree image + decoupage par seuil'
             PRINT *,'  > lecture du fichier de data : ',mesh
             CALL read_coordinates_data(mesh,data,coordmin,coordmax)
          ELSE
             PRINT *
             PRINT *,'format de donnees non reconnu !!!'
             CALL help
          ENDIF
          IF ((data%image==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
             !creation du tableau de correspondances pixel/coord
             PRINT *,'  > decodage du format image...'
             CALL tableau_image(data)
          ENDIF
       CASE('EPAISSEUR')
          ok=0
          READ(1,*) epsilon
          PRINT *,'  > epaisseur de la tranche :',epsilon
       CASE('NBLIMIT')
          ok=0
          READ(1,*) nblimit
          PRINT *,'  > nb maximal de clusters recherches :',nblimit
       CASE('NBCLUST')
          ok=0
          READ(1,*) listenbideal(:)
          PRINT *,'  > test pour nb de clusters=',listenbideal
       CASE('SIGMA')
          ok=0
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
          ok=0
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
             ALLOCATE(decoupe(data%dim))
          ELSEIF (data%image==1) THEN
             !decoupage par pixel
             ALLOCATE(decoupe(data%imgdim))
          ENDIF
          READ(1,*) decoupe(:)
          PRINT *, 'decoupe', decoupe
          IF (nbproc>1) THEN
             tot=1
             IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
                DO i=1,data%dim
                   tot=tot*decoupe(i)
                ENDDO
             ELSEIF (data%image==1) THEN
                !decoupage par pixel
                DO i=1,data%imgdim
                   tot=tot*decoupe(i)
                ENDDO
             ENDIF
             IF (tot/=nbproc-data%interface) THEN
                PRINT *,'decoupage non valide !'
                PRINT *,'le nombre de proc doit etre egal a',tot+data%interface
                CALL MPI_ABORT(ierr)
                STOP
             ENDIF
          ELSE
             !mode 1 proc
             IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
                DO i=1,data%dim
                   decoupe(i)=1
                ENDDO
                tot=1
             ELSEIF (data%image==1) THEN
                !decoupage par pixel
                DO i=1,data%imgdim
                   decoupe(i)=1
                ENDDO
                tot=1
             ENDIF
          ENDIF
          PRINT *,'  > decoupage :',decoupe
       CASE('END')
          ok=1
       CASE DEFAULT
          ok=0
          PRINT *,'mot cle inconnu :',mot
       END SELECT
    ENDDO
    !parametre de decoupage
    IF ((nbproc>1).AND.(decoupage==0)) THEN
       PRINT *
       PRINT *,'mot cle DECOUPAGE absent !'
       CALL help 
    ENDIF
    !cas monoproc
    IF (nbproc==1) THEN
       !initialisation par defaut a 1 de tous les parametres de decoupage
       IF (decoupage==1) DEALLOCATE(decoupe)
       IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
          ALLOCATE(decoupe(data%dim) )
       ELSEIF (data%image==1) THEN
          ALLOCATE(decoupe(data%imgdim))
       ENDIF
       decoupe(:)=1
       epsilon=1.0
    ENDIF   
    !validation des combinaisons de parametres d'entree
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


  SUBROUTINE read_coordinates_data(mesh, data, coordmin, coordmax)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    CHARACTER (LEN=30) :: mesh

    !=== IN/OUT ===
    TYPE(type_data) :: data

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmax
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmin

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    !lecture de donnees classiques
    OPEN(FILE=mesh,UNIT=2)
    READ(2,*) data%nb,data%dim
    data%nbclusters=0
    PRINT *,'    > nb de points :',data%nb
    PRINT *,'    > dimension :',data%dim
    ALLOCATE(data%point(data%nb))
    ALLOCATE(coordmax(data%dim))
    ALLOCATE(coordmin(data%dim))
    nb=0
    DO i=1,data%nb
       ALLOCATE(data%point(i)%coord(data%dim))
       READ(2,*,END=100) data%point(i)%coord(:)
       nb=nb+1
       data%point(i)%cluster=-1
       IF (i==1) THEN
          coordmax(:)=data%point(1)%coord(:)
          coordmin(:)=data%point(1)%coord(:)
       ELSE
          DO j=1,data%dim
             coordmin(j)=min(coordmin(j),data%point(i)%coord(j))
             coordmax(j)=max(coordmax(j),data%point(i)%coord(j))
          ENDDO
       ENDIF
    ENDDO
100 PRINT *,'nb de points',nb
    data%nb=nb
    CLOSE(2)
    PRINT *,'  > coordonnees min/max :'
    DO j=1,data%dim
       PRINT *,'    > ',j,':',coordmin(j),coordmax(j)
    ENDDO
    RETURN
  END SUBROUTINE read_coordinates_data


  SUBROUTINE read_picture_data(mesh, data, coordmin, coordmax)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    CHARACTER (LEN=30) :: mesh
    !=== IN/OUT ===
    TYPE(type_data) :: data
    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmax
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmin

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    OPEN(FILE=mesh,UNIT=2)
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
    ALLOCATE(coordmax(data%imgdim))
    ALLOCATE(coordmin(data%imgdim))
    coordmin(:)=0.9
    DO i=1,data%imgdim
       coordmax(i)=data%imgmap(i)+0.1
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
       PRINT *,'    > ',j,':',coordmin(j),coordmax(j)
    ENDDO
    RETURN
  END SUBROUTINE read_picture_data
  

  SUBROUTINE read_geometric_data(mesh, data, coordmin, coordmax)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    CHARACTER (LEN=30) :: mesh

    !=== IN/OUT ===
    TYPE(type_data) :: data

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmax
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmin

    !#### Variables  ####
    DOUBLE PRECISION :: pasmax
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    OPEN(FILE=mesh,UNIT=2)
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
    ALLOCATE(coordmax(data%dim))
    ALLOCATE(coordmin(data%dim))
    nb=0
    DO i=1,data%nb
       ALLOCATE(data%point(i)%coord(data%dim))
       data%point(i)%coord(:)=0.0
       READ(2,*,END=300) data%point(i)%coord(data%imgdim+1:data%imgdim+data%imgt)
       nb=nb+1
       data%point(i)%cluster=-1
       IF (i==1) THEN
          coordmax(:)=data%point(1)%coord(:)
          coordmin(:)=data%point(1)%coord(:)
       ELSE
          DO j=1,data%dim
             coordmin(j)=min(coordmin(j),data%point(i)%coord(j))
             coordmax(j)=max(coordmax(j),data%point(i)%coord(j))
          ENDDO
       ENDIF
    ENDDO
300 PRINT *,'    > nb de points lus',nb
    data%nb=nb
    CLOSE(2)
    PRINT *,'  > coordonnees min/max :'
    pasmax=1.e-13
    DO j=data%imgdim+1,data%imgdim+data%imgt
       pasmax=max(pasmax,coordmax(j)-coordmin(j))
       PRINT *,'    > ',j,':',coordmin(j),coordmax(j)
    ENDDO
    PRINT *,'  > pas max :',pasmax
    !recherche des pas par DIMENSION d'image
    DO j=1,data%imgdim
       data%pas(j)=pasmax/data%imgmap(j)
       PRINT *,'  > pas :',j,data%pas(j)
       coordmin(j)=0.9*data%pas(j)
       coordmax(j)=(data%imgmap(j)+1)*data%pas(j)
    ENDDO
    RETURN
  END SUBROUTINE read_geometric_data



  SUBROUTINE read_threshold_data(mesh, data, coordmin, coordmax)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    CHARACTER (LEN=30) :: mesh

    !=== IN/OUT ===
    TYPE(type_data) :: data

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmax
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmin

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    !lecture de donnees classiques
    OPEN(FILE=mesh,UNIT=2)
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
    ALLOCATE(coordmax(data%dim))
    ALLOCATE(coordmin(data%dim))
    nb=0
    DO i=1,data%nb
       ALLOCATE(data%point(i)%coord(data%dim))
       READ(2,*,END=400) data%point(i)%coord(:)
       nb=nb+1
       data%point(i)%cluster=-1
       IF (i==1) THEN
          coordmax(:)=data%point(1)%coord(:)
          coordmin(:)=data%point(1)%coord(:)
       ELSE
          DO j=1,data%dim
             coordmin(j)=min(coordmin(j),data%point(i)%coord(j))
             coordmax(j)=max(coordmax(j),data%point(i)%coord(j))
          ENDDO
       ENDIF
    ENDDO
400 PRINT *,'nb de points',nb
    data%nb=nb
    CLOSE(2)
    PRINT *,'  > coordonnees min/max :'
    DO j=1,data%dim
       PRINT *,'    > ',j,':',coordmin(j),coordmax(j)
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
    INTEGER :: ok !TODO utilise comme un booleen, modifier en LOGICAL ??

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    !creation du tableau de references points/coordonnes_images
    ALLOCATE(data%refimg(data%nb,data%imgdim))
    ALLOCATE(plan(data%imgdim)); plan(:)=1
    DO i=1,data%nb
       DO j=1,data%imgdim
          !index dans le tableau de reference points/pixel
          data%refimg(i,j)=plan(j)
          IF (data%geom==1) THEN
             !entree des coordonnees 1:imgdim pour le cluster geom
             data%point(i)%coord(j)=plan(j)*data%pas(j)
          ENDIF
       ENDDO
       ok=0
       k=data%imgdim
       DO WHILE(ok==0)
          IF (plan(k)<data%imgmap(k)) THEN
             plan(k)=plan(k)+1
             ok=1
          ELSE
             plan(k)=1
             k=k-1
          ENDIF
          IF (k==0) ok=1
       ENDDO
    ENDDO
    DEALLOCATE(plan)
    RETURN
  END SUBROUTINE tableau_image


END MODULE module_entree
