!>Contains methods enabling data and parameter file reading
MODULE module_entree
 USE module_structure
CONTAINS

!>Displays help on used formats and keywords for <em>param.in</em> file
  SUBROUTINE help
    IMPLICIT NONE
    !###########################################
    ! INSTRUCTIONS
    !###########################################
    PRINT *
    PRINT *,'Calling syntax : clusters input_file'
    PRINT *
    PRINT *,'Input file keywords : '
    PRINT *
    PRINT *,'DATA'
    PRINT *,'COORD (if data with coordinates)'
    PRINT *,'IMAGE (if mesh file = image + partitioning by pixel)'
    PRINT *,'GEOM  (if mesh file = image + geom partitioning)'
    PRINT *,'SEUIL (if mesh file = image + threshold partitioning)'
    PRINT *,'  mesh_file'
    PRINT *
    PRINT *,'EPAISSEUR'
    PRINT *,'  thickness_of_the_slice'
    PRINT *
    PRINT *,'NBLIMIT'
    PRINT *,'  max_nb_of_clusters'
    PRINT *
    PRINT *,'NBCLUST'
    PRINT *,'  (facultatif)'
    PRINT *,'  nb_of_clusters_by_subdomain'
    PRINT *
    PRINT *, 'CLUSTMETHID'
    PRINT *,'  (optional)'
    PRINT *,'  clustering_method_id : 1 Spectral, 2 MeanShift, 3 K-KMeans'
    PRINT *
    PRINT *, 'NBFINALCLUST'
    PRINT *,'  (optional)'
    PRINT *,'  nb_final_clust desired cluster number for KKMeans'
    PRINT *
    PRINT *, 'KERNEL'
    PRINT *,'  (optional)'
    PRINT *,'  ker_table 0 polynomial, 1 gaussian'
    PRINT *
    PRINT *,'SIGMA'
    PRINT *,'  (optional)'
    PRINT *,'  imposed_value_of_sigma'
    PRINT *
    PRINT *, 'GAMMA'
    PRINT *,'  (optional)'
    PRINT *,'  ker_table'
    PRINT *
    PRINT *, 'DELTA'
    PRINT *,'  (optional)'
    PRINT *,'  ker_table'
    PRINT *
    PRINT *, 'BANDWIDTH'
    PRINT *,'  (optional)'
    PRINT *,'  bandwidth_for_mean_shift'
    PRINT *
    PRINT *,'DECOUPAGE'
    PRINT *,'INTERFACE (partitioning with interface) '
    PRINT *,'RECOUVREMENT (partitioning with overlapping)'
    PRINT *,'  nb_of_subdomains_by_DIMENSION'
    PRINT *
    PRINT *,'END'
    PRINT *,'  (end of input file)'
    STOP
    RETURN
  END SUBROUTINE help


!>Reads a file in which there is the whole information on the data
!! @param[in] nb_proc the number of processors used
!! @param[in,out] data the entire data for computing
!! @param[out] input_file the name of the text file where input data is written
!! @param[out] coord_max the maxima along each dimension of the data (coordinates)
!! @param[out] coord_min the minima along each dimension of the data (coordinates)
!! @param[out] epsilon the slice thickness
!! @param[out] sigma the affinity parameter
!! @param[out] nb_clusters_max the maximum number of clusters
!! @param[out] list_nb_clusters the imposed numbers of clusters list for testing purpose (useless)
!! @param[out] partitioning the partitionning (number of processors along each dimension)
  SUBROUTINE read_params(data, epsilon, coord_min, coord_max, nb_proc, partitioning, &
       input_file, sigma, nb_clusters_max, list_nb_clusters,clust_param)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: nb_proc

    !=== IN/OUT ===
    TYPE(type_data) :: data
    TYPE(type_clustering_param) :: clust_param

    !====  OUT ====
    CHARACTER (LEN=30) :: input_file
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min
    DOUBLE PRECISION :: epsilon
    DOUBLE PRECISION :: sigma
    INTEGER, DIMENSION(:), POINTER :: partitioning
    INTEGER, DIMENSION(:), POINTER :: list_nb_clusters
    INTEGER :: nb_clusters_max
    

    !#### Variables  ####
    CHARACTER (LEN=30) :: word
    LOGICAL :: partitioning_bool
    INTEGER :: i
    INTEGER :: ierr
    INTEGER :: tot
    LOGICAL :: ok
    DOUBLE PRECISION :: gam
    DOUBLE PRECISION :: delta

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


    nb_clusters_max=4
    partitioning_bool=.FALSE.
    IF (nb_proc>1) THEN
       ALLOCATE(list_nb_clusters(0:nb_proc-1))
    ELSE
       ALLOCATE(list_nb_clusters(1))
    ENDIF
    list_nb_clusters(:)=0
    ! Reading
    ok=.FALSE.
    DO WHILE (.NOT. ok)
       ok=.TRUE.
       READ(1,*) word
       PRINT *, word
       SELECT CASE(word)
       CASE('DATA')
          ok=.FALSE.
          READ(1,*) input_file
          IF (input_file=='IMAGE') THEN
             data%image=1
             READ (1,*) input_file
             PRINT *, '> Input image format + partitioning by pixel'
             PRINT *, '> Reading input data file : ', input_file
             CALL read_picture_data(input_file,data,coord_min,coord_max)
          ELSEIF (input_file=='GEOM') THEN
             data%geom=1
             READ (1,*) input_file
             PRINT *, '> Input image format + geometric partitioning'
             PRINT *, '> Reading input data file : ', input_file
             CALL read_geometric_data(input_file,data,coord_min,coord_max)
          ELSEIF (input_file=='SEUIL') THEN
             data%seuil=1
             READ (1,*) input_file
             PRINT *, '> Input image format + partitioning by threshold'
             PRINT *, '> Reading input data file : ', input_file
             CALL read_threshold_data(input_file,data,coord_min,coord_max)
          ELSEIF (input_file=='COORD') THEN
             data%coord=1
             READ (1,*) input_file
             PRINT *, '> Input image format + partitioning by threshold'
             PRINT *, '> Reading input data file : ', input_file
             CALL read_coordinates_data(input_file,data,coord_min,coord_max)
          ELSE
             PRINT *
             PRINT *, 'Non-recognized data format !!!'
             CALL help
          ENDIF
          IF ((data%image==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
             ! Creation of array pixels/coordinates
             PRINT *, '> Decoding image format...'
             CALL assign_picture_array(data)
          ENDIF
       CASE('EPAISSEUR')
          ok=.FALSE.
          READ(1,*) epsilon
          PRINT *, '> Thickness of the slice :', epsilon
       CASE('NBLIMIT')
          ok=.FALSE.
          READ(1,*) nb_clusters_max
          PRINT *, '> Maximal number of searched clusters : ', nb_clusters_max
       CASE('NBCLUST')
          ok=.FALSE.
          READ(1,*) list_nb_clusters(:)
          PRINT *, '> Test for number of clusters=', list_nb_clusters
       CASE('KERNEL')
          ok=.FALSE.
          READ(1,*) clust_param%kernelfunindex
          PRINT *, '> Kernel fun index =', clust_param%kernelfunindex
       CASE('CLUSTMETHID')
          ok=.FALSE.
          READ(1,*) clust_param%clustering_method_id
          PRINT *, '> clustering method id =', clust_param%clustering_method_id
       CASE('GAMMA')
          ok=.FALSE.
          READ(1,*) clust_param%gam
          PRINT *, '> 	gamma =', clust_param%gam
       CASE('DELTA')
          ok=.FALSE.
          READ(1,*) clust_param%delta
          PRINT *, '> 	delta =', clust_param%delta
       CASE('SIGMA')
          ok=.FALSE.
          READ(1,*) sigma
          PRINT *, '> Imposed value of sigma : ', sigma
          clust_param%sigma=sigma
          IF (data%image==1) THEN
             IF (sigma<1.0) THEN
                PRINT *, 'Too small thickness for image mode !!!!'
                STOP
             ENDIF
          ENDIF
       CASE('BANDWIDTH')
          ok=.FALSE.
          READ(1,*) clust_param%bandwidth
          PRINT *, '> 	delta =', clust_param%bandwidth
       CASE('DECOUPAGE')
          partitioning_bool=.TRUE.
          ok=.FALSE.
          READ (1,*) word
          SELECT CASE(word)
          CASE('INTERFACE')
             data%interface=1
             PRINT *, '> Partitioning by interface activated.'
          CASE('RECOUVREMENT')
             data%recouvrement=1
             PRINT *, '> Partitioning by overlapping activated.'
          CASE DEFAULT
             PRINT *
             PRINT *, 'Bad partitioning format !!!'
             PRINT *
             CALL help
          END SELECT
          PRINT *, 'dim : ', data%dim
          IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
             ALLOCATE(partitioning(data%dim))
          ELSEIF (data%image==1) THEN
             ! Partitioning per pixel
             ALLOCATE(partitioning(data%imgdim))
          ENDIF
          READ(1,*) partitioning(:)
          PRINT *, 'partitioning', partitioning
          IF (nb_proc>1) THEN
             tot=1
             IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
                DO i=1,data%dim
                   tot=tot*partitioning(i)
                ENDDO
             ELSEIF (data%image==1) THEN
                ! Partitioning per pixel
                DO i=1,data%imgdim
                   tot=tot*partitioning(i)
                ENDDO
             ENDIF
             IF (tot/=nb_proc-data%interface) THEN
                PRINT *, 'Invalidated partitioning !'
                PRINT *, 'Number of process must be equal to ', tot+data%interface
                CALL MPI_ABORT(ierr)
                STOP
             ENDIF
          ELSE
             ! 1 proc
             IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
                DO i=1,data%dim
                   partitioning(i)=1
                ENDDO
                tot=1
             ELSEIF (data%image==1) THEN
                ! Partitioning per pixel
                DO i=1,data%imgdim
                   partitioning(i)=1
                ENDDO
                tot=1
             ENDIF
          ENDIF
          PRINT *, '> partitioning :',partitioning
       CASE('END')
          ok=.TRUE.
       CASE DEFAULT
          ok=.FALSE.
          PRINT *, 'Unknown keyword : ', word
       END SELECT
    ENDDO
    ! Partitioning parameter
    IF ((nb_proc>1).AND.(.NOT.partitioning_bool)) THEN
       PRINT *
       PRINT *, 'The keyword <<DECOUPAGE>> has not been found !'
       CALL help 
    ENDIF
    ! 1 proc
    IF (nb_proc==1) THEN
       ! Initialization to 1 by default of all the partitioning parameters
       IF (partitioning_bool) DEALLOCATE(partitioning)
       IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
          ALLOCATE(partitioning(data%dim) )
       ELSEIF (data%image==1) THEN
          ALLOCATE(partitioning(data%imgdim))
       ENDIF
       partitioning(:)=1
       epsilon=1.0
    ENDIF   
    ! Validation of the combinations of input parameters
    tot=data%geom+data%seuil+data%coord+data%image
    IF (tot/=1) THEN
       PRINT *
       PRINT *, 'Problem with data format !'
       CALL help
    ENDIF
    tot=data%interface+data%recouvrement
    IF (tot/=1) THEN
       PRINT *
       PRINT *, 'Problem with data format !'
       CALL help
    ENDIF
    RETURN
  END SUBROUTINE read_params


!>Reads data written in coordinates format
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
    PRINT *, '> Number of points : ', data%nb
    PRINT *, '> Dimension : ', data%dim
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
100 PRINT *, 'Number of points : ',nb
    data%nb=nb
    CLOSE(2)
    PRINT *, '> Min/max coordinates : '
    DO j=1,data%dim
       PRINT *, '> ', j, ' : ', coord_min(j), coord_max(j)
    ENDDO
    RETURN
  END SUBROUTINE read_coordinates_data


!>Reads data written in picture format
!!@details This type of file starts with two <em>integers</em>
!!separated by a blank space: the first one corresponds to the 
!!dimension of the image and the second one the number of attributes
!!(typically, the attributes could be the color channel intensities).
!!The following line is composed of two <em>integers</em> separated by 
!!a blank space: the first one corresponds to the number of pixels in
!!a column and the second one to the number of pixels in a row. 
!!Then, next lines are composed of <em>double</em> numbers separated 
!!by blank spaces. Each of these lines corresponds to the value of
!!the attributes for each pixel. The coordinates of the pixels are 
!!implicit as they are ordered by rows.
!!@note The attributes are stored in the field type_data::coord and so ONLY the color will be taken into account.
!!@note The outputs <em>coordmax</em> and <em>coordmin</em> are computed according to the position of the pixels.
!!@see assign_picture_array(), read_coordinates_data(), read_threshold_data(), read_geometric_data
!! @param[in] input_file the name of the text file where input data is written
!! @param[in,out] data the entire data for computing
!! @param[out] coord_max the maxima along each dimension of the data (coordinates)
!! @param[out] coord_min the minima along each dimension of the data (coordinates)
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
    PRINT *, '> Image dimension : ', data%imgdim
    PRINT *, '> Number of time : ', data%imgt
    ALLOCATE(data%imgmap(data%imgdim))
    READ(2,*) data%imgmap(:)
    PRINT *, '> Spatial partitioning : ', data%imgmap
    data%nb=1
    DO i=1,data%imgdim
       data%nb=data%nb*data%imgmap(i)
    ENDDO
    data%dim=data%imgt
    data%nbclusters=0
    PRINT *, '> Number of points to read : ', data%nb
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
200 PRINT *, '> Number of points read : ', nb       
    data%nb=nb
    CLOSE(2)
    PRINT *, '> Min/max coordinates :'
    DO j=1,data%dim
       PRINT *, '> ', j, ' : ', coord_min(j), coord_max(j)
    ENDDO
    RETURN
  END SUBROUTINE read_picture_data
  

!>Reads data written in geometric format
!!@details This type of file starts with two <em>integers</em>
!!separated by a blank space: the first one corresponds to the 
!!dimension of the image and the second one the number of attributes
!!(typically, the attributes could be the color channel intensities).
!!The following line is composed of two <em>integers</em> separated by 
!!a blank space: the first one corresponds to the number of pixels in
!!a column and the second one to the number of pixels in a row. 
!!Then, next lines are composed of <em>double</em> numbers separated 
!!by blank spaces. Each of these lines corresponds to the value of
!!the attributes for each pixel. The coordinates of the pixels are 
!!implicit as they are ordered by rows.
!!@note The attributes and the coordinates are stored in the field type_data::coord and so the position AND the color will be taken into account.
!!@note The partitioning is made according to the coordinates
!!@see assign_picture_array(), read_coordinates_data(), read_threshold_data(), read_picture_data
!! @param[in] input_file the name of the text file where input data is written
!! @param[in,out] data the entire data for computing
!! @param[out] coord_max the maxima along each dimension of the data (coordinates)
!! @param[out] coord_min the minima along each dimension of the data (coordinates)
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
    DOUBLE PRECISION :: max_step
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    OPEN(FILE=input_file,UNIT=2)
    READ(2,*) data%imgdim,data%imgt
    PRINT *, '> Image dimension : ', data%imgdim
    PRINT *, '> Number of time : ', data%imgt
    ALLOCATE(data%pas(data%imgdim))
    data%pas(:)=0.0
    ALLOCATE(data%imgmap(data%imgdim))
    READ(2,*) data%imgmap(:)
    PRINT *, '> Spatial partitioning : ', data%imgmap
    data%nb=1
    DO i=1,data%imgdim
       data%nb=data%nb*data%imgmap(i)
    ENDDO
    data%dim=data%imgdim+data%imgt
    data%nbclusters=0
    PRINT *,'> Number of points to read : ', data%nb
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
300 PRINT *, '> Number of points read : ', nb
    data%nb=nb
    CLOSE(2)
    PRINT *, '> Min/max coordinates : '
    max_step=1.e-13
    DO j=data%imgdim+1,data%imgdim+data%imgt
       max_step=max(max_step,coord_max(j)-coord_min(j))
       PRINT *, '> ', j, ' : ', coord_min(j), coord_max(j)
    ENDDO
    PRINT *,'> Maximal step : ', max_step
    ! Searching steps by picture dimension
    DO j=1,data%imgdim
       data%pas(j)=max_step/data%imgmap(j)
       PRINT *, '> Step : ', j, data%pas(j)
       coord_min(j)=0.9*data%pas(j)
       coord_max(j)=(data%imgmap(j)+1)*data%pas(j)
    ENDDO
    RETURN
  END SUBROUTINE read_geometric_data



!>Reads data written in threshold format
!!@details This type of file starts with two <em>integers</em>
!!separated by a blank space: the first one corresponds to the 
!!dimension of the image and the second one the number of attributes
!!(typically, the attributes could be the color channel intensities).
!!The following line is composed of two <em>integers</em> separated by 
!!a blank space: the first one corresponds to the number of pixels in
!!a column and the second one to the number of pixels in a row. 
!!Then, next lines are composed of <em>double</em> numbers separated 
!!by blank spaces. Each of these lines corresponds to the value of
!!the attributes for each pixel. The coordinates of the pixels are 
!!implicit as they are ordered by rows.
!!@note The attributes are stored in the field type_data::coord and so ONLY the color will be taken into account
!!@note The outputs <em>coordmax</em> and <em>coordmin</em> are computed according to the color of the pixels.
!!@see assign_picture_array(), read_picture_data(), read_coordinates_data(), read_geometric_data
!! @param[in] input_file the name of the text file where input data is written
!! @param[in,out] data the entire data for computing
!! @param[out] coord_max the maxima along each dimension of the data (coordinates)
!! @param[out] coord_min the minima along each dimension of the data (coordinates)
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
    PRINT *, '> Image dimension : ', data%imgdim
    PRINT *, '> Number of time : ', data%imgt
    ALLOCATE(data%imgmap(data%imgdim))
    READ(2,*) data%imgmap(:)
    PRINT *, '> Spatial dimension : ', data%imgmap
    data%nb=1
    DO i=1,data%imgdim
       data%nb=data%nb*data%imgmap(i)
    ENDDO
    data%dim=data%imgt
    data%nbclusters=0
    PRINT *, '> Number of points to read : ', data%nb
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
400 PRINT *, 'Number of points : ', nb
    data%nb=nb
    CLOSE(2)
    PRINT *, '> Min/max coordinates : '
    DO j=1,data%dim
       PRINT *, '> ', j, ' : ' , coord_min(j), coord_max(j)
    ENDDO
    RETURN
  END SUBROUTINE read_threshold_data



!>Puts the index number of the pixels into an array
!!@details Each point of the data set is mapped
!!to an array corresponding to an image. Thus, instead
!!of one index for accessing one point, two will be used
!!(row and column).
!! @param[in,out] data the entire data for computing
  SUBROUTINE assign_picture_array(data)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !=== IN/OUT ===
    TYPE(type_data) :: data

    !#### Variables  ####
    INTEGER, DIMENSION(:), POINTER :: plane
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    LOGICAL :: ok

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ! Creation of array points/image_coordinates
    ALLOCATE(data%refimg(data%nb,data%imgdim))
    ALLOCATE(plane(data%imgdim))
    plane(:)=1
    DO i=1,data%nb
       DO j=1,data%imgdim
          ! Index in the array points/pixel
          data%refimg(i,j)=plane(j)
          IF (data%geom==1) THEN
             ! Input of coordinates 1:imgdim for the geometric cluster
             data%point(i)%coord(j)=plane(j)*data%pas(j)
          ENDIF
       ENDDO
       ok=.FALSE.
       k=data%imgdim
       DO WHILE(.NOT. ok)
          IF (plane(k)<data%imgmap(k)) THEN
             plane(k)=plane(k)+1
             ok=.TRUE.
          ELSE
             plane(k)=1
             k=k-1
          ENDIF
          IF (k==0) ok=.TRUE.
       ENDDO
    ENDDO
    DEALLOCATE(plane)
    RETURN
  END SUBROUTINE assign_picture_array


END MODULE module_entree
