MODULE module_decoupe
  USE module_structure
  USE module_sortie
CONTAINS


  SUBROUTINE partition_data(data, epsilon, nbproc, coord_min, coord_max, partitionning,&
       points_by_domain, assignements, bounds)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION :: epsilon
    INTEGER, DIMENSION(:), POINTER :: partitionning
    INTEGER :: nbproc

    !=== IN/OUT ===
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bounds
    INTEGER, DIMENSION(:,:), POINTER :: assignements
    INTEGER, DIMENSION(:), POINTER :: points_by_domain

    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: domains

    !###########################################
    ! INSTRUCTIONS
    !########################################### 
    ! Bounds definition
    CALL define_bounds(data,coord_min,coord_max,bounds,partitionning,epsilon,nbproc)

    ! Subdomains definition
    CALL define_domains(nbproc,data,domains,bounds,partitionning)

    ! Writing of partionned subdomains
    CALL write_domains(data,nbproc,domains)

    ! Partitionning definition
    IF ((data%interface==1).OR.(nbproc==1)) THEN
       ! Partitionning by interfacing
       CALL partition_with_interface(nbproc,data,points_by_domain,assignements,domains,epsilon)
    ELSE
       ! Partitionning by overlapping
       CALL partition_with_overlappings(nbproc,data,points_by_domain,assignements,domains)
    ENDIF
    DEALLOCATE(domains)

    ! Saving partitionning
    CALL write_partitionning(nbproc,data,points_by_domain,assignements)

    RETURN
  END SUBROUTINE partition_data


  SUBROUTINE define_bounds(data, coord_min, coord_max, bounds, partitionning, epsilon, nbproc)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION :: epsilon
    INTEGER,DIMENSION(:), POINTER :: partitionning
    INTEGER :: nbproc

    !=== IN/OUT ===
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bounds

    !#### Variables  ####
    CHARACTER (LEN=30) :: num,files
    DOUBLE PRECISION :: prod
    DOUBLE PRECISION :: prod1
    DOUBLE PRECISION :: prod2
    DOUBLE PRECISION :: som1
    INTEGER :: i
    INTEGER :: j

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    prod1=1.0
    prod2=0.0
    som1=1.0
    prod=1.0
    ! Maximum volume
    DO i=1,data%dim
       prod=prod*(coord_max(i)-coord_min(i))
    ENDDO
    files='diminterface'
    WRITE(num,*),0
    num=adjustl(num)
    files=trim(files)//'.'//trim(num)
    OPEN(FILE=files,UNIT=20)
    DO i=1,data%dim
       som1=som1*(partitionning(i)-1)
       prod2=prod2+(partitionning(i)-1)*prod/(coord_max(i)-coord_min(i))    
    ENDDO    
    WRITE(20,*)  prod,epsilon*prod2-som1*(epsilon)**data%dim
    CLOSE(20)
 



 IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
       ! Processing : coordinates, coordinates picture or thresholded picture
       ALLOCATE(bounds(data%dim,max(nbproc-data%interface,1),2))
       bounds(:,:,:)=0.0
       DO i=1,data%dim
          coord_min(i)=coord_min(i)-epsilon*1.1
          coord_max(i)=coord_max(i)+epsilon*1.1
          DO j=1,partitionning(i)
             bounds(i,j,1)=coord_min(i)+(j-1)*(coord_max(i)-coord_min(i))/partitionning(i)
             bounds(i,j,2)=coord_min(i)+j*(coord_max(i)-coord_min(i))/partitionning(i)
          ENDDO
          IF (data%recouvrement==1) THEN
            ! Partitionning with interface mode
             DO j=1,partitionning(i)
                bounds(i,j,1)=bounds(i,j,1)-epsilon
                bounds(i,j,2)=bounds(i,j,2)+epsilon
             ENDDO
          ENDIF
          bounds(i,1,1)=coord_min(i)-0.01*abs(coord_min(i))
          bounds(i,partitionning(i),2)=coord_max(i)+0.01*abs(coord_max(i))
       ENDDO
    ELSEIF (data%image==1) THEN
       ! Processing for partionning pixels of picture
       ALLOCATE(bounds(data%imgdim,max(nbproc-1,1),2))
       bounds(:,:,:)=0.0
       IF ((data%imgdim/=2).AND.(data%imgdim/=3)) THEN
#if aff
          PRINT *
          PRINT *, 'DEBUG : Picture format /= 2D, 3D is not supported !!!!'
#endif
          STOP
       ENDIF
       DO i=1,data%imgdim
          coord_min(i)=1.0-epsilon*1.1
          coord_max(i)=data%imgmap(i)+epsilon*1.1
          DO j=1,partitionning(i)
             bounds(i,j,1)=coord_min(i)+(j-1)*(coord_max(i)-coord_min(i))/partitionning(i)
             bounds(i,j,2)=coord_min(i)+j*(coord_max(i)-coord_min(i))/partitionning(i)
          ENDDO
          IF (data%recouvrement==1) THEN
            ! Partitionning with interface mode
             DO j=1,partitionning(i)
                bounds(i,j,1)=bounds(i,j,1)-epsilon
                bounds(i,j,2)=bounds(i,j,2)+epsilon
             ENDDO
          ENDIF
          bounds(i,1,1)=coord_min(i)-0.01*abs(coord_min(i))
          bounds(i,partitionning(i),2)=coord_max(i)+0.01*abs(coord_max(i))
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE define_bounds


  SUBROUTINE define_domains(nbproc, data, domains, bounds, partitionning)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bounds
    INTEGER, DIMENSION(:), POINTER :: partitionning
    INTEGER :: nbproc

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: domains

    !#### Variables  ####
    INTEGER, DIMENSION(:), POINTER :: list
    INTEGER :: k
    INTEGER :: n
    LOGICAL :: ok


    !###########################################
    ! INSTRUCTIONS
    !###########################################
    IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
       ! Processing : coordinates, coordinates picture or thresholded picture
       ALLOCATE(domains(max(1,nbproc-data%interface),data%dim,2))
       domains(:,:,:)=0.0
       ALLOCATE(list(data%dim))
       list(:)=1
       IF (nbproc>1) THEN
          ! >1 proc
          DO n=1,nbproc-data%interface
             DO k=1,data%dim
                domains(n,k,:)=bounds(k,list(k),:)
             ENDDO
             ok=.TRUE.
             DO k=data%dim,1,-1
                IF (ok) THEN
                   list(k)=list(k)+1
                   IF (list(k)>partitionning(k)) THEN
                      list(k)=1
                   ELSE
                      ok=.FALSE.
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ELSE
          ! 1 proc
          DO k=1,data%dim
             domains(1,k,:)=bounds(k,1,:)
          ENDDO
       ENDIF
       DEALLOCATE(list)
    ELSEIF (data%image==1) THEN
       ! Processing for partitionning in pixels of picture
       ALLOCATE(domains(max(1,nbproc-data%interface),data%imgdim,2))
       domains(:,:,:)=0.0
       ALLOCATE(list(data%imgdim))
       list(:)=1
       IF (nbproc>1) THEN
          ! >1 proc
          DO n=1,nbproc-data%interface
             DO k=1,data%imgdim
                domains(n,k,:)=bounds(k,list(k),:)
             ENDDO
             ok=.TRUE.
             DO k=data%imgdim,1,-1
                IF (ok) THEN
                   list(k)=list(k)+1
                   IF (list(k)>partitionning(k)) THEN
                      list(k)=1
                   ELSE
                      ok=.FALSE.
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ELSE
          ! 1 proc
          DO k=1,data%imgdim
             domains(1,k,:)=bounds(k,1,:)
          ENDDO
       ENDIF
       DEALLOCATE(list)
    ENDIF
    RETURN
  END SUBROUTINE define_domains


  SUBROUTINE partition_with_interface(nbproc, data, points_by_domain, assignements, domains, epsilon)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: domains
    DOUBLE PRECISION :: epsilon
    INTEGER :: nbproc
    !====  OUT ====
    INTEGER, DIMENSION(:,:), POINTER :: assignements
    INTEGER, DIMENSION(:), POINTER :: points_by_domain

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: ierr
    INTEGER :: j
    INTEGER :: n
    LOGICAL :: ok

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ALLOCATE(points_by_domain(0:max(1,nbproc-1)))
    points_by_domain(:)=0
    ALLOCATE(assignements(0:max(1,nbproc-1),data%nb))
    assignements(:,:)=0
    DO i=1,data%nb
       ! Search of packages
       n=0
       ok=.FALSE.
       DO WHILE(.NOT. ok)
          n=n+1
          ok=.TRUE.
          IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
             ! Processing : coordinates, coordinates picture or thresholded picture
             DO j=1,data%dim
                IF ((data%point(i)%coord(j)>domains(n,j,2)).OR.&
                     (data%point(i)%coord(j)<domains(n,j,1))) ok=.FALSE.
             ENDDO
          ELSEIF (data%image==1) THEN
             ! Processing for partitionning in pixels of picture
             DO j=1,data%imgdim
                IF ((data%refimg(i,j)>domains(n,j,2)).OR.&
                     (data%refimg(i,j)<domains(n,j,1))) ok=.FALSE.
             ENDDO
          ENDIF
          IF ((n>nbproc-1).AND.(nbproc>1)) THEN
#if aff
             PRINT *, 'DEBUG : there is a bug in the partitioning ! n=', n, '. Number of process : ', nbproc-1
#endif
             IF (data%geom==0) THEN
#if aff
                PRINT *, 'DEBUG : ', data%point(i)%coord(:)
#endif
             ELSE
#if aff
                PRINT *, 'DEBUG :', i, data%refimg(i,:)
#endif
             ENDIF
             CALL MPI_ABORT(ierr)
             STOP
          ENDIF
       ENDDO
       points_by_domain(n)=points_by_domain(n)+1
       assignements(n,points_by_domain(n))=i
       IF (nbproc>1) THEN
          ! Search of interface if > 1 proc
          ok=.FALSE.
          IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
             ! Processing : coordinates, coordinates picture or thresholded picture
             DO j=1,data%dim
                IF ((abs(data%point(i)%coord(j)-domains(n,j,1))<epsilon).OR.&
                     (abs(data%point(i)%coord(j)-domains(n,j,2))<epsilon)) ok=.TRUE.
             ENDDO
          ELSEIF (data%image==1) THEN
             ! Processing for partitionning in pixels of picture
             DO j=1,data%imgdim
                IF ((abs(data%refimg(i,j)-domains(n,j,1))<epsilon).OR.&
                     (abs(data%refimg(i,j)-domains(n,j,2))<epsilon)) ok=.TRUE.
             ENDDO
          ENDIF
          IF (.NOT. ok) THEN
             points_by_domain(0)=points_by_domain(0)+1
             assignements(0,points_by_domain(0))=i
             WRITE(7,*) assignements(0,points_by_domain(0))
          ENDIF
       ENDIF
    ENDDO
    WRITE(7,*) points_by_domain(0)
    RETURN
  END SUBROUTINE partition_with_interface



  SUBROUTINE partition_with_overlappings(nbproc, data, points_by_domain, assignements, domains)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: domains
    INTEGER :: nbproc
    !====  OUT ====
    INTEGER, DIMENSION(:,:), POINTER :: assignements
    INTEGER, DIMENSION(:), POINTER :: points_by_domain

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: n
    LOGICAL :: ok

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ALLOCATE(points_by_domain(0:max(1,nbproc-1)))
    points_by_domain(:)=0
    ALLOCATE(assignements(0:max(1,nbproc-1),data%nb))
    assignements(:,:)=0
    DO i=1,data%nb
       ! Search of packages
       DO n=1,nbproc
          ok=.TRUE.
          IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
             ! Processing : coordinates, coordinates picture or thresholded picture
             DO j=1,data%dim
                IF ((data%point(i)%coord(j)>domains(n,j,2)).OR.&
                     (data%point(i)%coord(j)<domains(n,j,1))) ok=.FALSE.
             ENDDO
          ELSEIF (data%image==1) THEN
             ! Processing for partitionning in pixels of picture
             DO j=1,data%imgdim
                IF ((data%refimg(i,j)>domains(n,j,2)).OR.&
                     (data%refimg(i,j)<domains(n,j,1))) ok=.FALSE.
             ENDDO
          ENDIF
          IF (ok) THEN
             points_by_domain(n-1)=points_by_domain(n-1)+1
             assignements(n-1,points_by_domain(n-1))=i
          ENDIF
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE partition_with_overlappings


  SUBROUTINE group_clusters(nbclust, points_by_cluster, cluster_map, data)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: nbclust

    !=== IN/OUT ===
    INTEGER, DIMENSION(:,:), POINTER :: cluster_map
    INTEGER, DIMENSION(:), POINTER :: points_by_cluster

    !====  OUT ====
    TYPE(type_data) :: data

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: i2
    INTEGER :: j
    INTEGER :: j2
    INTEGER :: j3
    INTEGER :: k
    INTEGER :: n
    LOGICAL :: ok
    LOGICAL :: ok2
    LOGICAL :: ok3

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ok=.FALSE.
    i=1
    j=0
#if aff
    PRINT *, 'DEBUG : removing duplications...'
    PRINT *, 'DEBUG : grouping subcluster ', 1
#endif
    DO WHILE(.NOT.ok)
       j=j+1 
       IF (j>points_by_cluster(i)) THEN
          ! Line n°1 is entirely tested
#if aff
          PRINT *, 'DEBUG : number of elements after grouping :', points_by_cluster(i)
#endif
          i=i+1
          j=1
#if aff
          PRINT *, 'DEBUG : grouping cluster n', i
#endif
       ENDIF
       IF (i>nbclust-1) THEN
          ! No more points to test
          ok=.TRUE.
       ELSEIF (points_by_cluster(i)>0) THEN
          ! Storage of index
          data%point(cluster_map(i,j))%cluster=i
          ! Test of overlappings
          ok2=.FALSE.
          i2=i+1
          j2=1
          DO WHILE(.NOT. ok2)
             IF (j2>points_by_cluster(i2)) THEN
                ! Line n°i2 entirely tested for the point (i,j)
                i2=i2+1
                j2=1
             ENDIF
             IF (i2>nbclust) THEN
               ! End of test for the point (i,j)
                ok2=.TRUE.
             ELSE
                ! Intersections test
                IF (cluster_map(i,j)==cluster_map(i2,j2)) THEN
                   ! Intersection found : line n°i2 added to line n°i
                   n=0
                   DO k=1,points_by_cluster(i2)
                      ! Test of removal of duplications
                      ok3=.TRUE.
                      DO j3=1,points_by_cluster(i)
                         IF (cluster_map(i2,k)==cluster_map(i,j3)) ok3=.FALSE.
                      ENDDO
                      IF (ok3) THEN
                         n=n+1
                         cluster_map(i,points_by_cluster(i)+n)=cluster_map(i2,k)
                         cluster_map(i2,k)=0                         
                      ENDIF
                   ENDDO
                   points_by_cluster(i)=points_by_cluster(i)+n
                   points_by_cluster(i2)=0
                ELSE
                   ! Test of a new point
                   j2=j2+1
                ENDIF
             ENDIF
          ENDDO
       ENDIF
    ENDDO
#if aff
    PRINT *, 'DEBUG : number of elements after grouping : ', points_by_cluster(i)
#endif
    RETURN
  END SUBROUTINE group_clusters

END MODULE module_decoupe
