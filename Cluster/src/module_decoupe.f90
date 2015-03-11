!>Contains methods enabling the partitionning and the groupping of the data
MODULE module_decoupe
  USE module_structure
  USE module_sortie
CONTAINS


!>Partitions the data into subdomains for a latter processing by the slaves
!!@details The following is performed :
!!<ol>
!!<li> The bounds of each domain are defined (see define_bounds()).</li>
!!<li> The domains are defined using the bounds (see define_domains()).</li>
!!<li> The domains are written in a dedicated file (see write_domains()).</li>
!!<li> The data is partitioned using interface or overlapping (see partition_with_interface() and partition_with_overlapping()).</li>
!!<li> The partitioning is written in dedicated files (see write_partitioning()).</li>
!!</ol>
!! @param[in] data the entire data for computing
!! @param[in] epsilon the slice thickness
!! @param[in] nb_proc the number of processors used
!! @param[in] partitioning the partitionning (number of processors along each dimension)
!! @param[in,out] coord_max the maxima along each dimension of the data (coordinates)
!! @param[in,out] coord_min the minima along each dimension of the data (coordinates)
!! @param[out] bounds the intervals along each dimension representing the bounds of each partition
!! @param[out] assignments the assignement of each point in a partition
!! @param[out] points_by_domain the number of points in each partition
  SUBROUTINE partition_data(data, epsilon, nb_proc, coord_min, coord_max, partitioning,&
       points_by_domain, assignments, bounds)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION :: epsilon
    INTEGER, DIMENSION(:), POINTER :: partitioning
    INTEGER :: nb_proc

    !=== IN/OUT ===
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bounds
    INTEGER, DIMENSION(:,:), POINTER :: assignments
    INTEGER, DIMENSION(:), POINTER :: points_by_domain

    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: domains

    !###########################################
    ! INSTRUCTIONS
    !########################################### 
    ! Bounds definition
    CALL define_bounds(data,coord_min,coord_max,bounds,partitioning,epsilon,nb_proc)

    ! Subdomains definition
    CALL define_domains(nb_proc,data,domains,bounds,partitioning)

    ! Writing of partionned subdomains
    CALL write_domains(data,nb_proc,domains)

    ! Partitioning definition
    IF ((data%interface==1).OR.(nb_proc==1)) THEN
       ! Partitioning by interfacing
       CALL partition_with_interface(nb_proc,data,points_by_domain,assignments,domains,epsilon)
    ELSE
       ! Partitioning by overlapping
       CALL partition_with_overlapping(nb_proc,data,points_by_domain,assignments,domains)
    ENDIF
    DEALLOCATE(domains)

    ! Saving partitioning
    CALL write_partitioning(nb_proc,data,points_by_domain,assignments)

    RETURN
  END SUBROUTINE partition_data


!>Defines the bounds of each subdomain
!!@details The output <em>bounds</em> has to be interpreted as follows:
!!<ol>
!!<li> The first index corresponds to the dimensions </li>
!!<li> The second index corresponds to the partitioning along the dimension </li>
!!<li> The third index corresponds to the extrema of the bounds </li>
!!</ol>
!!@note The bounds are composed of two points.
!!@note In case of overlapping, the bounds overlapped each others.
!! @param[in] data the entire data for computing
!! @param[in] epsilon the slice thickness
!! @param[in] nb_proc the number of processors used
!! @param[in] partitioning the partitionning (number of processors along each dimension)
!! @param[in,out] coord_max the maxima along each dimension of the data (coordinates)
!! @param[in,out] coord_min the minima along each dimension of the data (coordinates)
!! @param[out] bounds the intervals along each dimension representing the bounds of each partition
  SUBROUTINE define_bounds(data, coord_min, coord_max, bounds, partitioning, epsilon, nb_proc)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION :: epsilon
    INTEGER,DIMENSION(:), POINTER :: partitioning
    INTEGER :: nb_proc

    !=== IN/OUT ===
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_max
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coord_min

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bounds

    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
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
       som1=som1*(partitioning(i)-1)
       prod2=prod2+(partitioning(i)-1)*prod/(coord_max(i)-coord_min(i))    
    ENDDO    
    WRITE(20,*)  prod,epsilon*prod2-som1*(epsilon)**data%dim
    CLOSE(20)
 



 IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
       ! Processing : coordinates, coordinates picture or thresholded picture
       ALLOCATE(bounds(data%dim,max(nb_proc-data%interface,1),2))
       bounds(:,:,:)=0.0
       DO i=1,data%dim
          coord_min(i)=coord_min(i)-epsilon*1.1
          coord_max(i)=coord_max(i)+epsilon*1.1
          DO j=1,partitioning(i)
             bounds(i,j,1)=coord_min(i)+(j-1)*(coord_max(i)-coord_min(i))/partitioning(i)
             bounds(i,j,2)=coord_min(i)+j*(coord_max(i)-coord_min(i))/partitioning(i)
          ENDDO
          IF (data%recouvrement==1) THEN
            ! Partitioning with interface mode
             DO j=1,partitioning(i)
                bounds(i,j,1)=bounds(i,j,1)-epsilon
                bounds(i,j,2)=bounds(i,j,2)+epsilon
             ENDDO
          ENDIF
          bounds(i,1,1)=coord_min(i)-0.01*abs(coord_min(i))
          bounds(i,partitioning(i),2)=coord_max(i)+0.01*abs(coord_max(i))
       ENDDO
    ELSEIF (data%image==1) THEN
       ! Processing for partionning pixels of picture
       ALLOCATE(bounds(data%imgdim,max(nb_proc-1,1),2))
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
          DO j=1,partitioning(i)
             bounds(i,j,1)=coord_min(i)+(j-1)*(coord_max(i)-coord_min(i))/partitioning(i)
             bounds(i,j,2)=coord_min(i)+j*(coord_max(i)-coord_min(i))/partitioning(i)
          ENDDO
          IF (data%recouvrement==1) THEN
            ! Partitioning with interface mode
             DO j=1,partitioning(i)
                bounds(i,j,1)=bounds(i,j,1)-epsilon
                bounds(i,j,2)=bounds(i,j,2)+epsilon
             ENDDO
          ENDIF
          bounds(i,1,1)=coord_min(i)-0.01*abs(coord_min(i))
          bounds(i,partitioning(i),2)=coord_max(i)+0.01*abs(coord_max(i))
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE define_bounds


!>Defines the different subdomains
!!@details The output <em>domains</em> has to be interpreted as follows:
!!<ol>
!!<li> The first index corresponds to the domain id</li>
!!<li> The second index corresponds to the dimensions </li>
!!<li> The third index corresponds to the extrema of the bounds </li>
!!</ol>
!!@sa define_bounds()
!! @param[in] data the entire data for computing
!! @param[in] bounds the intervals along each dimension representing the bounds of each partition
!! @param[in] nb_proc the number of processors used
!! @param[in] partitioning the partitionning (number of processors along each dimension)
!! @param[out] domains the domains constructed from the bounds
  SUBROUTINE define_domains(nb_proc, data, domains, bounds, partitioning)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bounds
    INTEGER, DIMENSION(:), POINTER :: partitioning
    INTEGER :: nb_proc

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
       ALLOCATE(domains(max(1,nb_proc-data%interface),data%dim,2))
       domains(:,:,:)=0.0
       ALLOCATE(list(data%dim))
       list(:)=1
       IF (nb_proc>1) THEN
          ! >1 proc
          DO n=1,nb_proc-data%interface
             DO k=1,data%dim
                domains(n,k,:)=bounds(k,list(k),:)
             ENDDO
             ok=.TRUE.
             DO k=data%dim,1,-1
                IF (ok) THEN
                   list(k)=list(k)+1
                   IF (list(k)>partitioning(k)) THEN
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
       ! Processing for partitioning in pixels of picture
       ALLOCATE(domains(max(1,nb_proc-data%interface),data%imgdim,2))
       domains(:,:,:)=0.0
       ALLOCATE(list(data%imgdim))
       list(:)=1
       IF (nb_proc>1) THEN
          ! >1 proc
          DO n=1,nb_proc-data%interface
             DO k=1,data%imgdim
                domains(n,k,:)=bounds(k,list(k),:)
             ENDDO
             ok=.TRUE.
             DO k=data%imgdim,1,-1
                IF (ok) THEN
                   list(k)=list(k)+1
                   IF (list(k)>partitioning(k)) THEN
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


!>Partitions the data using interface
!!@details The method partitions the entire data by defining
!!which point belongs to which domain. It consists of a loop
!!on all the point in data set. Then it "fills" the domains one
!!after another using a nested loop. When a point does not fit the 
!!bounds define in the input <em>domains</em>, the algorithm switch 
!!to the next domain. Finally, an extra domain is defined (the interface) 
!!which corresponds to the area around the bounds with a predefined 
!!slice thickness.
!!@see partition_with_overlapping()
!! @param[in] data the entire data for computing
!! @param[in] domains the domains constructed from the bounds
!! @param[in] epsilon the slice thickness
!! @param[in] nb_proc the number of processors used
!! @param[out] assignments the assignement of each point in a partition
!! @param[out] points_by_domain the number of points in each partition
  SUBROUTINE partition_with_interface(nb_proc, data, points_by_domain, assignments, domains, epsilon)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: domains
    DOUBLE PRECISION :: epsilon
    INTEGER :: nb_proc
    !====  OUT ====
    INTEGER, DIMENSION(:,:), POINTER :: assignments
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
    ALLOCATE(points_by_domain(0:max(1,nb_proc-1)))
    points_by_domain(:)=0
    ALLOCATE(assignments(0:max(1,nb_proc-1),data%nb_points))
    assignments(:,:)=0
    DO i=1,data%nb_points
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
             ! Processing for partitioning in pixels of picture
             DO j=1,data%imgdim
                IF ((data%refimg(i,j)>domains(n,j,2)).OR.&
                     (data%refimg(i,j)<domains(n,j,1))) ok=.FALSE.
             ENDDO
          ENDIF
          IF ((n>nb_proc-1).AND.(nb_proc>1)) THEN
#if aff
             PRINT *, 'DEBUG : there is a bug in the partitioning ! n=', n, '. Number of process : ', nb_proc-1
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
       assignments(n,points_by_domain(n))=i
       IF (nb_proc>1) THEN
          ! Search of interface if > 1 proc
          ok=.FALSE.
          IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
             ! Processing : coordinates, coordinates picture or thresholded picture
             DO j=1,data%dim
                IF ((abs(data%point(i)%coord(j)-domains(n,j,1))<epsilon).OR.&
                     (abs(data%point(i)%coord(j)-domains(n,j,2))<epsilon)) ok=.TRUE.
             ENDDO
          ELSEIF (data%image==1) THEN
             ! Processing for partitioning in pixels of picture
             DO j=1,data%imgdim
                IF ((abs(data%refimg(i,j)-domains(n,j,1))<epsilon).OR.&
                     (abs(data%refimg(i,j)-domains(n,j,2))<epsilon)) ok=.TRUE.
             ENDDO
          ENDIF
          IF (.NOT. ok) THEN
             points_by_domain(0)=points_by_domain(0)+1
             assignments(0,points_by_domain(0))=i
             WRITE(7,*) assignments(0,points_by_domain(0))
          ENDIF
       ENDIF
    ENDDO
    WRITE(7,*) points_by_domain(0)
    RETURN
  END SUBROUTINE partition_with_interface



!>Partitions the data using overlapping
!!@details The method partitions the entire data by defining
!!which point belongs to which domain. It consists of two nested
!!loop. The first one on the points and the second one on the
!!domains (ie the processes). For each point, it checks if it fits
!!the bounds defined in the input <em>domains</em> (and that for each
!!domain) and in that case add to it.
!!@note Some points will be present in different domains.
!!@see partition_with_interface()
!! @param[in] data the entire data for computing
!! @param[in] domains the domains constructed from the bounds
!! @param[in] nb_proc the number of processors used
!! @param[out] assignments the assignement of each point in a partition
!! @param[out] points_by_domain the number of points in each partition
  SUBROUTINE partition_with_overlapping(nb_proc, data, points_by_domain, assignments, domains)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: domains
    INTEGER :: nb_proc
    !====  OUT ====
    INTEGER, DIMENSION(:,:), POINTER :: assignments
    INTEGER, DIMENSION(:), POINTER :: points_by_domain

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: n
    LOGICAL :: ok

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ALLOCATE(points_by_domain(0:max(1,nb_proc-1)))
    points_by_domain(:)=0
    ALLOCATE(assignments(0:max(1,nb_proc-1),data%nb_points))
    assignments(:,:)=0
    DO i=1,data%nb_points
       ! Search of packages
       DO n=1,nb_proc
          ok=.TRUE.
          IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
             ! Processing : coordinates, coordinates picture or thresholded picture
             DO j=1,data%dim
                IF ((data%point(i)%coord(j)>domains(n,j,2)).OR.&
                     (data%point(i)%coord(j)<domains(n,j,1))) ok=.FALSE.
             ENDDO
          ELSEIF (data%image==1) THEN
             ! Processing for partitioning in pixels of picture
             DO j=1,data%imgdim
                IF ((data%refimg(i,j)>domains(n,j,2)).OR.&
                     (data%refimg(i,j)<domains(n,j,1))) ok=.FALSE.
             ENDDO
          ENDIF
          IF (ok) THEN
             points_by_domain(n-1)=points_by_domain(n-1)+1
             assignments(n-1,points_by_domain(n-1))=i
          ENDIF
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE partition_with_overlapping


!>Groups the clusters and removes duplicates from the set of found clusters
!!@details The method operates on all the computed clusters
!!and when an intersection is found between two of them (at least
!!one point in common), they are melted.
!! @param[in] nb_clusters the number of clusters
!! @param[in] nb_clusters the number of clusters
!! @param[in] nb_clusters the number of clusters
!! @param[in,out] cluster_map the cluster indices and the number of points in each cluster
!! @param[in,out] points_by_cluster the number of points in each cluster
!! @param[out] data the entire data for computing
  SUBROUTINE group_clusters(nb_clusters, points_by_cluster, cluster_map, data)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: nb_clusters

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
          ! Line nÂ°1 is entirely tested
#if aff
          PRINT *, 'DEBUG : number of elements after grouping :', points_by_cluster(i)
#endif
          i=i+1
          j=1
#if aff
          PRINT *, 'DEBUG : grouping cluster ', i
#endif
       ENDIF
       IF (i>nb_clusters-1) THEN
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
                ! Line nÂ°i2 entirely tested for the point (i,j)
                i2=i2+1
                j2=1
             ENDIF
             IF (i2>nb_clusters) THEN
               ! End of test for the point (i,j)
                ok2=.TRUE.
             ELSE
                ! Intersections test
                IF (cluster_map(i,j)==cluster_map(i2,j2)) THEN
                   ! Intersection found : line nÂ°i2 added to line nÂ°i
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
