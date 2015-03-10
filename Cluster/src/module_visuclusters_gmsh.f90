MODULE module_visuclusters_gmsh
  USE module_visuclusters_structure
CONTAINS



  SUBROUTINE write_partitioning_gmsh(params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params

    !#### Variables  ####
    CHARACTER (LEN=30) :: num
    DOUBLE PRECISION :: x0
    DOUBLE PRECISION :: x1
    DOUBLE PRECISION :: x_max
    DOUBLE PRECISION :: x_min
    DOUBLE PRECISION :: y0
    DOUBLE PRECISION :: y1
    DOUBLE PRECISION :: y_max
    DOUBLE PRECISION :: y_min
    DOUBLE PRECISION :: z0
    DOUBLE PRECISION :: z1
    DOUBLE PRECISION :: zmax
    DOUBLE PRECISION :: zmin
    INTEGER :: i

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    PRINT *,'-> decoupe.geo'
    OPEN(FILE='fort.2',UNIT=2)
    OPEN(FILE='decoupe.geo',UNIT=10)
    DO i=1,params%nbproc-params%interface
       IF (((params%coord==1).AND.(params%dim==2)) &
          .OR.((params%image==1).AND.(params%imgdim==2)) &
          .OR.((params%seuil==1).AND.(params%imgdim==2)) &
          .OR.((params%geom==1).AND.(params%imgdim+params%imgt==2))) THEN
          !2D
          IF (params%coord==1) THEN
             READ(2,*) x0,y0,num,x1,y1
             x_min=x0
             y_min=y0
             x_max=x1
             y_max=y1
          ELSEIF (params%seuil==1) THEN
             READ(2,*) x0,num,x1
             x_min=1
             y_min=-1
             x_max=params%imgmap(2)
             y_max=-params%imgmap(1)
             zmin=x0
             zmax=x1
          ELSEIF ((params%image==1).OR.(params%geom==1)) THEN
             READ(2,*) x0,y0,num,x1,y1
             x_min=y0
             y_min=-x0
             x_max=y1
             y_max=-x1
             zmin=x0
          ENDIF
          IF (params%seuil==1) THEN
             WRITE(10,*) 'Point(',8*(i-1)+1,')={',x_min,',',y_min,',',zmin,'};'
             WRITE(10,*) 'Point(',8*(i-1)+2,')={',x_max,',',y_min,',',zmin,'};'
             WRITE(10,*) 'Point(',8*(i-1)+3,')={',x_max,',',y_max,',',zmin,'};'
             WRITE(10,*) 'Point(',8*(i-1)+4,')={',x_min,',',y_max,',',zmin,'};'
             WRITE(10,*) 'Point(',8*(i-1)+5,')={',x_min,',',y_min,',',zmax,'};'
             WRITE(10,*) 'Point(',8*(i-1)+6,')={',x_max,',',y_min,',',zmax,'};'
             WRITE(10,*) 'Point(',8*(i-1)+7,')={',x_max,',',y_max,',',zmax,'};'
             WRITE(10,*) 'Point(',8*(i-1)+8,')={',x_min,',',y_max,',',zmax,'};'
             WRITE(10,*) 'Line(',12*(i-1)+1,')={',8*(i-1)+1,',',8*(i-1)+2,'};'
             WRITE(10,*) 'Line(',12*(i-1)+2,')={',8*(i-1)+2,',',8*(i-1)+3,'};'
             WRITE(10,*) 'Line(',12*(i-1)+3,')={',8*(i-1)+3,',',8*(i-1)+4,'};'
             WRITE(10,*) 'Line(',12*(i-1)+4,')={',8*(i-1)+4,',',8*(i-1)+1,'};'
             WRITE(10,*) 'Line(',12*(i-1)+5,')={',8*(i-1)+5,',',8*(i-1)+6,'};'
             WRITE(10,*) 'Line(',12*(i-1)+6,')={',8*(i-1)+6,',',8*(i-1)+7,'};'
             WRITE(10,*) 'Line(',12*(i-1)+7,')={',8*(i-1)+7,',',8*(i-1)+8,'};'
             WRITE(10,*) 'Line(',12*(i-1)+8,')={',8*(i-1)+8,',',8*(i-1)+5,'};'
             WRITE(10,*) 'Line(',12*(i-1)+9,')={',8*(i-1)+1,',',8*(i-1)+5,'};'
             WRITE(10,*) 'Line(',12*(i-1)+10,')={',8*(i-1)+2,',',8*(i-1)+6,'};'
             WRITE(10,*) 'Line(',12*(i-1)+11,')={',8*(i-1)+3,',',8*(i-1)+7,'};'
             WRITE(10,*) 'Line(',12*(i-1)+12,')={',8*(i-1)+4,',',8*(i-1)+8,'};'
          ELSE
             WRITE(10,*) 'Point(',4*(i-1)+1,')={',x_min,',',y_min,',0.};'
             WRITE(10,*) 'Point(',4*(i-1)+2,')={',x_max,',',y_min,',0.};'
             WRITE(10,*) 'Point(',4*(i-1)+3,')={',x_max,',',y_max,',0.};'
             WRITE(10,*) 'Point(',4*(i-1)+4,')={',x_min,',',y_max,',0.};'
             WRITE(10,*) 'Line(',4*(i-1)+1,')={',4*(i-1)+1,',',4*(i-1)+2,'};'
             WRITE(10,*) 'Line(',4*(i-1)+2,')={',4*(i-1)+2,',',4*(i-1)+3,'};'
             WRITE(10,*) 'Line(',4*(i-1)+3,')={',4*(i-1)+3,',',4*(i-1)+4,'};'
             WRITE(10,*) 'Line(',4*(i-1)+4,')={',4*(i-1)+4,',',4*(i-1)+1,'};'
          ENDIF
       ELSEIF (((params%coord==1).AND.(params%dim==3)) &
          .OR.((params%image==1).AND.(params%imgdim==3)) &
          .OR.((params%geom==1).AND.(params%imgdim+params%imgt==3))) THEN
          !3D
          IF (params%coord==1) THEN
             READ(2,*) x0,y0,z0,num,x1,y1,z1
             x_min=x0
             y_min=y0
             zmin=z0
             x_max=x1
             y_max=y1
             zmax=z1
          ELSEIF ((params%image==1).OR.(params%geom==1)) THEN
             READ(2,*) x0,y0,z0,num,x1,y1,z1
             x_min=y0
             y_min=-x0
             zmin=z0
             x_max=y1
             y_max=-x1
             zmax=z1
          ENDIF
          WRITE(10,*) 'Point(',8*(i-1)+1,')={',x_min,',',y_min,',',zmin,'};'
          WRITE(10,*) 'Point(',8*(i-1)+2,')={',x_max,',',y_min,',',zmin,'};'
          WRITE(10,*) 'Point(',8*(i-1)+3,')={',x_max,',',y_max,',',zmin,'};'
          WRITE(10,*) 'Point(',8*(i-1)+4,')={',x_min,',',y_max,',',zmin,'};'
          WRITE(10,*) 'Point(',8*(i-1)+5,')={',x_min,',',y_min,',',zmax,'};'
          WRITE(10,*) 'Point(',8*(i-1)+6,')={',x_max,',',y_min,',',zmax,'};'
          WRITE(10,*) 'Point(',8*(i-1)+7,')={',x_max,',',y_max,',',zmax,'};'
          WRITE(10,*) 'Point(',8*(i-1)+8,')={',x_min,',',y_max,',',zmax,'};'
          WRITE(10,*) 'Line(',12*(i-1)+1,')={',8*(i-1)+1,',',8*(i-1)+2,'};'
          WRITE(10,*) 'Line(',12*(i-1)+2,')={',8*(i-1)+2,',',8*(i-1)+3,'};'
          WRITE(10,*) 'Line(',12*(i-1)+3,')={',8*(i-1)+3,',',8*(i-1)+4,'};'
          WRITE(10,*) 'Line(',12*(i-1)+4,')={',8*(i-1)+4,',',8*(i-1)+1,'};'
          WRITE(10,*) 'Line(',12*(i-1)+5,')={',8*(i-1)+5,',',8*(i-1)+6,'};'
          WRITE(10,*) 'Line(',12*(i-1)+6,')={',8*(i-1)+6,',',8*(i-1)+7,'};'
          WRITE(10,*) 'Line(',12*(i-1)+7,')={',8*(i-1)+7,',',8*(i-1)+8,'};'
          WRITE(10,*) 'Line(',12*(i-1)+8,')={',8*(i-1)+8,',',8*(i-1)+5,'};'
          WRITE(10,*) 'Line(',12*(i-1)+9,')={',8*(i-1)+1,',',8*(i-1)+5,'};'
          WRITE(10,*) 'Line(',12*(i-1)+10,')={',8*(i-1)+2,',',8*(i-1)+6,'};'
          WRITE(10,*) 'Line(',12*(i-1)+11,')={',8*(i-1)+3,',',8*(i-1)+7,'};'
          WRITE(10,*) 'Line(',12*(i-1)+12,')={',8*(i-1)+4,',',8*(i-1)+8,'};'
       ENDIF
    ENDDO
    CLOSE(10)
    CLOSE(2)
    RETURN
  END SUBROUTINE write_partitioning_gmsh



  SUBROUTINE affectation_gmsh(params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params

    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coords
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: nb
    INTEGER :: offset
    INTEGER :: nb_slaves

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    OPEN(FILE='decoupe.visu',UNIT=1)
    WRITE(1,*) 'View "MPI" {'
    !Reading files
    IF (params%nbproc==1) THEN
       offset=1
       nb_slaves=1
    ELSE
       offset=0
       nb_slaves=params%nbproc-1
    ENDIF
    ALLOCATE(coords(1,params%dim))
    DO i=offset,nb_slaves
       ! File name
       WRITE(num,*) i
       files='decoupe.'//trim(adjustl(num))
       OPEN(FILE=files,UNIT=10)
       READ(10,*) nb
       PRINT *, '> ', i, ' : ', nb
       DO j=1,nb
          IF (params%coord==1) THEN
             ! Partitionning by coordinates
             coords(:,:)=0.
             READ(10,*) coords(1,:)
             CALL ecritpoint_gmsh(1,params%dim,coords,i,1)
          ELSE
             ! Partitionning 1D pictures
             READ (10,*) k
             CALL write_point_picture_format(1,params,i,k)
          ENDIF
       ENDDO
       CLOSE(10)
    ENDDO
    DEALLOCATE(coords)
    WRITE(1,*) '};'
    CLOSE(1)
    PRINT *, '-> decoupe.visu'
    RETURN
  END SUBROUTINE affectation_gmsh



  SUBROUTINE write_partial_clusters_gmsh(params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params

    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coords
    INTEGER, DIMENSION(:), POINTER :: matchings
    INTEGER :: i
    INTEGER :: ind ! TODO: renommer
    INTEGER :: j
    INTEGER :: k
    INTEGER :: length
    INTEGER :: nb

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    DO i=0,params%nbproc-1
       ! File name
       WRITE(num,*) i
       num=adjustl(num)
       length=len(trim(num))
       files='cluster.partiel.'//trim(num)
       OPEN(FILE=files,UNIT=10)
       ! Output
       files='cluster.partiel.'//num(1:length)//'.visu'
       OPEN(FILE=files,UNIT=1)
       WRITE(1,*) 'View "Clusters for part '//num(1:length)//'" {'
       READ(10,*) nb,k
       PRINT *, '> ', i, ' : ', nb, ' -> ', files
       ALLOCATE(coords(1,k))
       IF ((params%image==1).OR.(params%geom==1).OR.(params%seuil==1)) THEN
          ! Reading the matchings
          files='decoupe.'//trim(num)
          OPEN(FILE=files,UNIT=11)
          READ(11,*)
          ALLOCATE(matchings(nb))
          DO j=1,nb
             READ(11,*) matchings(j)
          ENDDO
          CLOSE(11)
       ENDIF
       DO j=1,nb
          IF (params%coord==1) THEN
             READ(10,*) coords(1,:),k
             CALL ecritpoint_gmsh(1,params%dim,coords,k,1)
          ELSE
             READ(10,*) k,ind
             CALL write_point_picture_format(1,params,ind,matchings(k))
          ENDIF
       ENDDO
       CLOSE(10) 
       WRITE(1,*) '};'
       CLOSE(1)
       DEALLOCATE(coords)
       IF ((params%image==1).OR.(params%geom==1).OR.(params%seuil==1)) THEN
          DEALLOCATE(matchings)
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE write_partial_clusters_gmsh



  SUBROUTINE write_final_clusters_gmsh(params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    !=== IN/OUT ===
    !====  OUT ====

    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coords
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: nb
    INTEGER :: nb0

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    OPEN(FILE=params%mesh,UNIT=1)
    IF (params%coord==1) THEN
       ! Reading with coordinates format
       READ(1,*) j,k
       PRINT *, 'Reading mesh file...', j, k
       ALLOCATE(coords(j,k))
       nb0=0
       DO i=1,j
          READ(1,*,END=100) coords(i,:)
          nb0=nb0+1
       ENDDO
100    PRINT *, 'Number of points : ', nb0
       j=nb0
       CLOSE(1)
    ENDIF

    nb0=k ! Points dimension
    OPEN(FILE='cluster.final.visu',UNIT=1)
    WRITE(1,*) 'View "Clusters" {'
    ! Reading files
    DO i=1,params%nbclusters
       ! File name
       WRITE(num,*) i
       files='cluster.final.'//trim(adjustl(num))
       OPEN(FILE=files,UNIT=10)
       READ(10,*) nb
       PRINT *,'  > ',i,' :',nb
       DO j=1,nb
          READ(10,*) k
          IF (params%coord==1) THEN
             ! Classic coordinates
             CALL ecritpoint_gmsh(1,nb0,coords,i,k)
          ELSE
             ! Pictures reassembly
             CALL write_point_picture_format(1,params,i,k)
          ENDIF
       ENDDO
       CLOSE(10)
    ENDDO
    WRITE(1,*) '};'
    CLOSE(1)
    PRINT *,'-> cluster.final.visu'
    IF (params%coord==1) DEALLOCATE(coords)
    RETURN
  END SUBROUTINE write_final_clusters_gmsh



  SUBROUTINE ecritpoint_gmsh(unit, dim, coords, ind, k)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coords
    INTEGER :: dim
    INTEGER :: ind
    INTEGER :: k
    INTEGER :: unit

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    IF (dim==2) THEN
       !2D
       WRITE(unit,*) 'SP(',coords(k,1),',',coords(k,2),',',0.,'){',ind,'};' 
    ELSEIF (dim==3) THEN
       !3D
       WRITE(unit,*) 'SP(',coords(k,1),',',coords(k,2),',',coords(k,3),'){',ind,'};'
    ENDIF
    RETURN
  END SUBROUTINE ecritpoint_gmsh



  SUBROUTINE write_point_picture_format(unit, params, ind, k)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    INTEGER :: ind
    INTEGER :: k
    INTEGER :: unit

    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:), POINTER :: data
    DOUBLE PRECISION :: kx
    DOUBLE PRECISION :: ky
    DOUBLE PRECISION :: kz
    INTEGER :: i
    INTEGER :: ix
    INTEGER :: iy

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    IF (((params%image==1).OR.(params%geom==1).OR.(params%seuil==1)) &
         .AND.(params%imgdim==2)) THEN
       ALLOCATE(data(params%nbp))
       OPEN(FILE=params%mesh,UNIT=50)
       READ(50,*)
       READ(50,*)
       DO i=1,params%nbp
          READ(50,*) data(i)
       ENDDO
       CLOSE(50)
    ENDIF
    ! Coordinates
    ix=params%refimg(k,1)
    iy=params%refimg(k,2)
    IF (params%imgdim==2) THEN
       ! 2D points
       IF (params%geom==1) THEN
          kx=iy*params%pas(2)
          ky=-ix*params%pas(1)
       ELSE
          kx=float(iy)
          ky=-float(ix)
       ENDIF
       IF (((params%image==1).OR.(params%geom==1).OR.(params%seuil==1)) &
            .AND.(params%imgdim==2)) THEN
          kz=data(k)
       ELSE
          kz=0.
       ENDIF
    ELSEIF (params%imgdim==3) THEN
       ! 3D points
       IF (params%geom==1) THEN
          kx=iy*params%pas(2)
          ky=-ix*params%pas(1)
          kz=float(params%refimg(k,3))*params%pas(3)
       ELSE
          kx=float(iy)
          ky=-float(ix)
          kz=float(params%refimg(k,3))
       ENDIF
    ENDIF
    ! Writing
    WRITE(unit,*) 'SP(',kx,',',ky,',',kz,'){',ind,'};'
    RETURN
  END SUBROUTINE write_point_picture_format



  SUBROUTINE list_commands_gmsh
    !###########################################
    ! INSTRUCTIONS
    !###########################################
    PRINT *,'gmsh decoupe.visu'
    PRINT *,'gmsh decoupe.geo'
    PRINT *,'gmsh decoupe.visu decoupe.geo'
    PRINT *,'gmsh cluster.partiel.*.visu'
    PRINT *,'gmsh cluster.final.visu'
    PRINT *,'gmsh decoupe.geo cluster.partiel.*.visu'
    PRINT *,'gmsh decoupe.geo cluster.final.visu'
    RETURN
  END SUBROUTINE list_commands_gmsh


END MODULE module_visuclusters_gmsh
