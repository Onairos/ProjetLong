MODULE module_visuclusters
  USE module_visuclusters_structure
  USE module_visuclusters_gmsh
  USE module_visuclusters_paraview
CONTAINS


  SUBROUTINE read_params(params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !=== IN/OUT ===
    TYPE(type_params) :: params
    
    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: n
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################     
    params%geom=0
    params%seuil=0
    PRINT *
    PRINT *, 'Reading info...'
    READ (3,*)
    READ(3,*) params%mesh
    PRINT *, '> Mesh file name : ', params%mesh
    READ (3,*)
    READ(3,*) params%nbp
    PRINT *, '> Number of points : ', params%nbp
    READ (3,*)
    READ (3,*) params%dim
    PRINT *, '> Dimension : ', params%dim
    READ (3,*)
    READ(3,*) params%nbproc
    PRINT *, '> Number of process : ', params%nbproc
    READ (3,*)
    READ(3,*) params%interface
    PRINT *, '> Partitioning by interfacing ? : ', params%interface
    READ (3,*)
    READ(3,*) params%recouvrement
    PRINT *, '> Partitioning by overlapping ? : ', params%recouvrement
    READ (3,*)
    READ(3,*) params%nbclusters  
    PRINT *, '> Number of clusters got : ', params%nbclusters
    READ (3,*)
    READ(3,*) params%coord
    PRINT *, '> Coordinates format ? : ', params%coord
    READ (3,*)
    READ(3,*) params%image
    PRINT *, '> Image format ? : ', params%image
    READ (3,*)
    READ(3,*) params%geom
    PRINT *, '> Geometric format ? : ', params%geom
    READ (3,*)
    READ(3,*) params%seuil
    PRINT *, '> Threshold format ? : ', params%seuil
    IF ((params%image==1).OR.(params%geom==1).OR.(params%seuil==1)) THEN
       READ (3,*)
       READ(3,*) params%imgdim
       PRINT *, '> Image dimension : ', params%imgdim
       ALLOCATE(params%imgmap(params%imgdim))
       READ (3,*)
       READ(3,*) params%imgmap(:)
       PRINT *, '> Image partitioning : ', params%imgmap
       READ (3,*)
       READ(3,*) params%imgt
       PRINT *, '> Number of time : ', params%imgt
       ! Points referencing
       ALLOCATE(params%refimg(params%nbp,params%imgdim))
       params%refimg(:,:)=0
       n=0
       DO i=1,params%imgmap(1)
          DO j=1,params%imgmap(2)
             IF (params%imgdim==2) THEN
                n=n+1
                params%refimg(n,1)=i
                params%refimg(n,2)=j
             ELSE
                DO k=1,params%imgmap(3)
                   n=n+1
                   params%refimg(n,1)=i
                   params%refimg(n,2)=j
                   params%refimg(n,3)=k
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       IF (params%geom==1) THEN
          ALLOCATE(params%pas(params%imgdim))
          READ (3,*)
          READ (3,*) params%pas(:)
          PRINT *, '> Mesh step : ', params%pas(:)
       ENDIF
    ENDIF
    RETURN
  END SUBROUTINE read_params


  SUBROUTINE write_partionning(format_output, params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    CHARACTER (LEN=30) :: format_output
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    PRINT *
    PRINT *, 'Writing partitioning geometry...'
    SELECT CASE(format_output)
    CASE('gmsh')
       CALL write_partionning_gmsh(params)
    CASE('paraview')
       CALL write_partitioning_paraview(params)
    END SELECT
    RETURN
  END SUBROUTINE write_partionning




  SUBROUTINE affectation(format_output, params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    CHARACTER (LEN=30) :: format_output
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    PRINT *
    PRINT *, 'Writing partitioning allocations...'
    SELECT CASE(format_output)
    CASE('gmsh')
       CALL affectation_gmsh(params)
    CASE('paraview')
       CALL affectation_paraview(params)
    END SELECT
    RETURN
  END SUBROUTINE affectation




  SUBROUTINE write_partial_clusters(format_output, params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    CHARACTER (LEN=30) :: format_output

    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    PRINT *
    PRINT *, 'Reading clusters before grouping...'
    SELECT CASE(format_output)
    CASE('gmsh')
       CALL write_partial_clusters_gmsh(params)
    CASE('paraview')
       CALL write_partial_clusters_paraview(params)
    END SELECT
    RETURN
  END SUBROUTINE write_partial_clusters




  SUBROUTINE write_final_clusters(format_output, params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    CHARACTER (LEN=30) :: format_output

    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
    PRINT *
    PRINT *, 'Reading clusters after grouping...'
    SELECT CASE(format_output)
    CASE('gmsh')
       CALL write_final_clusters_gmsh(params)
    CASE('paraview')
       CALL write_final_clusters_paraview(params)
    END SELECT
    RETURN
  END SUBROUTINE write_final_clusters




  SUBROUTINE list_commands(format_output)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    CHARACTER (LEN=30) :: format_output
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
    PRINT *
    PRINT *, '-------------------------------------'
    PRINT *, 'Command list of visualisation : '
    PRINT *
    SELECT CASE(format_output)
    CASE('gmsh')
       CALL list_commands_gmsh
    CASE('paraview')
       CALL list_commands_paraview
    END SELECT
    RETURN
  END SUBROUTINE list_commands


END MODULE module_visuclusters
