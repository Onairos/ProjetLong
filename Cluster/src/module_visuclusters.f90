!>Contains methods enabling writing results in a selected data file format (for now : Paraview or GMSH)
MODULE module_visuclusters
  USE module_visuclusters_structure
  USE module_visuclusters_gmsh
  USE module_visuclusters_paraview
CONTAINS


!>Reads a file containing metadata on the data and the computed clusters
!!@details This function extracts the following information from the input 
!!<em>fort.3</em> file :
!!<ol>
!!<li> The data file name </li>
!!<li> The number of the points in the entire data set </li>
!!<li> The number of processes used </li>
!!<li> The partitioning mode (by interface or overlapping) </li>
!!<li> The number of clusters found </li>
!!<li> The data file format </li>
!!</ol>
!!In the case of a picture format: this extra information is
!!written :
!!<ol>
!!<li> The image dimension </li>
!!<li> The image partitioning </li>
!!<li> The number of attributes </li>
!!<li> The number of steps (only in geometric format) </li>
!!</ol> 
!!@see write_metadata()
!! @param[in,out] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
  SUBROUTINE read_metadata(params)
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
    READ(3,*) params%nb_clusters  
    PRINT *, '> Number of clusters got : ', params%nb_clusters
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
  END SUBROUTINE read_metadata


!>Writes the geometry of the partitioning (Gmsh or Paraview) and calls the eponym function
!!@details This methods extracts the domain definitions
!!from the <em>fort.2</em> file.
!!@see module_calcul::write_partitioning()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
!! @param[in] format_output the file format for visualization
  SUBROUTINE write_partitioning(format_output, params)
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
       CALL write_partitioning_gmsh(params)
    CASE('paraview')
       CALL write_partitioning_paraview(params)
    END SELECT
    RETURN
  END SUBROUTINE write_partitioning




!>Initializes the file of the partitionning
!!@details This method extracts details on partitioning from the
!!<em>decoupe.x</em> files.
!!@see module_calcul::write_partial_clusters()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
!! @param[in] format_output the file format for visualization
  SUBROUTINE write_assignment(format_output, params)
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
       CALL write_assignment_gmsh(params)
    CASE('paraview')
       CALL write_assignment_paraview(params)
    END SELECT
    RETURN
  END SUBROUTINE write_assignment




!>Writes the clusters before grouping by calling the corresponding method (Gmsh or Paraview)
!!@details This methods extracts details on computed clusters
!!on each domain from <em>cluster.partiel.x</em> files.
!!@see module_calcul::write_partial_clusters()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
!! @param[in] format_output the file format for visualization
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




!>Writes the clusters after grouping by calling the corresponding method (Gmsh or Paraview)
!!@details This methods extracts details on computed clusters
!!from <em>cluster.final.x</em> files.
!!@see module_calcul::write_final_clusters()
!! @param[in] params the parameters defined in the \latexonly\textit{param.in}\endlatexonly\htmlonly<cite>param.in</cite>\endhtmlonly file
!! @param[in] format_output the file format for visualization
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




!>Lists the commands related to Gmsh or Paraview depending on the input parameter
!! @param[in] format_output the file format for visualization
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
