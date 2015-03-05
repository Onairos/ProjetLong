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
    PRINT *,'lecture des infos...'
    READ (3,*)
    READ(3,*) params%mesh
    PRINT *,'  > nom du fichier de maillage : ',params%mesh
    READ (3,*)
    READ(3,*) params%nbp
    PRINT *,'  > nb de points :',params%nbp
    READ (3,*)
    READ (3,*) params%dim
    PRINT *,'  > dimension : ',params%dim
    READ (3,*)
    READ(3,*) params%nbproc
    PRINT *,'  > nb de proc :',params%nbproc
    READ (3,*)
    READ(3,*) params%interface
    PRINT *,'  > decoupage par interface ?',params%interface
    READ (3,*)
    READ(3,*) params%recouvrement
    PRINT *,'  > decoupage par recouvrement ?',params%recouvrement
    READ (3,*)
    READ(3,*) params%nbclusters  
    PRINT *,'  > nb de clusters obtenus :',params%nbclusters
    READ (3,*)
    READ(3,*) params%coord
    PRINT *,'  > format coord ?',params%coord
    READ (3,*)
    READ(3,*) params%image
    PRINT *,'  > format image ?',params%image
    READ (3,*)
    READ(3,*) params%geom
    PRINT *,'  > format geom ?',params%geom
    READ (3,*)
    READ(3,*) params%seuil
    PRINT *,'  > format seuil ?',params%seuil
    IF ((params%image==1).OR.(params%geom==1).OR.(params%seuil==1)) THEN
       READ (3,*)
       READ(3,*) params%imgdim
       PRINT *,'  > dimension image :',params%imgdim
       ALLOCATE(params%imgmap(params%imgdim))
       READ (3,*)
       READ(3,*) params%imgmap(:)
       PRINT *,'  > decoupage image :',params%imgmap
       READ (3,*)
       READ(3,*) params%imgt
       PRINT *,'  > nb de temps :',params%imgt
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
          PRINT *,'    > pas de maillage :',params%pas(:)
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
    PRINT *,'ecriture de la geometrie du decoupage...'
    SELECT CASE(format_output)
    CASE('gmsh')
       CALL write_partionning_gmsh(params)
    CASE('paraview')
       CALL write_partionning_paraview(params)
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
    PRINT *,'ecriture des affectations du decoupage...'
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
    PRINT *,'lecture des clusters avant regroupement...'
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
    PRINT *,'lecture des clusters apres regroupement...'
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
    PRINT *,'-------------------------------------'
    PRINT *,'liste de commandes de visualisation :'
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
