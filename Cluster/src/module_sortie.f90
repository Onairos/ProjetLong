MODULE module_sortie
  USE module_structure
CONTAINS


  SUBROUTINE write_domains(data, nbproc, domaines)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: domaines  
    INTEGER :: nbproc
    
    !#### Variables  ####
    INTEGER :: i
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    DO i=1,nbproc-data%interface
       WRITE(2,*) domaines(i,:,1),'|', domaines(i,:,2)
    ENDDO
    CALL flush(2)
    CLOSE(2)
    RETURN
  END SUBROUTINE write_domains


  SUBROUTINE write_partitionning(nbproc, data, ldat, ddat)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    INTEGER,DIMENSION(:,:),POINTER :: ddat
    INTEGER,DIMENSION(:),POINTER :: ldat
    INTEGER :: nbproc
    
    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nbdom
    INTEGER :: offset
    
    !###########################################      
    ! INSTRUCTIONS
    !########################################### 
    PRINT *,'  > bilan decoupage :'
    offset=1
    nbdom=nbproc
    IF ((data%interface==1).AND.(nbproc>1)) THEN
       offset=0
       nbdom=nbproc-1
    ENDIF
    IF (data%recouvrement==1) THEN
       offset=0
       nbdom=nbproc-1
    ENDIF
    DO i=offset,nbdom
       PRINT *,'    > zone ',i,':',ldat(i)
       !nom du fichier
       WRITE(num,*) i
       files='decoupe.'//trim(adjustl(num))
       OPEN(FILE=files,UNIT=10)
       WRITE(10,*) ldat(i)
       DO j=1,ldat(i)
          IF (data%coord==1) THEN
             !ecriture en coordonnees
             WRITE(10,*) data%point(ddat(i,j))%coord(:)
          ELSEIF ((data%image==1).OR.(data%seuil==1).OR.(data%geom==1)) THEN
             !ecriture en image
             WRITE(10,*) ddat(i,j)
          ENDIF
       ENDDO
       CALL flush(10)
       CLOSE(10)
    ENDDO
    CALL flush(6)
    RETURN
  END SUBROUTINE write_partitionning


  SUBROUTINE write_partial_clusters(numproc, dataw)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====  
    TYPE(type_data) :: dataw
    INTEGER :: numproc
    
    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    INTEGER :: i
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################
    !nom du fichier
    WRITE(num,*),numproc
    num=adjustl(num)
    files='cluster.partiel.'//trim(num)
    PRINT *,numproc,'ecriture des clusters : ',files
    OPEN(FILE=files,UNIT=10)
    WRITE(10,*) dataw%nb,dataw%dim
    DO i=1,dataw%nb
       IF (dataw%coord==1) THEN
          WRITE(10,*) dataw%point(i)%coord(:),dataw%point(i)%cluster
       ELSE
          WRITE(10,*) i,dataw%point(i)%cluster
       ENDIF
    ENDDO
    CALL flush(10)
    CLOSE(10)
    RETURN
  END SUBROUTINE write_partial_clusters


  SUBROUTINE write_final_clusters(nbclust, iclust, clustermap)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER,DIMENSION(:,:),POINTER :: clustermap
    INTEGER,DIMENSION(:),POINTER :: iclust

    !=== IN/OUT === 
    INTEGER :: nbclust

    !#### Variables  ####
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    
    !###########################################      
    ! INSTRUCTIONS
    !########################################### 
    PRINT *,'  > Ecriture du resultat...'
    k=0
    DO i=1,nbclust
       IF (iclust(i)>0) THEN
          k=k+1
          WRITE(num,*) k
          files='cluster.final.'//trim(adjustl(num))
          OPEN(FILE=files,UNIT=20)     
          PRINT *,'    > cluster ',k,' :',iclust(i),' -> ',files
          WRITE(20,*) iclust(i)
          DO j=1,iclust(i)
             WRITE(20,*) clustermap(i,j)
          ENDDO
          CALL flush(20)
          CLOSE(20)
       ENDIF
    ENDDO
    nbclust=k
    RETURN
  END SUBROUTINE write_final_clusters


  SUBROUTINE write_metadata(mesh, data, nbproc, nbclust)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    CHARACTER (LEN=30) :: mesh
    INTEGER :: nbclust
    INTEGER :: nbproc
    
    !###########################################      
    ! INSTRUCTIONS
    !########################################### 
    WRITE(3,*) '# fichier de maillage :'
    WRITE(3,*) mesh
    WRITE(3,*) '#nb de points :'
    WRITE(3,*) data%nb
    WRITE(3,*) '# DIMENSION :'
    WRITE(3,*) data%dim
    WRITE(3,*) '# nb de proc :'
    WRITE(3,*) nbproc
    WRITE(3,*) '# decoupage par interface :'
    WRITE(3,*) data%interface
    WRITE(3,*) '# decoupage par recouvrement :'
    WRITE(3,*) data%recouvrement
    WRITE(3,*) '# nb de clusters :'
    WRITE(3,*) nbclust
    WRITE(3,*) '# format coord :'
    WRITE(3,*) data%coord
    WRITE(3,*) '# format image :'
    WRITE(3,*) data%image
    WRITE(3,*) '# format geom :'
    WRITE(3,*) data%geom
    WRITE(3,*) '# format seuil :'
    WRITE(3,*) data%seuil
    IF ((data%image==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
       WRITE(3,*) '# DIMENSION :'
       WRITE(3,*) data%imgdim
       WRITE(3,*) '# decoupage :'
       WRITE(3,*) data%imgmap(:)
       WRITE(3,*) '# nb de temps :'
       WRITE(3,*) data%imgt
       IF (data%geom==1) THEN
          WRITE(3,*) '## pas de maillage :'
          WRITE(3,*) data%pas(:)
       ENDIF
    ENDIF
    CALL flush(3)
    CLOSE(3)
    RETURN
  END SUBROUTINE write_metadata

END MODULE module_sortie
