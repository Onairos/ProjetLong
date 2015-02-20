MODULE module_sortie
  USE module_structure
CONTAINS

  !********************************
  !ecriture des domaines decoupes
  SUBROUTINE ecrit_domaines(data,nbproc,domaines)
    IMPLICIT NONE
    TYPE(type_data) :: data
    INTEGER :: nbproc
    DOUBLE PRECISION,DIMENSION(:,:,:),POINTER :: domaines
    INTEGER :: i
    DO i=1,nbproc-data%interface
       WRITE(2,*) domaines(i,:,1),'|', domaines(i,:,2)
    ENDDO
    CALL flush(2)
    CLOSE(2)
    RETURN
  END SUBROUTINE ecrit_domaines

  !********************************
  !ecriture des decoupages
  SUBROUTINE ecrit_decoupages(nbproc,data,ldat,ddat)
    IMPLICIT NONE
    INTEGER :: nbproc
    TYPE(type_data) :: data
    INTEGER,DIMENSION(:),POINTER :: ldat
    INTEGER,DIMENSION(:,:),POINTER :: ddat
    INTEGER :: i,j,offset,nbdom
    CHARACTER*30 :: files,num
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
  END SUBROUTINE ecrit_decoupages

  !********************************
  !ecriture des clusters regroupes
  SUBROUTINE ecritcluster(numproc,dataw)
    IMPLICIT NONE
    INTEGER :: numproc
    TYPE(type_data) :: dataw
    INTEGER :: i
    CHARACTER*30 :: files,num
    !nom du fichier
    WRITE(num,*),numproc
    num=adjustl(num)
    files='cluster.partiel.'//trim(num)
    !len=len(trim(num))
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
  END SUBROUTINE ecritcluster

  !****************************
  !ecriture de cluster.final.
  SUBROUTINE ecritclusterfinal(nbclust,iclust,clustermap)
    IMPLICIT NONE
    INTEGER :: nbclust
    INTEGER,DIMENSION(:),POINTER :: iclust
    INTEGER,DIMENSION(:,:),POINTER :: clustermap
    INTEGER :: i,j,k
    CHARACTER*30 :: files,num
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
  END SUBROUTINE ecritclusterfinal

  !***************************
  !ecriture des informations
  SUBROUTINE ecrit_info(mesh,data,nbproc,nbclust)
    IMPLICIT NONE
    CHARACTER*30 :: mesh
    INTEGER :: nbproc,nbclust
    TYPE(type_data) :: data
    !ecriture du fichier fort.3
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
  END SUBROUTINE ecrit_info

END MODULE module_sortie
