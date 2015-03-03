MODULE module_decoupe
  USE module_structure
  USE module_sortie
CONTAINS


  SUBROUTINE partition_data(data, epsilon, nbproc, coordmin, coordmax, decoupe,&
       ldat, ddat, bornes)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION :: epsilon
    INTEGER, DIMENSION(:), POINTER :: decoupe
    INTEGER :: nbproc

    !=== IN/OUT ===
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmax
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmin

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bornes
    INTEGER, DIMENSION(:,:), POINTER :: ddat
    INTEGER, DIMENSION(:), POINTER :: ldat

    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: domaines

    !###########################################
    ! INSTRUCTIONS
    !########################################### 
    !definition des bornes
    CALL define_bounds(data,coordmin,coordmax,bornes,decoupe,epsilon,nbproc)

    !definition des domaines
    CALL define_domains(nbproc,data,domaines,bornes,decoupe)

    !ecriture des domaines decoupes
    CALL write_domains(data,nbproc,domaines)

    !definition des decoupages
    IF ((data%interface==1).OR.(nbproc==1)) THEN
       !decoupage avec interface
       CALL partition_with_interfaces(nbproc,data,ldat,ddat,domaines,epsilon)
    ELSE
       !decoupage avec recouvrement
       CALL partition_with_overlappings(nbproc,data,ldat,ddat,domaines)
    ENDIF
    DEALLOCATE(domaines)

    !sauvegarde des decoupages
    CALL write_partitionning(nbproc,data,ldat,ddat)

    RETURN
  END SUBROUTINE partition_data


  SUBROUTINE define_bounds(data, coordmin, coordmax, bornes, decoupe, epsilon, nbproc)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION :: epsilon
    INTEGER,DIMENSION(:), POINTER :: decoupe
    INTEGER :: nbproc

    !=== IN/OUT ===
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmax
    DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmin

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bornes

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
    !volume max
    DO i=1,data%dim
       prod=prod*(coordmax(i)-coordmin(i))
    ENDDO
    files='diminterface'
    WRITE(num,*),0
    num=adjustl(num)
    files=trim(files)//'.'//trim(num)
    OPEN(FILE=files,UNIT=20)
    DO i=1,data%dim
       som1=som1*(decoupe(i)-1)
       prod2=prod2+(decoupe(i)-1)*prod/(coordmax(i)-coordmin(i))    
    ENDDO    
    WRITE(20,*)  prod,epsilon*prod2-som1*(epsilon)**data%dim
    CLOSE(20)
 



 IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
       !traitement en coordonnees ou image en coodonnees ou image seuillee
       ALLOCATE(bornes(data%dim,max(nbproc-data%interface,1),2))
       bornes(:,:,:)=0.0
       DO i=1,data%dim
          coordmin(i)=coordmin(i)-epsilon*1.1
          coordmax(i)=coordmax(i)+epsilon*1.1
          DO j=1,decoupe(i)
             bornes(i,j,1)=coordmin(i)+(j-1)*(coordmax(i)-coordmin(i))/decoupe(i)
             bornes(i,j,2)=coordmin(i)+j*(coordmax(i)-coordmin(i))/decoupe(i)
          ENDDO
          IF (data%recouvrement==1) THEN
             !decoupage en mode interface
             DO j=1,decoupe(i)
                bornes(i,j,1)=bornes(i,j,1)-epsilon
                bornes(i,j,2)=bornes(i,j,2)+epsilon
             ENDDO
          ENDIF
          bornes(i,1,1)=coordmin(i)-0.01*abs(coordmin(i))
          bornes(i,decoupe(i),2)=coordmax(i)+0.01*abs(coordmax(i))
       ENDDO
    ELSEIF (data%image==1) THEN
       !traitement pour decoupage en pixel d'images
       ALLOCATE(bornes(data%imgdim,max(nbproc-1,1),2));bornes(:,:,:)=0.0
       IF ((data%imgdim/=2).AND.(data%imgdim/=3)) THEN
#if aff
          PRINT *
          PRINT *,'format d images /= 2d, 3d pas supporte !!!!'
#endif
          STOP
       ENDIF
       DO i=1,data%imgdim
          coordmin(i)=1.0-epsilon*1.1
          coordmax(i)=data%imgmap(i)+epsilon*1.1
          DO j=1,decoupe(i)
             bornes(i,j,1)=coordmin(i)+(j-1)*(coordmax(i)-coordmin(i))/decoupe(i)
             bornes(i,j,2)=coordmin(i)+j*(coordmax(i)-coordmin(i))/decoupe(i)
          ENDDO
          IF (data%recouvrement==1) THEN
             !decoupage en mode interface
             DO j=1,decoupe(i)
                bornes(i,j,1)=bornes(i,j,1)-epsilon
                bornes(i,j,2)=bornes(i,j,2)+epsilon
             ENDDO
          ENDIF
          bornes(i,1,1)=coordmin(i)-0.01*abs(coordmin(i))
          bornes(i,decoupe(i),2)=coordmax(i)+0.01*abs(coordmax(i))
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE define_bounds


  SUBROUTINE define_domains(nbproc, data, domaines, bornes, decoupe)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bornes
    INTEGER, DIMENSION(:), POINTER :: decoupe
    INTEGER :: nbproc

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: domaines

    !#### Variables  ####
    INTEGER, DIMENSION(:), POINTER :: list
    INTEGER :: k
    INTEGER :: n
    INTEGER :: ok !TODO utilisÃÂÃÂ© comme un booleen, modifier en LOGICAL ??


    !###########################################
    ! INSTRUCTIONS
    !###########################################
    IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
       !traitement en coordonnees ou image en coodonnees ou image seuillee
       ALLOCATE(domaines(max(1,nbproc-data%interface),data%dim,2))
       domaines(:,:,:)=0.0
       ALLOCATE(list(data%dim)); list(:)=1
       IF (nbproc>1) THEN
          !** mode >1 proc
          DO n=1,nbproc-data%interface
             DO k=1,data%dim
                domaines(n,k,:)=bornes(k,list(k),:)
             ENDDO
             ok=1
             DO k=data%dim,1,-1
                IF (ok==1) THEN
                   list(k)=list(k)+1
                   IF (list(k)>decoupe(k)) THEN
                      list(k)=1
                   ELSE
                      ok=0
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ELSE
          !** mode 1 proc
          DO k=1,data%dim
             domaines(1,k,:)=bornes(k,1,:)
          ENDDO
       ENDIF
       DEALLOCATE(list)
    ELSEIF (data%image==1) THEN
       !traitement pour decoupage en pixel d'images
       ALLOCATE(domaines(max(1,nbproc-data%interface),data%imgdim,2))
       domaines(:,:,:)=0.0
       ALLOCATE(list(data%imgdim)); list(:)=1
       IF (nbproc>1) THEN
          !** mode >1 proc
          DO n=1,nbproc-data%interface
             DO k=1,data%imgdim
                domaines(n,k,:)=bornes(k,list(k),:)
             ENDDO
             ok=1
             DO k=data%imgdim,1,-1
                IF (ok==1) THEN
                   list(k)=list(k)+1
                   IF (list(k)>decoupe(k)) THEN
                      list(k)=1
                   ELSE
                      ok=0
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ELSE
          !** mode 1 proc
          DO k=1,data%imgdim
             domaines(1,k,:)=bornes(k,1,:)
          ENDDO
       ENDIF
       DEALLOCATE(list)
    ENDIF
    RETURN
  END SUBROUTINE define_domains


  SUBROUTINE partition_with_interfaces(nbproc, data, ldat, ddat, domaines, epsilon)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: domaines
    DOUBLE PRECISION :: epsilon
    INTEGER :: nbproc
    !====  OUT ====
    INTEGER, DIMENSION(:,:), POINTER :: ddat
    INTEGER, DIMENSION(:), POINTER :: ldat

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: ierr
    INTEGER :: j
    INTEGER :: ok
    INTEGER :: n

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ALLOCATE(ldat(0:max(1,nbproc-1)))
    ldat(:)=0
    ALLOCATE(ddat(0:max(1,nbproc-1),data%nb))
    ddat(:,:)=0
    DO i=1,data%nb
       !recherche des paquets
       n=0; ok=0
       DO WHILE(ok==0)
          n=n+1; ok=1
          IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
             !traitement en coordonnees ou image en coodonnees ou image seuillee
             DO j=1,data%dim
                IF ((data%point(i)%coord(j)>domaines(n,j,2)).OR.&
                     (data%point(i)%coord(j)<domaines(n,j,1))) ok=0
             ENDDO
          ELSEIF (data%image==1) THEN
             !traitement pour decoupage en pixel d'images
             DO j=1,data%imgdim
                IF ((data%refimg(i,j)>domaines(n,j,2)).OR.&
                     (data%refimg(i,j)<domaines(n,j,1))) ok=0
             ENDDO
          ENDIF
          IF ((n>nbproc-1).AND.(nbproc>1)) THEN
#if aff
             PRINT *,'bug dans le decoupage !',n,nbproc-1
#endif
             IF (data%geom==0) THEN
#if aff
                PRINT *,data%point(i)%coord(:)
#endif
             ELSE
#if aff
                PRINT *,i,data%refimg(i,:)
#endif
             ENDIF
             CALL MPI_ABORT(ierr)
             STOP
          ENDIF
       ENDDO
       ldat(n)=ldat(n)+1
       ddat(n,ldat(n))=i
       IF (nbproc>1) THEN
          !** recherche de l'interface si plus de 1 proc
          ok=0
          IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
             !traitement en coordonnees ou image en coodonnees ou image seuillee
             DO j=1,data%dim
                IF ((abs(data%point(i)%coord(j)-domaines(n,j,1))<epsilon).OR.&
                     (abs(data%point(i)%coord(j)-domaines(n,j,2))<epsilon)) ok=1
             ENDDO
          ELSEIF (data%image==1) THEN
             !traitement pour decoupage en pixel d'images
             DO j=1,data%imgdim
                IF ((abs(data%refimg(i,j)-domaines(n,j,1))<epsilon).OR.&
                     (abs(data%refimg(i,j)-domaines(n,j,2))<epsilon)) ok=1
             ENDDO
          ENDIF
          IF (ok==1) THEN
             ldat(0)=ldat(0)+1
             ddat(0,ldat(0))=i
             WRITE(7,*) ddat(0,ldat(0))
          ENDIF
       ENDIF
    ENDDO
    WRITE(7,*) ldat(0)
    RETURN
  END SUBROUTINE partition_with_interfaces



  SUBROUTINE partition_with_overlappings(nbproc, data, ldat, ddat, domaines)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: data
    DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: domaines
    INTEGER :: nbproc
    !====  OUT ====
    INTEGER, DIMENSION(:,:), POINTER :: ddat
    INTEGER, DIMENSION(:), POINTER :: ldat

    !#### Variables  ####
    INTEGER :: i
    INTEGER :: j
    INTEGER :: n
    INTEGER :: ok

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ALLOCATE(ldat(0:max(1,nbproc-1)))
    ldat(:)=0
    ALLOCATE(ddat(0:max(1,nbproc-1),data%nb))
    ddat(:,:)=0
    DO i=1,data%nb
       !recherche des paquets
       DO n=1,nbproc
          ok=1
          IF ((data%coord==1).OR.(data%geom==1).OR.(data%seuil==1)) THEN
             !traitement en coordonnees ou image en coodonnees ou image seuillee
             DO j=1,data%dim
                IF ((data%point(i)%coord(j)>domaines(n,j,2)).OR.&
                     (data%point(i)%coord(j)<domaines(n,j,1))) ok=0
             ENDDO
          ELSEIF (data%image==1) THEN
             !traitement pour decoupage en pixel d'images
             DO j=1,data%imgdim
                IF ((data%refimg(i,j)>domaines(n,j,2)).OR.&
                     (data%refimg(i,j)<domaines(n,j,1))) ok=0
             ENDDO
          ENDIF
          IF (ok==1) THEN
             ldat(n-1)=ldat(n-1)+1
             ddat(n-1,ldat(n-1))=i
          ENDIF
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE partition_with_overlappings


  SUBROUTINE group_clusters(nbclust, iclust, clustermap, data)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: nbclust

    !=== IN/OUT ===
    INTEGER, DIMENSION(:,:), POINTER :: clustermap
    INTEGER, DIMENSION(:), POINTER :: iclust

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
    INTEGER :: ok !TODO utilisÃÂÃÂ© comme un booleen, modifier en LOGICAL ??
    INTEGER :: ok2 !TODO utilisÃÂÃÂ© comme un booleen, modifier en LOGICAL ??
    INTEGER :: ok3 !TODO utilisÃÂÃÂ© comme un booleen, modifier en LOGICAL ??

    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ok=0; i=1; j=0
#if aff
    PRINT *,'  > elimination des doublons...'
    PRINT *,'    > regroupement du sous-cluster ',1
#endif
    DO WHILE(ok==0)
       j=j+1 
       IF (j>iclust(i)) THEN
          !la ligne i est completement testee
#if aff
          PRINT *,'      > nb d elements apres regroupement :',iclust(i)
#endif
          i=i+1; j=1
#if aff
          PRINT *,'    > regroupement du cluster ',i
#endif
       ENDIF
       IF (i>nbclust-1) THEN
          !plus de points a tester
          ok=1
       ELSEIF (iclust(i)>0) THEN
          !stockage de l'indice
          data%point(clustermap(i,j))%cluster=i
          !test des recouvrement
          ok2=0; i2=i+1; j2=1
          DO WHILE(ok2==0)
             IF (j2>iclust(i2)) THEN
                !ligne i2 completement teste pour le point (i,j)
                i2=i2+1; j2=1
             ENDIF
             IF (i2>nbclust) THEN
                !fin des tests pour le point (i,j)
                ok2=1
             ELSE
                !test d'intersections :
                IF (clustermap(i,j)==clustermap(i2,j2)) THEN
                   !intersection trouvee :
                   !ligne i2 ajoutee a la ligne i
                   n=0
                   DO k=1,iclust(i2)
                      !test d'elimination de doublon
                      ok3=1
                      DO j3=1,iclust(i)
                         IF (clustermap(i2,k)==clustermap(i,j3)) ok3=0
                      ENDDO
                      IF (ok3==1) THEN
                         n=n+1
                         clustermap(i,iclust(i)+n)=clustermap(i2,k)
                         clustermap(i2,k)=0                         
                      ENDIF
                   ENDDO
                   iclust(i)=iclust(i)+n
                   iclust(i2)=0
                ELSE
                   !test d'un nouveau point
                   j2=j2+1
                ENDIF
             ENDIF
          ENDDO
       ENDIF
    ENDDO
#if aff
    PRINT *,'      > nb d elements apres regroupement :',iclust(i)
#endif
    RETURN
  END SUBROUTINE group_clusters

END MODULE module_decoupe
