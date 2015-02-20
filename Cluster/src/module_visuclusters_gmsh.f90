MODULE module_visuclusters_gmsh
  USE module_visuclusters_structure
CONTAINS

  !************************
  !ecriture de la geometrie du decoupage
  SUBROUTINE ecrit_decoupage_gmsh(params)
    IMPLICIT NONE
    TYPE(type_params) :: params
    INTEGER :: i
    CHARACTER*30 :: num
    DOUBLE PRECISION :: xmin,ymin,zmin,xmax,ymax,zmax
    DOUBLE PRECISION :: x0,y0,z0,x1,y1,z1
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
             xmin=x0
             ymin=y0
             xmax=x1
             ymax=y1
          ELSEIF (params%seuil==1) THEN
             READ(2,*) x0,num,x1
             xmin=1
             ymin=-1
             xmax=params%imgmap(2)
             ymax=-params%imgmap(1)
             zmin=x0
             zmax=x1
          ELSEIF ((params%image==1).OR.(params%geom==1)) THEN
             READ(2,*) x0,y0,num,x1,y1
             xmin=y0
             ymin=-x0
             xmax=y1
             ymax=-x1
             zmin=x0
          ENDIF
          IF (params%seuil==1) THEN
             WRITE(10,*) 'Point(',8*(i-1)+1,')={',xmin,',',ymin,',',zmin,'};'
             WRITE(10,*) 'Point(',8*(i-1)+2,')={',xmax,',',ymin,',',zmin,'};'
             WRITE(10,*) 'Point(',8*(i-1)+3,')={',xmax,',',ymax,',',zmin,'};'
             WRITE(10,*) 'Point(',8*(i-1)+4,')={',xmin,',',ymax,',',zmin,'};'
             WRITE(10,*) 'Point(',8*(i-1)+5,')={',xmin,',',ymin,',',zmax,'};'
             WRITE(10,*) 'Point(',8*(i-1)+6,')={',xmax,',',ymin,',',zmax,'};'
             WRITE(10,*) 'Point(',8*(i-1)+7,')={',xmax,',',ymax,',',zmax,'};'
             WRITE(10,*) 'Point(',8*(i-1)+8,')={',xmin,',',ymax,',',zmax,'};'
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
             WRITE(10,*) 'Point(',4*(i-1)+1,')={',xmin,',',ymin,',0.};'
             WRITE(10,*) 'Point(',4*(i-1)+2,')={',xmax,',',ymin,',0.};'
             WRITE(10,*) 'Point(',4*(i-1)+3,')={',xmax,',',ymax,',0.};'
             WRITE(10,*) 'Point(',4*(i-1)+4,')={',xmin,',',ymax,',0.};'
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
             xmin=x0
             ymin=y0
             zmin=z0
             xmax=x1
             ymax=y1
             zmax=z1
          ELSEIF ((params%image==1).OR.(params%geom==1)) THEN
             READ(2,*) x0,y0,z0,num,x1,y1,z1
             !xmin=-x0
             !ymin=y0
             !zmin=z0
             !xmax=-x1
             !ymax=y1
             !zmax=z1
             xmin=y0
             ymin=-x0
             zmin=z0
             xmax=y1
             ymax=-x1
             zmax=z1
          ENDIF
          WRITE(10,*) 'Point(',8*(i-1)+1,')={',xmin,',',ymin,',',zmin,'};'
          WRITE(10,*) 'Point(',8*(i-1)+2,')={',xmax,',',ymin,',',zmin,'};'
          WRITE(10,*) 'Point(',8*(i-1)+3,')={',xmax,',',ymax,',',zmin,'};'
          WRITE(10,*) 'Point(',8*(i-1)+4,')={',xmin,',',ymax,',',zmin,'};'
          WRITE(10,*) 'Point(',8*(i-1)+5,')={',xmin,',',ymin,',',zmax,'};'
          WRITE(10,*) 'Point(',8*(i-1)+6,')={',xmax,',',ymin,',',zmax,'};'
          WRITE(10,*) 'Point(',8*(i-1)+7,')={',xmax,',',ymax,',',zmax,'};'
          WRITE(10,*) 'Point(',8*(i-1)+8,')={',xmin,',',ymax,',',zmax,'};'
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
  END SUBROUTINE ecrit_decoupage_gmsh

  !***********************
  !initialisation du fichier de decoupage
  SUBROUTINE affectation_gmsh(params)
    IMPLICIT NONE
    TYPE(type_params) :: params
    DOUBLE PRECISION,DIMENSION(:,:),POINTER :: coord
    CHARACTER*30 :: files,num
    INTEGER :: i,j,k,nb,offset,totnum
    OPEN(FILE='decoupe.visu',UNIT=1)
    WRITE(1,*) 'View "MPI" {'
    !lecture des fichiers
    IF (params%nbproc==1) THEN
       offset=1
       totnum=1
    ELSE
       offset=0
       totnum=params%nbproc-1
    ENDIF
    ALLOCATE(coord(1,params%dim))
    DO i=offset,totnum
       !nom du fichier
       WRITE(num,*) i
       files='decoupe.'//trim(adjustl(num))
       OPEN(FILE=files,UNIT=10)
       READ(10,*) nb
       PRINT *,'  > ',i,' :',nb
       DO j=1,nb
          IF (params%coord==1) THEN
             !decoupage par coordonnees
             coord(:,:)=0.
             READ(10,*) coord(1,:)
             CALL ecritpoint_gmsh(1,params%dim,coord,i,1)
          ELSE
             !decoupage d'image 1D
             READ (10,*) k
             CALL ecritpointimage_gmsh(1,params,i,k)
          ENDIF
       ENDDO
       CLOSE(10)
    ENDDO
    DEALLOCATE(coord)
    WRITE(1,*) '};'
    CLOSE(1)
    PRINT *,'-> decoupe.visu'
    RETURN
  END SUBROUTINE affectation_gmsh

  !***********************
  !ecriture des clusters avant regroupement
  SUBROUTINE sous_clusters_gmsh(params)
    IMPLICIT NONE
    TYPE(type_params) :: params
    INTEGER :: i,j,k,nb,ind,lenn
    CHARACTER*30 :: num,files
    DOUBLE PRECISION,DIMENSION(:,:),POINTER :: coord
    INTEGER,DIMENSION(:),POINTER :: corresp
    DO i=0,params%nbproc-1
       !nom du fichier
       !datas
       WRITE(num,*) i
       num=adjustl(num)
       lenn=len(trim(num))
       files='cluster.partiel.'//trim(num)
       OPEN(FILE=files,UNIT=10)
       !sortie
       files='cluster.partiel.'//num(1:lenn)//'.visu'
       OPEN(FILE=files,UNIT=1)
       WRITE(1,*) 'View "Clusters for part '//num(1:lenn)//'" {'
       READ(10,*) nb,k
       PRINT *,'  > ',i,' :',nb,' -> ',files
       ALLOCATE(coord(1,k))
       IF ((params%image==1).OR.(params%geom==1).OR.(params%seuil==1)) THEN
          !lecture des correspondances
          files='decoupe.'//trim(num)
          OPEN(FILE=files,UNIT=11)
          READ(11,*)
          ALLOCATE(corresp(nb))
          DO j=1,nb
             READ(11,*) corresp(j)
          ENDDO
          CLOSE(11)
       ENDIF
       DO j=1,nb
          IF (params%coord==1) THEN
             READ(10,*) coord(1,:),k
             CALL ecritpoint_gmsh(1,params%dim,coord,k,1)
          ELSE
             READ(10,*) k,ind
             CALL ecritpointimage_gmsh(1,params,ind,corresp(k))
          ENDIF
       ENDDO
       CLOSE(10); 
       WRITE(1,*) '};'
       CLOSE(1)
       DEALLOCATE(coord)
       IF ((params%image==1).OR.(params%geom==1).OR.(params%seuil==1)) THEN
          DEALLOCATE(corresp)
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE sous_clusters_gmsh

  !***********************
  !ecriture des clusters apres regroupement
  SUBROUTINE cluster_final_gmsh(params)
    IMPLICIT NONE
    TYPE(type_params) :: params
    INTEGER :: i,j,k,nb,nb0
    CHARACTER*30 :: num,files
    DOUBLE PRECISION,DIMENSION(:,:),POINTER :: coord
    OPEN(FILE=params%mesh,UNIT=1)
    IF (params%coord==1) THEN
       !lecture au format coordonnes
       READ(1,*) j,k
       PRINT *,'lecture du fichier de maillage...',j,k
       ALLOCATE(coord(j,k))
       nb0=0
       DO i=1,j
          READ(1,*,END=100) coord(i,:)
          nb0=nb0+1
       ENDDO
100    PRINT *,'nb de points ',nb0
       j=nb0
       CLOSE(1)
    ENDIF

    nb0=k !dim des points
    OPEN(FILE='cluster.final.visu',UNIT=1)
    WRITE(1,*) 'View "Clusters" {'
    !lecture des fichiers
    DO i=1,params%nbclusters
       !nom du fichier
       WRITE(num,*) i
       files='cluster.final.'//trim(adjustl(num))
       OPEN(FILE=files,UNIT=10)
       READ(10,*) nb
       PRINT *,'  > ',i,' :',nb
       DO j=1,nb
          READ(10,*) k
          IF (params%coord==1) THEN
             !coordonnees classiques
             CALL ecritpoint_gmsh(1,nb0,coord,i,k)
          ELSE
             !reassemblage des images
             CALL ecritpointimage_gmsh(1,params,i,k)
          ENDIF
       ENDDO
       CLOSE(10)
    ENDDO
    WRITE(1,*) '};'
    CLOSE(1)
    PRINT *,'-> cluster.final.visu'
    IF (params%coord==1) DEALLOCATE(coord)
    RETURN
  END SUBROUTINE cluster_final_gmsh

  !*************************
  !SUBROUTINE ecriture de point
  SUBROUTINE ecritpoint_gmsh(unit,dim,coord,ind,k)
    IMPLICIT NONE
    INTEGER :: unit,dim,ind,k
    DOUBLE PRECISION,DIMENSION(:,:),POINTER :: coord
    IF (dim==2) THEN
       !2D
       WRITE(unit,*) 'SP(',coord(k,1),',',coord(k,2),',',0.,'){',ind,'};' 
    ELSEIF (dim==3) THEN
       !3D
       WRITE(unit,*) 'SP(',coord(k,1),',',coord(k,2),',',coord(k,3),'){',ind,'};'
    ENDIF
    RETURN
  END SUBROUTINE ecritpoint_gmsh

  !*************************
  !SUBROUTINE ecriture de point en format image
  SUBROUTINE ecritpointimage_gmsh(unit,params,ind,k)
    IMPLICIT NONE
    TYPE(type_params) :: params
    INTEGER :: unit,ind,k,ix,iy,i
    DOUBLE PRECISION :: kx,ky,kz
    DOUBLE PRECISION,DIMENSION(:),POINTER :: data
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
    !coordonnees
    ix=params%refimg(k,1)
    iy=params%refimg(k,2)
    IF (params%imgdim==2) THEN
       !points en 2D
       IF (params%geom==1) THEN
          !coordonnees redimensionnees
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
       !points en 3D
       IF (params%geom==1) THEN
          kx=iy*params%pas(2)
          ky=-ix*params%pas(1)
          kz=float(params%refimg(k,3))*params%pas(3)
       ELSE
          !kx=-float(ix)
          !ky=float(iy)
          !kz=float(params%refimg(k,3))
          kx=float(iy)
          ky=-float(ix)
          kz=float(params%refimg(k,3))
       ENDIF
    ENDIF
    !ecriture
    WRITE(unit,*) 'SP(',kx,',',ky,',',kz,'){',ind,'};'
    RETURN
  END SUBROUTINE ecritpointimage_gmsh

  !************************
  !liste des commandes
  SUBROUTINE commandes_gmsh
    PRINT *,'gmsh decoupe.visu'
    PRINT *,'gmsh decoupe.geo'
    PRINT *,'gmsh decoupe.visu decoupe.geo'
    PRINT *,'gmsh cluster.partiel.*.visu'
    PRINT *,'gmsh cluster.final.visu'
    PRINT *,'gmsh decoupe.geo cluster.partiel.*.visu'
    PRINT *,'gmsh decoupe.geo cluster.final.visu'
    RETURN
  END SUBROUTINE commandes_gmsh


END MODULE module_visuclusters_gmsh
