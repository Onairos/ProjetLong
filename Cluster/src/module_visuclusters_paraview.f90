MODULE module_visuclusters_paraview
  USE module_visuclusters_structure
CONTAINS



  SUBROUTINE write_partionning_paraview(params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====  
    TYPE(type_params) :: params
    
    !#### Variables  ####
    CHARACTER (LEN=30) :: num
    DOUBLE PRECISION, DIMENSION(:), POINTER :: xmax
    DOUBLE PRECISION, DIMENSION(:), POINTER :: xmin
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ymax
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ymin
    DOUBLE PRECISION, DIMENSION(:), POINTER :: zmax
    DOUBLE PRECISION, DIMENSION(:), POINTER :: zmin
    DOUBLE PRECISION :: x0
    DOUBLE PRECISION :: x1
    DOUBLE PRECISION :: y0
    DOUBLE PRECISION :: y1
    DOUBLE PRECISION :: z0
    DOUBLE PRECISION :: z1
    INTEGER :: i
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################    
    ! Reading
    OPEN(FILE='fort.2',UNIT=2)
    ALLOCATE(xmin(params%nbproc))
    ALLOCATE(xmax(params%nbproc))
    ALLOCATE(ymin(params%nbproc))
    ALLOCATE(ymax(params%nbproc))
    ALLOCATE(zmin(params%nbproc))
    ALLOCATE(zmax(params%nbproc))
    DO i=1,params%nbproc-params%interface
       IF (((params%coord==1).AND.(params%dim==2)) &
          .OR.((params%image==1).AND.(params%imgdim==2)) &
          .OR.((params%seuil==1).AND.(params%imgdim==2)) &
          .OR.((params%geom==1).AND.(params%imgdim+params%imgt==2))) THEN
          !2D
          IF (params%coord==1) THEN
             READ(2,*) x0,y0,num,x1,y1
             xmin(i)=x0
             ymin(i)=y0
             xmax(i)=x1
             ymax(i)=y1
          ELSEIF (params%seuil==1) THEN
             READ(2,*) x0,num,x1
             xmin(i)=1
             ymin(i)=-1
             xmax(i)=params%imgmap(2)
             ymax(i)=-params%imgmap(1)
             zmin(i)=x0
             zmax(i)=x1
          ELSEIF ((params%image==1).OR.(params%geom==1)) THEN
             READ(2,*) x0,y0,num,x1,y1
             xmin(i)=y0
             ymin(i)=-x0
             xmax(i)=y1
             ymax(i)=-x1
             zmin(i)=x0
          ENDIF
       ELSEIF (((params%coord==1).AND.(params%dim==3)) &
          .OR.((params%image==1).AND.(params%imgdim==3)) &
          .OR.((params%geom==1).AND.(params%imgdim+params%imgt==3))) THEN
          !3D
          IF (params%coord==1) THEN
             READ(2,*) x0,y0,z0,num,x1,y1,z1
             xmin(i)=x0
             ymin(i)=y0
             zmin(i)=z0
             xmax(i)=x1
             ymax(i)=y1
             zmax(i)=z1
          ELSEIF ((params%image==1).OR.(params%geom==1)) THEN
             READ(2,*) x0,y0,z0,num,x1,y1,z1
             xmin(i)=y0
             ymin(i)=-x0
             zmin(i)=z0
             xmax(i)=y1
             ymax(i)=-x1
             zmax(i)=z1
          ENDIF
       ENDIF
    ENDDO
    CLOSE(2)
    ! Writing
    PRINT *,'-> visu/decoupe.geo'
    PRINT *,'-> visu/decoupe.indices'
    OPEN(FILE='visu/decoupe.geo',UNIT=10)
    OPEN(FILE='visu/decoupe.indices',UNIT=11)
    WRITE(10,*) '** sortie de visuclusters **'
    WRITE(10,*) '** decoupage des sous-clusters **'
    WRITE(10,'(a)') 'node id assign'
    WRITE(10,'(a)') 'element id assign'
    WRITE(11,*) '** indices des process **'
    WRITE(10,'(a)') 'part'
    WRITE(10,*) 1
    WRITE(10,*) '** decoupages **'
    WRITE(10,'(a)') 'coordinates'
    IF (((params%coord==1).AND.(params%dim==2)) &
         .OR.((params%image==1).AND.(params%imgdim==2)) &
         .OR.((params%seuil==1).AND.(params%imgdim==2)) &
         .OR.((params%geom==1).AND.(params%imgdim+params%imgt==2))) THEN
       !2D
       IF ((params%coord==1).OR.(params%geom==1).OR.(params%image==1)) THEN
          WRITE(10,*) 4*(params%nbproc-params%interface)
          DO i=1,params%nbproc-params%interface
             WRITE(10,*) xmin(i)
             WRITE(10,*) xmax(i)
             WRITE(10,*) xmax(i)
             WRITE(10,*) xmin(i)
          ENDDO
          DO i=1,params%nbproc-params%interface
             WRITE(10,*) ymin(i)
             WRITE(10,*) ymin(i)
             WRITE(10,*) ymax(i)
             WRITE(10,*) ymax(i)
          ENDDO
          DO i=1,params%nbproc-params%interface
             WRITE(10,*) 0.
             WRITE(10,*) 0.
             WRITE(10,*) 0.
             WRITE(10,*) 0.
          ENDDO
          WRITE(10,'(a)') 'quad4'
          WRITE(10,*) params%nbproc-params%interface
          WRITE(11,'(a)') 'part'
          WRITE(11,*) 1
          WRITE(11,'(a)') 'quad4'
          DO i=1,params%nbproc-params%interface
             WRITE(10,*) 4*(i-1)+1,4*(i-1)+2,4*(i-1)+3,4*(i-1)+4
             WRITE(11,*) i
          ENDDO
       ELSEIF (params%seuil==1) THEN
          WRITE(10,*) 8*(params%nbproc-params%interface)
          DO i=1,params%nbproc-params%interface
             WRITE(10,*) xmin(i)
             WRITE(10,*) xmax(i)
             WRITE(10,*) xmax(i)
             WRITE(10,*) xmin(i)
             WRITE(10,*) xmin(i)
             WRITE(10,*) xmax(i)
             WRITE(10,*) xmax(i)
             WRITE(10,*) xmin(i)
          ENDDO
          DO i=1,params%nbproc-params%interface
             WRITE(10,*) ymin(i)
             WRITE(10,*) ymin(i)
             WRITE(10,*) ymax(i)
             WRITE(10,*) ymax(i)
             WRITE(10,*) ymin(i)
             WRITE(10,*) ymin(i)
             WRITE(10,*) ymax(i)
             WRITE(10,*) ymax(i)
          ENDDO
          DO i=1,params%nbproc-params%interface
             WRITE(10,*) zmin(i)
             WRITE(10,*) zmin(i)
             WRITE(10,*) zmin(i)
             WRITE(10,*) zmin(i)
             WRITE(10,*) zmax(i)
             WRITE(10,*) zmax(i)
             WRITE(10,*) zmax(i)
             WRITE(10,*) zmax(i)
          ENDDO
          WRITE(10,'(a)') 'hexa8'
          WRITE(10,*) params%nbproc-params%interface
          WRITE(11,'(a)') 'part'
          WRITE(11,*) 1
          WRITE(11,'(a)') 'hexa8'
          DO i=1,params%nbproc-params%interface
             WRITE(10,*) 8*(i-1)+1,8*(i-1)+2,8*(i-1)+3,8*(i-1)+4,&
                  8*(i-1)+5,8*(i-1)+6,8*(i-1)+7,8*(i-1)+8
             WRITE(11,*) i
          ENDDO
       ENDIF
    ELSEIF (((params%coord==1).AND.(params%dim==3)) &
         .OR.((params%image==1).AND.(params%imgdim==3)) &
         .OR.((params%geom==1).AND.(params%imgdim+params%imgt==3))) THEN
       !3D
       IF ((params%coord==1).OR.(params%geom==1).OR.(params%image==1)) THEN
          WRITE(10,*) 8*(params%nbproc-params%interface)
          DO i=1,params%nbproc-params%interface
             WRITE(10,*) xmin(i)
             WRITE(10,*) xmax(i)
             WRITE(10,*) xmax(i)
             WRITE(10,*) xmin(i)
             WRITE(10,*) xmin(i)
             WRITE(10,*) xmax(i)
             WRITE(10,*) xmax(i)
             WRITE(10,*) xmin(i)
          ENDDO
          DO i=1,params%nbproc-params%interface
             WRITE(10,*) ymin(i)
             WRITE(10,*) ymin(i)
             WRITE(10,*) ymax(i)
             WRITE(10,*) ymax(i)
             WRITE(10,*) ymin(i)
             WRITE(10,*) ymin(i)
             WRITE(10,*) ymax(i)
             WRITE(10,*) ymax(i)
          ENDDO
          DO i=1,params%nbproc-params%interface
             WRITE(10,*) zmin(i)
             WRITE(10,*) zmin(i)
             WRITE(10,*) zmin(i)
             WRITE(10,*) zmin(i)
             WRITE(10,*) zmax(i)
             WRITE(10,*) zmax(i)
             WRITE(10,*) zmax(i)
             WRITE(10,*) zmax(i)
          ENDDO
          WRITE(10,'(a)') 'hexa8'
          WRITE(10,*) params%nbproc-params%interface
          WRITE(11,'(a)') 'part'
          WRITE(11,*) 1
          WRITE(11,'(a)') 'hexa8'
          DO i=1,params%nbproc-params%interface
             WRITE(10,*) 8*(i-1)+1,8*(i-1)+2,8*(i-1)+3,8*(i-1)+4,&
                  8*(i-1)+5,8*(i-1)+6,8*(i-1)+7,8*(i-1)+8
             WRITE(11,*) i
          ENDDO
       ENDIF
    ENDIF
    CLOSE(10)
    CLOSE(11)
    DEALLOCATE(xmin)
    DEALLOCATE(xmax)
    DEALLOCATE(ymin)
    DEALLOCATE(ymax)
    DEALLOCATE(zmin)
    DEALLOCATE(zmax)
    ! Master file
    PRINT *,'-> decoupe.CASE'
    OPEN(FILE='decoupe.CASE',UNIT=12)
    WRITE(12,'(a)') 'FORMAT'
    WRITE(12,'(a)') 'type: ensight gold'
    WRITE(12,*) 
    WRITE(12,'(a)') 'GEOMETRY'
    WRITE(12,'(a)') 'model:   visu/decoupe.geo'
    WRITE(12,*) 
    WRITE(12,'(a)') 'VARIABLE'
    WRITE(12,'(a)') 'scalar per element:   process   visu/decoupe.indices'
    CLOSE(12)
    RETURN
  END SUBROUTINE write_partionning_paraview



  SUBROUTINE affectation_paraview(params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    
    !#### Variables  ####
    CHARACTER (LEN=30) :: num
    CHARACTER (LEN=30) :: files
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coord
    INTEGER, DIMENSION(:), POINTER :: ind
    INTEGER, DIMENSION(:), POINTER :: indp
    INTEGER :: i
    INTEGER :: j
    INTEGER :: nb
    INTEGER :: offset
    INTEGER :: totnum
    
    !###########################################      
    ! INSTRUCTIONS
    !########################################### 
    ! Reading of files
    IF (params%nbproc==1) THEN
       offset=1
       totnum=1
    ELSE
       offset=0
       totnum=params%nbproc-1
    ENDIF
    DO i=offset,totnum
       IF ((i==0).AND.(params%interface==1).AND.(params%nbproc>1)) THEN
          ! File aside for interfacing 
          PRINT *,'-> visu/affectation-interface.geo'
          PRINT *,'-> visu/affectation-interface.indices'
          OPEN(FILE='visu/affectation-interface.geo',UNIT=10)
          OPEN(FILE='visu/affectation-interface.indices',UNIT=11)
          WRITE(10,*) '** sortie de visuclusters **'
          WRITE(10,*) '** points sur l interface **'
          WRITE(10,'(a)') 'node id assign'
          WRITE(10,'(a)') 'element id assign'
          WRITE(11,*) '** indices des process **'
       ELSEIF (((i==1).AND.(params%interface==1)).OR. &
            ((i==0).AND.(params%interface==0))) THEN
          ! General file for others subdomains
          PRINT *,'-> visu/affectation.geo'
          PRINT *,'-> visu/affectation.indices'
          OPEN(FILE='visu/affectation.geo',UNIT=10)
          OPEN(FILE='visu/affectation.indices',UNIT=11)
          WRITE(10,*) '** sortie de visuclusters **'
          WRITE(10,*) '** decoupage des sous-clusters **'
          WRITE(10,'(a)') 'node id assign'
          WRITE(10,'(a)') 'element id assign'
          WRITE(11,*) '** indices des process **'
       ENDIF
       ! File name
       WRITE(num,*) i
       files='decoupe.'//trim(adjustl(num))
       OPEN(FILE=files,UNIT=20)
       READ(20,*) nb
       PRINT *,'  > ',i,' :',nb
        ALLOCATE(coord(nb,params%dim))
       coord(:,:)=0.
       ALLOCATE(ind(nb))
       ind(:)=0
       ALLOCATE(indp(nb))
       indp(:)=0
       IF (nb>0) THEN
          DO j=1,nb
             IF (params%coord==1) THEN
                ! Partitionning by coordinates
                READ(20,*) coord(j,:)
                ind(j)=i
             ELSE
                ! Partitionning 1D picture
                READ (20,*) indp(j)
                ind(j)=i
             ENDIF
          ENDDO
          ! Writing
          IF (params%coord==1) THEN
             ! Partitionning by coordinates
             CALL ecritpoint_paraview(10,11,nb,params%dim,coord,ind,1)
          ELSE
             ! Partitionning 1D picture
             CALL write_picture_to_paraview(10,11,nb,params,ind,indp)
          ENDIF
       ENDIF
       DEALLOCATE(coord)
       DEALLOCATE(ind)
       CLOSE(20)
       IF ((i==0).AND.(params%interface==1)) THEN
          ! Closes interfacing files
          CLOSE(10)
          CLOSE(11)
       ENDIF
    ENDDO
    CLOSE(10)
    CLOSE(11)
    ! Writes interfacing data
    IF ((params%interface==1).AND.(params%nbproc>1)) THEN
       PRINT *,'-> affectation-interface.CASE'
       OPEN(FILE='affectation-interface.CASE',UNIT=12)
       WRITE(12,'(a)') 'FORMAT'
       WRITE(12,'(a)') 'type: ensight gold'
       WRITE(12,*) 
       WRITE(12,'(a)') 'GEOMETRY'
       WRITE(12,'(a)') 'model:   visu/affectation-interface.geo'
       WRITE(12,*) 
       WRITE(12,'(a)') 'VARIABLE'
       WRITE(12,'(a)') 'scalar per node:   process   visu/affectation-interface.indices'
       CLOSE(12)
    ENDIF
    ! writes others subdomains
    PRINT *,'-> affectation.CASE'
    OPEN(FILE='affectation.CASE',UNIT=12)
    WRITE(12,'(a)') 'FORMAT'
    WRITE(12,'(a)') 'type: ensight gold'
    WRITE(12,*) 
    WRITE(12,'(a)') 'GEOMETRY'
    WRITE(12,'(a)') 'model:   visu/affectation.geo'
    WRITE(12,*) 
    WRITE(12,'(a)') 'VARIABLE'
    WRITE(12,'(a)') 'scalar per node:   process   visu/affectation.indices'
    CLOSE(12)
    RETURN
  END SUBROUTINE affectation_paraview



  SUBROUTINE write_partial_clusters_paraview(params)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    
    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coord
    CHARACTER (LEN=30) :: files
    CHARACTER (LEN=30) :: num
    CHARACTER (LEN=30) :: star
    INTEGER, DIMENSION(:), POINTER :: corresp
    INTEGER, DIMENSION(:), POINTER :: ind
    INTEGER, DIMENSION(:), POINTER :: indp
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: lenn
    INTEGER :: nb
    INTEGER :: nbstar
    
    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
    ! Number of generic characters to be used
    nbstar=floor(log(real(params%nbproc-1))/log(real(10)))+1
    DO i=1,nbstar
       star(i:i)='0'
    ENDDO
    DO i=0,params%nbproc-1
       ! File name
       WRITE(num,*) i
       num=adjustl(num)
       lenn=len(trim(num))
       files='cluster.partiel.'//trim(num)
       OPEN(FILE=files,UNIT=20)
       ! Output
       ! General file for others subdomains
       files='cluster.partiel.'//star(1:nbstar-lenn)//num(1:lenn)
       PRINT *,'-> visu/'//trim(files)//'.geo'
       PRINT *,'-> visu/'//trim(files)//'.indices'
       OPEN(FILE='visu/'//trim(files)//'.geo',UNIT=10)
       OPEN(FILE='visu/'//trim(files)//'.indices',UNIT=11)
       WRITE(10,*) '** sortie de visuclusters **'
       WRITE(10,*) '** decoupage du sous-clusters '//trim(num)//' **'
       WRITE(10,'(a)') 'node id assign'
       WRITE(10,'(a)') 'element id assign'
       WRITE(11,*) '** indices des elements clusterises **'
       ! Reading
       READ(20,*) nb,k
       PRINT *,'  > ',i,' :',nb,' -> ',files
       ALLOCATE(coord(nb,k))
       ALLOCATE(ind(max(1,nb)))
       ind(:)=0
       ALLOCATE(indp(max(1,nb)))
       indp(:)=0
       IF ((params%image==1).OR.(params%geom==1).OR.(params%seuil==1)) THEN
          ! Reading the matchings
          files='decoupe.'//trim(num)
          OPEN(FILE=files,UNIT=21)
          READ(21,*)
          ALLOCATE(corresp(nb))
          DO j=1,nb
             READ(21,*) corresp(j)
          ENDDO
          CLOSE(21)
       ENDIF
       DO j=1,nb
          IF (params%coord==1) THEN
             READ(20,*) coord(j,:),ind(j)
          ELSE
             READ(20,*) indp(j),ind(j)
             indp(j)=corresp(indp(j))
          ENDIF
       ENDDO
       CLOSE(20) 
       ! Writing
       IF (params%coord==1) THEN
          ! partitionning by coordinates
          CALL ecritpoint_paraview(10,11,nb,params%dim,coord,ind,1)
       ELSE
          ! partitionning 1D picture
          CALL write_picture_to_paraview(10,11,nb,params,ind,indp)
       ENDIF
       DEALLOCATE(coord)
       DEALLOCATE(ind)
       DEALLOCATE(indp)
       CLOSE(10)
       CLOSE(11)
       IF ((params%image==1).OR.(params%geom==1).OR.(params%seuil==1)) THEN
          DEALLOCATE(corresp)
       ENDIF
    ENDDO
    PRINT *,'-> cluster.partiel.CASE'
    DO i=1,nbstar
       star(i:i)='*'
    ENDDO
    OPEN(FILE='cluster.partiel.CASE',UNIT=12)
    WRITE(12,'(a)') 'FORMAT'
    WRITE(12,'(a)') 'type: ensight gold'
    WRITE(12,*) 
    WRITE(12,'(a)') 'GEOMETRY'
    WRITE(12,'(a)') 'model: 1 visu/cluster.partiel.'//star(1:nbstar)//'.geo'
    WRITE(12,*) 
    WRITE(12,'(a)') 'VARIABLE'
    WRITE(12,'(a)') 'scalar per node: 1 cluster visu/cluster.partiel.'//&
         star(1:nbstar)//'.indices'
    WRITE(12,*) 
    WRITE(12,'(a)') 'TIME'
    WRITE(12,'(a)') 'time set:         1'
    WRITE(num,*) params%nbproc
    WRITE(12,'(a)') 'number of steps: '//trim(adjustl(num))
    WRITE(12,'(a)') 'filename start number: 0'
    WRITE(12,'(a)') 'filename increment: 1'
    WRITE(12,'(a)') 'time values:'
    DO i=0,params%nbproc-1
       WRITE(12,*) i
    ENDDO
    CLOSE(12)
    RETURN
  END SUBROUTINE write_partial_clusters_paraview



  SUBROUTINE write_final_clusters_paraview(params)
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
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coord
    INTEGER, DIMENSION(:), POINTER :: ind
    INTEGER, DIMENSION(:), POINTER :: indp
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: nb
    INTEGER :: nb0
    
    !###########################################
    ! INSTRUCTIONS
    !###########################################
    IF (params%coord==1) THEN
       OPEN(FILE=params%mesh,UNIT=1)
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

    nb0=k ! Points dimension
    ! Output
    PRINT *,'-> visu/cluster.final.geo'
    PRINT *,'-> visu/cluster.final.indices'
    OPEN(FILE='visu/cluster.final.geo',UNIT=10)
    OPEN(FILE='visu/cluster.final.indices',UNIT=11)
    WRITE(10,*) '** sortie de visuclusters **'
    WRITE(10,*) '** decoupage cluster final **'
    WRITE(10,'(a)') 'node id assign'
    WRITE(10,'(a)') 'element id assign'
    WRITE(11,*) '** clusters **'
    ALLOCATE(ind(max(1,params%nbp)))
    ind(:)=0
    ALLOCATE(indp(max(1,params%nbp)))
    indp(:)=0
    ! Reading files
    DO i=1,params%nbclusters
       ! File name
       WRITE(num,*) i
       files='cluster.final.'//trim(adjustl(num))
       OPEN(FILE=files,UNIT=20)
       READ(20,*) nb
       PRINT *,'  > ',i,' :',nb
       DO j=1,nb
          READ(20,*) k
          ind(k)=i
          indp(k)=k
       ENDDO
       CLOSE(20)
    ENDDO
    IF (params%coord==1) THEN
       ! classical coordinates
       CALL ecritpoint_paraview(10,11,params%nbp,params%dim,coord,ind,1)
    ELSE
       ! Pictures reassembly
       CALL write_picture_to_paraview(10,11,params%nbp,params,ind,indp)
    ENDIF
    CLOSE(10)
    CLOSE(11)
    DEALLOCATE(ind)
    DEALLOCATE(indp)
    IF (params%coord==1) DEALLOCATE(coord)
    ! Master file
    PRINT *,'-> cluster.final.CASE'
    OPEN(FILE='cluster.final.CASE',UNIT=12)
    WRITE(12,'(a)') 'FORMAT'
    WRITE(12,'(a)') 'type: ensight gold'
    WRITE(12,*) 
    WRITE(12,'(a)') 'GEOMETRY'
    WRITE(12,'(a)') 'model:   visu/cluster.final.geo'
    WRITE(12,*) 
    WRITE(12,'(a)') 'VARIABLE'
    WRITE(12,'(a)') 'scalar per element:   cluster   visu/cluster.final.indices'
    CLOSE(12)
    RETURN
  END SUBROUTINE write_final_clusters_paraview



  SUBROUTINE ecritpoint_paraview(unitgeo, unitind, nb, dim, coord, ind, k)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: coord
    INTEGER, DIMENSION(:), POINTER :: ind
    INTEGER :: unitgeo
    INTEGER :: unitind
    INTEGER :: k
    INTEGER :: nb
    INTEGER :: i
    INTEGER :: dim
    
    !###########################################
    ! INSTRUCTIONS
    !###########################################
    WRITE(unitgeo,'(a)') 'part'
    WRITE(unitgeo,*) ind(k)
    WRITE(unitgeo,*) '** decoupages **'
    WRITE(unitgeo,'(a)') 'coordinates'
    WRITE(unitgeo,*) nb
    WRITE(unitind,'(a)') 'part'
    WRITE(unitind,*) ind(k)
    WRITE(unitind,'(a)') 'point'
    DO i=1,nb
       WRITE(unitgeo,*) coord(i,1)
       WRITE(unitind,*) ind(i)
    ENDDO
    DO i=1,nb
       WRITE(unitgeo,*) coord(i,2)
    ENDDO
    DO i=1,nb
       IF (dim==2) THEN
          WRITE(unitgeo,*) 0.
       ELSE
          WRITE(unitgeo,*) coord(i,3)
       ENDIF
    ENDDO
    WRITE(unitgeo,'(a)') 'point'
    WRITE(unitgeo,*)nb
    DO i=1,nb
       WRITE(unitgeo,*) i
    ENDDO
    RETURN
  END SUBROUTINE ecritpoint_paraview



  SUBROUTINE write_picture_to_paraview(unitgeo, unitind, nbp, params, ind, indp)
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_params) :: params
    INTEGER, DIMENSION(:), POINTER :: ind
    INTEGER, DIMENSION(:), POINTER :: indp
    INTEGER :: nbp
    INTEGER :: unitgeo
    INTEGER :: unitind
    !=== IN/OUT ===
    !====  OUT ====
    
    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:), POINTER :: data
    DOUBLE PRECISION, DIMENSION(:), POINTER :: kx
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ky
    DOUBLE PRECISION, DIMENSION(:), POINTER :: kz
    INTEGER :: i
    INTEGER :: ix
    INTEGER :: iy
    INTEGER :: k
    
    !###########################################
    ! INSTRUCTIONS
    !###########################################
    ALLOCATE(kx(nbp))
    kx(:)=0
    ALLOCATE(ky(nbp))
    ky(:)=0
    ALLOCATE(kz(nbp))
    kz(:)=0
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
    ! Search the points
    DO i=1,nbp
       k=indp(i)
       ix=params%refimg(k,1)
       iy=params%refimg(k,2)
       IF (params%imgdim==2) THEN
          ! 2D points
          IF (params%geom==1) THEN
             kx(i)=iy*params%pas(2)
             ky(i)=-ix*params%pas(1)
          ELSE
             kx(i)=float(iy)
             ky(i)=-float(ix)
          ENDIF
          IF (((params%image==1).OR.(params%geom==1).OR.(params%seuil==1)) &
               .AND.(params%imgdim==2)) THEN
             kz(i)=data(k)
          ELSE
             kz(i)=0.
          ENDIF
       ELSEIF (params%imgdim==3) THEN
          ! 3D points
          IF (params%geom==1) THEN
             kx(i)=iy*params%pas(2)
             ky(i)=-ix*params%pas(1)
             kz(i)=float(params%refimg(k,3))*params%pas(3)
          ELSE

             kx(i)=float(iy)
             ky(i)=-float(ix)
             kz(i)=float(params%refimg(k,3))
          ENDIF
       ENDIF
    ENDDO
    ! Writing
    WRITE(unitgeo,'(a)') 'part'
    WRITE(unitgeo,*) ind(1)
    WRITE(unitgeo,*) '** decoupages **'
    WRITE(unitgeo,'(a)') 'coordinates'
    WRITE(unitgeo,*) nbp
    WRITE(unitind,'(a)') 'part'
    WRITE(unitind,*) ind(1)
    WRITE(unitind,'(a)') 'point'
    DO i=1,nbp
       WRITE(unitgeo,*) kx(i)
       WRITE(unitind,*) ind(i)
    ENDDO
    DO i=1,nbp
       WRITE(unitgeo,*) ky(i)
    ENDDO
    DO i=1,nbp
       WRITE(unitgeo,*) kz(i)
    ENDDO
    WRITE(unitgeo,'(a)') 'point'
    WRITE(unitgeo,*)nbp
    DO i=1,nbp
       WRITE(unitgeo,*) i
    ENDDO
    DEALLOCATE(kx)
    DEALLOCATE(ky)
    DEALLOCATE(kz)
    IF (((params%image==1).OR.(params%geom==1).OR.(params%seuil==1)) &
         .AND.(params%imgdim==2)) DEALLOCATE(data)
    RETURN
  END SUBROUTINE write_picture_to_paraview



  SUBROUTINE list_commands_paraview
    !###########################################
    ! INSTRUCTIONS
    !###########################################
    PRINT *,'paraview --data=decoupe.CASE'
    PRINT *,'paraview --data=affectation.CASE'
    PRINT *,'paraview --data=affectation-interface.CASE'
    PRINT *,'paraview --data=cluster.partiel.CASE'
    PRINT *,'paraview --data=cluster.final.CASE'
    RETURN
  END SUBROUTINE list_commands_paraview
  
END MODULE module_visuclusters_paraview
