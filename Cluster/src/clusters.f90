PROGRAM clusters
  USE module_structure
  USE module_entree
  USE module_decoupe
  USE module_MPI
  USE module_calcul
  USE module_sortie

  IMPLICIT NONE
 
  ! librairie MPI
  INCLUDE 'mpif.h'
  !###########################################
  ! DECLARATIONS
  !###########################################
  !#### Variables  ####
  TYPE(type_clusters), DIMENSION(:), POINTER :: nclust
  TYPE(type_data) :: data
  TYPE(type_data) :: dataw
  CHARACTER*80 :: procname ! variable MPI
  CHARACTER*30 :: entree
  CHARACTER*30 :: mesh
  DOUBLE PRECISION, DIMENSION(:,:,:), POINTER :: bornes
  DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmax
  DOUBLE PRECISION, DIMENSION(:), POINTER :: coordmin
  DOUBLE PRECISION :: elapsed(2) ! For receiving user and system time
  DOUBLE PRECISION :: endtime
  DOUBLE PRECISION :: epsilon
  DOUBLE PRECISION :: sigma
  DOUBLE PRECISION :: starttime
  DOUBLE PRECISION :: t_parall
  DOUBLE PRECISION :: t_parallg
  DOUBLE PRECISION :: t1
  DOUBLE PRECISION :: t2
  DOUBLE PRECISION :: temps
  INTEGER,DIMENSION(:,:), POINTER :: clustermap
  INTEGER,DIMENSION(:,:) ,POINTER :: ddat
  INTEGER,DIMENSION(:), POINTER :: decoupe
  INTEGER,DIMENSION(:), POINTER :: iclust
  INTEGER,DIMENSION(:), POINTER :: ldat
  INTEGER,DIMENSION(:), POINTER :: listenbideal
  INTEGER :: i
  INTEGER :: ierr ! variable MPI
  INTEGER :: j
  INTEGER :: len ! variable MPI
  INTEGER :: nbclust
  INTEGER :: nbideal
  INTEGER :: nblimit
  INTEGER :: nbproc ! variable MPI
  INTEGER :: nmax
  INTEGER :: numproc ! variable MPI
  INTEGER :: status(MPI_STATUS_SIZE) ! variable MPI
  INTEGER :: tag ! variable MPI
  INTEGER :: test
  LOGICAL :: existe

  !###########################################
  ! INSTRUCTIONS
  !###########################################
  !initialisation MPI
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,numproc,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nbproc,ierr)
  CALL MPI_GET_PROCESSOR_NAME(procname,len,ierr)
  PRINT*, 'rank', numproc, ' procname',procname

  IF(numproc==0) THEN
    starttime = MPI_WTIME();
  ENDIF
#if aff
  PRINT *,'lancement du proc',numproc,' de ',nbproc
#endif
  PRINT *,'----------------------------------------'

  !lecture des datas
  IF (numproc==0) THEN
     t1 = MPI_WTIME()
#if aff
     PRINT *
     PRINT *,'------------------------------'
     PRINT *,'spectral clustering de datas'
     PRINT *,'------------------------------'
     PRINT *
#endif
     IF (iargc()>0) THEN
        !recupere le nom du fichier d'entree
        CALL getarg(1,entree)
        INQUIRE(FILE=entree,EXIST=existe)
        IF (.NOT.existe) CALL help
     ELSE
        !si pas fichier d'entree defaut=fort.1
        CALL help
     ENDIF
     !lecture du fichier
#if aff
     PRINT *
     PRINT *,'lecture des data... ',entree
#endif
     OPEN(FILE=entree,UNIT=1)
     CALL read_file(data,epsilon,coordmin,coordmax,nbproc,decoupe,mesh,&
          sigma,nblimit,listenbideal)
     t2 = MPI_WTIME()
     PRINT *, 'temps lecture data ', t2-t1
     CLOSE(1)
  ENDIF

  !decoupage & calcul du sigma
  IF (numproc==0) THEN
     !decoupage
#if aff
     PRINT *
     PRINT *,'decoupage des datas...'
#endif
    t1 = MPI_WTIME();
    CALL partition_data(data,epsilon,nbproc,coordmin,coordmax,decoupe,&
         ldat,ddat,bornes)
    t2 = MPI_WTIME();
    PRINT *,'temps decoupage des datas...', t2-t1
  ENDIF
     
  IF(numproc==0) THEN
    t1 = MPI_WTIME();
  ENDIF

  !partie echanges
  IF (nbproc>1) THEN
     !** cas >1 proc

     !calcul du sigma si mode auto global
     CALL MPI_BCAST(sigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     IF ((sigma==0.0).AND.(numproc==0)) THEN
        CALL get_sigma_interface(numproc,data,sigma,bornes,decoupe,epsilon)
     ENDIF

     !envoi du sigma
     IF (sigma>=0.0) THEN
          CALL MPI_BCAST(sigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     ENDIF

     !envoi du nblimit
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     CALL MPI_BCAST(nblimit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

     !envoi du nbideal
     IF (numproc==0) THEN
        DO i=1,nbproc-1
           tag=i
           CALL MPI_SEND(listenbideal(i),1,MPI_INTEGER,i,tag,MPI_COMM_WORLD,ierr)
        ENDDO
        nbideal=listenbideal(0)
     ELSE
        tag=numproc
        CALL MPI_RECV(nbideal,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
     ENDIF
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

     !transferts des datas
     IF (numproc==0) THEN
        !envoi des datas
#if aff
        PRINT *
        PRINT *,'transfert des datas decoupees...'
#endif
        CALL send_partitionning(nbproc,data,ldat,ddat,dataw)
#if aff
        PRINT *
        PRINT *,'calcul des clusters...'
#endif
     ELSE
        !reception des datas
        CALL receive_partitionning(numproc,dataw)
     ENDIF

  ELSE
     !** cas 1 proc
     dataw%nb=data%nb; dataw%dim=data%dim; dataw%nbclusters=0
     ALLOCATE(dataw%point(data%nb))
     DO i=1,data%nb
        ALLOCATE(dataw%point(i)%coord(data%dim))
        dataw%point(i)%coord=data%point(i)%coord
        dataw%point(i)%cluster=0
     ENDDO
     nbideal=listenbideal(1)
  ENDIF

  !calcul du sigma si mode auto individuel
  IF (sigma<0.0) THEN
#if aff
     PRINT *
#endif
     IF (numproc==0) CALL get_sigma(data,sigma)
     CALL MPI_BCAST(sigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     IF (numproc==0) THEN
        PRINT *,numproc,'calcule le sigma global :',sigma
#if aff
        PRINT *,numproc,'calcule le sigma global :',sigma
#endif
        IF (data%interface==1) THEN 
           CALL get_sigma_interface(numproc,dataw,sigma,bornes,decoupe,epsilon) 
           PRINT *,numproc,'calcule le sigma interface :',sigma
#if aff
           PRINT *,numproc,'calcule le sigma interface :',sigma
#endif
        ENDIF
     ENDIF
  ENDIF

  IF(numproc==0) THEN
    t2 = MPI_WTIME();
    PRINT *,'temps envoi donnees et calcul sigma', t2-t1
  ENDIF

  !calcul des clusters
  ! partie //
  t1 = MPI_WTIME();
  IF (dataw%nb>0) THEN
#if aff
     PRINT *,numproc,'calcul des clusters...'
#endif
     CALL apply_spectral_clustering(numproc,nblimit,nbideal,dataw,sigma)
  ENDIF
  t2 = MPI_WTIME();
  t_parall = t2 - t1
  PRINT *, numproc, 'calcul cluster //', t_parall

  CALL MPI_REDUCE(t_parall, t_parallg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  IF(numproc==0) PRINT *, 'temps calcul cluster global //', t_parallg

  IF(numproc==0) THEN
    t1 = MPI_WTIME();
  ENDIF

  !sauvegarde des clusters regroupes
  CALL write_partial_clusters(numproc,dataw)

  !partie echanges
  IF (nbproc>1) THEN
     !regroupement des clusters
     IF (numproc==0) THEN
#if aff
        PRINT *
        PRINT *,'regroupement des clusters...'
#endif
        !creation de la liste des clusters avec doublon
        CALL receive_number_clusters(nbproc,nbclust,ldat,dataw,nclust)
#if aff
        PRINT *,'  > nb de clusters avec doublons obtenus :',nbclust
#endif
        !reception des infos de clusters
        ALLOCATE(clustermap(nbclust,data%nb))
        CALL receive_clusters(nbproc,nbclust,ldat,ddat,dataw,&
             clustermap,nclust,iclust)
     ELSE
        !creation de la liste des clusters avec doublon
        CALL send_number_clusters(numproc,dataw)
        !envoi des infos de clusters
        CALL send_clusters(numproc,dataw)
     ENDIF

     !fin du postprocess
     IF (numproc==0) THEN
        !regroupement des clusters et ecriture du resultat
        CALL group_clusters(nbclust,iclust,clustermap,data)
     ENDIF

  ELSE
     !** cas 1 proc
     nbclust=dataw%nbclusters
     ALLOCATE(iclust(nbclust)); iclust(:)=0; nmax=0
     DO i=1,dataw%nb
        j=dataw%point(i)%cluster
        iclust(j)=iclust(j)+1
        nmax=max(nmax,iclust(j))
     ENDDO
     ALLOCATE(clustermap(nbclust,nmax)); clustermap(:,:)=0
     iclust(:)=0
     DO i=1,dataw%nb
        j=dataw%point(i)%cluster
        iclust(j)=iclust(j)+1
        clustermap(j,iclust(j))=i
     ENDDO
  ENDIF

  IF(numproc==0) THEN
    t2 = MPI_WTIME();
    PRINT *,'temps regroupement des clusters', t2-t1
  ENDIF

  IF(numproc==0) THEN
    t1 = MPI_WTIME();
  ENDIF
  !sorties
  IF (numproc==0) THEN
     !ecriture des cluster.final.
     CALL write_final_clusters(nbclust,iclust,clustermap)

     !ecriture des informations
     CALL write_metadata(mesh,data,nbproc,nbclust)
  ENDIF
  IF(numproc==0) THEN
    t2 = MPI_WTIME();
    PRINT *,'temps ecriture des clusters', t2-t1
  ENDIF
  
  !fin du MPI
  IF (numproc==0) THEN
    endtime = MPI_WTIME();
    PRINT *,' Fin du calcul'
    PRINT *,'  > temps total=', endtime-starttime
  ENDIF
  
  CALL MPI_FINALIZE(ierr)
  STOP

END PROGRAM clusters
