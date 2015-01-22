program clusters
  use module_structure
  use module_entree
  use module_decoupe
  use module_MPI
  use module_calcul
  use module_sparse
  use module_sortie

  implicit none
 
  ! librairie MPI
  include 'mpif.h'

  ! variables MPI
  integer :: ierr,tag
  integer :: numproc,nbproc, len
  character*80 :: procname
  
  integer status(MPI_STATUS_SIZE)

 ! modif sandrine
  ! variables
  real :: elapsed(2)     ! For receiving user and system time
  real :: temps
  double precision :: epsilon
  type(type_data) :: data,dataw
  double precision,dimension(:),pointer :: coordmax,coordmin
  integer,dimension(:),pointer :: decoupe
  integer,dimension(:),pointer :: ldat,iclust,listenbideal
  type(type_clusters),dimension(:),pointer :: nclust
  integer,dimension(:,:),pointer :: ddat
  integer :: nblimit,nbideal
  logical :: existe
  integer :: test,nbclust,i,j,nmax
  double precision :: sigma
  integer,dimension(:,:),pointer :: clustermap
  real*8,dimension(:,:,:),pointer :: bornes
  character*30 :: mesh,entree
 
  double precision :: starttime, endtime
  double precision :: t1, t2
  double precision :: t_parall, t_parallg

  !initialisation MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,numproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nbproc,ierr)
  call MPI_GET_PROCESSOR_NAME(procname,len,ierr)
  print*, 'rank', numproc, ' procname',procname

  if(numproc==0) then
    starttime = MPI_WTIME();
  end if
#if aff
  print *,'lancement du proc',numproc,' de ',nbproc
#endif
  print *,'----------------------------------------'

  !lecture des datas
  if (numproc==0) then
     t1 = MPI_WTIME()
#if aff
     print *
     print *,'------------------------------'
     print *,'spectral clustering de datas'
     print *,'------------------------------'
     print *
#endif
     if (iargc()>0) then
        !recupere le nom du fichier d'entree
        call getarg(1,entree)
        inquire(file=entree,exist=existe)
        if (.not.existe) call help
     else
        !si pas fichier d'entree defaut=fort.1
        call help
     end if
     !lecture du fichier
#if aff
     print *
     print *,'lecture des data... ',entree
#endif
     open(file=entree,unit=1)
     call lit(data,epsilon,coordmin,coordmax,nbproc,decoupe,mesh,&
          sigma,nblimit,listenbideal)
     t2 = MPI_WTIME()
     print *, 'temps lecture data ', t2-t1
     close(1)
  endif

  !decoupage & calcul du sigma
  if (numproc==0) then
     !decoupage
#if aff
     print *
     print *,'decoupage des datas...'
#endif
    t1 = MPI_WTIME();
    call decoupedata(data,epsilon,nbproc,coordmin,coordmax,decoupe,&
         ldat,ddat,bornes)
    t2 = MPI_WTIME();
    print *,'temps decoupage des datas...', t2-t1
  end if
     
  if(numproc==0) then
    t1 = MPI_WTIME();
  end if

  !CALL MPI_FINALIZE(ierr)
  !STOP
  !partie echanges
  if (nbproc>1) then
     !** cas >1 proc

     !calcul du sigma si mode auto global
     call MPI_BCAST(sigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     if ((sigma==0.0).and.(numproc==0)) then
        call calculsigmainterface(numproc,data,sigma,bornes,decoupe,epsilon)
     end if

     !envoi du sigma
     if (sigma>=0.0) then
          call MPI_BCAST(sigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     end if

     !envoi du nblimit
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     call MPI_BCAST(nblimit,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

     !envoi du nbideal
     if (numproc==0) then
        do i=1,nbproc-1
           tag=i
           call MPI_SEND(listenbideal(i),1,MPI_INTEGER,i,tag,MPI_COMM_WORLD,ierr)
        end do
        nbideal=listenbideal(0)
     else
        tag=numproc
        call MPI_RECV(nbideal,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status,ierr)
     endif
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     !transferts des datas
     if (numproc==0) then
        !envoi des datas
#if aff
        print *
        print *,'transfert des datas decoupees...'
#endif
        call  envoidecoupes(nbproc,data,ldat,ddat,dataw)
#if aff
        print *
        print *,'calcul des clusters...'
#endif
     else
        !reception des datas
        call recoitdecoupes(numproc,dataw)
     endif

  else
     !** cas 1 proc
     dataw%nb=data%nb; dataw%dim=data%dim; dataw%nbclusters=0
     allocate(dataw%point(data%nb))
     do i=1,data%nb
        allocate(dataw%point(i)%coord(data%dim))
        dataw%point(i)%coord=data%point(i)%coord
        dataw%point(i)%cluster=0
     enddo
     nbideal=listenbideal(1)
  endif

  !calcul du sigma si mode auto individuel
  if (sigma<0.0) then
#if aff
     print *
#endif
     if (numproc==0) call calculsigma(data,sigma)
     call MPI_BCAST(sigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     if (numproc==0) then
        print *,numproc,'calcule le sigma global :',sigma
#if aff
        print *,numproc,'calcule le sigma global :',sigma
#endif
        if (data%interface==1) then 
           call calculsigmainterface(numproc,dataw,sigma,bornes,decoupe,epsilon) 
           print *,numproc,'calcule le sigma interface :',sigma
#if aff
           print *,numproc,'calcule le sigma interface :',sigma
#endif
        end if
     end if
  end if

  if(numproc==0) then
    t2 = MPI_WTIME();
    print *,'temps envoi données et calcul sigma', t2-t1
  end if

  !calcul des clusters
  ! partie //
  t1 = MPI_WTIME();
  if (dataw%nb>0) then
#if aff
     print *,numproc,'calcul des clusters...'
#endif
     call calculclusters(numproc,nblimit,nbideal,dataw,sigma)
     !call sp_calculclusters(numproc,nblimit,nbideal,dataw,sigma)
  endif
  t2 = MPI_WTIME();
  t_parall = t2 - t1
  print *, numproc, 'calcul cluster //', t_parall

  call MPI_REDUCE(t_parall, t_parallg, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  if(numproc==0) print *, 'temps calcul cluster global //', t_parallg

  if(numproc==0) then
    t1 = MPI_WTIME();
  end if

  !sauvegarde des clusters regroupes
  call ecritcluster(numproc,dataw)

  !partie echanges
  if (nbproc>1) then
     !regroupement des clusters
     if (numproc==0) then
#if aff
        print *
        print *,'regroupement des clusters...'
#endif
        !creation de la liste des clusters avec doublon
        call preparecpclusters(nbproc,nbclust,ldat,dataw,nclust)
#if aff
        print *,'  > nb de clusters avec doublons obtenus :',nbclust
#endif
        !reception des infos de clusters
        allocate(clustermap(nbclust,data%nb))
        call recepclusters(nbproc,nbclust,ldat,ddat,dataw,&
             clustermap,nclust,iclust)
     else
        !creation de la liste des clusters avec doublon
        call prepaenvclusters(numproc,dataw)
        !envoi des infos de clusters
        call envoiclusters(numproc,dataw)
     endif

     !fin du postprocess
     if (numproc==0) then
        !regroupement des clusters et ecriture du resultat
        call regroupe(nbclust,iclust,clustermap,data)
     endif

  else
     !** cas 1 proc
     nbclust=dataw%nbclusters
     allocate(iclust(nbclust)); iclust(:)=0; nmax=0
     do i=1,dataw%nb
        j=dataw%point(i)%cluster
        iclust(j)=iclust(j)+1
        nmax=max(nmax,iclust(j))
     enddo
     allocate(clustermap(nbclust,nmax)); clustermap(:,:)=0
     iclust(:)=0
     do i=1,dataw%nb
        j=dataw%point(i)%cluster
        iclust(j)=iclust(j)+1
        clustermap(j,iclust(j))=i
     enddo
  endif

  if(numproc==0) then
    t2 = MPI_WTIME();
    print *,'temps regroupement des clusters', t2-t1
  end if

  if(numproc==0) then
    t1 = MPI_WTIME();
  end if
  !sorties
  if (numproc==0) then
     !ecriture des cluster.final.
     call ecritclusterfinal(nbclust,iclust,clustermap)

     !ecriture des informations
     call ecrit_info(mesh,data,nbproc,nbclust)
  endif
  if(numproc==0) then
    t2 = MPI_WTIME();
    print *,'temps écriture des clusters', t2-t1
  end if
  
  !fin du MPI
  if (numproc==0) then
    endtime = MPI_WTIME();
    print *,' Fin du calcul'
    print *,'  > temps total=', endtime-starttime
  end if
  
  call MPI_FINALIZE(ierr)
  stop

end program clusters
