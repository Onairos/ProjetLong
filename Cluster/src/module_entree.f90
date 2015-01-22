module module_entree
 use module_structure
contains

  !*********************************************
  !fichier d'aide
  subroutine help
    implicit none
    !integer :: ierr
    print *
    print *,'syntaxe d appel : clusters fichier_d_entree'
    print *
    print *,'mots cles du fichier d entree :'
    print *
    print *,'DATA'
    print *,'COORD (si donnees par coordonnees)'
    print *,'IMAGE (si fichier de maillage sous forme image + decoupage par pixel)'
    print *,'GEOM  (si fichier de maillage sous forme image + decoupage geom)'
    print *,'SEUIL  (si fichier de maillage sous forme image + decoupage par seuil)'
    print *,'  fichier_de_maillage'
    print *
    print *,'EPAISSEUR'
    print *,'  epaisseur_de_la_tranche'
    print *
    print *,'NBLIMIT'
    print *,'  nb_max_de_clusters'
    print *
    print *,'NBCLUST'
    print *,'  (facultatif)'
    print *,'  nb_de_clusters_par_sous-domaine'
    print *
    print *,'SIGMA'
    print *,'  (facultatif)'
    print *,'  valeur_de_sigma_imposee'
    print *
    print *,'DECOUPAGE'
    print *,'INTERFACE (decoupage + interface) '
    print *,'RECOUVREMENT (decoupage avec recouvrement)'
    print *,'  nb_de_sous-domaines_par_dimension'
    print *
    print *,'END'
    print *,'  (fin du fichier d entree)'
    stop
    return
  end subroutine help

  !**********************************************
  !lecture du fichier d'entree
  subroutine lit(data,epsilon,coordmin,coordmax,nbproc,decoupe,&
       mesh,sigma,nblimit,listenbideal)
    implicit none
    integer :: nbproc,nblimit,decoupage
    integer,dimension(:),pointer :: decoupe,listenbideal
    type(type_data) :: data
    double precision :: epsilon,sigma
    double precision,dimension(:),pointer :: coordmax,coordmin
    character*30 :: mot,mesh
    integer :: ok,i,ierr,tot
    !initialisation
    epsilon=0.0
    sigma=-1.0
    data%coord=0
    data%image=0
    data%geom=0
    data%seuil=0
    data%interface=0
    data%recouvrement=0
    nblimit=4
    decoupage=0
    if (nbproc>1) then
       allocate(listenbideal(0:nbproc-1))
    else
       allocate(listenbideal(1))
    end if
    listenbideal(:)=0
    !lecture
    ok=0
    do while (ok/=1)
       ok=1
       read(1,*) mot
       print *,mot
       select case(mot)
       case('DATA')
          ok=0
          read(1,*) mesh
          if (mesh=='IMAGE') then
             data%image=1
             read (1,*) mesh
             print *,'  > format d entree image + decoupage par pixel'
             print *,'  > lecture du fichier de data : ',mesh
             call lit_mesh_image(mesh,data,coordmin,coordmax)
          elseif (mesh=='GEOM') then
             data%geom=1
             read (1,*) mesh
             print *,'  > format d entree image + decoupage geometrique'
             print *,'  > lecture du fichier de data : ',mesh
             call lit_mesh_geom(mesh,data,coordmin,coordmax)
          elseif (mesh=='SEUIL') then
             data%seuil=1
             read (1,*) mesh
             print *,'  > format d entree image + decoupage par seuil'
             print *,'  > lecture du fichier de data : ',mesh
             call lit_mesh_seuil(mesh,data,coordmin,coordmax)
          elseif (mesh=='COORD') then
             data%coord=1
             read (1,*) mesh
             print *,'  > format d entree image + decoupage par seuil'
             print *,'  > lecture du fichier de data : ',mesh
             call lit_mesh_coord(mesh,data,coordmin,coordmax)
          else
             print *
             print *,'format de donnees non reconnu !!!'
             call help
          end if
          if ((data%image==1).or.(data%geom==1).or.(data%seuil==1)) then
             !creation du tableau de correspondances pixel/coord
             print *,'  > decodage du format image...'
             call tableau_image(data)
          end if
       case('EPAISSEUR')
          ok=0
          read(1,*) epsilon
          print *,'  > epaisseur de la tranche :',epsilon
       case('NBLIMIT')
          ok=0
          read(1,*) nblimit
          print *,'  > nb maximal de clusters recherches :',nblimit
       case('NBCLUST')
          ok=0
          read(1,*) listenbideal(:)
          print *,'  > test pour nb de clusters=',listenbideal
       case('SIGMA')
          ok=0
          read(1,*) sigma
          print *,'  > valeur de sigma imposee :',sigma
          if (data%image==1) then
             if (sigma<1.0) then
                print *,'epaisseur trop petite pour le mode image !!!!'
                stop
             end if
          end if
       case('DECOUPAGE')
          decoupage=1
          ok=0
          read (1,*) mot
          select case(mot)
          case('INTERFACE')
             data%interface=1
             print *,'  > decoupage par interface active.'
          case('RECOUVREMENT')
             data%recouvrement=1
             print *,'  > decoupage par recouvrement active.'
          case default
             print *
             print *,'mauvais format de decoupage !!!'
             print *
             call help
          end select
          print *, 'dim', data%dim
          if ((data%coord==1).or.(data%geom==1).or.(data%seuil==1)) then
             allocate(decoupe(data%dim))
          elseif (data%image==1) then
             !decoupage par pixel
             allocate(decoupe(data%imgdim))
          end if
          read(1,*) decoupe(:)
          print *, 'decoupe', decoupe
          if (nbproc>1) then
             tot=1
             if ((data%coord==1).or.(data%geom==1).or.(data%seuil==1)) then
                do i=1,data%dim
                   tot=tot*decoupe(i)
                enddo
             elseif (data%image==1) then
                !decoupage par pixel
                do i=1,data%imgdim
                   tot=tot*decoupe(i)
                enddo
             end if
             if (tot/=nbproc-data%interface) then
                print *,'decoupage non valide !'
                print *,'le nombre de proc doit etre egal a',tot+data%interface
                call MPI_ABORT(ierr)
                stop
             end if
          else
             !mode 1 proc
             if ((data%coord==1).or.(data%geom==1).or.(data%seuil==1)) then
                do i=1,data%dim
                   decoupe(i)=1
                enddo
                tot=1
             elseif (data%image==1) then
                !decoupage par pixel
                do i=1,data%imgdim
                   decoupe(i)=1
                enddo
                tot=1
             end if
          end if
          print *,'  > decoupage :',decoupe
       case('END')
          ok=1
       case default
          ok=0
          print *,'mot clé inconnu :',mot
       end select
    end do
    !parametre de decoupage
    if ((nbproc>1).and.(decoupage==0)) then
       print *
       print *,'mot cle DECOUPAGE absent !'
       call help 
    end if
    !cas monoproc
    if (nbproc==1) then
       !initialisation par defaut a 1 de tous les parametres de decoupage
       if (decoupage==1) deallocate(decoupe)
       if ((data%coord==1).or.(data%geom==1).or.(data%seuil==1)) then
          allocate(decoupe(data%dim) )
       elseif (data%image==1) then
          allocate(decoupe(data%imgdim))
       end if
       decoupe(:)=1
       epsilon=1.0
    end if   
    !validation des combinaisons de parametres d'entree
    tot=data%geom+data%seuil+data%coord+data%image
    if (tot/=1) then
       print *
       print *,'pb dans les formats de data entres !'
       call help
    end if
    tot=data%interface+data%recouvrement
    if (tot/=1) then
       print *
       print *,'pb dans les formats de decoupage entres !'
       call help
    end if
    return
  end subroutine lit

  !**********************************
  !lecture des datas en format coord
  subroutine lit_mesh_coord(mesh,data,coordmin,coordmax)
    implicit none
    character*30 :: mesh
    type(type_data) :: data
    double precision,dimension(:),pointer :: coordmax,coordmin
    integer :: i,j,nb
    !lecture de donnees classiques
    open(file=mesh,unit=2)
    read(2,*) data%nb,data%dim
    data%nbclusters=0
    print *,'    > nb de points :',data%nb
    print *,'    > dimension :',data%dim
    allocate(data%point(data%nb))
    allocate(coordmax(data%dim))
    allocate(coordmin(data%dim))
    nb=0
    do i=1,data%nb
       allocate(data%point(i)%coord(data%dim))
       read(2,*,end=100) data%point(i)%coord(:)
       nb=nb+1
       data%point(i)%cluster=-1
       if (i==1) then
          coordmax(:)=data%point(1)%coord(:)
          coordmin(:)=data%point(1)%coord(:)
       else
          do j=1,data%dim
             coordmin(j)=min(coordmin(j),data%point(i)%coord(j))
             coordmax(j)=max(coordmax(j),data%point(i)%coord(j))
          enddo
       endif
    enddo
100 print *,'nb de points',nb
    data%nb=nb
    close(2)
    print *,'  > coordonnees min/max :'
    do j=1,data%dim
       print *,'    > ',j,':',coordmin(j),coordmax(j)
    enddo
    return
  end subroutine lit_mesh_coord

  !**********************************
  !lecture d'image
  subroutine lit_mesh_image(mesh,data,coordmin,coordmax)
    implicit none
    character*30 :: mesh
    type(type_data) :: data
    double precision,dimension(:),pointer :: coordmax,coordmin
    integer :: i,j,nb
    open(file=mesh,unit=2)
    read(2,*) data%imgdim,data%imgt
    print *,'    >dimension de l image:',data%imgdim
    print *,'    >nb de temps:',data%imgt
    allocate(data%imgmap(data%imgdim))
    read(2,*) data%imgmap(:)
    print *,'    >decoupage spatial:',data%imgmap
    data%nb=1
    do i=1,data%imgdim
       data%nb=data%nb*data%imgmap(i)
    end do
    data%dim=data%imgt
    data%nbclusters=0
    print *,'    > nb de points a lire :',data%nb
    allocate(data%point(data%nb))
    allocate(coordmax(data%imgdim))
    allocate(coordmin(data%imgdim))
    coordmin(:)=0.9
    do i=1,data%imgdim
       coordmax(i)=data%imgmap(i)+0.1
    end do
    nb=0
    do i=1,data%nb
       allocate(data%point(i)%coord(data%dim))
       read(2,*,end=200) data%point(i)%coord(:)
       nb=nb+1
       data%point(i)%cluster=-1
    enddo
200 print *,'    > nb de points lus',nb       
    data%nb=nb
    close(2)
    print *,'  > coordonnees min/max :'
    do j=1,data%dim
       print *,'    > ',j,':',coordmin(j),coordmax(j)
    enddo
    return
  end subroutine lit_mesh_image
  
  !**********************************
  !lecture image en mode geom
  subroutine lit_mesh_geom(mesh,data,coordmin,coordmax)
    implicit none
    character*30 :: mesh
    type(type_data) :: data
    double precision,dimension(:),pointer :: coordmax,coordmin
    double precision :: pasmax
    integer :: i,j,nb
    open(file=mesh,unit=2)
    read(2,*) data%imgdim,data%imgt
    print *,'    >dimension de l image:',data%imgdim
    print *,'    >nb de temps:',data%imgt
    allocate(data%pas(data%imgdim))
    data%pas(:)=0.0
    allocate(data%imgmap(data%imgdim))
    read(2,*) data%imgmap(:)
    print *,'    >decoupage spatial:',data%imgmap
    data%nb=1
    do i=1,data%imgdim
       data%nb=data%nb*data%imgmap(i)
    end do
    data%dim=data%imgdim+data%imgt
    data%nbclusters=0
    print *,'    > nb de points a lire :',data%nb
    allocate(data%point(data%nb))
    allocate(coordmax(data%dim))
    allocate(coordmin(data%dim))
    nb=0
    do i=1,data%nb
       allocate(data%point(i)%coord(data%dim))
       data%point(i)%coord(:)=0.0
       read(2,*,end=300) data%point(i)%coord(data%imgdim+1:data%imgdim+data%imgt)
       nb=nb+1
       data%point(i)%cluster=-1
       if (i==1) then
          coordmax(:)=data%point(1)%coord(:)
          coordmin(:)=data%point(1)%coord(:)
       else
          do j=1,data%dim
             coordmin(j)=min(coordmin(j),data%point(i)%coord(j))
             coordmax(j)=max(coordmax(j),data%point(i)%coord(j))
          enddo
       endif
    enddo
300 print *,'    > nb de points lus',nb
    data%nb=nb
    close(2)
    print *,'  > coordonnees min/max :'
    pasmax=1.e-13
    do j=data%imgdim+1,data%imgdim+data%imgt
       pasmax=max(pasmax,coordmax(j)-coordmin(j))
       print *,'    > ',j,':',coordmin(j),coordmax(j)
    enddo
    print *,'  > pas max :',pasmax
    !recherche des pas par dimension d'image
    do j=1,data%imgdim
       data%pas(j)=pasmax/data%imgmap(j)
       print *,'  > pas :',j,data%pas(j)
       coordmin(j)=0.9*data%pas(j)
       coordmax(j)=(data%imgmap(j)+1)*data%pas(j)
    end do
    return
  end subroutine lit_mesh_geom

  !**********************************
  !lecture des datas en format seuil
  subroutine lit_mesh_seuil(mesh,data,coordmin,coordmax)
    implicit none
    character*30 :: mesh
    type(type_data) :: data
    double precision,dimension(:),pointer :: coordmax,coordmin
    integer :: i,j,nb
    !lecture de donnees classiques
    open(file=mesh,unit=2)
    read(2,*) data%imgdim,data%imgt
    print *,'    >dimension de l image:',data%imgdim
    print *,'    >nb de temps:',data%imgt
    allocate(data%imgmap(data%imgdim))
    read(2,*) data%imgmap(:)
    print *,'    >decoupage spatial:',data%imgmap
    data%nb=1
    do i=1,data%imgdim
       data%nb=data%nb*data%imgmap(i)
    end do
    data%dim=data%imgt
    data%nbclusters=0
    print *,'    > nb de points a lire :',data%nb
    allocate(data%point(data%nb))
    allocate(coordmax(data%dim))
    allocate(coordmin(data%dim))
    nb=0
    do i=1,data%nb
       allocate(data%point(i)%coord(data%dim))
       read(2,*,end=400) data%point(i)%coord(:)
       nb=nb+1
       data%point(i)%cluster=-1
       if (i==1) then
          coordmax(:)=data%point(1)%coord(:)
          coordmin(:)=data%point(1)%coord(:)
       else
          do j=1,data%dim
             coordmin(j)=min(coordmin(j),data%point(i)%coord(j))
             coordmax(j)=max(coordmax(j),data%point(i)%coord(j))
          enddo
       endif
    enddo
400 print *,'nb de points',nb
    data%nb=nb
    close(2)
    print *,'  > coordonnees min/max :'
    do j=1,data%dim
       print *,'    > ',j,':',coordmin(j),coordmax(j)
    enddo
    return
  end subroutine lit_mesh_seuil

  !**********************************
  !mise en tableau des indices de points pour les formats image
  subroutine tableau_image(data)
    implicit none
    type(type_data) :: data
    integer :: i,j,k,ok
    integer,dimension(:),pointer :: plan
    !creation du tableau de references points/coordonnes_images
    allocate(data%refimg(data%nb,data%imgdim))
    allocate(plan(data%imgdim)); plan(:)=1
    do i=1,data%nb
       do j=1,data%imgdim
          !index dans le tableau de reference points/pixel
          data%refimg(i,j)=plan(j)
          if (data%geom==1) then
             !entree des coordonnees 1:imgdim pour le cluster geom
             data%point(i)%coord(j)=plan(j)*data%pas(j)
          end if
       end do
       ok=0
       k=data%imgdim
       do while(ok==0)
          if (plan(k)<data%imgmap(k)) then
             plan(k)=plan(k)+1
             ok=1
          else
             plan(k)=1
             k=k-1
          end if
          if (k==0) ok=1
       end do
    end do
    deallocate(plan)
    return
  end subroutine tableau_image


end module module_entree
