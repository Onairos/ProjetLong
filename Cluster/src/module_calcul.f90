module module_calcul
  use module_structure
  use module_solve
  use module_embed
contains

  !*****************************************
  !calcul du sigma
  subroutine calculsigma(dataw,sigma)
    implicit none
    type(type_data) :: dataw
    double precision :: sigma,norme,norme1,sigma1
    integer :: i1,j1,k1

    sigma=0.0
    sigma1=0.0
    do i1=1,dataw%nb
       do j1=i1+1,dataw%nb
          norme=0.0
          do k1=1,dataw%dim
             norme=norme+&
                  (dataw%point(i1)%coord(k1)-dataw%point(j1)%coord(k1))**2
          end do
          sigma=max(sigma,sqrt(norme))
       end do
    end do
    sigma=sigma/(2*exp(log(float(dataw%nb))*(1.0/float(dataw%dim))))
    !securite
    if (sigma==0.0) sigma=1.0

    !do i1=1,dataw%nb
    !   norme1=10**8
    !   do j1=1,dataw%nb
    !      norme=0.0
    !      if (j1/=i1) then
    !         do k1=1,dataw%dim
    !            norme=norme+&
    !                 (dataw%point(i1)%coord(k1)-dataw%point(j1)%coord(k1))**2
    !         end do
    !         norme1=min(norme1,sqrt(norme))
    !      end if
    !   enddo
    !   sigma1=sigma1+sqrt(norme1)/dataw%nb
    !!end do
    !sigma1=sqrt(sigma1)


    !print *
    !print *,numproc,'valeur de sigma calculee :',sigma
    !print *,numproc,'2eme valeur de sigma calculee :',sigma1
    return
  end subroutine calculsigma

  !*****************************************
  !calcul du sigma pour l'interface
  subroutine calculsigmainterface(numproc,dataw,sigma,bornes,decoupe,epsilon)
    implicit none
    type(type_data) :: dataw
    double precision :: sigma
    integer :: i,j,k,numproc,nb,ok
    double precision,dimension(:,:,:),pointer :: bornes
    integer,dimension(:),pointer :: decoupe,decoupe0
    integer,dimension(:,:),pointer :: tableau
    double precision :: epsilon
    double precision :: long,sigma0
    double precision :: volext,volint
    !nb de decoupes
    nb=1
    do i=1,dataw%dim
       nb=nb*decoupe(i)
    end do
    !creation des decoupes
    allocate(tableau(nb,0:dataw%dim))
    allocate(decoupe0(dataw%dim)); decoupe0(:)=1
    do i=1,nb
       do j=1,dataw%dim
          tableau(i,j)=decoupe0(j)
       enddo
       decoupe0(1)=decoupe0(1)+1
       ok=0; k=1
       do while(ok==0)
          if (decoupe0(k)>decoupe(k)) then
             decoupe0(k)=1
             if (k<dataw%dim) decoupe0(k+1)=decoupe0(k+1)+1
          else
             ok=1
          end if
       end do
    end do
    deallocate(decoupe0)
    !valeur de sigma
    sigma0=0.0
    do i=1,nb
       volext=1.0; volint=1.0
       do j=1,dataw%dim
          k=tableau(i,j)
          long=bornes(j,k,2)-bornes(j,k,1)
          volext=volext*long
          volint=volint*max(0.0,long-2.0*epsilon)
       enddo
       !print *,'volext',volext,'volint=',volint
       sigma0=sigma0+volext-volint
    end do
    deallocate(tableau)
    !calcul du grandeur equivalente
    sigma0=exp(1.0/float(dataw%dim)*log(sigma0))
    !calcul du sigma
    sigma0=sigma0/(2.0*exp(log(float(dataw%nb))*(1.0/float(dataw%dim))))
    !calcul du sigma formule globale
    call calculsigma(dataw,sigma)
#if aff
    print *,numproc,'valeur de sigma calculee pour interface:',sigma0
#endif
    !sigma=min(sigma,sigma0)
#if aff
    print *,numproc,'valeur sigma interface',sigma
#endif
    return
  end subroutine calculsigmainterface


  !*****************************************
  !calcul des clusters
  subroutine calculclusters(numproc,nblimit,nbideal,dataw,sigma)
    implicit integer(i,j,q)

    include 'mpif.h'
    type(type_data) :: dataw
    integer :: numproc,nbproc
    double precision :: sigma
    double precision,dimension(:,:),pointer :: A,Z,A2,cluster_center
    double precision,dimension(:),pointer :: W
    integer :: n,k,nbcluster
    double precision, dimension(:),pointer :: D
    double precision,dimension(:),pointer :: ratiomax,cluster_energy,&
         ratiomin,ratiomoy,ratiorii,ratiorij
    integer,dimension(:),pointer ::cluster,cluster_population,nbinfo
    integer :: nblimit, nbideal, nb, nbvp
    double precision :: norme,ratio,ratio1,ratio2,seuilrij
    character*30 :: num,files
    double precision :: t1, t2, t_cons_vp

    ! deux valeurs qui quand elles ne sont pas déclarées
    ! et donc implicitement des real font que ça marche mieux
    real :: val, value

    ! solveur au valeur propre => paramètre de contrôle
    integer :: solver

    !creation de la matrice
    print *,numproc,'valeur du sigma',sigma
#if aff
    print *,numproc,'valeur du sigma',sigma
#endif
    n=dataw%nb
    allocate(A(n,n));A(:,:)=0.0

    ! A(i,i) = 0
    do i=1,n-1
       do j=i+1,n
          norme=0.0
          do k=1,dataw%dim
             norme=norme+(dataw%point(i)%coord(k)-dataw%point(j)%coord(k))**2
          end do
          value=exp(-norme/sigma)
          ! partie triangulaire sup
          A(i,j) = value
          ! partie triangulaire inf
          A(j,i)=A(i,j)
       end do
    end do

    !Normalisation de la matrice affinite
    allocate(D(n)); D(:)=0.0
    do i=1,n
       do j=1,i-1
          D(i)=D(i)+A(i,j)
       end do
       do j=i+1, n
          D(i)=D(i)+A(j,i)
       end do
    end do

    do i=1,n
       do j=1,n
          ! la matrice A n'est plus symétique
          A(i,j)=A(i,j)/D(i)
       end do
    end do
    deallocate(D)

    ! solver 0 = lapack - 1 = arpack
    solver = 0

    if(solver == 0) then
      print *, numproc, 'solveur lapack'

      nbvp = n

      allocate(A2(n,n)); A2(:,:)=0.0
      do i=1,n
        do j=1,n
          A2(i,j)=A(i,j)
        end do
      end do

      ! lapack
      ! n dimension in
      ! A2 matrice in
      ! Z(N,N) vecteurs propres out
      ! W(N) valeurs propres out

      t1 = MPI_WTIME()
      call solve_dgeev(n,A2,Z,W)
    else
      print *, numproc, 'solveur arpack'

      ! arpack
      nb = 2*nblimit
      nbvp = nb
      call solve_arpack_full(A, n, nb, W, Z)

    end if

    t2 = MPI_WTIME()

    t_cons_vp = t2 - t1
    print *, numproc, 'cout construction vp', t_cons_vp

    !print *,numproc,'reordonne les vp...'
    do i=1,nbvp-1
       do j=i+1,nbvp
          if (W(i)<W(j)) then
             val=W(i); W(i)=W(j); W(j)=val
             do k=1,n
                val=Z(k,i); Z(k,i)=Z(k,j); Z(k,j)=val
             end do
          end if
       end do
       !print *,'vp ',i,W(i)
    end do

    !Test spectral embedding avec different nbcluster   
    !***********************
    ! Spectral embedding
    !print *,numproc,'Spectral Embedding'

    if ((nbideal==0).and.(n>2)) then
       !** recherche du meilleur decoupage
       allocate(ratiomax(nblimit)); ratiomax(:)=0
       allocate(ratiomin(nblimit)); ratiomin(:)=0
       allocate(ratiomoy(nblimit)); ratiomoy(:)=0
       allocate(ratiorii(nblimit)); ratiorii(:)=0
       allocate(ratiorij(nblimit)); ratiorij(:)=0

       allocate(nbinfo(nblimit)); nbinfo(:)=0
       do nbcluster=2,min(n,nblimit)
          !print *,numproc,'teste avec nbcluster=',nbcluster
          !call flush(6)

          allocate(cluster(n));cluster(:)=0.0
          allocate(cluster_center(nbcluster,nbcluster)); cluster_center(:,:)=0.0
          allocate(cluster_population(nbcluster));cluster_population(:)=0.0
          allocate(cluster_energy(nbcluster));cluster_energy(:)=0.0

          call spectral_embedding(nbcluster,n,Z,A,&
               ratiomax(nbcluster),cluster,cluster_center,cluster_population,&
               cluster_energy,nbinfo(nbcluster),numproc,ratiomoy(nbcluster), &
               ratiorij(nbcluster),ratiorii(nbcluster))

          !write(num,*),numproc
          !num=adjustl(num)
          !files='clusterk.partiel'//trim(num)
          !write(num,*),nbcluster
          !num=adjustl(num)
          !files=trim(files)//'.'//trim(num)
          !len=len(trim(num))
          !print *,numproc,'ecriture des clusters : ',files
          !open(file=files,unit=20)
          !write(20,*) dataw%nb,dataw%dim
          !do i=1,n
          !   write(20,*) dataw%point(i)%coord(:),cluster(i)
          !enddo
          !close(20)

          deallocate(cluster);deallocate(cluster_center);
          deallocate(cluster_energy)
          deallocate(cluster_population);
       end do

       !write(num,*),numproc
       !num=adjustl(num)
!       files='ratio'
!          write(num,*),numproc
!          num=adjustl(num)
!          files=trim(files)//'.'//trim(num)
         ! len=len(trim(num))
         ! print *,numproc,'ecriture des ratio : '
!          open(file=files,unit=20)
!          do i=1,min(n,nblimit)
!             write(20,*) ratiorij(i),ratiorii(i),ratiomoy(i)
!          enddo
!          close(20)



#if aff
print *, 'ratio de frobenius'
#endif
       !*******************************
       ! Ratio de norme de frobenius
       !print *,numproc,'Ratio max par cluster',ratiomax(2:nblimit)
       !print *,numproc,'Ratio min par cluster',ratiomin(2:nblimit)
       !print *,numproc,'Ratio Rii par cluster',ratiorii(2:nblimit)
       !print *,numproc,'Ratio Rij par cluster',ratiorij(2:nblimit)
       ratio=ratiomax(nblimit); dataw%nbclusters=nblimit
       ratio1=0.0;ratio2=1e+10
       do i=2,nblimit
          !if  ((nbinfo(i)==i).and.(ratio(i)<ratio)) then
          !if  ((nbinfo(i)==i).and.(ratiomoy(i)<ratio)) then
          if ((numproc==0).and.(nbproc>1)) then 
             seuilrij=1e-1
          else
             seuilrij=1e-4
          end if
          if ((ratiorii(i)>=0.95*ratio1).and.(ratiorij(i)-ratio2<=seuilrij)) then  
             !if (ratiomoy(i)-ratio1<=1e-4) then
             !2eme critère
             !(ratiorij(i)/ratiorii(i)<=1e-4)
             dataw%nbclusters=i
             ! ratio=ratiomax(i)
             ratio1=ratiorii(i)
             ratio2=ratiorij(i)
             !ratio1=ratiomoy(i)
          end if
       end do
      ! print *,numproc,'nb de clusters final :',dataw%nbclusters

    elseif ((nbideal==1).and.(n>nbideal)) then
       !** test avec un cluster impose
       allocate(nbinfo(nbideal)); nbinfo(:)=0
       allocate(ratiomin(1)); ratiomin(:)=0.0
       dataw%nbclusters=nbideal
    else
       !** cas d'un domaine avec moins de points que nbideal ou 1 seul point
       allocate(nbinfo(n)); nbinfo(:)=0
       allocate(ratiomin(1)); ratiomin(:)=0.0
       dataw%nbclusters=n
       allocate(ratiomax(n)); ratiomax(:)=0
       allocate(ratiomoy(n)); ratiomoy(:)=0
       allocate(ratiomin(n)); ratiomin(:)=0
       allocate(ratiorii(n)); ratiorii(:)=0
       allocate(ratiorij(n)); ratiorij(:)=0
    endif
    ! cas où nbcluster==1
    if (dataw%nbclusters==2) then
       print *, 'difference ratio',ratiorij(2)/ratiorii(2)
       if (ratiomax(2)>=0.6) then 
          dataw%nbclusters=1
       else 
          dataw%nbclusters=2
       end if
    end if
#if aff
    print *,numproc,'cluster final obtenu : ',dataw%nbclusters
#endif

    !** calcul du clustering final
    if (dataw%nbclusters>1) then
       call spectral_embedding(dataw%nbclusters,n,Z,A,ratio,cluster,&
            cluster_center,cluster_population,cluster_energy,&
            nbinfo(dataw%nbclusters),numproc,ratiomin(1),ratiorij(1),ratiorii(1))
       do i=1,dataw%nb
          dataw%point(i)%cluster=cluster(i)
       enddo
       deallocate(cluster)
       deallocate(cluster_population)
       deallocate(ratiomax)
       deallocate(cluster_energy)
       deallocate(ratiomin)
       deallocate(ratiomoy)
       deallocate(ratiorii)
       deallocate(ratiorij)
       deallocate(A)
       deallocate(Z)
       if(solver == 0) deallocate(A2)
       deallocate(cluster_center)
       deallocate(W)
    else 
       !cluster_population(1)=dataw%nb
#if aff
       print *, numproc, 'ok'
#endif
       do i=1,dataw%nb
          dataw%point(i)%cluster=1
       end do
#if aff
       print *,numproc,'cluster'
#endif
    end if


    !affichage
    !do j=1,dataw%nbclusters
    !   print *,numproc,'centres des clusters',cluster_center(:,j)
    !end do
    !maj du cluster
  

    !deallocations
   
!    call flush()
    return
  end subroutine calculclusters

end module module_calcul
