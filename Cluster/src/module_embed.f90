module module_embed
  use module_structure
contains



  subroutine spectral_embedding(nbcluster,n,Z,A,ratio,cluster,&
       cluster_center,cluster_population,cluster_energy,nbinfo,numproc,&
       ratiomoy,ratiorij,ratiorii)

    !*****************************************
    ! spectral embedding
    !
    ! nbcluster = nbre de cluster
    ! dataw : points
    ! Z : matrice des vecteurs propres
    ! M : nbre de vp trouvéesx
    ! ratio : max des ration de frob sur matrice aff réordonnancée suivant
    ! les clusters
    ! cluster : appartenance des clusters
    ! cluster_center : centre des nbclusters clusters
    ! cluster_population : nbre de points par cluster
    ! cluster_energy : somme des énergies par cluster
    !

    implicit none
    !type(type_data) :: dataw
    double precision,dimension(:,:),pointer:: Z,A,cluster_center
    integer ::nbcluster,n,nbinfo,numproc
    double precision ::ratio,test,ratiomin,ratiorii,ratiorij, ratiomoy
    double precision,dimension(:),pointer :: cluster_energy,Z3
    integer,dimension(:),pointer ::cluster,cluster_population
    !integer,dimension(:),pointer::ordaffperclus
    double precision, dimension(:,:),pointer :: Frob
    double precision,dimension(:,:),pointer::Z1,Z2
    integer :: it_max,it_num,i,j,k
    integer,dimension(:,:),pointer :: clustercorresp
    integer :: ki,kj,ni,nj,ok,nbmax
    allocate(cluster(n));
    allocate(cluster_center(nbcluster,nbcluster));
    allocate(cluster_population(nbcluster));
    allocate(cluster_energy(nbcluster));
    allocate(Z1(n,nbcluster));  allocate(Z2(nbcluster,n));
    allocate(Z3(n));Z3(:)=0.0

    do i=1,n
       do j=1,nbcluster
          Z1(i,j)=Z(i,j)
          !if (i==1) print *,numproc,'matrice vp',j,W(j)
          Z3(i)=Z3(i)+Z1(i,j)**2
       end do
    end do

    do i=1,n
       test=0.0
       do j=1,nbcluster
          Z2(j,i)=Z1(i,j)/(sqrt(Z3(i)))
          test=test+Z2(j,i)**2
       end do
    end do

!    print *, numproc,'methode kmeans'

    it_max=n*n !1000.

    call kmeans_01 ( nbcluster, n, nbcluster, it_max, it_num,Z2,&
         cluster, cluster_center, cluster_population, cluster_energy, &
         numproc)
    !print *,numproc,'fin de kmeans. nb d iterations effectuees : ',it_num

    !print *,numproc,'Nombre points par cluster', cluster_population
    ! print *,'vecteur cluster',cluster(1:5) 

    !*****************************
    ! Mesure de qualité
    !print *,'Indexation'
    nbmax=0
    do i=1,nbcluster
       nbmax=max(nbmax,cluster_population(i))
    end do
    allocate(clustercorresp(nbcluster,nbmax)); clustercorresp(:,:)=0
    do i=1,n
       j=cluster(i)
       ok=0;k=1
       do while(ok==0)
          if (clustercorresp(j,k)==0) then
             ok=1
          else
             k=k+1
          end if
       end do
       clustercorresp(j,k)=i
    end do
    allocate(Frob(nbcluster,nbcluster)); Frob(:,:)=0.0 
    do i=1,nbcluster
       do j=1,nbcluster
          do ki=1,cluster_population(i)
             ni=clustercorresp(i,ki)
             do kj=1,cluster_population(j)
                nj=clustercorresp(j,kj)
                Frob(i,j)=Frob(i,j)+A(ni,nj)**2
             end do
          end do
       end do
    end do
    deallocate(clustercorresp)
    ratio=0.0; ratiomin=1.D+16;ratiorii=0.0;ratiorij=0.0
    nbinfo=nbcluster
    do i=1,nbcluster
       if ((cluster_population(i)/=0).and.(Frob(i,i)/=0)) then
          do j=1,nbcluster
             if (i/=j) then
                ratio=ratio+Frob(i,j)/Frob(i,i)
                !ratio=max(ratio,Frob(i,j)/Frob(i,i))
                ratiomoy=ratiomoy+Frob(i,j)/Frob(i,i)
                ratiorij=ratiorij+Frob(i,j)
                ratiorii=ratiorii+Frob(i,i)
                ratiomin=min(ratiomin,Frob(i,j)/Frob(i,i))
             endif
          end do
       else
          nbinfo=nbinfo-1
       end if
       ratiorij=ratiorij*2/(nbcluster*(nbcluster-1))
       ratiorii=ratiorii!/nbcluster
    end do
    deallocate(Frob)

!!$    p2=1
!!$    dataw%point(:)%cluster=cluster(:)
!!$    allocate(P(n,n)); P(:,:)=0.0
!!$    !print *,'petit test', dataw%point(1:5)%cluster
!!$    do i=1,nbcluster
!!$       do j=1,n
!!$          ! P(j,:)=0.0
!!$          if (dataw%point(j)%cluster==i) then 
!!$             P(p2,j)=1.0
!!$             p2=p2+1
!!$          end if
!!$       end do
!!$    end do
!!$    !print *,
!!$    !print *, 'matrice permut',P
!!$    allocate(P1(n,n)); P1(:,:)=0.0
!!$    P1=MATMUL(P,A)
!!$    print *,numproc,'Mesure de qualite : Frobenius norm'
!!$    allocate(ordaffperclus(nbcluster+1));ordaffperclus(:)=0.0
!!$    ordaffperclus(1)=0.0
!!$    do i=2,nbcluster+1
!!$       ordaffperclus(i)=ordaffperclus(i-1)+cluster_population(i-1)
!!$    end do
!!$    ! print *, 'ordreaffperclus',  ordaffperclus
!!$    allocate(Frob(nbcluster,nbcluster)); Frob(:,:)=0.0 
!!$    nbinfo=nbcluster
!!$    do i=1,nbcluster
!!$       if (cluster_population(i)/=0) then
!!$          do p3=i,nbcluster
!!$             normfro=0.0
!!$             do j=ordaffperclus(i)+1,ordaffperclus(i+1)
!!$                do q=ordaffperclus(p3)+1,ordaffperclus(p3+1)
!!$                   normfro=normfro+P1(j,q)**2
!!$                end do
!!$             end do
!!$             if (i/=p3) then                
!!$                Frob(i,p3)=normfro/Frob(i,i)
!!$                ratio=ratio+Frob(i,p3)/nbcluster
!!$             else
!!$                Frob(i,i)=normfro
!!$             end if
!!$          end do
!!$       else
!!$          nbinfo=nbinfo-1
!!$       end if
!!$    end do
!!$    ! print *, 
!!$    ! print *,'norme de frobenius ',Frob
#if aff
    print *,numproc,'nbinfo=', nbinfo,' nbcluster=',nbcluster
#endif
    !deallocate(ordaffperclus);deallocate(Frob)
    !deallocate(Z1), deallocate(Z2),deallocate(P);deallocate(P1)

    return 
  end subroutine spectral_embedding


  !********************************************************
  !K-means
  subroutine kmeans_01 (dim_num,point_num,cluster_num,it_max, it_num, point, &
       cluster, cluster_center, cluster_population, cluster_energy, numproc)

    !*****************************************************************************80
    !
    !! KMEANS_01 applies the K-Means algorithm.
    !
    !  Discussion:
    !
    !    Given a matrix of POINT_NUM observations on DIM_NUM variables, the
    !    observations are to be allocated to CLUSTER_NUM clusters in such 
    !    a way that the within-cluster sum of squares is minimized.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    18 November 2004
    !
    !  Author:
    !
    !    FORTRAN77 original version by David Sparks
    !    FORTRAN90 version by John Burkardt
    !
    !  Reference:
    !
    !    David Sparks,
    !    Algorithm AS 58: 
    !    Euclidean Cluster Analysis,
    !    Applied Statistics,
    !    Volume 22, Number 1, 1973, pages 126-130.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
    !
    !    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
    !
    !    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
    !
    !    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
    !
    !    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
    !
    !    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the points.
    !
    !    Output, integer ( kind = 4 ) CLUSTER(POINT_NUM), indicates which cluster
    !    each point belongs to.
    !
    !    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
    !    the cluster centers.
    !
    !    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number 
    !    of points in each cluster.
    !
    !    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the 
    !    cluster energies.
    !
    implicit none
    integer ( kind = 4 ) cluster_num
    integer ( kind = 4 ) dim_num
    integer ( kind = 4 ) point_num
    integer ( kind = 4 ) cluster(point_num)
    real    ( kind = 8 ) cluster_center(dim_num,cluster_num)
    real    ( kind = 8 ) stockcenter(dim_num,cluster_num)
    !real    ( kind = 8 ) diffcenter(dim_num)
    real    ( kind = 8 ) listnorm(point_num,cluster_num)
    real    ( kind = 8 ) cluster_energy(cluster_num),stockenergy(cluster_num)
    integer ( kind = 4 ) cluster_population(cluster_num)
    integer ( kind = 4 ) stockpopulation(cluster_num)
    !real    ( kind = 8 ) dc,de,epsilon
    !real    ( kind = 8 ) f(point_num)
    integer ( kind = 4 ) i
    !integer ( kind = 4 ) il
    !integer ( kind = 4 ) diffpopulation
    integer ( kind = 4 ) it_max
    integer ( kind = 4 ) it_num
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    !integer ( kind = 4 ) list(1)
    real    ( kind = 8 ) point(dim_num,point_num)
    integer ( kind = 4 ) swap
    integer  :: ok,p,ok2
    real (kind=8) :: val,valmax,seuil,norme !,diffenergy
    integer :: cluster_id(cluster_num)

    integer :: numproc
    character*30 :: files, num

    it_num = 0
    !
    !  Idiot checks.
    !
    if ( cluster_num < 1 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       write ( *, '(a)' ) '  CLUSTER_NUM < 1.0'
       stop
    end if

    if ( dim_num < 1 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       write ( *, '(a)' ) '  DIM_NUM < 1.0'
       stop
    end if

    if ( point_num < 1 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       write ( *, '(a)' ) '  POINT_NUM < 1.0'
       stop
    end if
    !
    !  For each observation, calculate the distance from each cluster
    !  center, and assign to the nearest.
    !

    !
    !  Assign one point to each cluster center.
    !
    cluster_center(:,1) = point(:,1)
    cluster_id(:)=0; cluster_id(1)=1
    p=2
    seuil=0.4
    !do i=2,cluster_num
    !   ok=0
    !   do while (ok==0)
    !      valmax=2.0*seuil
    !      do j=1,i-1
    !         val=0.0
    !         do k=1,dim_num
    !            val=max(val,abs(cluster_center(k,j)-point(k,p)))
    !         end do
    !         valmax=min(valmax,val)
    !      end do
    !      if (valmax>seuil) ok=1
    !      p=p+1
    !   end do
    !   cluster_center(:,i)=point(:,p-1)
    !end do
#if aff
print *, 'recherche des centres'
#endif
    do i = 2, cluster_num
       ok=0
       do while(ok==0)
          valmax=2.0*seuil
          !recherche si le point est deja utilise dans comme centre
          ok2=0
          do j=1,i-1
             if (cluster_id(j)==p) ok2=1
          end do
          !si point pas centre, teste par rapport au seuil
          if (ok2==0) then
             do j=1,i-1
                val=0.0; norme=0.0
                do k=1,dim_num
                   val=max(val,abs(cluster_center(k,j)-point(k,p)))
                  !norme=norme+abs(cluster_center(k,j)-point(k,p))**2
                enddo
               !norme=sqrt(norme)
                !valmax=min(valmax,val/norme)
                valmax=min(val,valmax)
               !if (valmax>seuil) ok=1
             enddo
             if (valmax>=seuil) ok=1
          end if
         p=p+1
!if (seuil<=1e-4) then
!do i1=p+1,cluster_num

         !abaisse le seuil si pas assez de centre sont trouves
         if ((p>point_num).and.(ok==0)) then 
            seuil=0.9*seuil
#if aff
            print *,'abaisse seuil :',seuil
#endif
            p=1
          end if
       end do
       p=p-1
       cluster_center(:,i)=point(:,p)
       cluster_id(i)=p
    end do
#if aff
   print *,'centres initiaux',p
#endif
!cluster_center

!!! boucle            
    it_num = 0
    swap=1
    cluster(:)=1
    do while ((it_num<it_max).and.(swap/=0))
       it_num = it_num + 1
       swap=0
       do i=1,cluster_num
          stockenergy(i)=cluster_energy(i);
          stockpopulation(i)=cluster_population(i);
          do j=1,dim_num
             stockcenter(j,i)=cluster_center(j,i)
          end do
       end do

       !! Calcul de toutes les distances
       cluster_population(1:cluster_num) = 1
       ! allocate(listnorm(point_num,cluster_num));
       listnorm(:,:)=0.0
       do i=1,point_num
          do j=1,cluster_num
             do k=1,dim_num
                listnorm(i,j)=listnorm(i,j)+(point(k,i)-cluster_center(k,j))**2
             end do
          end do
       end do

       !!assignation par rapport au min des distances
       cluster_population(:)=0.0
       do i=1,point_num
          do j=1,cluster_num
             if (listnorm(i,j)<listnorm(i,cluster(i))) then
                cluster(i)=j
                swap=swap+1
             end if
          end do
          cluster_energy(cluster(i))=cluster_energy(cluster(i))&
               +listnorm(i,cluster(i))
          cluster_population(cluster(i))=cluster_population(cluster(i))+1
       end do

       !write(num,*) numproc
       !files = 'kmeans.'//trim(adjustl(num))
       !OPEN(file=files, UNIT=15, ACCESS="APPEND")
       !write (15,*) numproc,'iteration',it_num,'test clustering',cluster(1:3)
       !CALL FLUSH(15)
       !CLOSE(15)

       !print *,numproc,'nbre de permut',swap

       !! mise à jour des centres
       cluster_center(:,:)=0.0
       do j=1,point_num
          i=cluster(j) 
          do k=1,dim_num
             cluster_center(k,i)=cluster_center(k,i)+point(k,j)
          end do
       end do
       do i=1,cluster_num
          cluster_center(:,i)=cluster_center(:,i)/cluster_population(i)
       end do



    end do










!!$
!!$    do i = 1, point_num
!!$       do j = 1, cluster_num
!!$          cluster_energy(j) = sum ( &
!!$               ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
!!$       end do
!!$       list = minloc ( cluster_energy(1:cluster_num) )
!!$       cluster(i) = list(1)
!!$    end do
!!$    !
!!$    !  Determine the cluster population counts.
!!$    !
!!$    cluster_population(1:cluster_num) = 0
!!$
!!$    do i = 1, point_num
!!$       j = cluster(i)
!!$       cluster_population(j) = cluster_population(j) + 1
!!$    end do
!!$    !
!!$    !  Calculate the mean and sum of squares for each cluster.
!!$    !
!!$
!!$    do i = 1, point_num
!!$       j = cluster(i)
!!$       cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
!!$            + point(1:dim_num,i)
!!$    end do
!!$
!!$    do i = 1, cluster_num
!!$       if ( 0 < cluster_population(i) ) then
!!$          cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
!!$               / real ( cluster_population(i), kind = 8 )
!!$       end if
!!$    end do
!!$    !
!!$    !  Set the point energies.
!!$    !
!!$    f(1:point_num) = 0.0
!!$
!!$    do i = 1, point_num
!!$       j = cluster(i)
!!$       f(i) = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
!!$    end do
!!$    !
!!$    !  Set the cluster energies.
!!$    !
!!$    cluster_energy(1:cluster_num) = 0.0
!!$
!!$    do i = 1, point_num
!!$       j = cluster(i)
!!$       cluster_energy(j) = cluster_energy(j) + f(i)
!!$    end do
!!$    !
!!$    !  Adjust the point energies by a weight factor.
!!$    !
!!$    do i = 1, point_num
!!$       j = cluster(i)
!!$       if ( 1 < cluster_population(j) ) then
!!$          f(i) = f(i) * real ( cluster_population(j), kind = 8 ) &
!!$               / real ( cluster_population(j) - 1, kind = 8 )
!!$       end if
!!$    end do
!!$    !
!!$    !  Examine each observation in turn to see if it should be
!!$    !  reassigned to a different cluster.
!!$    !
!!$    it_num = 0
!!$    do while ( it_num < it_max )
!!$       it_num = it_num + 1
!!$       swap = 0
!!$       do i = 1, point_num
!!$          il = cluster(i)
!!$          ir = il
!!$          if ( cluster_population(il) <= 1 ) then
!!$             cycle
!!$          end if
!!$          dc = f(i)
!!$          do j = 1, cluster_num
!!$             if ( j /= il ) then
!!$                de = sum ( &
!!$                     ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 ) &
!!$                     * real ( cluster_population(j), kind = 8 ) &
!!$                     / real ( cluster_population(j) + 1, kind = 8 )
!!$                if ( de < dc ) then
!!$                   dc = de
!!$                   ir = j
!!$                end if
!!$             end if
!!$          end do
!!$          !
!!$          !  If the lowest value was obtained by staying in the current cluster,
!!$          !  then cycle.
!!$          !
!!$          if ( ir == il ) then
!!$             cycle
!!$          end if
!!$          !
!!$          !  Reassign the point from cluster IL to cluster IR.
!!$          !
!!$          cluster_center(1:dim_num,il) = &
!!$               ( cluster_center(1:dim_num,il) &
!!$               * real ( cluster_population(il), kind = 8 ) &
!!$               - point(1:dim_num,i) ) / real ( cluster_population(il) - 1, kind = 8 )
!!$
!!$          cluster_center(1:dim_num,ir) = &
!!$               ( cluster_center(1:dim_num,ir) &
!!$               * real ( cluster_population(ir), kind = 8 ) &
!!$               + point(1:dim_num,i) ) / real ( cluster_population(ir) + 1, kind = 8 )
!!$
!!$          cluster_energy(il) = cluster_energy(il) - f(i)
!!$          cluster_energy(ir) = cluster_energy(ir) + dc
!!$          cluster_population(ir) = cluster_population(ir) + 1
!!$          cluster_population(il) = cluster_population(il) - 1
!!$
!!$          cluster(i) = ir
!!$          !
!!$          !  Adjust the value of F for points in clusters IL and IR.
!!$          !
!!$          do j = 1, point_num
!!$             k = cluster(j)
!!$             if ( k == il .or. k == ir ) then
!!$                f(j) = sum ( &
!!$                     ( point(1:dim_num,j) - cluster_center(1:dim_num,k) )**2 )
!!$                if ( 1 < cluster_population(k) ) then
!!$                   f(j) = f(j) * real ( cluster_population(k), kind = 8 ) &
!!$                        / ( real ( cluster_population(k) - 1, kind = 8 ) )
!!$                end if
!!$             end if
!!$          end do
!!$          swap = swap + 1
!!$       end do
!!$       !
!!$       !  Exit if no reassignments were made during this iteration.
!!$       !
!!$       if ( swap == 0 ) then
!!$          exit
!!$       end if
!!$    end do
!!$    !
!!$    !  Compute the cluster energies.
!!$    !
!!$    cluster_energy(1:cluster_num) = 0.0
!!$    do i = 1, point_num
!!$       j = cluster(i)
!!$       cluster_energy(j) = cluster_energy(j) + sum ( &
!!$            ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
!!$    end do
!!$

    return
  end subroutine kmeans_01



end module module_embed
