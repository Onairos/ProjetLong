MODULE module_embed
  USE module_structure
CONTAINS



  SUBROUTINE spectral_embedding(nbcluster,n,Z,A,ratio,cluster,&
       cluster_center,cluster_population,cluster_energy,nbinfo,numproc,&
       ratiomoy,ratiorij,ratiorii)

    !*****************************************
    ! spectral embedding
    !
    ! nbcluster = nbre de cluster
    ! dataw : points
    ! Z : matrice des vecteurs propres
    ! M : nbre de vp trouveesx
    ! ratio : max des ration de frob sur matrice aff reordonnancee suivant
    ! les clusters
    ! cluster : appartenance des clusters
    ! cluster_center : centre des nbclusters clusters
    ! cluster_population : nbre de points par cluster
    ! cluster_energy : somme des energies par cluster
    !

    IMPLICIT NONE
    DOUBLE PRECISION,DIMENSION(:,:),POINTER:: Z,A,cluster_center
    INTEGER ::nbcluster,n,nbinfo,numproc
    DOUBLE PRECISION ::ratio,test,ratiomin,ratiorii,ratiorij, ratiomoy
    DOUBLE PRECISION,DIMENSION(:),POINTER :: cluster_energy,Z3
    INTEGER,DIMENSION(:),POINTER ::cluster,cluster_population
    DOUBLE PRECISION, DIMENSION(:,:),POINTER :: Frob
    DOUBLE PRECISION,DIMENSION(:,:),POINTER::Z1,Z2
    INTEGER :: it_max,it_num,i,j,k
    INTEGER,DIMENSION(:,:),POINTER :: clustercorresp
    INTEGER :: ki,kj,ni,nj,ok,nbmax
    ALLOCATE(cluster(n));
    ALLOCATE(cluster_center(nbcluster,nbcluster));
    ALLOCATE(cluster_population(nbcluster));
    ALLOCATE(cluster_energy(nbcluster));
    ALLOCATE(Z1(n,nbcluster))
    ALLOCATE(Z2(nbcluster,n));
    ALLOCATE(Z3(n));Z3(:)=0.0

    DO i=1,n
       DO j=1,nbcluster
          Z1(i,j)=Z(i,j)
          Z3(i)=Z3(i)+Z1(i,j)**2
       ENDDO
    ENDDO

    DO i=1,n
       test=0.0
       DO j=1,nbcluster
          Z2(j,i)=Z1(i,j)/(sqrt(Z3(i)))
          test=test+Z2(j,i)**2
       ENDDO
    ENDDO

    it_max=n*n

    CALL kmeans_01 ( nbcluster, n, nbcluster, it_max, it_num,Z2,&
         cluster, cluster_center, cluster_population, cluster_energy, &
         numproc)

    !*****************************
    ! Mesure de qualite
    nbmax=0
    DO i=1,nbcluster
       nbmax=max(nbmax,cluster_population(i))
    ENDDO
    ALLOCATE(clustercorresp(nbcluster,nbmax)); clustercorresp(:,:)=0
    DO i=1,n
       j=cluster(i)
       ok=0;k=1
       DO WHILE(ok==0)
          IF (clustercorresp(j,k)==0) THEN
             ok=1
          ELSE
             k=k+1
          ENDIF
       ENDDO
       clustercorresp(j,k)=i
    ENDDO
    ALLOCATE(Frob(nbcluster,nbcluster)); Frob(:,:)=0.0 
    DO i=1,nbcluster
       DO j=1,nbcluster
          DO ki=1,cluster_population(i)
             ni=clustercorresp(i,ki)
             DO kj=1,cluster_population(j)
                nj=clustercorresp(j,kj)
                Frob(i,j)=Frob(i,j)+A(ni,nj)**2
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(clustercorresp)
    ratio=0.0
    ratiomin=1.D+16
    ratiorii=0.0
    ratiorij=0.0
    nbinfo=nbcluster
    DO i=1,nbcluster
       IF ((cluster_population(i)/=0).AND.(Frob(i,i)/=0)) THEN
          DO j=1,nbcluster
             IF (i/=j) THEN
                ratio=ratio+Frob(i,j)/Frob(i,i)
                ratiomoy=ratiomoy+Frob(i,j)/Frob(i,i)
                ratiorij=ratiorij+Frob(i,j)
                ratiorii=ratiorii+Frob(i,i)
                ratiomin=min(ratiomin,Frob(i,j)/Frob(i,i))
             ENDIF
          ENDDO
       ELSE
          nbinfo=nbinfo-1
       ENDIF
       ratiorij=ratiorij*2/(nbcluster*(nbcluster-1))
       ratiorii=ratiorii!/nbcluster
    ENDDO
    DEALLOCATE(Frob)

#if aff
    PRINT *,numproc,'nbinfo=', nbinfo,' nbcluster=',nbcluster
#endif

    RETURN 
  END SUBROUTINE spectral_embedding


  !********************************************************
  !K-means
  SUBROUTINE kmeans_01 (dim_num,point_num,cluster_num,it_max, it_num, point, &
       cluster, cluster_center, cluster_population, cluster_energy, numproc)

    !*****************************************************************************80
    !
    !! KMEANS_01 applies the K-Means algorithm.
    !
    !  Discussion:
    !
    !    Given a matrix of POINT_NUM observations on DIM_NUM variables, the
    !    observations are to be ALLOCATEd to CLUSTER_NUM clusters in such 
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
    !    Input, INTEGER DIM_NUM, the number of spatial dimensions.
    !
    !    Input, INTEGER POINT_NUM, the number of points.
    !
    !    Input, INTEGER CLUSTER_NUM, the number of clusters.
    !
    !    Input, INTEGER IT_MAX, the maximum number of iterations.
    !
    !    Output, INTEGER IT_NUM, the number of iterations taken.
    !
    !    Input, DOUBLE PRECISION POINT(DIM_NUM,POINT_NUM), the points.
    !
    !    Output, INTEGER CLUSTER(POINT_NUM), indicates which cluster
    !    each point belongs to.
    !
    !    Input/output, DOUBLE PRECISION CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
    !    the cluster centers.
    !
    !    Output, INTEGER CLUSTER_POPULATION(CLUSTER_NUM), the number 
    !    of points in each cluster.
    !
    !    Output, DOUBLE PRECISION CLUSTER_ENERGY(CLUSTER_NUM), the 
    !    cluster energies.
    !
    IMPLICIT NONE
    INTEGER cluster_num
    INTEGER dim_num
    INTEGER point_num
    INTEGER cluster(point_num)
    DOUBLE PRECISION cluster_center(dim_num,cluster_num)
    DOUBLE PRECISION stockcenter(dim_num,cluster_num)
    DOUBLE PRECISION listnorm(point_num,cluster_num)
    DOUBLE PRECISION cluster_energy(cluster_num),stockenergy(cluster_num)
    INTEGER cluster_population(cluster_num)
    INTEGER stockpopulation(cluster_num)
    INTEGER i
    INTEGER it_max
    INTEGER it_num
    INTEGER j
    INTEGER k
    DOUBLE PRECISION point(dim_num,point_num)
    INTEGER swap
    INTEGER  :: ok,p,ok2
    DOUBLE PRECISION :: val,valmax,seuil,norme !,diffenergy
    INTEGER :: cluster_id(cluster_num)

    INTEGER :: numproc
    CHARACTER*30 :: files, num

    it_num = 0
    !
    !  Idiot checks.
    !
    IF ( cluster_num < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  CLUSTER_NUM < 1.0'
       STOP
    ENDIF

    IF ( dim_num < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  DIM_NUM < 1.0'
       STOP
    ENDIF

    IF ( point_num < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  POINT_NUM < 1.0'
       STOP
    ENDIF
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
#if aff
PRINT *, 'recherche des centres'
#endif
    DO i = 2, cluster_num
       ok=0
       DO WHILE(ok==0)
          valmax=2.0*seuil
          !recherche si le point est deja utilise dans comme centre
          ok2=0
          DO j=1,i-1
             IF (cluster_id(j)==p) ok2=1
          ENDDO
          !si point pas centre, teste par rapport au seuil
          IF (ok2==0) THEN
             DO j=1,i-1
                val=0.0; norme=0.0
                DO k=1,dim_num
                   val=max(val,abs(cluster_center(k,j)-point(k,p)))
                ENDDO
                valmax=min(val,valmax)
             ENDDO
             IF (valmax>=seuil) ok=1
          ENDIF
         p=p+1

         !abaisse le seuil si pas assez de centre sont trouves
         IF ((p>point_num).AND.(ok==0)) THEN 
            seuil=0.9*seuil
#if aff
            PRINT *,'abaisse seuil :',seuil
#endif
            p=1
          ENDIF
       ENDDO
       p=p-1
       cluster_center(:,i)=point(:,p)
       cluster_id(i)=p
    ENDDO
#if aff
   PRINT *,'centres initiaux',p
#endif
!cluster_center

!!! boucle            
    it_num = 0
    swap=1
    cluster(:)=1
    DO WHILE ((it_num<it_max).AND.(swap/=0))
       it_num = it_num + 1
       swap=0
       DO i=1,cluster_num
          stockenergy(i)=cluster_energy(i);
          stockpopulation(i)=cluster_population(i);
          DO j=1,dim_num
             stockcenter(j,i)=cluster_center(j,i)
          ENDDO
       ENDDO

       !! Calcul de toutes les distances
       cluster_population(1:cluster_num) = 1
       listnorm(:,:)=0.0
       DO i=1,point_num
          DO j=1,cluster_num
             DO k=1,dim_num
                listnorm(i,j)=listnorm(i,j)+(point(k,i)-cluster_center(k,j))**2
             ENDDO
          ENDDO
       ENDDO

       !!assignation par rapport au min des distances
       cluster_population(:)=0.0
       DO i=1,point_num
          DO j=1,cluster_num
             IF (listnorm(i,j)<listnorm(i,cluster(i))) THEN
                cluster(i)=j
                swap=swap+1
             ENDIF
          ENDDO
          cluster_energy(cluster(i))=cluster_energy(cluster(i))&
               +listnorm(i,cluster(i))
          cluster_population(cluster(i))=cluster_population(cluster(i))+1
       ENDDO

       !! mise a jour des centres
       cluster_center(:,:)=0.0
       DO j=1,point_num
          i=cluster(j) 
          DO k=1,dim_num
             cluster_center(k,i)=cluster_center(k,i)+point(k,j)
          ENDDO
       ENDDO
       DO i=1,cluster_num
          cluster_center(:,i)=cluster_center(:,i)/cluster_population(i)
       ENDDO



    ENDDO

    RETURN
  END SUBROUTINE kmeans_01



END MODULE module_embed
