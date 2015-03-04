MODULE clustering_method
CONTAINS

  DOUBLE PRECISION, DIMENSION(:,:) function poly_kernel( dataw, gam, delta )
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: dataw
    DOUBLE PRECISION gam 
    DOUBLE PRECISION delta

    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: K
    INTEGER n
    INTEGER i
    INTEGER j
    INTEGER k

    n=dataw%nb
    ALLOCATE(K(n,n))
    K(:,:)=0.0

    DO i=1,n-1
      DO j=1,n-1
        DO k=1,dataw%dim
        K(i,j)=K(i,j)+dataw%point(i)%coord(k)*dataw%point(j)%coord(k)
        ENDDO 
        K(i,j)=(K(i,j)+gam)^delta
      ENDDO
    ENDDO
    poly_kernel=K
    RETURN
  END

  DOUBLE PRECISION, DIMENSION(:,:) function gaussian_kernel( dataw, sigma )
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: dataw
    DOUBLE PRECISION sigma

    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: K
    INTEGER n
    INTEGER i
    INTEGER j
    INTEGER k

    n=dataw%nb
    ALLOCATE(K(n,n))
    K(:,:)=0.0

    DO i=1,n-1
      DO j=i+1,n
        DO k=1,dataw%dim
        K(i,j)=K(i,j)+(dataw%point(i)%coord(k)-dataw%point(j)%coord(k))**2
        ENDDO
        K(i,j)=exp(- K(i,j)/2*sigma^2)
        ! Symetry
        K(j,i)=K(i,j)
      ENDDO
    ENDDO
    gaussian_kernel=K
    RETURN
  END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



||phi(xi)-mk||² = Kii - (2*sum(Indicatrice(xj E Ck)*Kij)/ sum (Indicatrice(xj E Ck))
+ sum (sum(Indicatrice(xj E Ck)*Indicatrice(xl E Ck)*Kjl))/sum (sum(Indicatrice(xj E Ck)*Indicatrice(xl E Ck)

For i=1:nbCluster
For j=1:n
c*(xj)=argmin_i(||phi(xj)-mi||² renvoit l'indice du centre de cluster le plus proche
Cluster[c*(xj)] << xj

Stop when converged compute E = sum_N(sum_M( Indicatrice (xi E Ck)*||phi(xi)-mk||²))



SUBROUTINE apply_kernel_k_means(numproc,nblimit,nbideal,dataw,hash,kernelFunIndex)

   INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
  
    INTEGER :: nbideal
    INTEGER :: nblimit
    INTEGER :: numproc
    INTEGER ::kernelFunIndex !Kernel matrix selection value

    !=== IN/OUT ===
    TYPE(type_data) :: dataw

    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################      
    INTEGER ::point_num=dataw%nb
    INTEGER ::dim_num=dataw%dim

    !#### Variables  ####
   !!!!!!!!!!!!!!!!!!!!! DOUBLE PRECISION :: point (dim_num, point_num) ! the points
    INTEGER :: cluster_num ! the number of clusters
    cluster_num=dataw%nbclusters

    INTEGER :: point_num ! the number of points
    point_num=dataw%nb

    INTEGER :: dim_num ! the number of spatial dimensions
    dim_num=dataw%dim
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: K
    INTEGER :: it_max ! the maximum number of iterations
    DOUBLE PRECISION :: cluster_center (dim_num, cluster_num) ! the cluster centers
    DOUBLE PRECISION :: cluster_energy (cluster_num) ! the cluster energies
    INTEGER :: it_num ! the number of iterations taken
    INTEGER :: cluster (point_num) ! indicates which cluster each point belongs to
    INTEGER :: cluster_population (cluster_num) ! the number of points in each cluster
    DOUBLE PRECISION :: listnorm (point_num, cluster_num)
    DOUBLE PRECISION :: stockcenter (dim_num, cluster_num)
    DOUBLE PRECISION :: stockenergy (cluster_num)
    DOUBLE PRECISION :: norme
    DOUBLE PRECISION :: seuil
    DOUBLE PRECISION :: val
    DOUBLE PRECISION :: valmax
    INTEGER :: cluster_id (cluster_num)
    INTEGER :: stockpopulation (cluster_num)
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: l
    INTEGER :: k
    DOUBLE PRECISION :: num1
    DOUBLE PRECISION :: den1
    DOUBLE PRECISION :: num2
    DOUBLE PRECISION :: den2
    INTEGER :: ok ! TODO: booleen ?
    INTEGER :: ok2 ! TODO: booleen ?
    INTEGER :: swap
    INTEGER :: p
    

    ALLOCATE(cluster(n));
    ALLOCATE(cluster_center(nbcluster,nbcluster));
    ALLOCATE(cluster_population(nbcluster));
    ALLOCATE(cluster_energy(nbcluster));

    ALLOCATE(K(n,n))
    K(:,:)=0.0


    !###########################################      
    ! INSTRUCTIONS
    !###########################################   
    it_num = 0
    !
    !  Idiot checks.
    !
    IF ( cluster_num < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  CLUSTER_NUM < 1.0'
       STOP
    ENDIF

    IF ( dim_num < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  DIM_NUM < 1.0'
       STOP
    ENDIF

    IF ( point_num < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  POINT_NUM < 1.0'
       STOP
    ENDIF

    IF (kernelfunindex==0 AND ( hash("gamma")=NULL OR hash("delta")=NULL) ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  GAMMA AND DELTA NOT INITIALIZED IN POLYNOMIAL KERNEL'
       STOP  
    

    ELSEIF (kernelfunindex==1 AND hash("sigma")=NULL) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  SIGMA NOT INITIALIZED IN GAUSSIAN KERNEL'
       STOP 
    ENDIF


    IF (kernelfunindex==0) THEN
        K=poly_kernel( dataw, hash("gamma"), hash("delta"))
    ELSEIF (kernelfunindex==1) THEN
        K=gaussian_kernel(dataw, hash("sigma"))
    END



    !
    !  For each observation, calculate the distance from each cluster
    !  center, and assign to the nearest.
    !

    !
    !  Assign one point to each cluster center.
    !
    cluster_center(:,1) = dataw%point(1)%coord(:) !point(:,1) %%
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
             IF (dataw%point(j)%cluster==p) ok2=1
          ENDDO
          !si point pas centre, teste par rapport au seuil
          IF (ok2==0) THEN
             DO j=1,i-1
                val=0.0; norme=0.0
                DO k=1,dataw%dim
                   val=max(val,abs(cluster_center(k,j)-point(k,p))) 
!VOIR SI CELA DOIT ÊTRE MODIFIE EN FONCTION DES KERNEL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
       cluster_center(:,i)= dataw%point(P)%coord(:) !point(:,p) 
       cluster_id(i)=p
    ENDDO
#if aff
   PRINT *,'centres initiaux',p
#endif
!cluster_center

!!! boucle            
    it_num = 0
    swap=1
    dataw%point(:)%cluster=1 !  cluster(:)=1
    DO WHILE ((it_num<it_max).AND.(swap/=0))
       it_num = it_num + 1
       swap=0
       DO i=1,dataw%nbclusters
          stockenergy(i)=cluster_energy(i);
          stockpopulation(i)=cluster_population(i);
          DO j=1,dataw%dim
             stockcenter(j,i)=cluster_center(j,i)
          ENDDO
       ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !! Calcul de toutes les distances
       cluster_population(1:cluster_num) = 1
       listnorm(:,:)=0.0
       num1=0.0
       den1=0.0
       num2=0.0
       den2=0.0
       DO k=1,dataw%nbclusters
           DO i=1,dataw%nb
               DO j=1,dataw%nb
                   IF ( dataw%point(j)%cluster.EQ.k) THEN
                   num1=num1 + 2*(K(i,j))
                   den1=den1+1 
                   END
                   DO l=1,dataw%nb
                       IF ( dataw%point(j)%cluster.EQ.k AND dataw%point(l)%cluster.EQ.k) THEN
                       num2=num2 + K(j,l)
                       den2=den2+1
                       END
                   ENDDO
               ENDDO
           listnorm(i,k)= K(i,i)
               IF (den1.NE.0.0) THEN
               listnorm(i,k)=  listnorm(i,k) - (num1/den1)
               END
               IF (den2.NE.0.0) THEN
               listnorm(i,k)=  listnorm(i,k)  + (num2/den2)
               END
           ENDDO
       ENDDO
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!assignation par rapport au min des distances
       cluster_population(:)=0
       DO i=1,dataw%nb
          DO j=1,dataw%nbclusters
             IF (listnorm(i,j)<listnorm(i,dataw%point(i)%cluster)) THEN
                dataw%point(i)%cluster=j
                swap=swap+1
             ENDIF
          ENDDO
          cluster_energy(dataw%point(i)%cluster)=cluster_energy(dataw%point(i)%cluster)&
               +listnorm(i,dataw%point(i)%cluster)
          cluster_population(dataw%point(i)%cluster)=cluster_population(dataw%point(i)%cluster)+1
       ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !! mise a jour des centres
       cluster_center(:,:)=0.0
       DO j=1,dataw%nb
          i=dataw%point(j)%cluster 
          DO k=1,dataw%dim
             cluster_center(k,i)=cluster_center(k,i)+dataw%point(j)%coord(k)
          ENDDO
       ENDDO
       DO i=1,dataw%nbclusters
          cluster_center(:,i)=cluster_center(:,i)/cluster_population(i)
       ENDDO



    ENDDO

    RETURN
  END SUBROUTINE apply_kmeans
