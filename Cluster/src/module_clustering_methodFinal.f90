MODULE module_clustering_methodFinal
  USE module_look_up_table
  USE module_structure

  
CONTAINS

    FUNCTION poly_kernel( dataw, gam, delta )
    IMPLICIT NONE
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: dataw
    DOUBLE PRECISION ::gam 
    DOUBLE PRECISION :: delta

    !====  OUT  ====
    DOUBLE PRECISION, DIMENSION(dataw%nb,dataw%nb) :: poly_kernel
    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: K
    INTEGER :: n
    INTEGER :: i
    INTEGER :: j
    INTEGER :: d

    n=dataw%nb
    ALLOCATE(K(n,n))
    !ALLOCATE(polyKernel)
    K(:,:)=0.0

    DO i=1,n-1
      DO j=1,n-1
        DO d=1,dataw%dim
        K(i,j)=K(i,j)+dataw%point(i)%coord(d)*dataw%point(j)%coord(d)
        ENDDO 
        K(i,j)=(K(i,j)+gam)**delta
      ENDDO
    ENDDO
    poly_kernel=K
    RETURN
  END

    FUNCTION gaussian_kernel( dataw, sigma )
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
    TYPE(type_data) :: dataw
    DOUBLE PRECISION sigma


    !====  OUT  ====
    DOUBLE PRECISION, DIMENSION(dataw%nb,dataw%nb) :: gaussian_kernel

    !#### Variables  ####
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: K
    INTEGER n
    INTEGER i
    INTEGER j
    INTEGER d

    n=dataw%nb
    ALLOCATE(K(n,n))
    !ALLOCATE(gaussian_kernel)
    K(:,:)=0.0

    DO i=1,n-1
      DO j=i+1,n
        DO d=1,dataw%dim
        K(i,j)=K(i,j)+(dataw%point(i)%coord(d)-dataw%point(j)%coord(d))**2
        ENDDO
        K(i,j)=exp(- K(i,j)/(2*sigma**2))
        ! Symetry
        K(j,i)=K(i,j)
      ENDDO
    ENDDO
    gaussian_kernel=K
    RETURN
  END


!Stop when converged compute E = sum_N(sum_M( Indicatrice (xi E Ck)*||phi(xi)-mk||²))



SUBROUTINE apply_kernel_k_means(numproc,nblimit,nbideal,dataw,lu_table,kernelFunIndex)
    IMPLICIT NONE

  ! INCLUDE 'mpif.h'
    !###########################################
    ! DECLARATIONS
    !###########################################
    !#### Parameters ####
    !====  IN  ====
  
    INTEGER :: nbideal
    INTEGER :: nblimit
    INTEGER :: numproc
    INTEGER ::kernelFunIndex !Kernel matrix selection value
    TYPE(look_up_table) :: lu_table

    !=== IN/OUT ===
    TYPE(type_data) :: dataw


    !###########################################
    ! DECLARATIONS
    !###########################################      
   ! INTEGER ::point_num
   ! point_num=dataw%nb

    !#### Variables  ####

    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Ker
    INTEGER :: it_max ! the maximum number of iterations
    DOUBLE PRECISION :: cluster_center (dataw%dim, dataw%nbclusters) ! the cluster centers
    DOUBLE PRECISION :: cluster_energy (dataw%nbclusters) ! the cluster energies
    INTEGER :: it_num ! the number of iterations taken
    INTEGER :: cluster (dataw%nb) ! indicates which cluster each point belongs to
    INTEGER :: cluster_population (dataw%nbclusters) ! the number of points in each cluster
    DOUBLE PRECISION :: listnorm (dataw%nb, dataw%nbclusters)
    DOUBLE PRECISION :: stockcenter (dataw%dim, dataw%nbclusters)
    DOUBLE PRECISION :: stockenergy (dataw%nbclusters)
    DOUBLE PRECISION :: norme
    DOUBLE PRECISION :: seuil
    DOUBLE PRECISION :: val
    DOUBLE PRECISION :: valmax
    INTEGER :: cluster_id (dataw%nbclusters)
    INTEGER :: stockpopulation (dataw%nbclusters)
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: l

   ! LOGICAL :: found
    DOUBLE PRECISION :: num1
    DOUBLE PRECISION :: den1
    DOUBLE PRECISION :: num2
    DOUBLE PRECISION :: den2

    INTEGER :: ok ! TODO: booleen ?
    INTEGER :: ok2 ! TODO: booleen ?
    INTEGER :: swap
    INTEGER :: p
    

    ALLOCATE(Ker(dataw%nb,dataw%nb))
    Ker(:,:)=0.0


    !###########################################      
    ! INSTRUCTIONS
    !###########################################   
    it_num = 0
    !
    !  Idiot checks.
    !
    IF ( dataw%nbclusters < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  CLUSTER_NUM < 1.0'
       STOP
    ENDIF

    IF ( dataw%dim < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  DIM_NUM < 1.0'
       STOP
    ENDIF

    IF ( dataw%nb < 1 ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  POINT_NUM < 1.0'
       STOP
    ENDIF

 
    IF (kernelfunindex==0 .AND. ( .NOT. look_up_exists(lu_table,"gamma") .OR. .NOT. look_up_exists(lu_table,"delta") ) ) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  GAMMA AND DELTA NOT INITIALIZED IN POLYNOMIAL KERNEL'
       STOP  
    

    ELSEIF (kernelfunindex==1 .AND. .NOT. look_up_exists(lu_table,"sigma")) THEN
       WRITE ( *, '(a)' ) ' '
       WRITE ( *, '(a)' ) 'KERNELKMEANS_01 - Fatal error!'
       WRITE ( *, '(a)' ) '  SIGMA NOT INITIALIZED IN GAUSSIAN KERNEL'
       STOP 
    ENDIF


    IF (kernelfunindex==0) THEN
        Ker=poly_kernel( dataw, look_up_get(lu_table,"gamma"), look_up_get(lu_table,"delta"))
    ELSEIF (kernelfunindex==1) THEN
        Ker=gaussian_kernel(dataw, look_up_get(lu_table,"sigma"))
    ENDIF



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
!#if aff
PRINT *, 'recherche des centres'
!#endif
    DO i = 2, dataw%nbclusters
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
                   val=max(val,abs(cluster_center(k,j)-dataw%point(p)%coord(k))) 
!VOIR SI CELA DOIT ÊTRE MODIFIE EN FONCTION DES KERNEL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ENDDO
                valmax=min(val,valmax)
             ENDDO
             IF (valmax>=seuil) ok=1
          ENDIF
         p=p+1

         !abaisse le seuil si pas assez de centre sont trouves
         IF ((p>dataw%nb).AND.(ok==0)) THEN 
            seuil=0.9*seuil
!#if aff
            PRINT *,'abaisse seuil :',seuil
!#endif
            p=1
          ENDIF
       ENDDO
       p=p-1
       cluster_center(:,i)= dataw%point(P)%coord(:) !point(:,p) 
       cluster_id(i)=p
    ENDDO
!#if aff
   PRINT *,'centres initiaux',p
!#endif
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

       !! Calcul de toutes les distances
       cluster_population(1:dataw%nbclusters) = 1
       listnorm(:,:)=0.0
       num1=0.0
       den1=0.0
       num2=0.0
       den2=0.0
       DO k=1,dataw%nbclusters
           DO i=1,dataw%nb
               DO j=1,dataw%nb
                   IF ( dataw%point(j)%cluster.EQ.k) THEN
                   num1=num1 + 2*(Ker(i,j))
                   den1=den1+1 
                   ENDIF
                   DO l=1,dataw%nb
                       IF ( dataw%point(j)%cluster.EQ.k .AND. dataw%point(l)%cluster.EQ.k) THEN
                       num2=num2 + Ker(j,l)
                       den2=den2+1
                       ENDIF
                   ENDDO
               ENDDO
           listnorm(i,k)= Ker(i,i)
               IF (den1.NE.0.0) THEN
               listnorm(i,k)=  listnorm(i,k) - (num1/den1)
               ENDIF
               IF (den2.NE.0.0) THEN
               listnorm(i,k)=  listnorm(i,k)  + (num2/den2)
               ENDIF
           ENDDO
       ENDDO
     

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
  END SUBROUTINE apply_kernel_k_means
END MODULE module_clustering_methodFinal
