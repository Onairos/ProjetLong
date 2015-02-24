
SUBROUTINE (numproc,nblimit,nbideal,dataw,sigma,kernelFunIndex)

   INCLUDE 'mpif.h'
    TYPE(type_data) :: dataw
    INTEGER :: numproc,nbproc,kernelFunIndex
    DOUBLE PRECISION :: sigma
    DOUBLE PRECISION,DIMENSION(:,:),POINTER :: A,Z,A2,cluster_center
    DOUBLE PRECISION,DIMENSION(:),POINTER :: W
    INTEGER :: n,k,nbcluster
    DOUBLE PRECISION, DIMENSION(:),POINTER :: D
    DOUBLE PRECISION,DIMENSION(:),POINTER :: ratiomax,cluster_energy,&
         ratiomin,ratiomoy,ratiorii,ratiorij
    INTEGER,DIMENSION(:),POINTER ::cluster,cluster_population,nbinfo
    INTEGER :: nblimit, nbideal, nb, nbvp
    DOUBLE PRECISION :: norme,ratio,ratio1,ratio2,seuilrij
    CHARACTER*30 :: num,files
    DOUBLE PRECISION :: t1, t2, t_cons_vp

    ! deux valeurs qui quand elles ne sont pas declarees
    ! et donc implicitement des REAL font que ca marche mieux
    DOUBLE PRECISION :: val, value

! Poly Kernel K(xi,xj)=(xi'xj+gamma)^delta
!Gaussian Kernel K(xi,xj)=exp(-||xi-xj||^2/2sigma^2)

||phi(xi)-mk||² = Kii - (2*sum(Indicatrice(xj E Ck)*Kij)/ sum (Indicatrice(xj E Ck))
+ sum (sum(Indicatrice(xj E Ck)*Indicatrice(xl E Ck)*Kjl))/sum (sum(Indicatrice(xj E Ck)*Indicatrice(xl E Ck)

For i=1:nbCluster
For j=1:n
c*(xj)=argmin_i(||phi(xj)-mi||² renvoit l'indice du centre de cluster le plus proche
Cluster[c*(xj)] << xj

Stop when converged compute E = sum_N(sum_M( Indicatrice (xi E Ck)*||phi(xi)-mk||²))
