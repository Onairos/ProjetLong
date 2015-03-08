MODULE module_sparse
  USE module_structure
  USE module_solve
  USE module_embed
CONTAINS

  !*****************************************
  !calcul des clusters
  SUBROUTINE sp_calculclusters(numproc, nblimit, nbideal, partitioned_data, sigma)

    IMPLICIT INTEGER(i, j, q)
    INCLUDE 'mpif.h'
    TYPE(type_data) :: partitioned_data
    INTEGER :: numproc, nbproc
    DOUBLE PRECISION :: sigma
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: clusters_centers
    INTEGER :: n, k, nbcluster
    DOUBLE PRECISION, DIMENSION(:), POINTER :: ratiomax, clusters_energies, &
         ratiomin, ratiomoy, ratiorii, ratiorij
    INTEGER, DIMENSION(:), POINTER ::clusters, points_by_clusters, nbinfo
    INTEGER :: nblimit, nbideal
    DOUBLE PRECISION :: norme, ratio, ratio1, ratio2, seuilrij
    CHARACTER (LEN=30) :: num, files

! sparsification debut
    DOUBLE PRECISION :: t1, t2, t_cons_a, t_cons_vp
    INTEGER :: nnz, nnz2, nb
    DOUBLE PRECISION :: facteur
    INTEGER :: l
    DOUBLE PRECISION :: treshold
    DOUBLE PRECISION, DIMENSION(:), POINTER :: AS
    INTEGER, DIMENSION(:), POINTER :: IAS, JAS
    DOUBLE PRECISION, DIMENSION(:), POINTER :: D
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z
    DOUBLE PRECISION, DIMENSION(:), POINTER :: W
! sparsification fin

    !creation de la matrice
#if aff
    PRINT *,numproc,'valeur du sigma',sigma
#endif
    n=partitioned_data%nb

! sparsification debut
    nnz = 0
    ! valeur de treshold arbitraire -> parametre du sp ou calcul interne
    !                                  (voir avec S.)
    ! TODO : mettre la valeur du facteur dans le fichier param
    facteur = 3.0
    treshold = facteur*sigma

    t1 = MPI_WTIME()
    DO i=1,n-1  ! borne ?
       DO j=i+1,n ! borne ?

          norme=0.0

          DO k=1,partitioned_data%dim
             norme=norme+(partitioned_data%point(i)%coord(k)-partitioned_data%point(j)%coord(k))**2
          ENDDO

          IF(sqrt(norme) <= treshold) THEN
            nnz = nnz + 1
          ENDIF

       ENDDO
    ENDDO

    t2 = MPI_WTIME()
    t_cons_a = t2 - t1
    PRINT *, numproc, 'surcout A', t_cons_a

    t1 = MPI_WTIME()
    nnz2 = nnz*2

    ALLOCATE(AS(nnz2))
    ALLOCATE(IAS(nnz2))
    ALLOCATE(JAS(nnz2))
    l = 1
    DO i=1,n-1
       DO j=i+1,n
          norme=0.0
          DO k=1,partitioned_data%dim
             norme=norme+(partitioned_data%point(i)%coord(k)-partitioned_data%point(j)%coord(k))**2
          ENDDO
          value=exp(-norme/sigma)
          ! on garde si value <= treshold
          ! (si on veut tout garder, commenter ligne IF, ENDIF)
          IF(sqrt(norme) <= treshold) THEN
            AS(l) = value
            IAS(l) = i
            JAS(l) = j
            l = l+1
            !------
            AS(l) = value
            IAS(l) = j
            JAS(l) = i
            l = l+1
          ENDIF
       ENDDO
    ENDDO
    WRITE(*,*) '========== facteur, n*n nnz2 = ', facteur, n*n, nnz2

  ALLOCATE(D(n))
  D(:)=0.0
  DO l=1, nnz2
    D(IAS(l)) = D(IAS(l)) + AS(l)
  ENDDO

  DO l=1, nnz2
    AS(l)=AS(l)/D(IAS(l))
  ENDDO

  DEALLOCATE(D)

    ! nb et nblimit meme valeur ?
    nb = 2*nblimit

    t1 = MPI_WTIME()
    CALL solve_arpack(AS, IAS, JAS, n, nnz2, nb, W, Z)
    PRINT *, "---------- W -------------"
    DO i=1,nb
       PRINT *,'valeurs propres arpack brutes',i, W(i)
    ENDDO
    

    ! reordonne les vp... QUESTION: necessaire avec arpack ?
    DO i=1,nb-1
       DO j=i+1,nb
          IF (W(i)<W(j)) THEN
             value=W(i)
             W(i)=W(j)
             W(j)=value
             DO k=1,n
                value=Z(k,i)
                Z(k,i)=Z(k,j)
                Z(k,j)=value
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    DO i=1,nb
       PRINT *,'valeurs propres arpack reordonnees',i, W(i)
    ENDDO

    !Test spectral embedding avec different nbcluster   
    !***********************
    ! Spectral embedding

    IF ((nbideal==0).AND.(n>2)) THEN
       !** recherche du meilleur decoupage
       ALLOCATE(ratiomax(nblimit))
       ratiomax(:)=0
       ALLOCATE(ratiomin(nblimit))
       ratiomin(:)=0
       ALLOCATE(ratiomoy(nblimit))
       ratiomoy(:)=0
       ALLOCATE(ratiorii(nblimit))
       ratiorii(:)=0
       ALLOCATE(ratiorij(nblimit))
       ratiorij(:)=0

       ALLOCATE(nbinfo(nblimit))
       nbinfo(:)=0

       DO nbcluster = 2 ,min(n,nblimit)

          ALLOCATE(clusters(n))
          clusters(:)=0.0
          ALLOCATE(clusters_centers(nbcluster,nbcluster))
          clusters_centers(:,:)=0.0
          ALLOCATE(points_by_clusters(nbcluster))
          points_by_clusters(:)=0.0
          ALLOCATE(clusters_energies(nbcluster))
          clusters_energies(:)=0.0

          CALL sp_spectral_embedding(nbcluster, n, Z, nnz2, AS, IAS, JAS, &
               ratiomax(nbcluster),clusters,clusters_centers,points_by_clusters, &
               clusters_energies,nbinfo(nbcluster),numproc,ratiomoy(nbcluster), &
               ratiorij(nbcluster),ratiorii(nbcluster))

          DEALLOCATE(clusters)
          DEALLOCATE(clusters_centers)
          DEALLOCATE(clusters_energies)
          DEALLOCATE(points_by_clusters)
       ENDDO


#if aff
PRINT *, 'ratio de frobenius'
#endif
       !*******************************
       ! Ratio de norme de frobenius
       ratio=ratiomax(nblimit)
       partitioned_data%nbclusters=nblimit
       ratio1=0.0
       ratio2=1e+10

       DO i=2,nblimit
          IF ((numproc==0).AND.(nbproc>1)) THEN 
             seuilrij=1e-1
          ELSE
             seuilrij=1e-4
          ENDIF

          IF ((ratiorii(i)>=0.95*ratio1).AND.(ratiorij(i)-ratio2<=seuilrij)) THEN  
             partitioned_data%nbclusters=i
             ratio1=ratiorii(i)
             ratio2=ratiorij(i)
          ENDIF
       ENDDO

    ELSEIF ((nbideal==1).AND.(n>nbideal)) THEN
       !** test avec un cluster impose
       ALLOCATE(nbinfo(nbideal))
       nbinfo(:) = 0
       ALLOCATE(ratiomin(1))
       ratiomin(:) = 0.0
       partitioned_data%nbclusters = nbideal
    ELSE
       !** cas d'un domaine avec moins de points que nbideal ou 1 seul point
       ALLOCATE(nbinfo(n))
       nbinfo(:)=0
       ALLOCATE(ratiomin(1))
       ratiomin(:)=0.0
       partitioned_data%nbclusters=n
       ALLOCATE(ratiomax(n))
       ratiomax(:)=0
       ALLOCATE(ratiomoy(n))
       ratiomoy(:)=0
       ALLOCATE(ratiomin(n))
       ratiomin(:)=0
       ALLOCATE(ratiorii(n))
       ratiorii(:)=0
       ALLOCATE(ratiorij(n))
       ratiorij(:)=0
    ENDIF
    ! cas avec nbcluster==1
    IF (partitioned_data%nbclusters==2) THEN
       PRINT *, 'difference ratio',ratiorij(2)/ratiorii(2)
       IF (ratiomax(2)>=0.6) THEN 
          partitioned_data%nbclusters=1
       ELSE 
          partitioned_data%nbclusters=2
       ENDIF
    ENDIF
#if aff
    PRINT *,numproc,'cluster final obtenu : ',partitioned_data%nbclusters
#endif

    !** calcul du clustering final
    IF (partitioned_data%nbclusters>1) THEN

       CALL sp_spectral_embedding(partitioned_data%nbclusters, n, Z, nnz2, AS, IAS, JAS,ratio,clusters,&
            clusters_centers,points_by_clusters,clusters_energies,&
            nbinfo(partitioned_data%nbclusters),numproc,ratiomin(1),ratiorij(1),&
            ratiorii(1))

       DO i=1,partitioned_data%nb
          partitioned_data%point(i)%clusters=clusters(i)
       ENDDO

       DEALLOCATE(clusters)
       DEALLOCATE(points_by_clusters)
       DEALLOCATE(ratiomax)
       DEALLOCATE(clusters_energies)
       DEALLOCATE(ratiomin)
       DEALLOCATE(ratiomoy)
       DEALLOCATE(ratiorii)
       DEALLOCATE(ratiorij)
       DEALLOCATE(clusters_centers)

    ELSE 
#if aff
       PRINT *, numproc, 'ok'
#endif
       DO i=1,partitioned_data%nb
          partitioned_data%point(i)%clusters=1
       ENDDO
#if aff
       PRINT *,numproc,'cluster'
#endif
    ENDIF

    !deallocations
    DEALLOCATE(AS)
    DEALLOCATE(IAS)
    DEALLOCATE(JAS)
    DEALLOCATE(W)
    DEALLOCATE(Z)

    RETURN
  END SUBROUTINE sp_calculclusters

    SUBROUTINE sp_spectral_embedding(nbcluster, n, Z, nnz, AS, IAS, JAS, ratio, clusters, &
       clusters_centers, points_by_clusters, clusters_energies, nbinfo, numproc, &
       ratiomoy, ratiorij, ratiorii)

    !*****************************************
    ! spectral embedding
    !
    ! nbcluster = nbre de cluster
    ! Z : matrice des vecteurs propres
    ! M : nbre de vp trouvees
    ! ratio : max des ration de frob sur matrice aff reordonnancee suivant
    ! les clusters
    ! clusters : appartenance des clusters
    ! clusters_centers : centre des nbclusters clusters
    ! points_by_clusters : nbre de points par cluster
    ! clusters_energies : somme des energies par cluster
    !

    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z, clusters_centers
    INTEGER ::nbcluster, n, nbinfo, numproc
    DOUBLE PRECISION ::ratio, test, ratiomin, ratiorii, ratiorij, ratiomoy
    DOUBLE PRECISION, DIMENSION(:), POINTER :: clusters_energies, Z3
    INTEGER, DIMENSION(:), POINTER ::clusters, points_by_clusters
    !INTEGER,DIMENSION(:),POINTER::ordaffperclus
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Frob
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: Z1, Z2
    INTEGER :: it_max, it_num, i, j, k
    INTEGER, DIMENSION(:,:), POINTER :: clustercorresp
    INTEGER :: ki, kj, ni, nj, ok, nbmax

    DOUBLE PRECISION, DIMENSION(:), POINTER:: AS
    INTEGER, DIMENSION(:), POINTER :: IAS, JAS
    INTEGER :: nnz

    INTEGER :: l
    INTEGER :: num1, num2

    ALLOCATE(clusters(n))
    ALLOCATE(clusters_centers(nbcluster,nbcluster))
    ALLOCATE(points_by_clusters(nbcluster))
    ALLOCATE(clusters_energies(nbcluster))
    ALLOCATE(Z1(n,nbcluster))
    ALLOCATE(Z2(nbcluster,n))
    ALLOCATE(Z3(n))
    Z3(:)=0.0

    PRINT *, '************ sp_spectral_embedding *************'
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

    PRINT *, numproc,'methode kmeans'

    it_max=n*n !1000.0

    CALL apply_kmeans( nbcluster, n, nbcluster, it_max, it_num, Z2,&
         clusters, clusters_centers, points_by_clusters, clusters_energies, &
         numproc)

    !*****************************
    ! Mesure de qualite
    !PRINT *,'Indexation'

    nbmax=0
    DO i=1,nbcluster
       nbmax=max(nbmax,points_by_clusters(i))
    ENDDO
    PRINT *, points_by_clusters
    ALLOCATE(clustercorresp(nbcluster,nbmax))
    clustercorresp(:,:)=0
    DO i=1,n
       j=clusters(i)
       ok=0
       k=1
       DO WHILE(ok==0)
          IF (clustercorresp(j,k)==0) THEN
             ok=1
          ELSE
             k=k+1
          ENDIF
       ENDDO
       clustercorresp(j,k)=i
    ENDDO


! sparsification debut
    ALLOCATE(Frob(nbcluster,nbcluster))
    Frob(:,:)=0.0
    DO i=1, nnz
      num1 = clusters(IAS(i))
      num2 = clusters(JAS(i))
      Frob(num1, num2) = Frob(num1, num2) + AS(i)**2
    ENDDO
! sparsification fin


! sparsification debut
    ratio=0.0
    ratiomin=1.D+16
    ratiorii=0.0
    ratiorij=0.0
    ratiomoy = 0.0
    nbinfo=nbcluster
    DO i=1,nbcluster
       IF ((points_by_clusters(i)/=0).AND.(Frob(i,i)/=0)) THEN
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
       ratiomoy=ratiomoy*2/(nbcluster*(nbcluster-1))
       ratiorii=ratiorii!/nbcluster
    ENDDO

    PRINT *, "============= ratio ================", ratiomoy, ratiorij

    DEALLOCATE(Frob)
! sparsification fin

#if aff
    PRINT *,numproc,'nbinfo=', nbinfo,' nbcluster=',nbcluster
#endif

    RETURN 
  END SUBROUTINE sp_spectral_embedding

  SUBROUTINE sp_matvec(A, IA, JA, X, Y, n, nnz)
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN), DIMENSION(nnz) :: A
  INTEGER, INTENT(IN), DIMENSION(nnz) :: IA, JA
  DOUBLE PRECISION, INTENT(IN), DIMENSION(n) :: X
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(n) :: Y
  INTEGER, INTENT(IN) :: n, nnz

  INTEGER :: l

  Y(:) = dfloat(0)

  DO l = 1, nnz
    Y(IA(l)) = Y(IA(l)) + A(l)*X(JA(l))
  ENDDO

  RETURN

  END SUBROUTINE sp_matvec

  SUBROUTINE solve_arpack(A, IA, JA, ndim, nnz, nblimit, W, Z)

  DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: A
  INTEGER, INTENT(IN), DIMENSION(:) :: IA, JA
  INTEGER, INTENT(IN) :: ndim, nnz, nblimit

  DOUBLE PRECISION, INTENT(OUT), POINTER :: W(:)
  DOUBLE PRECISION, INTENT(OUT), POINTER :: Z(:,:)
  
  INTEGER :: maxn, maxnev, maxncv, ldv, i, nbite
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
  INTEGER           iparam(11), ipntr(14)
  LOGICAL, DIMENSION(:), ALLOCATABLE :: SELECT

!     d valeurs propres
!     v vecteurs propres

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ax, resid, workd, workev, workl
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: d, v
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
  CHARACTER         bmat*1, which*2
  INTEGER           ido, n, nx, nev, ncv, lworkl, info, ierr, &
                    j, ishfts, maxitr, mode1, nconv
  DOUBLE PRECISION  tol, sigmar, sigmai
  LOGICAL           first, rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
  DOUBLE PRECISION  zero
  PARAMETER         (zero = 0.0D+0)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      DOUBLE PRECISION  dlapy2, dnrm2
      EXTERNAL          dlapy2, dnrm2, daxpy 
!
!     %--------------------%
!     | Intrinsic FUNCTION |
!     %--------------------%
!
      INTRINSIC         abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------------------------%
!     | The following INCLUDE statement and assignments |
!     | initiate trace output from the internal         |
!     | actions of ARPACK.  See debug.doc in the        |
!     | DOCUMENTS directory for usage.  Initially, the  |
!     | most useful information will be a breakdown of  |
!     | time spent in the various stages of computation |
!     | given by setting mnaupd = 1.                    |
!     %-------------------------------------------------%
!
      INCLUDE 'debug.h'

! lien entre les tailles
      maxn = ndim
      ldv = maxn
      maxnev = nblimit
      maxncv = 2*maxnev + 1

! allocation memoire
      ALLOCATE(SELECT(maxn))
      ALLOCATE(ax(maxn), resid(maxn), workd(3*maxn), &
               workev(3*maxncv), workl(3*maxncv*maxncv+6*maxncv))
      ALLOCATE(d(maxncv, 3), v(ldv, maxncv))

      ndigit = -3
      logfil = 6
      mnaitr = 0
      mnapps = 0
      mnaupd = 1
      mnaup2 = 0
      mneigh = 0
      mneupd = 0
!
!     %-------------------------------------------------%
!     | The following sets dimensions for this problem. |
!     %-------------------------------------------------%
!
      n     = ndim

!     %-----------------------------------------------%
!     |                                               |
!     | Specifications for ARPACK usage are set       |
!     | below:                                        |
!     |                                               |
!     |    1) NEV = 4  asks for 4 eigenvalues to be   |
!     |       computed.                               |
!     |                                               |
!     |    2) NCV = 20 sets the length of the Arnoldi |
!     |       factorization.                          |
!     |                                               |
!     |    3) This is a standard problem.             |
!     |         (indicated by bmat  = 'I')            |
!     |                                               |
!     |    4) Ask for the NEV eigenvalues of          |
!     |       largest magnitude.                      |
!     |         (indicated by which = 'LM')           |
!     |       See DOcumentation in DNAUPD for the     |
!     |       other options SM, LR, SR, LI, SI.       |
!     |                                               |
!     | Note: NEV and NCV must satisfy the following  |
!     | conditions:                                   |
!     |              NEV <= MAXNEV                    |
!     |          NEV + 2 <= NCV <= MAXNCV             |
!     |                                               |
!     %-----------------------------------------------%
!
      nev   = maxnev
      ncv   = maxncv
      bmat  = 'I'
      which = 'LR'
!
      IF ( n .GT. maxn ) THEN
         PRINT *, ' ERROR with _NSIMP: N is greater than MAXN '
         GOTO 9000
      ELSEIF ( nev .GT. maxnev ) THEN
         PRINT *, ' ERROR with _NSIMP: NEV is greater than MAXNEV '
         GOTO 9000
      ELSEIF ( ncv .GT. maxncv ) THEN
         PRINT *, ' ERROR with _NSIMP: NCV is greater than MAXNCV '
         GOTO 9000
      ENDIF
!
!     %-----------------------------------------------------%
!     |                                                     |
!     | Specification of stopping rules and initial         |
!     | conditions before calling DNAUPD                    |
!     |                                                     |
!     | TOL  determines the stopping criterion.             |
!     |                                                     |
!     |      Expect                                         |
!     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
!     |               computed   true                       |
!     |                                                     |
!     |      If TOL .LE. 0,  THEN TOL <- macheps            |
!     |           (machine precision) is used.              |
!     |                                                     |
!     | IDO  is the REVERSE COMMUNICATION PARAMETER         |
!     |      used to specify actions to be taken on return  |
!     |      from DNAUPD. (see usage below)                 |
!     |                                                     |
!     |      It MUST initially be set to 0 before the first |
!     |      call to DNAUPD.                                |
!     |                                                     |
!     | INFO on entry specifies starting vector information |
!     |      and on return indicates error codes            |
!     |                                                     |
!     |      Initially, setting INFO=0 indicates that a     |
!     |      random starting vector is requested to         |
!     |      start the ARNOLDI iteration.  Setting INFO to  |
!     |      a nonzero value on the initial call is used    |
!     |      if you want to specify your own starting       |
!     |      vector (This vector must be placed in RESID).  |
!     |                                                     |
!     | The work array WORKL is used in DNAUPD as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.                                  |
!     |                                                     |
!     %-----------------------------------------------------%
!
      lworkl  = 3*ncv**2+6*ncv 
      tol    = 1.D-6
      ido    = 0
      info   = 0
!
!     %---------------------------------------------------%
!     | Specification of Algorithm Mode:                  |
!     |                                                   |
!     | This program uses the exact shift strategy        |
!     | (indicated by setting IPARAM(1) = 1).             |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the DOcumentation in |
!     | DNAUPD.                                           |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode1 = 1
!
      iparam(1) = ishfts
!
      iparam(3) = maxitr
!
      iparam(7) = mode1
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) | 
!     %-------------------------------------------%
!
      nbite = 1
 10   CONTINUE
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         CALL dnaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
                       v, ldv, iparam, ipntr, workd, workl, lworkl, & 
                       info )
!
         IF (ido .EQ. -1 .OR. ido .EQ. 1) THEN
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- Op*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and RETURN the matrix vector      |
!           | product to workd(ipntr(2)).               | 
!           %-------------------------------------------%
!
            CALL sp_matvec(A, IA, JA, workd(ipntr(1)), workd(ipntr(2)), &
                           ndim, nnz)

            nbite = nbite + 1
!
!           %-----------------------------------------%
!           | L O O P   B A C K to CALL DNAUPD again. |
!           %-----------------------------------------%
!
            GOTO 10
!
         ENDIF
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      IF ( info .LT. 0 ) THEN
!
!        %--------------------------%
!        | Error message, check the |
!        | DOcumentation in DNAUPD. |
!        %--------------------------%
!
         PRINT *, ' '
         PRINT *, ' Error with _naupd, info = ',info
         PRINT *, ' Check the DOcumentation of _naupd'
         PRINT *, ' '
!
      ELSE 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may be also computed now if  |
!        | desired.  (indicated by rvec = .TRUE.)    |
!        |                                           |
!        | The routine DNEUPD now CALLed to DO this  |
!        | post processing (Other modes may require  |
!        | more complicated post processing than     |
!        | mode1,)                                   |
!        |                                           |
!        %-------------------------------------------%
!
         rvec = .TRUE.
!
         CALL dneupd ( rvec, 'A', SELECT, d, d(1,2), v, ldv, &
              sigmar, sigmai, workev, bmat, n, which, nev, tol, & 
              resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
              lworkl, ierr )
!
!        %------------------------------------------------%
!        | The REAL parts of the eigenvalues are returned |
!        | in the first column of the two dimensional     |
!        | array D, and the IMAGINARY part are returned   |
!        | in the second column of D.  The corresponding  |
!        | eigenvectors are returned in the first         |
!        | NCONV (= IPARAM(5)) columns of the two         |
!        | dimensional array V if requested.  Otherwise,  |
!        | an orthogonal basis for the invariant subspace |
!        | corresponding to the eigenvalues in D is       |
!        | returned in V.                                 |
!        %------------------------------------------------%
!
         IF ( ierr .NE. 0) THEN
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DNEUPD. |
!           %------------------------------------%
!
            PRINT *, ' '
            PRINT *, ' Error with _neupd, info = ', ierr
            PRINT *, ' Check the documentation of _neupd. '
            PRINT *, ' '
!
         ELSE
!
            first = .TRUE.
            nconv =  iparam(5)
            DO 20 j=1, nconv
!
!              %---------------------------%
!              | Compute the residual norm |
!              |                           |
!              |   ||  A*x - lambda*x ||   |
!              |                           |
!              | for the NCONV accurately  |
!              | computed eigenvalues and  |
!              | eigenvectors.  (IPARAM(5) |
!              | indicates how many are    |
!              | accurate to the requested |
!              | tolerance)                |
!              %---------------------------%
!
               IF (d(j,2) .EQ. zero)  THEN
!
!                 %--------------------%
!                 | Ritz value is REAL |
!                 %--------------------%
!
                  !CALL av(nx, v(1,j), ax)
                  CALL sp_matvec(A, IA, JA, v(1,j), ax, ndim, nnz)
                  CALL daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                  d(j,3) = d(j,3) / abs(d(j,1))
!
               ELSEIF (first) THEN
!
!                 %------------------------%
!                 | Ritz value is COMPLEX. |
!                 | Residual of one Ritz   |
!                 | value of the conjugate |
!                 | pair is computed.      |
!                 %------------------------%
!
                  !CALL av(nx, v(1,j), ax)
                  CALL sp_matvec(A, IA, JA, v(1,j), ax, ndim, nnz)
                  CALL daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  CALL daxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                  !CALL av(nx, v(1,j+1), ax)
                  CALL sp_matvec(A, IA, JA, v(1,j+1), ax, ndim, nnz)
                  CALL daxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                  CALL daxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                  d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
                  d(j,3) = d(j,3) / dlapy2(d(j,1),d(j,2))
                  d(j+1,3) = d(j,3)
                  first = .FALSE.
               ELSE
                  first = .TRUE.
               ENDIF
!
 20         CONTINUE
!
!           %-----------------------------%
!           | Display computed residuals. |
!           %-----------------------------%
!
            CALL dmout(6, nconv, 3, d, maxncv, -6, &
                 'Ritz values (Real, Imag) and residual residuals')
         ENDIF
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
         IF ( info .EQ. 1) THEN
             PRINT *, ' '
             PRINT *, ' Maximum number of iterations reached.'
             PRINT *, ' '
         ELSEIF ( info .EQ. 3) THEN
             PRINT *, ' ' 
             PRINT *, ' No shifts could be applied during IMPLICIT &
                        Arnoldi update, try increasing NCV.'
             PRINT *, ' '
         ENDIF      
!
         PRINT *, ' '
         PRINT *, ' _NSIMP '
         PRINT *, ' ====== '
         PRINT *, ' '
         PRINT *, ' Size of the matrix is ', n
         PRINT *, ' The number of Ritz values requested is ', nev
         PRINT *, ' The number of Arnoldi vectors generated', &
                  ' (NCV) is ', ncv
         PRINT *, ' What portion of the spectrum: ', which
         PRINT *, ' The number of converged Ritz values is ', &
                    nconv 
         PRINT *, ' The number of Implicit Arnoldi update', &
                  ' iterations taken is ', iparam(3)
         PRINT *, ' The number of OP*x is ', iparam(9)
         PRINT *, ' The convergence criterion is ', tol
         PRINT *, ' '
!
      ENDIF
!
!     %---------------------------%
!     | Done with program dnsimp. |
!     %---------------------------%
!
 9000 CONTINUE

      ALLOCATE(W(nblimit))
      ALLOCATE(Z(n, nblimit))

      W(1:nblimit) = d(1:nblimit, 1)

      DO i = 1, nblimit
        Z(:,i) = v(:,i)
      ENDDO

      DEALLOCATE(SELECT)
      DEALLOCATE(ax, resid, workd, workev, workl)
      DEALLOCATE(d, v)

  END SUBROUTINE solve_arpack

  END MODULE module_sparse
