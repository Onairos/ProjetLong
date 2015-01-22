 module module_solve
contains

  !************************************************
  !DGEEVX
  subroutine solve_dgeevx(N,A,VR,WR)
    implicit none
    external DGEEVX
    integer :: N,LWORK,INFO,LDA,LDVL,LDVR,ILO,IHI
    integer,dimension(:),pointer :: IWORK
    DOUBLE PRECISION,dimension(:,:),pointer :: A,VL,VR
    double precision,dimension(:),pointer :: WORK,WR,WI,SCALE,RCONDE,RCONDV
    DOUBLE PRECISION :: ABNRM
    CHARACTER :: JOBVL,JOBVR,BALANC,SENSE
    LDA = N
    LDVL=N
    LDVR=N
    LWORK=4*N
    INFO=0
    JOBVL='N'
    JOBVR='V'
    BALANC='N'
    SENSE='B'
    allocate(WORK(LWORK)); WORK(:)=0.0
    allocate(WR(N)); WR(:)=0.0
    allocate(WI(N)); WI(:)=0.0
    allocate(VL(LDVL,N)); VL(:,:)=0.0
    allocate(VR(LDVR,N)); VR(:,:)=0.0
    call DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI,&
         VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM,&
         RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )

    deallocate(WORK)
    deallocate(WI)
    deallocate(VL)
    if (INFO/=0) print *,'erreur dans DSYERVR ? INFO=',INFO
    return
  end subroutine solve_dgeevx


  !************************************************
  !DGEEV
  subroutine solve_dgeev(N,A,VR,WR)
    implicit none
    external DGEEV
    integer :: N,LWORK,INFO,LDA,LDVL,LDVR
    DOUBLE PRECISION,dimension(:,:),pointer :: A,VL,VR
    double precision,dimension(:),pointer :: WORK,WR,WI
    CHARACTER :: JOBVL,JOBVR
    LDA = N
    LDVL=N
    LDVR=N
    LWORK=4*N
    INFO=0
    JOBVL='N'
    JOBVR='V'
    allocate(WORK(LWORK)); WORK(:)=0.0
    allocate(WR(N)); WR(:)=0.0
    allocate(WI(N)); WI(:)=0.0
    allocate(VL(LDVL,N)); VL(:,:)=0.0
    allocate(VR(LDVR,N)); VR(:,:)=0.0
    call DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,&
         LDVR, WORK, LWORK, INFO )
    deallocate(WORK)
    deallocate(WI)
    deallocate(VL)
    if (INFO/=0) print *,'erreur dans DSYERVR ? INFO=',INFO
    return
  end subroutine solve_dgeev

  !************************************************
  !DSYEV
  subroutine solve_dsyev(N,A,W,LWORK)
    implicit none
    external DSYEV
    integer :: N,LWORK,INFO,LDA
    DOUBLE PRECISION :: A(N,N),WORK(LWORK),W(N)
    CHARACTER :: JOBZ,UPLO
    JOBZ = 'V'
    INFO=0
    UPLO = 'U'
    LDA = N
    call DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
    if (INFO/=0) print *,'erreur dans DSYERVR ? INFO=',INFO
    return
  end subroutine solve_dsyev

  !**********************************************
  !DSYEVR
  subroutine solve_dsyevr(k,N,A,Z,LWORK,LIWORK,W,M)
    implicit none
    external DSYEVR
    integer :: N,LWORK,LIWORK
    double precision :: A(N,N), W(N),Z(N,N)
    INTEGER :: LDA,IL,IU,M,INFO,LDZ,ISUPPZ(2*N),IWORK(LIWORK)
    DOUBLE PRECISION :: ABSTOL,VL,VU,WORK(LWORK)
    CHARACTER :: JOBZ,RANGE,UPLO
    integer :: i,j,k
    JOBZ = 'V'
    INFO=0
    RANGE = 'A'
    UPLO = 'U'
    LDA = N
    VL = -100
    VU = 5
    ABSTOL = 0.D+0
    IL = N-k+1
    IU = N
    LDZ = N
    !on enleve la partie triangulaire inferieure
    do i=1,N
       do j=1,i-1
          A(i,j)=0.0
       end do
    end do
    call DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL,&
         M, W, Z, LDZ, ISUPPZ, WORK, LWORK,IWORK, LIWORK, INFO )
    LWORK = WORK(1)
    LIWORK = IWORK(1)
    if (INFO/=0) print *,'erreur dans DSYERVR ? INFO=',INFO
   ! print *,'nb de vp trouvees :',M
    !do i=1,M
    !   print *,'vp trouvees :', W(i)!,' vect ass :',Z(:,i)
    !enddo
  end subroutine solve_dsyevr

  !***********************************************
  !DSYEVX
  subroutine solve_dsyevx(k,N,A,Z,LWORK,LIWORK,W,M)
    implicit none
    external DSYEVX
    integer :: N,LWORK,LIWORK,IFAIL
    double precision :: A(N,N), W(N),Z(N,N)
    INTEGER :: LDA,IL,IU,M,INFO,LDZ,IWORK(LIWORK)
    DOUBLE PRECISION :: ABSTOL,VL,VU,WORK(LWORK)
    CHARACTER :: JOBZ,RANGE,UPLO
    integer :: i,j,k
    JOBZ = 'V'
    INFO=0
    RANGE = 'A'
    UPLO = 'U'
    LDA = N
    VL = -100
    VU = 5
    ABSTOL = 1.e-12 !0.0
    IL = N-k+1
    IU = N
    LDZ = N
    IFAIL=0
    !on enleve la partie triangulaire inferieure
    do i=1,N
       do j=1,i-1
          A(i,j)=0.0
       end do
    end do
    call DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,&
         ABSTOL, M, W, Z, LDZ, WORK, LWORK,IWORK, IFAIL, INFO )
    LWORK = WORK(1)
    LIWORK = IWORK(1)
    if (INFO/=0) print *,'erreur dans DSYERVR ? INFO=',INFO
   ! print *,'nb de vp trouvees :',M
    !do i=1,M
    !   print *,'vp trouvees :', W(i)!,' vect ass :',Z(:,i)
    !enddo
  end subroutine solve_dsyevx

  !***********************************************
  !DGEHRD
!!$  subroutine solve_dgehrd(k,N,A,Z,LWORK,M)
!!$    implicit none
!!$    external DGEHRD
!!$    integer :: N, ILO, IHI
!!$    double precision :: A(N,N),TAU(N-1),WORK(LWORK)
!!$
!!$
!!$    INTEGER :: LDA,IL,IU,M,INFO,LDZ,ISUPPZ(2*N),IWORK(LIWORK)
!!$    DOUBLE PRECISION :: ABSTOL,VL,VU
!!$    CHARACTER :: JOBZ,RANGE,UPLO
!!$    integer :: i,j,k
!!$    INFO=0
!!$    ILO=1
!!$    IHI=N
!!$    LDA=N
!!$    TAU(:)=0.0
!!$    !on enleve la partie triangulaire inferieure
!!$    do i=1,N
!!$       do j=1,i-1
!!$          A(i,j)=0.0
!!$       end do
!!$    end do
!!$    call DGEHRD( N, ILO,IHI,A,LDA,TAU,WORK,LWORK,INFO)
!!$    LWORK = WORK(1)
!!$    LIWORK = IWORK(1)
!!$    if (INFO/=0) print *,'erreur dans DSYERVR ? INFO=',INFO
!!$   ! print *,'nb de vp trouvees :',M
!!$    !do i=1,M
!!$    !   print *,'vp trouvees :', W(i)!,' vect ass :',Z(:,i)
!!$    !enddo
!!$  end subroutine solve_dgehrd

! Subroutine ARPACK



  subroutine solve_arpack_full(A, ndim, nblimit, W, Z)

  double precision, intent(in), dimension(:,:) :: A
  integer, intent(in) :: ndim, nblimit

  double precision, intent(out), pointer :: W(:)
  double precision, intent(out), pointer :: Z(:, :)

  integer :: maxn, maxnev, maxncv, ldv, i, nbite
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
  integer           iparam(11), ipntr(14)
  logical, dimension(:), allocatable :: select

!     d valeurs propres
!     v vecteurs propres

  Double precision, dimension(:), allocatable :: ax, resid, workd, workev, workl
  Double precision, dimension(:,:), allocatable :: d, v
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
  character         bmat*1, which*2
  integer           ido, n, nx, nev, ncv, lworkl, info, ierr, &
                    j, ishfts, maxitr, mode1, nconv
  Double precision  tol, sigmar, sigmai
  logical           first, rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
  Double precision  zero
  parameter         (zero = 0.0D+0)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Double precision  dlapy2, dnrm2
      external          dlapy2, dnrm2, daxpy
!
!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
      intrinsic         abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------------------------%
!     | The following include statement and assignments |
!     | initiate trace output from the internal         |
!     | actions of ARPACK.  See debug.doc in the        |
!     | DOCUMENTS directory for usage.  Initially, the  |
!     | most useful information will be a breakdown of  |
!     | time spent in the various stages of computation |
!     | given by setting mnaupd = 1.0                    |
!     %-------------------------------------------------%
!
      include 'debug.h'

! lien entre les tailles
      maxn = ndim
      ldv = maxn
      maxnev = nblimit
      maxncv = 2*maxnev + 1

! allocation memoire
      allocate(select(maxn))
      allocate(ax(maxn), resid(maxn), workd(3*maxn), &
               workev(3*maxncv), workl(3*maxncv*maxncv+6*maxncv))
      allocate(d(maxncv, 3), v(ldv, maxncv))

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
!     |       See documentation in DNAUPD for the     |
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
      print *, "nev ***********", nev, ncv, n
      bmat  = 'I'
      which = 'LR'
!
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NSIMP: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NSIMP: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NSIMP: NCV is greater than MAXNCV '
         go to 9000
      end if
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
!     |      If TOL .le. 0,  then TOL <- macheps            |
!     |           (machine precision) is used.              |
!     |                                                     |
!     | IDO  is the REVERSE COMMUNICATION parameter         |
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
!     | by the user. For details see the documentation in |
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
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call dnaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
                       v, ldv, iparam, ipntr, workd, workl, lworkl, &
                       info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- Op*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and return the matrix vector      |
!           | product to workd(ipntr(2)).               |
!           %-------------------------------------------%
!
            !call av (nx, workd(ipntr(1)), workd(ipntr(2)))

            call dgemv('N', ndim, ndim, 1.D0, A, ndim, workd(ipntr(1)), 1, 0.D0, workd(ipntr(2)), 1)

            nbite = nbite + 1
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DNAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         endif
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in DNAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _naupd, info = ',info
         print *, ' Check the documentation of _naupd'
         print *, ' '
!
      else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may be also computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        |                                           |
!        | The routine DNEUPD now called to do this  |
!        | post processing (Other modes may require  |
!        | more complicated post processing than     |
!        | mode1,)                                   |
!        |                                           |
!        %-------------------------------------------%
!
         rvec = .true.
!
         call dneupd ( rvec, 'A', select, d, d(1,2), v, ldv, &
              sigmar, sigmai, workev, bmat, n, which, nev, tol, &
              resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
              lworkl, ierr )
!
!        %------------------------------------------------%
!        | The real parts of the eigenvalues are returned |
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
         if ( ierr .ne. 0) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DNEUPD. |
!           %------------------------------------%
!
            print *, ' '
            print *, ' Error with _neupd, info = ', ierr
            print *, ' Check the documentation of _neupd. '
            print *, ' '
!
         else
!
            first = .true.
            nconv =  iparam(5)
            do 20 j=1, nconv
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
               if (d(j,2) .eq. zero)  then
!
!                 %--------------------%
!                 | Ritz value is real |
!                 %--------------------%
!
                  !call av(nx, v(1,j), ax)

                  call dgemv('N', ndim, ndim, 1.D0, A, ndim, v(1,j), 1, 0.D0, ax, 1)

                  call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                  d(j,3) = d(j,3) / abs(d(j,1))
!
               else if (first) then
!
!                 %------------------------%
!                 | Ritz value is complex. |
!                 | Residual of one Ritz   |
!                 | value of the conjugate |
!                 | pair is computed.      |
!                 %------------------------%
!
                  !call av(nx, v(1,j), ax)

                  call dgemv('N', ndim, ndim, 1.D0, A, ndim, v(1,j), 1, 0.D0, ax, 1)

                  call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  call daxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                  !call av(nx, v(1,j+1), ax)

                  call dgemv('N', ndim, ndim, 1.D0, A, ndim, v(1,j+1), 1, 0.D0, ax, 1)
                  call daxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                  call daxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                  d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
                  d(j,3) = d(j,3) / dlapy2(d(j,1),d(j,2))
                  d(j+1,3) = d(j,3)
                  first = .false.
               else
                  first = .true.
               end if
!
 20         continue
!
!           %-----------------------------%
!           | Display computed residuals. |
!           %-----------------------------%
!
            call dmout(6, nconv, 3, d, maxncv, -6, &
                 'Ritz values (Real, Imag) and residual residuals')
         end if
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' '
             print *, ' No shifts could be applied during implicit &
                        Arnoldi update, try increasing NCV.'
             print *, ' '
         end if
!
         print *, ' '
         print *, ' _NSIMP '
         print *, ' ====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated', &
                  ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', &
                    nconv
         print *, ' The number of Implicit Arnoldi update', &
                  ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
!
      end if
!
!     %---------------------------%
!     | Done with program dnsimp. |
!     %---------------------------%
!
 9000 continue

      allocate(W(nblimit))
      allocate(Z(n, nblimit))

      W(1:nblimit) = d(1:nblimit, 1)

      do i = 1, nblimit
        Z(:,i) = v(:,i)
      end do

      deallocate(select)
      deallocate(ax, resid, workd, workev, workl)
      deallocate(d, v)

  end subroutine solve_arpack_full

end module module_solve
