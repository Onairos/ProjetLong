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





end module module_solve
