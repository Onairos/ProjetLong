 MODULE module_solve
CONTAINS

  !************************************************
  !DGEEVX
  SUBROUTINE solve_dgeevx(N, A, VR, WR)
    IMPLICIT NONE
    EXTERNAL DGEEVX
    INTEGER :: N, LWORK, INFO, LDA, LDVL, LDVR, ILO, IHI
    INTEGER, DIMENSION(:), POINTER :: IWORK
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A, VL, VR
    DOUBLE PRECISION, DIMENSION(:), POINTER :: WORK, WR, WI, SCALE, RCONDE, RCONDV
    DOUBLE PRECISION :: ABNRM
    CHARACTER :: JOBVL, JOBVR, BALANC, SENSE

    LDA = N
    LDVL=N
    LDVR=N
    LWORK=4*N
    INFO=0
    JOBVL='N'
    JOBVR='V'
    BALANC='N'
    SENSE='B'
    ALLOCATE(WORK(LWORK)); WORK(:)=0.0
    ALLOCATE(WR(N)); WR(:)=0.0
    ALLOCATE(WI(N)); WI(:)=0.0
    ALLOCATE(VL(LDVL,N)); VL(:,:)=0.0
    ALLOCATE(VR(LDVR,N)); VR(:,:)=0.0
    CALL DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI,&
         VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM,&
         RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )

    DEALLOCATE(WORK)
    DEALLOCATE(WI)
    DEALLOCATE(VL)
    IF (INFO/=0) PRINT *,'erreur dans DSYERVR ? INFO=',INFO
    RETURN
  END SUBROUTINE solve_dgeevx


  !************************************************
  !DGEEV
  SUBROUTINE solve_dgeev(N, A, VR, WR)
    IMPLICIT NONE
    EXTERNAL DGEEV
    INTEGER :: N, LWORK, INFO, LDA, LDVL, LDVR
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A, VL, VR
    DOUBLE PRECISION, DIMENSION(:), POINTER :: WORK, WR, WI
    CHARACTER :: JOBVL, JOBVR

    LDA = N
    LDVL=N
    LDVR=N
    LWORK=4*N
    INFO=0
    JOBVL='N'
    JOBVR='V'
    ALLOCATE(WORK(LWORK)); WORK(:)=0.0
    ALLOCATE(WR(N)); WR(:)=0.0
    ALLOCATE(WI(N)); WI(:)=0.0
    ALLOCATE(VL(LDVL,N)); VL(:,:)=0.0
    ALLOCATE(VR(LDVR,N)); VR(:,:)=0.0
    CALL DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,&
         LDVR, WORK, LWORK, INFO )
    DEALLOCATE(WORK)
    DEALLOCATE(WI)
    DEALLOCATE(VL)
    IF (INFO/=0) PRINT *,'erreur dans DSYERVR ? INFO=',INFO
    RETURN
  END SUBROUTINE solve_dgeev

  !************************************************
  !DSYEV
  SUBROUTINE solve_dsyev(N, A, W, LWORK)
    IMPLICIT NONE
    EXTERNAL DSYEV
    INTEGER :: N, LWORK, INFO, LDA
    DOUBLE PRECISION :: A(N,N), WORK(LWORK), W(N)
    CHARACTER :: JOBZ, UPLO

    JOBZ = 'V'
    INFO=0
    UPLO = 'U'
    LDA = N
    CALL DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
    IF (INFO/=0) PRINT *,'erreur dans DSYERVR ? INFO=',INFO
    RETURN
  END SUBROUTINE solve_dsyev

  !**********************************************
  !DSYEVR
  SUBROUTINE solve_dsyevr(k, N, A, Z, LWORK, LIWORK, W, M)
    IMPLICIT NONE
    EXTERNAL DSYEVR
    INTEGER :: N, LWORK, LIWORK
    DOUBLE PRECISION :: A(N,N), W(N), Z(N,N)
    INTEGER :: LDA, IL, IU, M, INFO, LDZ, ISUPPZ(2*N), IWORK(LIWORK)
    DOUBLE PRECISION :: ABSTOL, VL, VU, WORK(LWORK)
    CHARACTER :: JOBZ, RANGE, UPLO
    INTEGER :: i, j, k

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
    DO i=1,N
       DO j=1,i-1
          A(i,j)=0.0
       ENDDO
    ENDDO
    CALL DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL,&
         M, W, Z, LDZ, ISUPPZ, WORK, LWORK,IWORK, LIWORK, INFO )
    LWORK = WORK(1) !TODO : Warning: Possible change of value in conversion from REAL(8) to INTEGER(4) at (1)
    LIWORK = IWORK(1)
    IF (INFO/=0) PRINT *,'erreur dans DSYERVR ? INFO=',INFO
  END SUBROUTINE solve_dsyevr

  !***********************************************
  !DSYEVX
  SUBROUTINE solve_dsyevx(k, N, A, Z, LWORK, LIWORK, W, M)
    IMPLICIT NONE
    EXTERNAL DSYEVX
    INTEGER :: N, LWORK, LIWORK, IFAIL
    DOUBLE PRECISION :: A(N,N), W(N), Z(N,N)
    INTEGER :: LDA, IL, IU, M, INFO, LDZ, IWORK(LIWORK)
    DOUBLE PRECISION :: ABSTOL, VL, VU, WORK(LWORK)
    CHARACTER :: JOBZ, RANGE, UPLO
    INTEGER :: i, j, k

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
    DO i=1,N
       DO j=1,i-1
          A(i,j)=0.0
       ENDDO
    ENDDO
    CALL DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,&
         ABSTOL, M, W, Z, LDZ, WORK, LWORK,IWORK, IFAIL, INFO )
    LWORK = WORK(1) !TODO : Warning: Possible change of value in conversion from REAL(8) to INTEGER(4) at (1)
    LIWORK = IWORK(1)
    IF (INFO/=0) PRINT *,'erreur dans DSYERVR ? INFO=',INFO
  END SUBROUTINE solve_dsyevx

END MODULE module_solve
