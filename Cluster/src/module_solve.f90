 MODULE module_solve
CONTAINS



  SUBROUTINE solve_dgeevx(N, A, VR, WR)
    IMPLICIT NONE
    EXTERNAL DGEEVX

    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A
    INTEGER :: N

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: VR
    DOUBLE PRECISION, DIMENSION(:), POINTER :: WR

    !#### Variables  ####
    CHARACTER :: BALANC
    CHARACTER :: JOBVL
    CHARACTER :: JOBVR
    CHARACTER :: SENSE
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: VL
    DOUBLE PRECISION, DIMENSION(:), POINTER :: RCONDE
    DOUBLE PRECISION, DIMENSION(:), POINTER :: RCONDV
    DOUBLE PRECISION, DIMENSION(:), POINTER :: SCALE
    DOUBLE PRECISION, DIMENSION(:), POINTER :: WI
    DOUBLE PRECISION, DIMENSION(:), POINTER :: WORK
    DOUBLE PRECISION :: ABNRM
    INTEGER, DIMENSION(:), POINTER :: IWORK
    INTEGER :: IHI
    INTEGER :: ILO
    INTEGER :: INFO
    INTEGER :: LDA
    INTEGER :: LDVL
    INTEGER :: LDVR
    INTEGER :: LWORK

    !###########################################      
    ! INSTRUCTIONS
    !###########################################      
    LDA = N
    LDVL=N
    LDVR=N
    LWORK=4*N
    INFO=0
    JOBVL='N'
    JOBVR='V'
    BALANC='N'
    SENSE='B'
    ALLOCATE(WORK(LWORK))
    WORK(:)=0.0
    ALLOCATE(WR(N))
    WR(:)=0.0
    ALLOCATE(WI(N))
    WI(:)=0.0
    ALLOCATE(VL(LDVL,N))
    VL(:,:)=0.0
    ALLOCATE(VR(LDVR,N))
    VR(:,:)=0.0
    CALL DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI,&
         VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM,&
         RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )

    DEALLOCATE(WORK)
    DEALLOCATE(WI)
    DEALLOCATE(VL)
    IF (INFO/=0) PRINT *, 'Error in DSYERVR ? INFO=', INFO
    RETURN
  END SUBROUTINE solve_dgeevx




  SUBROUTINE solve_dgeev(N, A, VR, WR)
    IMPLICIT NONE
    EXTERNAL DGEEV
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: N
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: VR
    DOUBLE PRECISION, DIMENSION(:), POINTER :: WR

    !#### Variables  ####
    CHARACTER :: JOBVL
    CHARACTER :: JOBVR
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: VL
    DOUBLE PRECISION, DIMENSION(:), POINTER :: WI
    DOUBLE PRECISION, DIMENSION(:), POINTER :: WORK
    INTEGER :: INFO
    INTEGER :: LDA
    INTEGER :: LDVL
    INTEGER :: LDVR
    INTEGER :: LWORK

    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
    LDA = N
    LDVL=N
    LDVR=N
    LWORK=4*N
    INFO=0
    JOBVL='N'
    JOBVR='V'
    ALLOCATE(WORK(LWORK))
    WORK(:)=0.0
    ALLOCATE(WR(N))
    WR(:)=0.0
    ALLOCATE(WI(N))
    WI(:)=0.0
    ALLOCATE(VL(LDVL,N))
    VL(:,:)=0.0
    ALLOCATE(VR(LDVR,N))
    VR(:,:)=0.0
    CALL DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,&
         LDVR, WORK, LWORK, INFO )
    DEALLOCATE(WORK)
    DEALLOCATE(WI)
    DEALLOCATE(VL)
    IF (INFO/=0) PRINT *, 'Error in DSYERVR ? INFO=',INFO
    RETURN
  END SUBROUTINE solve_dgeev



  SUBROUTINE solve_dsyev(N, A, W, LWORK)
    IMPLICIT NONE
    EXTERNAL DSYEV
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    DOUBLE PRECISION :: A(N,N)
    INTEGER :: N

    !====  OUT ====
    DOUBLE PRECISION :: W(N)
    INTEGER :: LWORK

    !#### Variables  ####
    CHARACTER :: JOBZ
    CHARACTER :: UPLO
    DOUBLE PRECISION :: WORK(LWORK)
    INTEGER :: INFO
    INTEGER :: LDA

    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
    JOBZ = 'V'
    INFO=0
    UPLO = 'U'
    LDA = N
    CALL DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
    IF (INFO/=0) PRINT *, 'Error in DSYERVR ? INFO=',INFO
    RETURN
  END SUBROUTINE solve_dsyev



  SUBROUTINE solve_dsyevr(k, N, A, Z, LWORK, LIWORK, W, M)
    IMPLICIT NONE
    EXTERNAL DSYEVR
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    DOUBLE PRECISION :: A(N,N)
    DOUBLE PRECISION :: Z(N,N)
    DOUBLE PRECISION :: W(N)
    INTEGER :: k
    INTEGER :: M
    INTEGER :: N
    INTEGER :: LIWORK
    INTEGER :: LWORK

    !#### Variables  ####
    CHARACTER :: JOBZ
    CHARACTER :: RANGE
    CHARACTER :: UPLO
    DOUBLE PRECISION :: WORK(LWORK)
    DOUBLE PRECISION :: ABSTOL
    DOUBLE PRECISION :: VL
    DOUBLE PRECISION :: VU
    INTEGER :: ISUPPZ(2*N)
    INTEGER :: IWORK(LIWORK)
    INTEGER :: INFO
    INTEGER :: i
    INTEGER :: IL
    INTEGER :: IU
    INTEGER :: j
    INTEGER :: LDA
    INTEGER :: LDZ

    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
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
    ! Removing lower triangular part
    DO i=1,N
       DO j=1,i-1
          A(i,j)=0.0
       ENDDO
    ENDDO
    CALL DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL,&
         M, W, Z, LDZ, ISUPPZ, WORK, LWORK,IWORK, LIWORK, INFO )
    LWORK = WORK(1)
    LIWORK = IWORK(1)
    IF (INFO/=0) PRINT *, 'Error in DSYERVR ? INFO=',INFO
  END SUBROUTINE solve_dsyevr



  SUBROUTINE solve_dsyevx(k, N, A, Z, LWORK, LIWORK, W, M)
    IMPLICIT NONE
    EXTERNAL DSYEVX
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !=== IN/OUT ===
    DOUBLE PRECISION :: A(N,N)
    DOUBLE PRECISION :: Z(N,N)
    DOUBLE PRECISION :: W(N)
    INTEGER :: LIWORK
    INTEGER :: LWORK
    INTEGER :: k
    INTEGER :: M
    INTEGER :: N

    !#### Variables  ####
    CHARACTER :: JOBZ
    CHARACTER :: UPLO
    CHARACTER :: RANGE
    DOUBLE PRECISION :: WORK(LWORK)
    DOUBLE PRECISION :: ABSTOL
    DOUBLE PRECISION :: VL
    DOUBLE PRECISION :: VU
    INTEGER :: i
    INTEGER :: IFAIL
    INTEGER :: IL
    INTEGER :: INFO
    INTEGER :: IU
    INTEGER :: IWORK(LIWORK)
    INTEGER :: j
    INTEGER :: LDA
    INTEGER :: LDZ

    !###########################################      
    ! INSTRUCTIONS
    !###########################################  
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
    ! removing upper triangular part
    DO i=1,N
       DO j=1,i-1
          A(i,j)=0.0
       ENDDO
    ENDDO
    CALL DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,&
         ABSTOL, M, W, Z, LDZ, WORK, LWORK,IWORK, IFAIL, INFO )
    LWORK = WORK(1)
    LIWORK = IWORK(1)
    IF (INFO/=0) PRINT *, 'Error in DSYERVR ? INFO=',INFO
  END SUBROUTINE solve_dsyevx


END MODULE module_solve
