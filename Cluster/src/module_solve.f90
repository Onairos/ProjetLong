!>Contains methods from Lapack library dealing with eigen values computing
 MODULE module_solve
CONTAINS



!>
!! @param A the affinity matrix
!! @param VR 
!! @param WR 
!! @param N 
  SUBROUTINE solve_dgeevx(N, A, VR, WR)
    IMPLICIT NONE
    EXTERNAL DGEEVX

    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: N
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: A

    !====  OUT ====
    DOUBLE PRECISION, DIMENSION(:), POINTER :: WR
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: VR

    !#### Variables  ####
    CHARACTER :: BALANC
    CHARACTER :: JOBVL
    CHARACTER :: JOBVR
    CHARACTER :: SENSE
    INTEGER :: IHI
    INTEGER :: ILO
    INTEGER :: INFO
    INTEGER :: LDA
    INTEGER :: LDVL
    INTEGER :: LDVR
    INTEGER :: LWORK
    INTEGER, DIMENSION(:), POINTER :: IWORK
    DOUBLE PRECISION :: ABNRM
    DOUBLE PRECISION, DIMENSION(:), POINTER :: RCONDE
    DOUBLE PRECISION, DIMENSION(:), POINTER :: RCONDV
    DOUBLE PRECISION, DIMENSION(:), POINTER :: SCALE
    DOUBLE PRECISION, DIMENSION(:), POINTER :: WI
    DOUBLE PRECISION, DIMENSION(:), POINTER :: WORK
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: VL

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




!>
!! @param[in] A the affinity matrix
!! @param[in] N 
!! @param[out] VR 
!! @param[out] WR 
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
    DOUBLE PRECISION, DIMENSION(:), POINTER :: WR
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: VR

    !#### Variables  ####
    CHARACTER :: JOBVL
    CHARACTER :: JOBVR
    INTEGER :: INFO
    INTEGER :: LDA
    INTEGER :: LDVL
    INTEGER :: LDVR
    INTEGER :: LWORK
    DOUBLE PRECISION, DIMENSION(:), POINTER :: WI
    DOUBLE PRECISION, DIMENSION(:), POINTER :: WORK
    DOUBLE PRECISION, DIMENSION(:,:), POINTER :: VL

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



!>
!! @param A the affinity matrix
!! @param W 
!! @param LWORK 
!! @param N 
  SUBROUTINE solve_dsyev(N, A, W, LWORK)
    IMPLICIT NONE
    EXTERNAL DSYEV
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: N
    DOUBLE PRECISION :: A(N,N)

    !====  OUT ====
    INTEGER :: LWORK
    DOUBLE PRECISION :: W(N)

    !#### Variables  ####
    CHARACTER :: JOBZ
    CHARACTER :: UPLO
    INTEGER :: INFO
    INTEGER :: LDA
    DOUBLE PRECISION :: WORK(LWORK)

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



!>
!! @param A the affinity matrix
!! @param Z the matrix of eigen vectors
!! @param W 
!! @param k 
!! @param LIWORK 
!! @param LWORK 
!! @param M 
!! @param N 
  SUBROUTINE solve_dsyevr(k, N, A, Z, LWORK, LIWORK, W, M)
    IMPLICIT NONE
    EXTERNAL DSYEVR
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !====  IN  ====
    INTEGER :: k
    INTEGER :: M
    INTEGER :: N
    INTEGER :: LIWORK
    INTEGER :: LWORK
    DOUBLE PRECISION :: W(N)
    DOUBLE PRECISION :: A(N,N)
    DOUBLE PRECISION :: Z(N,N)

    !#### Variables  ####
    CHARACTER :: JOBZ
    CHARACTER :: RANGE
    CHARACTER :: UPLO
    INTEGER :: INFO
    INTEGER :: i
    INTEGER :: IL
    INTEGER :: IU
    INTEGER :: j
    INTEGER :: LDA
    INTEGER :: LDZ
    INTEGER :: ISUPPZ(2*N)
    INTEGER :: IWORK(LIWORK)
    DOUBLE PRECISION :: WORK(LWORK)
    DOUBLE PRECISION :: ABSTOL
    DOUBLE PRECISION :: VL
    DOUBLE PRECISION :: VU

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



!>
!! @param A the affinity matrix
!! @param Z the matrix of eigen vectors
!! @param W 
!! @param k 
!! @param LIWORK 
!! @param LWORK 
!! @param M 
!! @param N 
  SUBROUTINE solve_dsyevx(k, N, A, Z, LWORK, LIWORK, W, M)
    IMPLICIT NONE
    EXTERNAL DSYEVX
    !###########################################
    ! DECLARATIONS
    !###########################################      
    !#### Parameters ####
    !=== IN/OUT ===
    INTEGER :: LIWORK
    INTEGER :: LWORK
    INTEGER :: k
    INTEGER :: M
    INTEGER :: N
    DOUBLE PRECISION :: W(N)
    DOUBLE PRECISION :: A(N,N)
    DOUBLE PRECISION :: Z(N,N)

    !#### Variables  ####
    CHARACTER :: JOBZ
    CHARACTER :: UPLO
    CHARACTER :: RANGE
    INTEGER :: i
    INTEGER :: IFAIL
    INTEGER :: IL
    INTEGER :: INFO
    INTEGER :: IU
    INTEGER :: j
    INTEGER :: LDA
    INTEGER :: LDZ
    INTEGER :: IWORK(LIWORK)
    DOUBLE PRECISION :: ABSTOL
    DOUBLE PRECISION :: VL
    DOUBLE PRECISION :: VU
    DOUBLE PRECISION :: WORK(LWORK)

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
