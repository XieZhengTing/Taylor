    !************************* MEGA *************************
    !********************************************************
    !*                                                      *
    !********* Meshfree Explicit Galerkin Analysis **********
    !*                                                      *
    !*                                                      *
    !*                                                      *
    !*     Copyright 2016 Michael C Hillman                 *
    !*                                                      *
    !********************************************************

MODULE HYPERELASTIC_MOD
    USE FINT_FUNCTIONS
    USE DETERMINANT_MOD
    USE INVERSE_MOD
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: HYPERELASTIC

CONTAINS
    SUBROUTINE HYPERELASTIC(LPROP,LSTRESS,FTEN,LSTRAIN, ierr)
    !$ACC ROUTINE SEQ

    ! USE THIS SUBROUTINE TO GET THE 2ND PK STRESS BY DEFORMATION GRADIENT FOR COMPRESSIBLE NEO-HOOKEAN MATERAIL
!    USE FINT_FUNCTIONS
!    IMPLICIT NONE
    DOUBLE PRECISION, INTENT (IN):: LPROP(30)
    DOUBLE PRECISION, INTENT (OUT):: LSTRESS(6)
    DOUBLE PRECISION, INTENT (OUT):: LSTRAIN(6)
    DOUBLE PRECISION, INTENT (IN):: FTEN(3,3)
    INTEGER, INTENT(OUT) :: ierr ! Error flag (0 for success, non-zero for error)

!F2PY INTENT(IN)::LPROP(30),FTEN(3,3)
!F2PY INTENT(OUT)::LSTRESS(6),LSTRAIN(6)
    DOUBLE PRECISION:: LSTRESSTEN(3,3)
    DOUBLE PRECISION:: C(3,3),FT(3,3),CINV(3,3),DET
    DOUBLE PRECISION:: I1,I2,I3,I3T1,I3T2,JT
    DOUBLE PRECISION:: Q1,PTEMP
    INTEGER:: I,J
    INTEGER :: inv_err_flag, det_err_flag
    DOUBLE PRECISION:: IDENT(3,3)
    DATA IDENT/ 1.0d0, 0.0d0, 0.0d0, &
        0.0d0, 1.0d0, 0.0d0, &
        0.0d0, 0.0d0, 1.0d0/
    ierr = 0 ! Initialize error flag

    FT=TRANSPOSE(FTEN)

    C=MATMUL(FT,FTEN)   !RIGHT CAUCHY GREEN TENSOR
    
    LSTRAIN=TENSOR_2_VTENSOR(C)

    CALL DETERMINANT(FTEN,DET)
    IF (ABS(DET) < 1.0D-12) THEN
        ierr = 2 ! Error: Deformation gradient is singular or near-singular
        LSTRESS = HUGE(0.0D0)
        LSTRAIN = HUGE(0.0D0)
        RETURN
    END IF

    I1= C(1,1)+C(2,2)+C(3,3)

    CALL DETERMINANT(C,I3)
    IF (I3 < 0.0D0 .OR. ABS(I3 - DET**2) > 1.0D-9 * DET**2 ) THEN ! Check for consistency and positivity
        ierr = 3 ! Error: Determinant of C is invalid or inconsistent
        LSTRESS = HUGE(0.0D0)
        LSTRAIN = HUGE(0.0D0)
        RETURN
    END IF    
    JT=SQRT(I3)

    CALL INVERSE(C,3,CINV, inv_err_flag)
    IF (inv_err_flag /= 0) THEN
        ierr = 4 ! Error: Failed to invert C
        LSTRESS = HUGE(0.0D0)
        LSTRAIN = HUGE(0.0D0)
        RETURN
    END IF

    I2=C(1,1)*C(2,2)+C(2,2)*C(3,3)+C(1,1)*C(3,3)-C(1,2)*C(2,1)-C(2,3)*C(3,2)-C(1,3)*C(3,1)

    SELECT CASE (INT(LPROP(1)))

!    CASE(1) !NEO_HOOKEAN(LPROP,LSTRESS,FMAT)
!
!        Q1=LPROP(4)
!        I3T1=Q1*I3**(-1.D0/3.D0)
!        PTEMP=LPROP(2)*(JT-1.0D0)*JT
!
!        !CASE(2)
!
!        !CALL MOONEY_RIVILING(LPROP,LSTRESS,DET,FMAT)

    CASE(3) !J.S. CHEN 96 SECOND EXAMPLE

        Q1=LPROP(4)+2*LPROP(5)*(I1-3)+3*LPROP(6)*(I1-3)**2
        I3T1=Q1*I3**(-1.D0/3.D0)
        PTEMP=LPROP(2)*(JT-1.0D0)*JT

    CASE DEFAULT

        ! CALL EXIT_PROGRAM('NOT A HYPERELASTIC MATERIAL',1)
        ierr = 1 ! Set error flag for invalid material type
        LSTRESS = 0.0D0 ! Or some other indicator
        LSTRAIN = 0.0D0 ! Or some other indicator

    END SELECT

    DO I=1,3
        DO J=1,3
            LSTRESSTEN(I,J)=2.0D0*(I3T1*(IDENT(I,J)-(1.0D0/3.0D0)*I1*(CINV(I,J))))+PTEMP*CINV(I,J)
        END DO
    END DO

    CALL SPK2CAUCHY(LSTRESSTEN,DET,FTEN,FT, det_err_flag) ! Use a different error flag variable
    IF (det_err_flag /= 0) THEN
        ierr = 5 ! Error: SPK2CAUCHY failed
        LSTRESS = HUGE(0.0D0)
        LSTRAIN = HUGE(0.0D0) ! LSTRAIN is already set, but good to be consistent
        RETURN
    END IF

    LSTRESS=TENSOR_2_VTENSOR(LSTRESSTEN)



    END SUBROUTINE

    !*********************************************************

    SUBROUTINE SPK2CAUCHY(LSTRESSTEN_IO,DET_IN,FTEN_IN,FT_IN, ierr_out)
    !$ACC ROUTINE SEQ
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT):: LSTRESSTEN_IO(3,3)
    DOUBLE PRECISION, INTENT(IN):: FTEN_IN(3,3),FT_IN(3,3),DET_IN
    INTEGER, INTENT(OUT) :: ierr_out ! Error flag (0 for success, non-zero for error)
!    DOUBLE PRECISION:: LSTEMP(3,3)
    DOUBLE PRECISION:: LSTEMP_LOCAL(3,3)
!    LSTEMP=MATMUL(FTEN,LSTRESSTEN)

!    LSTRESSTEN=(1.0D0/DET)*MATMUL(LSTEMP,FT)
    ierr_out = 0 ! Initialize error flag

    IF (ABS(DET_IN) < 1.0E-12) THEN ! Check for determinant close to zero
        LSTRESSTEN_IO = 0.0D0 ! Or some other error indicator like NaN
        ierr_out = 2 ! Specific error code for singular transformation
        RETURN
    END IF

!    LSTEMP=MATMUL(FTEN_IN,LSTRESSTEN_IO)
!    LSTRESSTEN_IO=(1.0D0/DET_IN)*MATMUL(LSTEMP,FT_IN)
    ! Cauchy_Stress = (1/J) * F * Second_PK_Stress * F_Transpose
    LSTEMP_LOCAL=MATMUL(FTEN_IN,LSTRESSTEN_IO)      ! LSTEMP_LOCAL = F * Second_PK_Stress
    LSTRESSTEN_IO=(1.0D0/DET_IN)*MATMUL(LSTEMP_LOCAL,FT_IN) ! LSTRESSTEN_IO = (1/J) * (LSTEMP_LOCAL * F_T)

    END SUBROUTINE
    END MODULE HYPERELASTIC_MOD