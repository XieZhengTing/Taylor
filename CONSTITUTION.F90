
	  !********************************************************
	  !*                                                      *
	  !********************************************************
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
	  
MODULE CONSTITUTION_MOD
    USE FINT_FUNCTIONS
    USE CONSTITUTE_VONMISES_MOD    ! For VON_MISES
    USE CONSTITUTE_DRUCKPRAG_MOD   ! For DRUCK_PRAG
    USE CONSTITUTE_VISCOELASTIC_MOD! For VISCO_ELASTIC
    USE CONSTITUTE_VONMISES_DAM_MOD! For VON_MISES_DAM
    USE HYPERELASTIC_MOD           ! For HYPERELASTIC
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: CONSTITUTION

CONTAINS	  
	  
	  SUBROUTINE CONSTITUTION(LSTRESS_PREDICTOR, LMAT_TYPE, LSTRAIN, STRAIN, LPROP, DLT, FMAT, & !IN
		                  LSTATE, LSTRESS, S_STRESS, H_STRESS, ierr) !IN/OUT, OUT
	  !$ACC ROUTINE SEQ
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! COMPUTE THE CAUCHY STRESS GIVEN THE ELASTIC PREDICTOR STRESS AND MATERIAL TYPE
	  !
!      USE FINT_FUNCTIONS             ! Assuming this is needed by some logic within or for general functions
!      USE CONSTITUTE_VONMISES_MOD    ! For VON_MISES
!      USE CONSTITUTE_DRUCKPRAG_MOD   ! For DRUCK_PRAG (assuming module name)
!      USE CONSTITUTE_VISCOELASTIC_MOD! For VISCO_ELASTIC (assuming module name)
!      USE CONSTITUTE_VONMISES_DAM_MOD! For VON_MISES_DAM (assuming module name)
      ! USE HYPERELASTIC_MOD           ! If HYPERELASTIC is called directly and is in a module

!	  IMPLICIT NONE
	  !
	  !GLOBAL IN-OUT
	  !
	  DOUBLE PRECISION, INTENT(IN):: LSTRESS_PREDICTOR(6)
	  DOUBLE PRECISION, INTENT (IN):: LPROP(30)
	  INTEGER, INTENT(IN):: LMAT_TYPE
	  DOUBLE PRECISION, INTENT (OUT):: LSTRESS(6)
	  DOUBLE PRECISION, INTENT (INOUT):: LSTATE(20)
	  DOUBLE PRECISION, INTENT (IN):: LSTRAIN(6), STRAIN(6)
	  DOUBLE PRECISION, INTENT(IN):: DLT !GC 
	  DOUBLE PRECISION, INTENT(IN):: FMAT(3,3) !GC
	  DOUBLE PRECISION, INTENT (INOUT)::S_STRESS(6), H_STRESS(6) !GC 
	  INTEGER, INTENT(OUT) :: ierr ! Error flag: 0 = success, non-zero = error	  
	  !
	  ! LOCAL
	  !
	  !
	  INTEGER :: constitution_sub_error_flag ! To capture error from called routines

	  !*********************************************************
	  !******************** EXECUTABLE CODE ********************
	  !*********************************************************
	  ierr = 0
      SELECT CASE (LMAT_TYPE)
	  
	    CASE(1)
		  !
		  !ELASTIC
		  !
		  LSTRESS = LSTRESS_PREDICTOR
		  !
        CASE(2)
		  !
		  !VON MISES PLASTICITY WITH EXPOENTIAL HARDENING
		  !
          CALL VON_MISES(LSTRESS, LSTRAIN, LSTRESS_PREDICTOR, LSTATE, LPROP, constitution_sub_error_flag)
          IF (constitution_sub_error_flag /= 0) ierr = 20 + constitution_sub_error_flag ! Propagate error with offset

		  !
        CASE(3)
          !
          ! DRUCKER-PRAGER WITH ISOTROPIC DAMAGE AND TENSION CUT-OFF
          !
          CALL DRUCK_PRAG(LSTRESS, LSTRAIN, STRAIN, LSTRESS_PREDICTOR, LSTATE, LPROP, constitution_sub_error_flag)
          IF (constitution_sub_error_flag /= 0) ierr = 30 + constitution_sub_error_flag ! Propagate error with offset

          !
		CASE(4)
		  !
		  ! VISCOELASTIC USING THE NEO-HOOKEAN HYPERELASTICITY LAW
		  !
		  CALL VISCO_ELASTIC(LPROP, DLT, FMAT, LSTRESS, H_STRESS, S_STRESS, constitution_sub_error_flag)
		  IF (constitution_sub_error_flag /= 0) ierr = 40 + constitution_sub_error_flag ! Propagate error with offset

		  !
        CASE(6)
		  !
		  !VON MISES PLASTICITY WITH EXPOENTIAL HARDENING AND ISOTROPIC DAMAGE
		  !
          CALL VON_MISES_DAM(LSTRESS, LSTRAIN, STRAIN, LSTRESS_PREDICTOR, LSTATE, LPROP, constitution_sub_error_flag)
          IF (constitution_sub_error_flag /= 0) ierr = 60 + constitution_sub_error_flag ! Propagate error with offset

		  !
		  !
		CASE DEFAULT
		  !
		  !SOMETHING WENT WRONG
		  !
!		  CALL EXIT_PROGRAM('INVALID MATERIAL TYPE IN SUBROUTINE CONSTITUTION',1)
		  !
		  ierr = 1 ! Set error flag for invalid material type
		  LSTRESS = 0.0D0 ! Or some other indicator for stress in case of error		  
      END SELECT
		
      RETURN
		
      END SUBROUTINE
END MODULE CONSTITUTION_MOD
		