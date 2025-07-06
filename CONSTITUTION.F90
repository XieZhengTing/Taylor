
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
	  
	  
	  
	  SUBROUTINE CONSTITUTION(LSTRESS_PREDICTOR, LMAT_TYPE, LSTRAIN, STRAIN, LPROP, DLT, FMAT, & !IN
		                  LSTATE, LSTRESS, S_STRESS, H_STRESS) !IN/OUT, OUT
	  !
	  USE GPU_ERROR
	  !$ACC ROUTINE SEQ
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! COMPUTE THE CAUCHY STRESS GIVEN THE ELASTIC PREDICTOR STRESS AND MATERIAL TYPE
	  !
	  IMPLICIT NONE
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
	  !
	  ! LOCAL
	  !
	  !
	  !*********************************************************
	  !******************** EXECUTABLE CODE ********************
	  !*********************************************************
	  !
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
          CALL VON_MISES(LSTRESS, LSTRAIN, LSTRESS_PREDICTOR, LSTATE, LPROP)
		  !
        CASE(3)
          !
          ! DRUCKER-PRAGER WITH ISOTROPIC DAMAGE AND TENSION CUT-OFF
          !
          CALL DRUCK_PRAG(LSTRESS, LSTRAIN, STRAIN, LSTRESS_PREDICTOR, LSTATE, LPROP)
          !
		CASE(4)
		  !
		  ! VISCOELASTIC USING THE NEO-HOOKEAN HYPERELASTICITY LAW
		  !
		  CALL VISCO_ELASTIC(LPROP, DLT, FMAT, LSTRESS, H_STRESS, S_STRESS)
		  !
        CASE(6)
		  !
		  !VON MISES PLASTICITY WITH EXPOENTIAL HARDENING AND ISOTROPIC DAMAGE
		  !
          CALL VON_MISES_DAM(LSTRESS, LSTRAIN, STRAIN, LSTRESS_PREDICTOR, LSTATE, LPROP)
		  !
		  !
		CASE DEFAULT
		  !
		  ! GPU version: set error flag
		  !
         CALL SET_GPU_ERROR('INVALID_MAT_TYPE', 1)
         GPU_ERROR_FLAG = 2
         LSTRESS = LSTRESS_PREDICTOR  ! Default to elastic
		  !
      END SELECT
		
      RETURN
		
      END SUBROUTINE
		