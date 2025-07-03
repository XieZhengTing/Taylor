
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
                                  LSTATE, LSTRESS, S_STRESS, H_STRESS, FAILED) !IN/OUT, OUT

		!$ACC ROUTINE SEQ
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! COMPUTE THE CAUCHY STRESS GIVEN THE ELASTIC PREDICTOR STRESS AND MATERIAL TYPE
	  !
	  IMPLICIT NONE
            IMPLICIT NONE
            INTERFACE
                SUBROUTINE VON_MISES(STRESS, STRAIN, STRESS_PREDICT, STATE, PROPS, FAILED)
                    !$ACC ROUTINE SEQ
                    DOUBLE PRECISION, INTENT(INOUT) :: STRESS(6)
                    DOUBLE PRECISION, INTENT(INOUT) :: STATE(20)
                    DOUBLE PRECISION, INTENT(IN) :: STRAIN(6)
                    DOUBLE PRECISION, INTENT(IN) :: STRESS_PREDICT(6)
                    DOUBLE PRECISION, INTENT(IN) :: PROPS(30)
                    LOGICAL, INTENT(OUT) :: FAILED
                END SUBROUTINE VON_MISES
                SUBROUTINE DRUCK_PRAG(STRESS, TOTAL_STRAIN, INC_STRAIN, STRESS_PREDICT, STATE, PROPS)
                    !$ACC ROUTINE SEQ
                    DOUBLE PRECISION, INTENT(INOUT) :: STRESS(6)
                    DOUBLE PRECISION, INTENT(INOUT) :: STATE(20)
                    DOUBLE PRECISION, INTENT(IN) :: INC_STRAIN(6)
                    DOUBLE PRECISION, INTENT(IN) :: TOTAL_STRAIN(6)
                    DOUBLE PRECISION, INTENT(IN) :: STRESS_PREDICT(6)
                    DOUBLE PRECISION, INTENT(IN) :: PROPS(30)
                END SUBROUTINE DRUCK_PRAG
                SUBROUTINE VISCO_ELASTIC(PROPS, DLT, FMAT, LSTRESS, H_STRESS, S_STRESS)
                    !$ACC ROUTINE SEQ
                    DOUBLE PRECISION, INTENT(IN) :: PROPS(30), DLT, FMAT(3,3)
                    DOUBLE PRECISION, INTENT(OUT) :: LSTRESS(6), H_STRESS(6), S_STRESS(6)
                END SUBROUTINE VISCO_ELASTIC
                SUBROUTINE VON_MISES_DAM(STRESS, TOTAL_STRAIN, INC_STRAIN, STRESS_PREDICT, STATE, PROPS, FAILED)
                    !$ACC ROUTINE SEQ
                    DOUBLE PRECISION, INTENT(INOUT) :: STRESS(6)
                    DOUBLE PRECISION, INTENT(INOUT) :: STATE(20)
                    DOUBLE PRECISION, INTENT(IN) :: INC_STRAIN(6)
                    DOUBLE PRECISION, INTENT(IN) :: TOTAL_STRAIN(6)
                    DOUBLE PRECISION, INTENT(IN) :: STRESS_PREDICT(6)
                    DOUBLE PRECISION, INTENT(IN) :: PROPS(30)
                    LOGICAL, INTENT(OUT) :: FAILED
                END SUBROUTINE VON_MISES_DAM
            END INTERFACE
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
	  LOGICAL, INTENT(OUT) :: FAILED
	  !
	  ! LOCAL
	  !
	  !
	  !*********************************************************
	  !******************** EXECUTABLE CODE ********************
	  !*********************************************************
	  !
      FAILED = .FALSE.
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
            CALL VON_MISES(LSTRESS, LSTRAIN, LSTRESS_PREDICTOR, LSTATE, LPROP, FAILED)
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
            CALL VON_MISES_DAM(LSTRESS, LSTRAIN, STRAIN, LSTRESS_PREDICTOR, LSTATE, LPROP, FAILED)

		  !
		  !
                  CASE DEFAULT
                    !
                    !SOMETHING WENT WRONG
                    !
                    FAILED = .TRUE.
		  !
      END SELECT
		
      RETURN
		
      END SUBROUTINE
		