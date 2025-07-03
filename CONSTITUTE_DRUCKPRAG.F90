
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
	  
	  
	  

	  

	SUBROUTINE DRUCK_PRAG(STRESS, TOTAL_STRAIN, INC_STRAIN, STRESS_PREDICT, STATE, PROPS)
    !$ACC ROUTINE SEQ	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! ROTATE A (VOIGT NOTATION) TENSOR USING GIVEN ROTATION MATRIX
	  !
      USE FINT_FUNCTIONS
      !
	  IMPLICIT NONE
	  !
	  !GLOBAL IN-OUT
	  DOUBLE PRECISION, INTENT(INOUT):: STRESS(6)
	  DOUBLE PRECISION, INTENT(INOUT):: STATE(20)
	  DOUBLE PRECISION, INTENT(IN)::    INC_STRAIN(6) !INCREMENTAL OBJECTIVE STRAIN
	  DOUBLE PRECISION, INTENT(IN)::	TOTAL_STRAIN(6) !TOTAL STRAIN
	  DOUBLE PRECISION, INTENT(IN)::    STRESS_PREDICT(6)
	  DOUBLE PRECISION, INTENT(IN)::    PROPS(30)
	  !
	  !LOOP INDEX VARIABLES
	  INTEGER:: I, J, K
	  !PROPERTY VARIABLES
	  DOUBLE PRECISION::  YOUNG,POISS,Q_PHI,K_PHI,Q_PSI,T_CUT,K_I,K_C,DAM_MAX, K_I_CHECK
	  !STATE VARIABLES
	  DOUBLE PRECISION:: EPS, EPS_VOL, EPS_SHEAR, DAM, DAM_SV
	  !TEMPORARY VARIABLES
	  DOUBLE PRECISION:: BULK, DIFF, TAU_P, ALPHA_P, SPH_DIFF, H
	  DOUBLE PRECISION:: F_SHEAR_YEILD, MU
	  LOGICAL:: SHEAR_YEILD, TENSION_YEILD, ELASTIC
	  DOUBLE PRECISION:: TAU_PRE,SPH_PRE,ETA,TAU, SHEAR, GAMMA, SPH
      DOUBLE PRECISION:: EFF_STRESS_PREDICT(6), DEV_STRESS_PRE(6),DEV_STRESS(6)
	  DOUBLE PRECISION:: ELAS_MAT(6,6)
	  !
	  !
	  !
	  ! ASSIGN PROPERTY VARIABLES
	  !
	  POISS = PROPS(1)
	  YOUNG = PROPS(2)
	  Q_PHI = PROPS(25)
	  K_PHI = PROPS(26)
	  Q_PSI = PROPS(27)
	  T_CUT = PROPS(8)
	  K_I = PROPS(9) !DAM PAR #1
	  K_C = PROPS(10) !DAM PAR #2
	  DAM_MAX = PROPS(11)
	  !
	  ! ASSIGN STATE VARIABLES
	  !
	  EPS_SHEAR = STATE(1)
	  EPS_VOL = STATE(2)
	  EPS = STATE(3)
	  DAM = STATE(4)
	  DAM_SV = STATE(5) !STATE VARIABLE FOR CURRENT DAMAGE (FOR IRREVERSABILITY)
	  !
	  K_I_CHECK = DAM_SV + K_I
      !
      ELAS_MAT = FORM_CMAT(PROPS)
	  !
	  ! SET THE STRESS TO BE IN EFFECTIVE SPACE
	  !
	  EFF_STRESS_PREDICT = STRESS / (1.0d0-DAM)
      EFF_STRESS_PREDICT = EFF_STRESS_PREDICT + MATMUL(ELAS_MAT,INC_STRAIN)
	  !
	  ! COMPUTE AND ACCUMULATE DAMAGE
	  !
	  ETA = TENSOR_NORM(TOTAL_STRAIN)
	  !
	  IF (ETA.LT.K_I_CHECK) THEN
	  !
	  ! DO NOT INDUCE DAMAGE
	  !
	  !DAM = 0.0d0
	  !
	  ELSEIF ((ETA.GT.K_I).AND.(ETA.LT.K_C)) THEN
	  !
	  ! DAMAGE HAS BEEN INICIATED, COMPUTE IT
	  !
	  DAM = K_C*(ETA-K_I)/(ETA*(K_C-K_I))
	  !
	  DAM_SV = ETA
	  DAM_SV = DAM_SV - K_I
	  !
	  !ELSEIF (ETA.GT.K_C) THEN
	  !
	  ! DAMAGE IS "ONE"
	  !
	  !DAM = DAM_MAX
	  !
	  END IF
	  !
	  ! LIMIT THE DAMAGE
	  !
	  IF (DAM.GT.(DAM_MAX)) DAM = DAM_MAX
	  !
	  ! DEVIATORIC PREDICTOR STRESS
	  !
	  DEV_STRESS_PRE = DEV_PROJ(EFF_STRESS_PREDICT)
	  !
	  ! EFFECTIVE SHEAR PREDICTOR STRESS
	  !
	  TAU_PRE = TENSOR_NORM(DEV_STRESS_PRE)
	  !
	  ! SPHERICAL PREDICTOR STRESS
	  !
	  SPH_PRE = (EFF_STRESS_PREDICT(1) + EFF_STRESS_PREDICT(2) + EFF_STRESS_PREDICT(3))/3.0d0
	  !
	  ! CHECK THE TENSION CUT-OFF AND OTHER SURFACES
	  !
	  TAU_P = K_PHI - Q_PHI*T_CUT
	  !
	  ALPHA_P = SQRT(1.0d0 + Q_PSI**2) - Q_PSI
	  !
	  SPH_DIFF = SPH_PRE - T_CUT
	  !
	  H = TAU_PRE - TAU_P - ALPHA_P*SPH_DIFF
	  !
      F_SHEAR_YEILD = TAU_PRE + Q_PHI*SPH_PRE - K_PHI
	  !
		  IF (F_SHEAR_YEILD.LT.(0.0d0)) THEN
			  !
			  ! ELASTIC
			  !
			  SHEAR_YEILD = .FALSE.
			  ELASTIC = .TRUE.
			  TENSION_YEILD = .FALSE.
			  !
		  ELSE 
			  !
			  ! SHEAR YEILD
			  !
			  SHEAR_YEILD = .TRUE.
			  ELASTIC = .FALSE.
			  TENSION_YEILD = .FALSE.
			  !
		  END IF
          
          
     
     
		  IF (ELASTIC) THEN
			!
			! NO YEILDING, THE TRIAL STRESS IS THE FINAL STRESS
			!
			STRESS = EFF_STRESS_PREDICT
			!
		  END IF
	  
		  IF (TENSION_YEILD) THEN
			!
			!MAP SPHERICAL STRESS TO YEILD SURFACE 
			!
			STRESS = EFF_STRESS_PREDICT
			STRESS(1) = STRESS(1) - SPH_DIFF
			STRESS(2) = STRESS(2) - SPH_DIFF
			STRESS(3) = STRESS(3) - SPH_DIFF
            
            !ONLY ALLOW T_CUT TENSION SPERICAL STRESS
            SPH_PRE = T_CUT
            
			!
			! ACCUMULATE PLASTIC INC_STRAIN
			!
			BULK = BULK_MOD(YOUNG,POISS)
			EPS_VOL = EPS_VOL + DIFF/BULK
			EPS = EPS + EPS_VOL
			!
		  END IF
		  
		  IF (SHEAR_YEILD) THEN
			!
			!MAP SHEAR BACK TO YEILD SURFACE
			!
			TAU = K_PHI - Q_PHI*SPH_PRE
			!
			BULK = BULK_MOD(YOUNG,POISS)
			SHEAR = SHEAR_MOD(YOUNG,POISS)
			!
			!MAP SPHERICAL STRESS TO YEILD SURFACE 
			!(IF NON-ASSOCIATIVE, Q_PSI= 0, PRESSURE WILL STAY THE SAME: NO DILATION)
			!

            !SET ZERO TO ONLY ALLOW NON-ASSOCIATE FLOW RULE FIRST
            !Q_PSI = 0.D0 

            
			GAMMA = F_SHEAR_YEILD / (2.D0*SHEAR + BULK*Q_PHI*Q_PSI)
			!
			SPH = SPH_PRE - BULK*Q_PSI*GAMMA
            
			!
			! CONSTRUCT THE STRESS
			!
			DEV_STRESS = DEV_STRESS_PRE * TAU/TAU_PRE
            
            STRESS = DEV_STRESS
            
			STRESS(1) = DEV_STRESS(1) + SPH
			STRESS(2) = DEV_STRESS(2) + SPH
			STRESS(3) = DEV_STRESS(3) + SPH
			!
			! ACCUMULATE PLASTIC INC_STRAIN
			!
			EPS_SHEAR = EPS_SHEAR + GAMMA*SQRT(2.0d0/3.0d0+2.0d0/9.0d0*Q_PSI**2)
			EPS = EPS + EPS_SHEAR
			!
		  END IF
		  !
		  ! RETURN TO THE REAL STRESS SPACE
		  !
	      STRESS = STRESS * (1.0d0-DAM)
		  !
		  ! GRAB THE FINAL VALUES
		  STATE(1) = EPS_SHEAR
		  STATE(2) = EPS_VOL
		  STATE(3) = EPS
		  STATE(4) = DAM
	      STATE(5) = DAM_SV
	      !
		  !
	  RETURN
	  END SUBROUTINE
		