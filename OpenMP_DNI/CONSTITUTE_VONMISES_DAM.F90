
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
	  
	  
	  

	  

      SUBROUTINE VON_MISES_DAM(STRESS, TOTAL_STRAIN, INC_STRAIN, STRESS_PREDICT, STATE, PROPS)
	  !
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
	  DOUBLE PRECISION:: TAU_PRE,SPH_PRE(6),ETA,TAU, SHEAR, GAMMA, SPH
      DOUBLE PRECISION:: EFF_STRESS_PREDICT(6), DEV_STRESS_PRE(6)
	  DOUBLE PRECISION:: DELTA_GAMMA, SIG_NOT, SIG23, TNORM_DEV_STRESS_PRE, YLD_FN
	  DOUBLE PRECISION:: DYLD_FN,XK,XKP,EPST,B,CE
	  DOUBLE PRECISION:: ELAS_MAT(6,6)
	  !
	  !
	  !
	  ! ASSIGN PROPERTY VARIABLES
	  !
	  POISS = PROPS(1)
	  YOUNG = PROPS(2)
	  SIG_NOT = PROPS(4)
	  MU = SHEAR_MOD(YOUNG,POISS)
      B = PROPS(5) 
      CE = PROPS(6) 
	  K_I = PROPS(9) !DAM PAR #1
	  K_C = PROPS(10) !DAM PAR #2
	  DAM_MAX = PROPS(11)
	  !
	  ! ASSIGN STATE VARIABLES
	  !
	  EPS = STATE(3)
	  DAM = STATE(4)
	  DAM_SV = STATE(5) !STATE VARIABLE FOR CURRENT DAMAGE (FOR IRREVERSABILITY)
	  !
	  K_I_CHECK = DAM_SV + K_I
	  !
	  !FORM ELASTIC STIFFNESS FOR PREDICTOR
	  !
      ELAS_MAT = FORM_CMAT(PROPS)
	  !
	  ! SET THE STRESS TO BE IN EFFECTIVE SPACE
	  !
	  EFF_STRESS_PREDICT = STRESS / (1.0d0-DAM)
      EFF_STRESS_PREDICT = EFF_STRESS_PREDICT + MATMUL(ELAS_MAT,INC_STRAIN)
	  !
	  ! DEVIATORIC PREDICTOR FOR EFFECTIVE STRESS
	  !
	  DEV_STRESS_PRE = DEV_PROJ(EFF_STRESS_PREDICT)
	  !
	  ! COMPUTE AND ACCUMULATE DAMAGE
	  !
	  ETA = TENSOR_NORM(TOTAL_STRAIN)
	  !
	  !ETA = (TOTAL_STRAIN(1) + TOTAL_STRAIN(2)+TOTAL_STRAIN(3))*3.0d0
	  !ETA = MAX(ETA,0.0d0)
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
	  ! SPHERICAL PREDICTOR FOR EFFECTIVE STRESS
	  !
	  SPH_PRE = SPHR_PROJ(EFF_STRESS_PREDICT)
	  !
	  ! GET THE J2 INVARIENT
	  !
	  !J2_PREDICT = J2(STRESS_PREDICT)
	  !
	  ! GET THE TENSOR NORM OF DEVIATOR OF STRESS
	  !
	  TNORM_DEV_STRESS_PRE = TENSOR_NORM(DEV_STRESS_PRE)
	  !
	  ! CHECK THE YEILD FUNCTION
	  !
	  SIG23 = DSQRT(2.0d0/3.0d0)*SIG_NOT
      ! HARDENING GIVEN IN RKPM 1996
      !XK = SIG23*(1.D0+125.D0*EPS)**(0.1D0)
      XK = SIG23*(1.D0 + B*EPS)**(CE)     
	  !
	  YLD_FN =  TNORM_DEV_STRESS_PRE - XK
	  !
	  ! MAP BACK TO THE YEILD SURFACE (PURFECT PLASTICITY, FOR NOW)
	  !
	  IF (YLD_FN.GT.(0.0d0)) THEN
		
      !PERFORM NEWTON ITERATION
      DELTA_GAMMA = 0.D0
      EPST = EPS
      DO I = 1,20  !MAX ITERATE 20 TIMES
       !XK = SIG23*(1.D0+125.D0*EPS)**(0.1D0)
       XK = SIG23*(1.D0 + B*EPS)**(CE)       
       YLD_FN = TNORM_DEV_STRESS_PRE - XK - (2.D0*MU*DELTA_GAMMA)
      IF(ABS(YLD_FN).GT. SIG23*1.0E-05) THEN
       !XKP = 12.5D0*SIG_NOT*(1.D0+125.D0*EPS)**(-0.9D0)
       XKP = B*CE*SIG_NOT*(1.D0 + B*EPS)**(CE - 1.D0)       
       DYLD_FN = -2.D0*MU*(1.D0 + XKP/3.D0/MU)         
       DELTA_GAMMA = DELTA_GAMMA - YLD_FN/DYLD_FN
       !WRITE(*,*)"ITEGRATION:", I, YLD_FN
          IF(I.EQ.20) THEN
             !!WRITE(*,*)"TOO MANY ITERATIONS IN PlASTICITY, EQIT"
             !!PAUSE
             CALL EXIT_PROGRAM('PLASTICITY MATERIAL MODEL DID NOT CONVERGE',0)
             STOP
          ENDIF
      ELSE  !CONVERGED
            !WRITE(*,*)"CONVERGED"
            EPS = EPST + DSQRT(2.0d0/3.0d0)*DELTA_GAMMA
          EXIT
      ENDIF
      
      
      EPS = EPST + DSQRT(2.0d0/3.0d0)*DELTA_GAMMA
      END DO
      DEV_STRESS_PRE = DEV_STRESS_PRE - DELTA_GAMMA*2.D0*MU*DEV_STRESS_PRE/TNORM_DEV_STRESS_PRE
      
	  END IF
	 
		  !
		  ! GRAB THE FINAL VALUES
		  STATE(3) = EPS
		  STATE(4) = DAM
	      STATE(5) = DAM_SV
	      !
	      ! ASSEMBLE THE FINAL STRESS AFTER RADIAL RETURN
	      STRESS = -SPH_PRE + DEV_STRESS_PRE
          !
		  !
		  ! RETURN TO THE REAL STRESS SPACE
		  !
	      STRESS = STRESS * (1.0d0-DAM)
		  !
	  RETURN
	  END SUBROUTINE
		