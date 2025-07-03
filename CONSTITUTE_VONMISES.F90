
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
	  
	  
	  

	  

      SUBROUTINE VON_MISES(STRESS, STRAIN, STRESS_PREDICT, STATE, PROPS)
          !$ACC ROUTINE SEQ
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
	  DOUBLE PRECISION, INTENT(IN)::    STRAIN(6)
	  DOUBLE PRECISION, INTENT(IN)::    STRESS_PREDICT(6)
	  DOUBLE PRECISION, INTENT(IN)::    PROPS(30)
	  !
	  !LOOP INDEX VARIABLES
	  INTEGER:: I, J, K
	  !PROPERTY VARIABLES
	  DOUBLE PRECISION::  YOUNG,POISS,Q_PHI,K_PHI,Q_PSI,T_CUT,K_I,K_C,DAM_MAX
	  !STATE VARIABLES
	  DOUBLE PRECISION:: EPS, EPS_VOL, EPS_SHEAR, DAM
	  !TEMPORARY VARIABLES
	  DOUBLE PRECISION:: BULK, DIFF, TAU_P, ALPHA_P, SPH_DIFF, H
	  DOUBLE PRECISION:: F_SHEAR_YEILD, MU
	  LOGICAL:: SHEAR_YEILD, TENSION_YEILD, ELASTIC
	  DOUBLE PRECISION:: TAU_PRE,SPH_PRE(6),ETA,TAU, SHEAR, GAMMA, SPH
      DOUBLE PRECISION:: EFF_STRESS_PREDICT(6), DEV_STRESS_PRE(6)
	  DOUBLE PRECISION:: DELTA_GAMMA, SIG_NOT, SIG23, TNORM_DEV_STRESS_PRE, YLD_FN
	  DOUBLE PRECISION:: DYLD_FN,XK,XKP,EPST,B,CE
	  !
	  !
	  !
	  ! ASSIGN PROPERTY VARIABLES
	  !
	  POISS = PROPS(1)
	  YOUNG = PROPS(2)
	  SIG_NOT = PROPS(4)
	  MU = SHEAR_MOD(YOUNG,POISS)
      B = PROPS(5)   !125.D0
      CE = PROPS(6)  !0.1D0
      !TESTING: HARDCODE SOME PARAMETERS
      !SIG_NOT=2.70000E+02
      

      
	  !
	  ! ASSIGN STATE VARIABLES
	  !
	  EPS = STATE(3)
	  !
	  ! DEVIATORIC PREDICTOR STRESS
	  !
	  DEV_STRESS_PRE = DEV_PROJ(STRESS_PREDICT)
	  !
	  ! SPHERICAL PREDICTOR STRESS
	  !
	  SPH_PRE = SPHR_PROJ(STRESS_PREDICT)
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
      
    !!!!!      ! PERFECT PLASTICITY
    !!!!!      !
		  !!!!!! MAP BACK TO THE YEILD SURFACE
		  !!!!!!
		  !!!!!DEV_STRESS_PRE = DEV_STRESS_PRE * SIG23 / TNORM_DEV_STRESS_PRE
		  !!!!!!
		  !!!!!! ACCUMULATE THE PLASTIC STAIN
		  !!!!!! 
		  !!!!!DELTA_GAMMA = YLD_FN / (2.0d0 * MU)
		  !!!!!!
          
          
          
		  !!!!EPS = EPS + DSQRT(2.0d0/3.0d0)*DELTA_GAMMA
		  !
	  END IF
	 
		  !
		  ! GRAB THE FINAL VALUES
		  STATE(3) = EPS
	      !
		  !
     ! ASSEMBLE THE FINAL STRESS AFTER RADIAL RETURN
     STRESS = -SPH_PRE + DEV_STRESS_PRE
          
	  RETURN
	  END SUBROUTINE
		