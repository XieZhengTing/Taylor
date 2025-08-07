
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

      SUBROUTINE ESTIMATE_MODULI(STRESS_INC,STRAIN_INC,SHEAR_TRIAL,BULK_TRIAL, &
	                               SHEAR, BULK,NSNI_FLAG)
	!$ACC ROUTINE SEQ
      USE FINT_FUNCTIONS
	  !
	  ! FUNCTION OF THIS LOGIC:
	  !
	  ! GENERAL ALGORITHM FOR PROBING BULK AND SHEAR MODULI BASED ON STRESS
	  ! AND STRAIN INCREMENTS
	  !
	  
      !
	  IMPLICIT NONE
      !
      !LOCAL VARIABLES
	  DOUBLE PRECISION:: STRESS_INC(6), STRAIN_INC(6), POISS, YOUNG, BULK, SHEAR,  &
	                     STRESS_INC_DEV(6), STRESS_INC_SPHR(6), STRAIN_INC_DEV(6), STRAIN_INC_SPHR(6),  &
	                     NORM_STRESS_INC_DEV, NORM_STRESS_INC_SPHR, NORM_STRAIN_INC_DEV, NORM_STRAIN_INC_SPHR,   &
	                     SHEAR_TRIAL, BULK_TRIAL, PMOD, DENSITY, MAXMOD, MAX_VEL,   &
	                     DIST, XJ, YJ, ZJ, CHAR_DIST, DLT_TEMP, XI, YI, ZI 
	  LOGICAL:: NOEST_SHEAR, NOEST_BULK, NSNI_FLAG
	                     
	    STRESS_INC_DEV = DEV_PROJ(STRESS_INC)
	    STRESS_INC_SPHR = SPHR_PROJ(STRESS_INC)
	    
	    STRAIN_INC_DEV = DEV_PROJ(STRAIN_INC)
	    STRAIN_INC_SPHR = SPHR_PROJ(STRAIN_INC)
	    
	    NORM_STRESS_INC_DEV = TENSOR_NORM(STRESS_INC_DEV)
	    NORM_STRESS_INC_SPHR = TENSOR_NORM(STRESS_INC_SPHR)
	    
	    NORM_STRAIN_INC_DEV = TENSOR_NORM_STRAIN(STRAIN_INC_DEV)
	    NORM_STRAIN_INC_SPHR = TENSOR_NORM_STRAIN(STRAIN_INC_SPHR)
	    
	    NOEST_SHEAR = .FALSE.
	    NOEST_BULK = .FALSE.
	          
	        IF ((NORM_STRAIN_INC_DEV.GT.(1.0D-8)).AND. &
	            (NORM_STRESS_INC_DEV.GT.(1.0D-8))) THEN
	          SHEAR_TRIAL = NORM_STRESS_INC_DEV / NORM_STRAIN_INC_DEV / 2.0d0
	        ELSE
              SHEAR_TRIAL = SHEAR
              END IF
			  
			
			IF (.NOT.NSNI_FLAG) THEN
	        SHEAR_TRIAL = MAX(SHEAR_TRIAL, SHEAR*0.5d0)
	        SHEAR_TRIAL = MIN(SHEAR_TRIAL, SHEAR)
	        END IF
			
			SHEAR = SHEAR_TRIAL
			
	        IF ((NORM_STRAIN_INC_SPHR.GT.(1.0D-8)).AND. &
	            (NORM_STRESS_INC_SPHR.GT.(1.0D-8))) THEN
	          BULK_TRIAL = NORM_STRESS_INC_SPHR / NORM_STRAIN_INC_SPHR / 3.0d0
              ELSE
              BULK_TRIAL = BULK
	        END IF
			
              !BULK_TRIAL = BULK
			
			IF (.NOT.NSNI_FLAG) THEN
	        BULK_TRIAL = MAX(BULK_TRIAL, BULK*0.5d0)
	        BULK_TRIAL = MIN(BULK_TRIAL, BULK)
	        END IF
			
			BULK = BULK_TRIAL
	        
	        
	    RETURN
	        
	    END SUBROUTINE