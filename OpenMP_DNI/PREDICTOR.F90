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
	  
	  
      SUBROUTINE PREDICTOR(TOTAL_LOCAL_SIZE,ACL,VEL,DSP,DLT)
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! GET THE PREDICTOR VALUES
	  !
      INTEGER,INTENT(IN):: TOTAL_LOCAL_SIZE
	  DOUBLE PRECISION, INTENT(IN):: ACL(TOTAL_LOCAL_SIZE*3)
	  DOUBLE PRECISION, INTENT(INOUT):: DSP(TOTAL_LOCAL_SIZE*3)
	  DOUBLE PRECISION, INTENT(INOUT):: VEL(TOTAL_LOCAL_SIZE*3)
	  DOUBLE PRECISION, INTENT(IN):: DLT
	  
	  !PREDICT THE DISPLACEMENT INCREMENT AND VELOCITY FROM THE PREVIOUS ACCELERATION
	  DSP = DLT * VEL + DLT**2 * 0.5d0 * ACL
	  
	  VEL = VEL + DLT * 0.5d0 * ACL
	  
	  
	  RETURN
	  END SUBROUTINE