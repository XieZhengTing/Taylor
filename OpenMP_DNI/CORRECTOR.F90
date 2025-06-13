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
	  
	  
      SUBROUTINE CORRECTOR(NUMP,VEL,ACL,DLT)
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! GET THE CORRECTED VALUES (JUST VELOCITY)
	  !
      INTEGER,INTENT(IN):: NUMP
	  DOUBLE PRECISION, INTENT(IN):: ACL(NUMP*3)
	  DOUBLE PRECISION, INTENT(INOUT):: VEL(NUMP*3)
	  DOUBLE PRECISION, INTENT(IN):: DLT
	  
	  VEL = VEL + 0.5d0 * DLT * ACL
	  
	  RETURN
	  END SUBROUTINE