
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
	  
	  
	  

	  

      SUBROUTINE ROTATE_TENSOR(ROT,VTENSOR)
      !$ACC ROUTINE SEQ
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! ROTATE A (VOIGT NOTATION) TENSOR USING GIVEN ROTATION MATRIX
	  !
      USE FINT_FUNCTIONS
      !
	  IMPLICIT NONE
	  !
	  !GLOBAL IN-OUT
	  DOUBLE PRECISION, INTENT(IN):: ROT(3,3)
	  DOUBLE PRECISION, INTENT(INOUT):: VTENSOR(6)
	  !
	  !LOCAL VARIABLES
	  DOUBLE PRECISION:: ROT_T(3,3), TENSOR(3,3)
	  DOUBLE PRECISION:: TEMP(3,3)
	  DOUBLE PRECISION:: TENSOR_ROTATED(3,3)
	  !
	  !
	  ! TRANSPOSE ROTATION
	  !
	  ROT_T = TRANSPOSE(ROT)
	  !
	  ! ASSIGN FULL TENSOR VALUES
	  !
	  !VOIGT ORDERING:
      !XX, YY, ZZ, XY, YZ, XZ
	  !
	  TENSOR = VTENSOR_2_TENSOR(VTENSOR)
	  !
	  ! THE ABOVE CODE REPLACES THE FOLLOWING:
	  !
	  ! TENSOR(1,1) = VTENSOR(1) !XX
	  ! TENSOR(2,2) = VTENSOR(2) !YY
	  ! TENSOR(3,3) = VTENSOR(3) !ZZ
	  ! TENSOR(1,2) = VTENSOR(4) !XY
	  ! TENSOR(2,1) = VTENSOR(4) !YX
	  ! TENSOR(2,3) = VTENSOR(5) !YZ
	  ! TENSOR(3,2) = VTENSOR(5) !ZY
	  ! TENSOR(1,3) = VTENSOR(6) !XZ
	  ! TENSOR(3,1) = VTENSOR(6) !ZX
	  !
	  TEMP = MATMUL(ROT,TENSOR)
	  TENSOR_ROTATED = MATMUL(TEMP,ROT_T)
	  ! 
	  ! ASSIGN VOIGT TENSOR VALUES
	  !
	  !VOIGT ORDERING:
      !XX, YY, ZZ, XY, YZ, XZ
	  !
	  VTENSOR = TENSOR_2_VTENSOR(TENSOR_ROTATED)
	  !
	  ! THE ABOVE CODE REPLACES THE FOLLOWING:
	  !
	  ! VTENSOR(1) = TENSOR_ROTATED(1,1) !XX
	  ! VTENSOR(2) = TENSOR_ROTATED(2,2) !YY
	  ! VTENSOR(3) = TENSOR_ROTATED(3,3) !ZZ
	  ! VTENSOR(4) = TENSOR_ROTATED(1,2) !XY
	  ! VTENSOR(4) = TENSOR_ROTATED(2,1) !YX
	  ! VTENSOR(5) = TENSOR_ROTATED(2,3) !YZ
	  ! VTENSOR(5) = TENSOR_ROTATED(3,2) !ZY
	  ! VTENSOR(6) = TENSOR_ROTATED(1,3) !XZ
	  ! VTENSOR(6) = TENSOR_ROTATED(3,1) !ZX
	  !
	  RETURN
	  END SUBROUTINE
		