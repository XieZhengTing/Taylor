
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

    SUBROUTINE GET_VOL1(DVOL,XYZEL)
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! RETRUN THE VOLUME ASSOCATED WITH A SINGLE TET (?)
	  !

    IMPLICIT NONE

    !GLOBAL
    DOUBLE PRECISION:: DVOL
    DOUBLE PRECISION:: XYZEL(3,9)
    
    !LOCAL
    DOUBLE PRECISION:: A(3,3)
    DOUBLE PRECISION:: AD
    INTEGER:: I

    
    DVOL = 0.0d0
    
    DO I=1,3
    A(1,I)=XYZEL(I,4)-XYZEL(I,3)
    A(2,I)=XYZEL(I,1)-XYZEL(I,3)
    A(3,I)=XYZEL(I,2)-XYZEL(I,3)
    END DO
      CALL DETERMINANT(A,AD)

    
      DVOL =AD/6.0d0

    RETURN
    
    
    END SUBROUTINE
    
    
    

