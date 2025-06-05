
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

    SUBROUTINE GET_VOL2(DVOL,XYZEL)
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! RETRUN THE VOLUME ASSOCATED WITH A SINGLE HEX (?)
	  !

    IMPLICIT NONE

    !GLOBAL
    DOUBLE PRECISION:: DVOL
    DOUBLE PRECISION:: XYZEL(3,9)
    
    !LOCAL
    DOUBLE PRECISION:: A(3,3),B(3,3),C(3,3)
    DOUBLE PRECISION:: AD,BD,CD

    
    DVOL = 0.0d0
    
    A(:,1)=XYZEL(:,7)-XYZEL(:,1)
    B(:,1)=XYZEL(:,7)-XYZEL(:,1)
    C(:,1)=XYZEL(:,7)-XYZEL(:,1)
    A(:,2)=XYZEL(:,2)-XYZEL(:,1)
    B(:,2)=XYZEL(:,5)-XYZEL(:,1)
    C(:,2)=XYZEL(:,4)-XYZEL(:,1)
    A(:,3)=XYZEL(:,3)-XYZEL(:,6)
    B(:,3)=XYZEL(:,6)-XYZEL(:,8)
    C(:,3)=XYZEL(:,8)-XYZEL(:,3)

    
    
      CALL DETERMINANT(A,AD)
      CALL DETERMINANT(B,BD)
      CALL DETERMINANT(C,CD)
    
      DVOL =(AD+BD+CD)/6.0d0

    RETURN
    
    
    END SUBROUTINE
    
    

