
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

    SUBROUTINE GET_VOL(DVOL,XYZEL,T)
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! RETRUN THE VOLUME ASSOCATED WITH A SINGLE HEX (?)
	  !

    IMPLICIT NONE

    !GLOBAL
    DOUBLE PRECISION:: DVOL
    INTEGER:: T(12,4)
    DOUBLE PRECISION:: XYZEL(3,9)
    
    !LOCAL
    INTEGER:: I,J
    INTEGER:: N1,N2,N3,N4
    DOUBLE PRECISION:: A,B,C,D
    DOUBLE PRECISION:: DVOL_TET
    DOUBLE PRECISION:: AMD(3), BMD(3), CMD(3)
    DOUBLE PRECISION:: M(3,3)
    
    DVOL = 0.0d0
    
    DO I=1,12
    
    
       N1=T(I,1)
       N2=T(I,2)
       N3=T(I,3)
       N4=T(I,4)
    
      DO J=1,3
        
        A = XYZEL(J,N1)
        B = XYZEL(J,N2)
        C = XYZEL(J,N3)
        D = XYZEL(J,N4)
        
        AMD(J) = A - D
        BMD(J) = B - D
        CMD(J) = C - D
        
      END DO
      
      M(1,1) = AMD(1)
      M(2,1) = AMD(2)
      M(3,1) = AMD(3)
      
      M(1,2) = BMD(1)
      M(2,2) = BMD(2)
      M(3,2) = BMD(3)
      
      M(1,3) = CMD(1)
      M(2,3) = CMD(2)
      M(3,3) = CMD(3)
    
      CALL DETERMINANT(M,DVOL_TET)
    
      DVOL_TET = DVOL_TET/6.0d0
    
      DVOL = DVOL + DVOL_TET
    
    END DO
    
    
    
    
    RETURN
    
    
    END SUBROUTINE
    
    

