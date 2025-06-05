
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
    
    MODULE DETERMINANT_MOD
      IMPLICIT NONE
    CONTAINS      
      
    SUBROUTINE DETERMINANT(A,DET)
    !$ACC ROUTINE SEQ         
    IMPLICIT NONE
    
    !GLOBAL
    DOUBLE PRECISION:: A(3,3), DET

    !LOCAL
    INTEGER:: I,J
  
    DET =       A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) &

                   + A(1,2)*(A(3,1)*A(2,3) - A(2,1)*A(3,3))  &
                     + A(1,3)*(A(2,1)*A(3,2) - A(3,1)*A(2,2))


    RETURN
    
    
    END SUBROUTINE
    END MODULE DETERMINANT_MOD
    
    
    
    