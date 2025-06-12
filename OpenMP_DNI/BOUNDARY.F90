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
	  !  
      SUBROUTINE BOUNDARY(NUMP,FRC,ACL,VEL,DSP,MODEL_EBC)
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! GET THE PREDICTOR VALUES
	  ! 
      INTEGER,INTENT(IN):: NUMP, MODEL_EBC(3,NUMP)
	  DOUBLE PRECISION, INTENT(INOUT):: DSP(NUMP*3)
	  DOUBLE PRECISION, INTENT(INOUT):: VEL(NUMP*3)
	  DOUBLE PRECISION, INTENT(INOUT):: ACL(NUMP*3)
	  DOUBLE PRECISION, INTENT(INOUT):: FRC(NUMP*3)
	  
	  INTEGER:: I
	  
	  DO I=1, NUMP
	  
	    IF (MODEL_EBC(1,I).EQ.1) THEN
	       DSP((I-1)*3+1) = 0.0d0
	       VEL((I-1)*3+1) = 0.0d0
	       ACL((I-1)*3+1) = 0.0d0
	       FRC((I-1)*3+1) = 0.0d0
	    END IF
	    
	    IF (MODEL_EBC(2,I).EQ.1) THEN
	       DSP((I-1)*3+2) = 0.0d0
	       VEL((I-1)*3+2) = 0.0d0
	       ACL((I-1)*3+2) = 0.0d0
	       FRC((I-1)*3+2) = 0.0d0
	    END IF
	    
	    IF (MODEL_EBC(3,I).EQ.1) THEN
	       DSP((I-1)*3+3) = 0.0d0
	       VEL((I-1)*3+3) = 0.0d0
	       ACL((I-1)*3+3) = 0.0d0
	       FRC((I-1)*3+3) = 0.0d0
	    END IF
	    
	  END DO
	  
	  RETURN
	  
	  END SUBROUTINE
            
            
            
            
            
            
            