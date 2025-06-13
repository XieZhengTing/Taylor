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
	  
	  
      SUBROUTINE EOM(NUMP,FINT,FEXT,MASS,ACL)
                                
      IMPLICIT NONE
                                
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! SOLVE FOR THE CURRENT ACCELERATION
	  !
      INTEGER,INTENT(IN):: NUMP
	  DOUBLE PRECISION, INTENT(IN):: FINT(NUMP*3)
	  DOUBLE PRECISION, INTENT(IN):: FEXT(NUMP*3)
	  DOUBLE PRECISION, INTENT(IN):: MASS(NUMP*3)
	  DOUBLE PRECISION, INTENT(OUT):: ACL(NUMP*3)
	  !
	  ! LOCAL
	  !
	  INTEGER:: I
	  
	  
	  DO I=1,NUMP*3
	  
	    ACL(I) = (FEXT(I) - FINT(I))/MASS(I)
	  
	  END DO
	  
	  RETURN
	  END SUBROUTINE