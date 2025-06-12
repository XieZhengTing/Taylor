
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
    
    

      
    SUBROUTINE GET_NUM_PROP(MAT_TYPE,NUM_PROP)
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! RETURN THE NUMBER OF PROPERTIES ASSOCIATED WITH THIS MATERIAL ID
	  !
	  IMPLICIT NONE
	  
    INTEGER,INTENT(IN)::MAT_TYPE
    INTEGER,INTENT(OUT)::NUM_PROP
    
    NUM_PROP = 10
    IF (MAT_TYPE.EQ.1) THEN
    NUM_PROP = 3
    ELSEIF(MAT_TYPE.EQ.2) THEN
    NUM_PROP = 8
    ! FIX LATER
    ENDIF
    
    RETURN
    END SUBROUTINE