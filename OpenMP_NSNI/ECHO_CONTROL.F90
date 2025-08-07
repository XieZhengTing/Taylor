
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
    
    

	  SUBROUTINE ECHO_CONTROL
      
      !
	  ! FUNCTION OF THIS LOGIC:
	  !
	  ! ECHO THE CONTROL PARAMETERS OF THE MODEL
	  !
      USE MODEL
      USE CONTROL
      
      IMPLICIT NONE
      !
      !LOCAL
      INTEGER:: I, J, K
      INTEGER:: NUM_PROP
      CHARACTER(50)::CTEMP
      
      OPEN(50,FILE='MEGA.log',STATUS='OLD',POSITION='APPEND')

      WRITE(50,*)
      WRITE(50,*) '*************** NODES - CONTROL PROPERTIES ***************'
      WRITE(50,*)
      !WRITE(50,'(2A9,2X,A10, 3A7,   7A15, A7,A15)') 'NODE ID','BODY ID','BODY NAME','X-FIX','Y-FIX','Z-FIX', &
      !'VX-INITIAL','VY-INITIAL','VZ-INITIAL', &
      !'NORM-WIN','X-PHYS-WIN','Y-PHYS-WIN','Z-PHYS-WIN','MATID','PROPS'
      
      WRITE(50,'(2A9,2X, 3A7,   7A15, A7,A15)') 'NODE ID','BODY ID','X-FIX','Y-FIX','Z-FIX', &
      'VX-INITIAL','VY-INITIAL','VZ-INITIAL', &
      'NORM-WIN','X-PHYS-WIN','Y-PHYS-WIN','Z-PHYS-WIN','MATID','PROPS'
     
      DO I=1,MODEL_NUMP
        CALL GET_NUM_PROP(MODEL_MAT_TYPE(I),NUM_PROP)
       
        !CTEMP = MODEL_SET_NAMES(MODEL_BODY_ID(I))
        !CTEMP= TRIM(CTEMP)
        
        !WRITE(50,'(2I9,2X,A10,3I7,   7E15.5,I7,30E15.5)') I, MODEL_BODY_ID(I),CTEMP, MODEL_EBC(:,I), MODEL_VINIT(:,I), MODEL_NORM_WIN(I), MODEL_WIN(:,I),MODEL_MAT_TYPE(I), &
        ! MODEL_PROP(1:NUM_PROP,I)

        WRITE(50,'(2I9,2X,3I7,   7E15.5,I7,30E15.5)') I, MODEL_BODY_ID(I), MODEL_EBC(:,I), MODEL_VINIT(:,I), MODEL_NORM_WIN(I), MODEL_WIN(:,I),MODEL_MAT_TYPE(I), &
         MODEL_PROP(1:NUM_PROP,I)         

      END DO
      WRITE(50,*)

      
      
      CLOSE(50)
      
      
      RETURN
      
      END SUBROUTINE
      
      