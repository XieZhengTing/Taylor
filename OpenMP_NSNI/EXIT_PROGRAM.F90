
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
	  
      SUBROUTINE EXIT_PROGRAM(REASON,ITYPE)
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! STOP THE PROGRAM, AND TELL THE USER WHY
	  !
	  IMPLICIT NONE
	  CHARACTER(*):: REASON
	  INTEGER:: ITYPE !(0 = dev error, 1 = user error)
	  !
	  !*********************************************************
	  !******************** EXECUTABLE CODE ********************
	  !*********************************************************
	  !
	  !
	  ! WRITE OUT THE REASON
	  !
	  WRITE(*,*) ' '
	  WRITE(*,*) ' '
	  WRITE(*,*) ' '
      WRITE(*,*) '**********************************************************************'
      WRITE(*,*) '**********************************************************************'
      WRITE(*,*) '********************    A FATAL ERROR OCCURED   **********************'
      WRITE(*,*) '**********************************************************************'
      WRITE(*,*) '**********************************************************************'
	  WRITE(*,*) ' '
	  WRITE(*,*) ' '
	  WRITE(*,*) ' '
	  WRITE(*,*) REASON
	  WRITE(*,*) ' '
      CALL LOG_APPEND(' ')
      CALL LOG_APPEND(' ')
      CALL LOG_APPEND(' ')
      CALL LOG_APPEND('******************************************************************************')
      CALL LOG_APPEND('******************************************************************************')
      CALL LOG_APPEND('*************************    A FATAL ERROR OCCURED   *************************')
      CALL LOG_APPEND('******************************************************************************')
      CALL LOG_APPEND('******************************************************************************')
      CALL LOG_APPEND(' ')
      CALL LOG_APPEND(' ')
      CALL LOG_APPEND(' ')
      CALL LOG_APPEND(REASON)
      CALL LOG_APPEND(' ')
	  !
	  ! SAY IF IT WAS USER ERROR OR CODE ERROR (CODE NOT BULLETPROOF)
	  !
	  
	  !OPTIONS:
	  ! 
	  ! -1 = NOT IMPLEMENTED YET
	  ! 
	  !  0 = THE FAULT OF DEVELOPER
	  !
	  !  1 = LIKELY THE FUALT OF THE USER
	  !
	  !  2 = THE FUALT OF THE USER
	  !
      SELECT CASE (ITYPE)
	     CASE(-1)
	       WRITE(*,*) ' '
	       WRITE(*,*) 'THIS WAS THE FAULT OF THE DEVELOPER, THIS'
	       WRITE(*,*) 'OPTION IS NOT IMPLEMENTED YET'
	       WRITE(*,*) ' '
	       CALL LOG_APPEND( ' ')
	       CALL LOG_APPEND( 'THIS WAS THE FAULT OF THE DEVELOPER, THIS')
	       CALL LOG_APPEND( 'OPTION IS NOT IMPLEMENTED YET')
	       CALL LOG_APPEND( ' ')
	     CASE(0)
	       WRITE(*,*) ' '
	       WRITE(*,*) 'THIS WAS LIKELY THE FAULT OF THE DEVELOPER,'
	       WRITE(*,*) 'BUT COULD ALSO BE DUE TO INVALID USER INPUT,'
	       WRITE(*,*) 'WHICH WAS NOT DETECTED IN THE PREPROCESSOR'
	       WRITE(*,*) ' '
	       CALL LOG_APPEND( ' ')
	       CALL LOG_APPEND( 'THIS WAS LIKELY THE FAULT OF THE DEVELOPER,')
	       CALL LOG_APPEND( 'BUT COULD ALSO BE DUE TO INVALID USER INPUT,')
	       CALL LOG_APPEND( 'WHICH WAS NOT DETECTED IN THE PREPROCESSOR')
	       CALL LOG_APPEND( ' ')
		 CASE(1)
	       WRITE(*,*) ' '
	       WRITE(*,*) 'THIS WAS LIKELY THE FAULT OF THE USER,'
	       WRITE(*,*) 'OR THIS MISTAKE WAS NOT DETECTED IN THE PREPROCESSOR'
	       WRITE(*,*) ' '
	       CALL LOG_APPEND( ' ')
	       CALL LOG_APPEND( 'THIS WAS LIKELY THE FAULT OF THE USER,')
	       CALL LOG_APPEND( 'OR THIS MISTAKE WAS NOT DETECTED IN THE PREPROCESSOR')
	       CALL LOG_APPEND( ' ')
		 CASE(2)
	       WRITE(*,*) ' '
	       WRITE(*,*) 'THIS WAS THE FAULT OF THE USER'
	       WRITE(*,*) ' '
	       CALL LOG_APPEND( ' ')
	       CALL LOG_APPEND( 'THIS WAS THE FAULT OF THE USER')
	       CALL LOG_APPEND( ' ')
	  END SELECT
	  !
	  ! WRITE THE FINAL STATEMENT REGARDING FORTRAN PAUSE AND EXIT
	  !
	  WRITE(*,*) 'THE EXECUTION OF MEGA WILL NOW TERMINATE, INSERTING FORTRAN PAUSE'
	  WRITE(*,*) ' '
	  CALL LOG_APPEND( 'THE EXECUTION OF MEGA WILL NOW TERMINATE, INSERTING FORTRAN PAUSE')
	 CALL LOG_APPEND( ' ')
	  PAUSE
	  STOP
	  RETURN
	  END SUBROUTINE
	  
      
      