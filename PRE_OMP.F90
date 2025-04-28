 
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
 
      SUBROUTINE PRE_OMP
	  
      !
      ! PURPOSE OF THIS LOGIC:
	  !
	  ! READ IN THE *.k FILE (FOR NOW USE LS-DYNA FILES) AND ASSIGN
	  ! THE INFO TO MEGA ARRAYS
	  !
      USE CONTROL
    use omp_lib
    implicit none
    
    INTEGER:: I
    CHARACTER*80 CTEMP
    CHARACTER*80 CHAR_VAR
	
	integer :: num_args,IERROR
	
	
    ! NOTE: In the Visual Studio environment, go to the Fortran property page Language > Process OpenMP Directives > Enable.
	
	
	  LUSER_OMP_CORES = .FALSE.
	
      num_args = command_argument_count()
 
      IF (num_args.GT.0) THEN
 
      call get_command_argument(1,CTEMP)
		
	    CTEMP=TRIM(CTEMP)
		
	    IF (CTEMP.eq.'-j') THEN
			
		  call get_command_argument(2,CTEMP)
		  
		  READ(CTEMP,'(I1000)',iostat=IERROR) NCORES_INPUT
			
			  IF (IERROR.EQ.0) THEN !A NUMBER WAS READ, TELL OPENMP TO USE X CORES
			  
				LUSER_OMP_CORES = .TRUE.
				
				MAX_CORES = OMP_get_max_threads()
				
                WRITE(CHAR_VAR,'(A21,I10)') 'USING OMP, MAX_CORES=' , MAX_CORES
                CALL WRITE_OUT(CHAR_VAR)
			
                WRITE(CHAR_VAR,'(A23,I10)') 'NCORES_INPUT=' , NCORES_INPUT
                CALL WRITE_OUT(CHAR_VAR)
				
				IF (NCORES_INPUT.GT.MAX_CORES) THEN
 
				WRITE(CTEMP,'(A57,I4,A1)') 'WARNING, (NCORES_INPUT.GT.MAX_CORES), using MAX_CORES (=',MAX_CORES,')' 
				CALL WARN(CTEMP)
 
				NCORES_INPUT = MAX_CORES
				END IF
			  
                WRITE(CHAR_VAR,'(A23,I10)') 'FINAL CORES=' , NCORES_INPUT
                CALL WRITE_OUT(CHAR_VAR)
				
				CALL OMP_set_num_threads(NCORES_INPUT)
 
			  END IF
	    END IF
	  END IF
	  
	  IF (.NOT.LUSER_OMP_CORES) THEN
	    !DONT USE OMP
		CALL OMP_set_num_threads(1)
        NCORES_INPUT = 1
	  END IF
	  
	  RETURN
	  
	  END SUBROUTINE
	  
      
      
      
      
      
      
      
      

      
