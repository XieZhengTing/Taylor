
          SUBROUTINE EST_TIME(NCORES_INPUT,TIMER_STEPS,SIM_TIME_1,SIM_TIME_2, &
                              REAL_TIME_1,REAL_TIME_2, &
                              SIM_TIME_LEFT,TIME,TIME_END, &
                              REAL_TIME_REMAINING,LINIT_TIME)
	      !
	      ! FUNCTION OF THIS LOGIC:
	      !
	      ! ESTIMATE THE TIME REMAINING IN SECONDS
	      !
          !
	      IMPLICIT NONE
          ! 
          INTEGER:: TIMER_STEPS,NCORES_INPUT
          REAL:: REAL_TIME_1,REAL_TIME_2
          DOUBLE PRECISION::SIM_TIME_1,SIM_TIME_2, TIME_END
          DOUBLE PRECISION::TIME, SIM_TIME_LEFT, REAL_TIME_REMAINING
          LOGICAL:: LINIT_TIME
          !
          !LOCAL
          !
          DOUBLE PRECISION:: REAL_TIME_TEMP, SIM_TIME_TEMP, SPEED
          INTEGER:: ICHECK
          
          
          
          IF (LINIT_TIME) THEN
          ICHECK = 0
          ELSE
          ICHECK = 50
          END IF
          
          IF (TIMER_STEPS.GT.ICHECK) THEN
            !
            TIMER_STEPS = 0
            !
            SIM_TIME_2 = TIME
            CALL CPU_TIME(REAL_TIME_2) !REAL_TIME_1,REAL_ TIME_2, REAL_TIME_TEMP
            !
            ! HOW LONG IN SIMULATION AND REAL TIME BETWEEN 100 STEPS
            !
            SIM_TIME_TEMP = SIM_TIME_2 - SIM_TIME_1
            REAL_TIME_TEMP = REAL_TIME_2 - REAL_TIME_1 
            
            SPEED = DBLE(REAL_TIME_TEMP) / SIM_TIME_TEMP
			
			SPEED = SPEED/DBLE(NCORES_INPUT)
            
            SIM_TIME_LEFT = TIME_END - TIME
            
            REAL_TIME_REMAINING = SIM_TIME_LEFT * SPEED
            
            SIM_TIME_1 = TIME
            CALL CPU_TIME(REAL_TIME_1) 
            
          END IF  
          
          RETURN
          
          END SUBROUTINE


            SUBROUTINE MAKE_CTIME(CTIME_ALL,REAL_TIME_REMAINING,TIMER_FLAG)
	      !
	      ! FUNCTION OF THIS LOGIC:
	      !
	      ! MAKE THE TIME REMAINING IN SECONDS HUMAN READABLE (HOURS ETC)
	      !
          !
	      IMPLICIT NONE
	      !
	      !GLOBAL
	      !
          CHARACTER*13:: CTIME_ALL
          LOGICAL::TIMER_FLAG
          DOUBLE PRECISION:: REAL_TIME_REMAINING
          !
          !LOCAL
          !
          INTEGER:: DAYS, HOURS, MINS, SECS
          CHARACTER*4:: CDAYS
          CHARACTER*2:: CHOURS, CMINS, CSECS
          
            DAYS  = FLOOR (REAL_TIME_REMAINING/(3600.0d0*24.0d0))
            HOURS = FLOOR (REAL_TIME_REMAINING/(3600.0d0)-DBLE(DAYS)*24.0d0)
            MINS  = FLOOR (REAL_TIME_REMAINING/60.0d0 - DBLE(HOURS)*60.0d0-DBLE(DAYS)*24.0d0*60.0d0)
            SECS  = FLOOR (REAL_TIME_REMAINING - DBLE(MINS)*60.0d0 - DBLE(HOURS)*3600.0d0-DBLE(DAYS)*24.0d0*3600.0d0)
            
          
              TIMER_FLAG = .FALSE.
              if (DAYS.lt.10) THEN
                WRITE(CDAYS,'(A3,I1)') '000',DAYS
              elseif (DAYS.lt.100) THEN
                WRITE(CDAYS,'(A2,I2)') '00',DAYS
              elseif (DAYS.lt.1000) THEN
                WRITE(CDAYS,'(A1,I3)') '0',DAYS
              elseif (DAYS.lt.10000) THEN
                WRITE(CDAYS,'(I4)') DAYS
              else
                TIMER_FLAG = .TRUE.
              end if
              
              if (HOURS.lt.10) THEN
                WRITE(CHOURS,'(A1,I1)') '0',HOURS
              elseif (HOURS.lt.100) THEN
                WRITE(CHOURS,'(I2)') HOURS
              end if
      
              if (MINS.lt.10) THEN
                WRITE(CMINS,'(A1,I1)') '0',MINS
              elseif (HOURS.lt.100) THEN
                WRITE(CMINS,'(I2)') MINS
              end if
            
              if (SECS.lt.10) THEN
                WRITE(CSECS,'(A1,I1)') '0',SECS
              elseif (HOURS.lt.100) THEN
                WRITE(CSECS,'(I2)') SECS
              end if
              
              WRITE(CTIME_ALL,'(A4,A1,A2,A1,A2,A1,A2)') CDAYS, ':', CHOURS, ':', CMINS, ':', CSECS
              
          RETURN
          
          END SUBROUTINE