
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
    
      SUBROUTINE COUNT_NODES(MODEL_FILE, & !IN
                       MODEL_NUMP) !OUT
      
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! COUNT THE NUMBER OF NODES IN A *.K FILE
	  !
	  IMPLICIT NONE
      ! 
      ! GLOBAL IN-OUT
      !
      CHARACTER(50), INTENT(IN):: MODEL_FILE
      INTEGER, INTENT(OUT):: MODEL_NUMP
      !
	  !
	  ! LOCAL I/O FLAGS ETC
	  !
	  CHARACTER(100):: CTEMP
	  INTEGER:: IERROR
      !
      INTEGER:: I_TEMP, I2_TEMP(2), I3_TEMP(3), I4_TEMP(4)
      INTEGER:: I5_TEMP(5), I6_TEMP(6), I7_TEMP(7), I8_TEMP(8)
      !
      DOUBLE PRECISION:: D_TEMP, D2_TEMP(2), D3_TEMP(3), D4_TEMP(4)
      DOUBLE PRECISION:: D5_TEMP(5), D6_TEMP(6), D7_TEMP(7), D8_TEMP(8)
      
      
	  OPEN(10,FILE=MODEL_FILE,STATUS='OLD',ACTION='READ',iostat=IERROR)
      !
      ! COUNT THE NUMBER OF NODES
      !
      MODEL_NUMP = 0
      !
	  DO 
	  
100     READ(10,'(A100)',iostat=IERROR) CTEMP
	  
	    IF (IERROR.EQ.(-1)) EXIT
	    CTEMP=TRIM(CTEMP)
	    IF (CTEMP.eq.'*NODE') THEN
	  
	      READ (10,*) 
          
	        DO 
            
              READ(10,'(I8,3F16.0,2I8)',IOSTAT=IERROR) I_TEMP, D3_TEMP, I2_TEMP
             
              IF (IERROR.EQ.0) THEN
                MODEL_NUMP = MODEL_NUMP + 1
              ELSE
                GOTO 100
              END IF
            
            END DO
            
	    END IF !*NODE LIST CARD
	  
	  END DO !READ FILE WITH ARBITRARY DO
      
      CLOSE(10)
      
      RETURN
      
	  END SUBROUTINE
      
      
      
      
      
      

      SUBROUTINE COUNT_ELEMENTS(MODEL_FILE, & !IN
                       MODEL_NUMEL, MODEL_NUMBLOCK, MODEL_NUMEL_BLOCKS) !OUT
      
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! COUNT THE NUMBER OF NODES IN A *.K FILE
	  !
	  IMPLICIT NONE
      ! 
      ! GLOBAL IN-OUT
      !
      CHARACTER(50), INTENT(IN):: MODEL_FILE
      INTEGER, INTENT(OUT):: MODEL_NUMEL, MODEL_NUMBLOCK,MODEL_NUMEL_BLOCKS(1000)
      !
	  !
	  ! LOCAL I/O FLAGS ETC
	  !
	  CHARACTER(100):: CTEMP
	  INTEGER:: IERROR
      !
      INTEGER:: I_TEMP, I2_TEMP(2), I3_TEMP(3), I4_TEMP(4), I5_TEMP(5)
      INTEGER:: I6_TEMP(6), I7_TEMP(7), I8_TEMP(8), I9_TEMP(9), I10_TEMP(8)
      !
      DOUBLE PRECISION:: D_TEMP, D2_TEMP(2), D3_TEMP(3), D4_TEMP(4)
      DOUBLE PRECISION:: D5_TEMP(5), D6_TEMP(6), D7_TEMP(7), D8_TEMP(8)
      
      
	  OPEN(10,FILE=MODEL_FILE,STATUS='OLD',ACTION='READ',iostat=IERROR)
      !
      ! COUNT THE NUMBER OF ELEMENTS
      !
      MODEL_NUMBLOCK = 0
      MODEL_NUMEL = 0
      !
	  DO 
	  
100     READ(10,'(A100)',iostat=IERROR) CTEMP
	  
	    IF (IERROR.EQ.(-1)) EXIT
	    CTEMP=TRIM(CTEMP)
	    IF (CTEMP.eq.'*ELEMENT_SOLID') THEN
        MODEL_NUMBLOCK= MODEL_NUMBLOCK + 1

	        DO 
200           READ(10,'(A100)',iostat=IERROR) CTEMP
	          IF (CTEMP.eq.'*ELEMENT_SOLID') THEN
                 MODEL_NUMBLOCK= MODEL_NUMBLOCK + 1
                 GOTO 200
              ELSE
                READ(CTEMP,'(10I8)',IOSTAT=IERROR) I10_TEMP
              END IF
              
              IF (IERROR.EQ.0) THEN
                MODEL_NUMEL = MODEL_NUMEL + 1
                MODEL_NUMEL_BLOCKS(MODEL_NUMBLOCK) = MODEL_NUMEL_BLOCKS(MODEL_NUMBLOCK) + 1
              ELSE
                GOTO 100
              END IF
            
            END DO
            
	    END IF !*ELEMENT_SOLID LIST CARD
	  
	  END DO !READ FILE WITH ARBITRARY DO
      
      CLOSE(10)
      
      RETURN
      
	  END SUBROUTINE
      
      
      
      


      SUBROUTINE COUNT_NODE_SETS(MODEL_FILE, & !IN
                       NUM_NODESET,MAX_NINODESET) !OUT
      
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! COUNT THE NUMBER OF NODE SETS IN A *.K FILE
	  !
	  IMPLICIT NONE
      ! 
      ! GLOBAL IN-OUT
      !
      CHARACTER(50), INTENT(IN):: MODEL_FILE
      INTEGER, INTENT(OUT):: NUM_NODESET, MAX_NINODESET
      !
	  !
	  ! LOCAL I/O FLAGS ETC
	  !
	  CHARACTER(100):: CTEMP
	  INTEGER:: IERROR
	  !
	  INTEGER:: TEMP_NINODESET
      !
      INTEGER:: I_TEMP, I2_TEMP(2), I3_TEMP(3), I4_TEMP(4), I5_TEMP(5)
      INTEGER:: I6_TEMP(6), I7_TEMP(7), I8_TEMP(8), I9_TEMP(9), I10_TEMP(8)
      !
      DOUBLE PRECISION:: D_TEMP, D2_TEMP(2), D3_TEMP(3), D4_TEMP(4)
      DOUBLE PRECISION:: D5_TEMP(5), D6_TEMP(6), D7_TEMP(7), D8_TEMP(8)
      !
      INTEGER:: NUM_READ
      !
      
      
	  OPEN(10,FILE=MODEL_FILE,STATUS='OLD',ACTION='READ',iostat=IERROR)
      !
      ! COUNT THE NUMBER OF ELEMENTS
      !
      NUM_NODESET = 0
      MAX_NINODESET = 0
      !
	  DO 
	  
100     READ(10,'(A50)',iostat=IERROR) CTEMP
	  
	    IF (IERROR.EQ.(-1)) EXIT
	    CTEMP=TRIM(CTEMP)
	    IF (CTEMP.eq.'*SET_NODE_LIST') THEN
        NUM_NODESET= NUM_NODESET + 1
        TEMP_NINODESET = 0
          READ(10,*) !SET ID

	        DO 
200           READ(10,'(A100)',iostat=IERROR) CTEMP
	          IF (CTEMP.eq.'*SET_NODE_LIST') THEN
                 READ(10,*) !SET ID
                 TEMP_NINODESET = 0
                 NUM_NODESET= NUM_NODESET + 1
                 GOTO 200
              ELSE
                I8_TEMP = 0
                READ(CTEMP,'(8I10)',IOSTAT=IERROR) I8_TEMP
              END IF
              
              IF (IERROR.EQ.0) THEN
                IF (I8_TEMP(8).NE.0) NUM_READ = 8
                IF (I8_TEMP(8).EQ.0) NUM_READ = 7
                IF (I8_TEMP(7).EQ.0) NUM_READ = 6
                IF (I8_TEMP(6).EQ.0) NUM_READ = 5
                IF (I8_TEMP(5).EQ.0) NUM_READ = 4
                IF (I8_TEMP(4).EQ.0) NUM_READ = 3
                IF (I8_TEMP(3).EQ.0) NUM_READ = 2
                IF (I8_TEMP(2).EQ.0) NUM_READ = 1
                IF (I8_TEMP(1).EQ.0) NUM_READ = 0
                !NO ERROR, COUNT THE NUMBER OF NODES
                TEMP_NINODESET = TEMP_NINODESET + NUM_READ
                
                IF (TEMP_NINODESET.GT.MAX_NINODESET) THEN
                  MAX_NINODESET = TEMP_NINODESET
                END IF
                
              ELSE
                GOTO 100
              END IF
            
            END DO
            
	    END IF !*ELEMENT_SOLID LIST CARD
	  
	  END DO !READ FILE WITH ARBITRARY DO
      
      CLOSE(10)
      
      RETURN
      
	  END SUBROUTINE
      

SUBROUTINE COUNT_FIXITY2_STEPS(MODEL_FILE, &! IN
                          FIXITY2_STEPS) !OUT

  !GC: FUNCTION OF THIS SUBROUTINE: READ IN THE NUMBER OF STEPS FOR FIXITY 2
  IMPLICIT NONE

  CHARACTER(50), INTENT(IN):: MODEL_FILE
  INTEGER,INTENT(OUT)::FIXITY2_STEPS

  CHARACTER(100):: CTEMP
  INTEGER:: IERROR
  INTEGER::FIXITY(3)

  FIXITY2_STEPS = 0

  OPEN(10,FILE='Control.dat',STATUS='OLD',ACTION='READ',iostat=IERROR)

    DO 
    
100     READ(10,'(A100)',iostat=IERROR) CTEMP
    
      IF (IERROR.EQ.(-1)) EXIT
      CTEMP=TRIM(CTEMP)
      IF (CTEMP.eq.'*BOUNDARY CONDITION') THEN
        READ(10,*) !COMMENT
        READ(10,*) !NODESET NAME
        READ(10,*) !FIXITY
        READ(10,*,iostat=IERROR) FIXITY
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE VARIABLE FIXITY IN THE CARD *ORDERED BOUNDARY CONDITION IN THE CONTROL FILE',2) !2 = USER ERROR
        
		IF(FIXITY(1).EQ.2 .OR. FIXITY(2).EQ.2 .OR. FIXITY(3).EQ.2) THEN
			READ(10,*) !STEPS
			READ(10,*) FIXITY2_STEPS
			IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE VARIABLE FIXITY STEPS IN THE CARD *ORDERED BOUNDARY CONDITION IN THE CONTROL FILE',2) !2 = USER ERROR
		END IF
      ELSE IF (CTEMP.eq.'*ORDERED BOUNDARY CONDITION') THEN
        READ(10,*) !COMMENT
        READ(10,*) !NODESET NAME
        READ(10,*) !FIXITY
        READ(10,*,iostat=IERROR) FIXITY
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE VARIABLE FIXITY IN THE CARD *ORDERED BOUNDARY CONDITION IN THE CONTROL FILE',2) !2 = USER ERROR
        
		IF(FIXITY(1).EQ.2 .OR. FIXITY(2).EQ.2 .OR. FIXITY(3).EQ.2) THEN
			READ(10,*) !STEPS
			READ(10,*) FIXITY2_STEPS
			IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE VARIABLE FIXITY STEPS IN THE CARD *ORDERED BOUNDARY CONDITION IN THE CONTROL FILE',2) !2 = USER ERROR
		END IF
      END IF

    
    END DO !READ FILE WITH ARBITRARY DO

END SUBROUTINE