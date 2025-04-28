
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
	  
      SUBROUTINE READ_CONTROL(NUMP,MODEL_NODE_SET_LIST,MODEL_NODE_SET_LENGTH, &
                               MODEL_NODE_SET_ID,NUM_NODESET,MAX_NINODESET, FIXITY2_STEPS,TIME_END,&
                               MODEL_EBC,MODEL_NONZERO_EBC,FIXITY2_TIME,FIXITY2_NONZERO_EBC,MODEL_VINIT, MODEL_MAT_TYPE,MODEL_PROP,MODEL_BODY_ID, &
                               MODEL_NORM_WIN,MODEL_SET_NAMES,MODEL_BODYFORCE)
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! SET THE BOUNDARY CONDITIONS
	  !
	  IMPLICIT NONE
      ! 
      ! GLOBAL IN-OUT
      INTEGER, INTENT(IN):: NUMP
      INTEGER, INTENT(IN):: NUM_NODESET
      INTEGER, INTENT(IN):: MAX_NINODESET
      INTEGER, INTENT(IN):: MODEL_NODE_SET_LIST(NUM_NODESET,MAX_NINODESET)
      INTEGER, INTENT(IN):: MODEL_NODE_SET_LENGTH(NUM_NODESET)
      INTEGER, INTENT(IN):: MODEL_NODE_SET_ID(NUM_NODESET)
      INTEGER, INTENT(IN):: FIXITY2_STEPS !GC
      DOUBLE PRECISION,INTENT(IN)::TIME_END
      
      
      DOUBLE PRECISION, INTENT(OUT):: MODEL_VINIT(3,NUMP)
      INTEGER, INTENT(OUT):: MODEL_MAT_TYPE(NUMP)
      DOUBLE PRECISION, INTENT(OUT):: MODEL_PROP(30,NUMP)
      INTEGER, INTENT(OUT):: MODEL_EBC(3,NUMP)
      DOUBLE PRECISION, INTENT(OUT)::MODEL_NONZERO_EBC(3,NUMP)
      INTEGER, INTENT(OUT):: MODEL_BODY_ID(NUMP)
      DOUBLE PRECISION, INTENT(OUT):: MODEL_NORM_WIN(NUMP)
	  CHARACTER(100), INTENT(OUT):: MODEL_SET_NAMES(1000)

      DOUBLE PRECISION, INTENT(OUT)::FIXITY2_TIME(FIXITY2_STEPS,2),FIXITY2_NONZERO_EBC(FIXITY2_STEPS,3) !GC
      !
	  !
	  ! LOCAL I/O FLAGS ETC
	  !
	  CHARACTER(100):: CTEMP,CTEMP2
	  INTEGER:: IERROR
	  CHARACTER(100):: SET_NAME
	  !
	  INTEGER:: SET_NO, FIXITY(3), SET_INDEX, NODE_ID
	  DOUBLE PRECISION:: VELOCITIES(3), NORM_WINDOW
	  INTEGER:: MAT_TYPE, NUM_PROP, BODY_ID
	  DOUBLE PRECISION:: PROPS(30)
      INTEGER::I, J
      !
      INTEGER:: I_TEMP, I2_TEMP(2), I3_TEMP(3), I4_TEMP(4), I5_TEMP(5)
      INTEGER:: I6_TEMP(6), I7_TEMP(7), I8_TEMP(8), I9_TEMP(9), I10_TEMP(8)
      !
      DOUBLE PRECISION:: D_TEMP, D2_TEMP(2), D3_TEMP(3), D4_TEMP(4)
      DOUBLE PRECISION:: D5_TEMP(5), D6_TEMP(6), D7_TEMP(7), D8_TEMP(8)
      LOGICAL:: READ_BC, READ_IC, READ_SET_NAMES, READ_MAT, READ_WIN,READ_BODY, &
                READ_ORDERED_BC
      
      CHARACTER(8):: SNUM_MISSING
      CHARACTER(100):: CHAR41
      INTEGER:: NUM_MISSING
      
      LOGICAL:: USE_STRING
	  
      DOUBLE PRECISION::PHI,C,PSI,PHI_RAD,PSI_RAD,SRT32,Q_PHI,K_PHI,Q_PSI
      
      INTEGER:: EBCTEMP,DP_TYPE
			  
      DOUBLE PRECISION:: EPOISS, EYOUNG, LAMDA, MU, LAMDA_PLUS_2MU
      
      
      !0329
      DOUBLE PRECISION:: BODYFORCE(3)
      !DOUBLE PRECISION, INTENT(OUT):: MODEL_BODYFORCE(3)
      !0702
      DOUBLE PRECISION, INTENT(OUT):: MODEL_BODYFORCE(3,NUMP)
      
      READ_SET_NAMES = .FALSE.
      FIXITY2_TIME = 0.d0 
      FIXITY2_NONZERO_EBC = 0.d0 
      
	  OPEN(10,FILE='Control.dat',STATUS='OLD',ACTION='READ',iostat=IERROR)
	  

	  DO 
	  
        READ(10,'(A50)',iostat=IERROR) CTEMP
	  
	    IF (IERROR.EQ.(-1)) EXIT
	    CTEMP=TRIM(CTEMP)
	    IF (CTEMP.eq.'*SET NAME') THEN
          
          !THIS PART IS PURELY FOR CONVINIENCE. CHANGED IT SO IT
          !IS JUST TO GIVE NODESETS NAMES. WILL HAVE TO CHECK
          !FOR UNDEFINED BEHAVOIR NOW FOR BODY IDs WHEN MULTIPLE
          !VALUES ARE ASSIGNED TO A SINGLE NODE
       
	  
      
          READ_SET_NAMES = .TRUE.

          READ(10,*) !COMMENT
          READ(10,* ,iostat=IERROR) SET_NO
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE VARIABLE SET_NO IN THE CARD *SET NAME IN THE CONTROL FILE',2) !2 = USER ERROR
          
          READ(10,*) !COMMENT
          READ(10,*,iostat=IERROR) SET_NAME
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE VARIABLE SET_NAME IN THE CARD *SET NAME IN THE CONTROL FILE',2) !2 = USER ERROR
          
	      SET_NAME=TRIM(SET_NAME)
          MODEL_SET_NAMES(SET_NO) = SET_NAME
          
          SET_INDEX = 0
          DO I=1,NUM_NODESET
            IF (MODEL_NODE_SET_ID(I).EQ.SET_NO) SET_INDEX = I
          END DO
          
          IF (SET_INDEX.EQ.0) THEN
            !THE SET WAS NOT IN THE LIST!
	        IF (IERROR.NE.0) CALL EXIT_PROGRAM('SET NUMBER NOT FOUND IN MODEL FOR SET IN CARD *SET NAME',2)
          END IF
          
          
        
        END IF
        
	  END DO !READ FILE WITH ARBITRARY DO
	  
	  CLOSE(10)
      
      IF (.NOT.READ_SET_NAMES) THEN
      
        CALL WARN('SET NAMES WERE NOT FOUND IN THE CONTROL FILE')
        
      END IF
        
      READ_BC = .FALSE.
      READ_ORDERED_BC = .FALSE.
      READ_IC = .FALSE.
      READ_MAT = .FALSE.
      READ_WIN = .FALSE.
      READ_BODY = .FALSE.
      
      MODEL_EBC = -1
      MODEL_NONZERO_EBC = 0.D0
      MODEL_VINIT = 0.0d0
      MODEL_MAT_TYPE = -1
      MODEL_PROP = 0.0d0
      MODEL_NORM_WIN = 0.0d0
      MODEL_BODY_ID = -1
      MODEL_BODYFORCE = 0.0d0
      
      
	  OPEN(10,FILE='Control.dat',STATUS='OLD',ACTION='READ',iostat=IERROR)
	  

	  DO 
	  
100     READ(10,'(A50)',iostat=IERROR) CTEMP
	  
	    IF (IERROR.EQ.(-1)) EXIT
	    CTEMP=TRIM(CTEMP)
	    IF (CTEMP.eq.'*BOUNDARY CONDITION') THEN
          READ_BC = .TRUE.
          USE_STRING = .FALSE.
          
          
          READ(10,*) !COMMENT
          !
          ! NOW, DETERMINE IF A NUMBER OF NAME WAS GIVEN
          !
          CALL READ_SET_NO(10,MODEL_SET_NAMES,SET_INDEX, NUM_NODESET,MODEL_NODE_SET_ID,'*BOUNDARY CONDITION')
          
          READ(10,*) !COMMENT
          READ(10,*,iostat=IERROR) FIXITY
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE VARIABLE FIXITY IN THE CARD *BOUNDARY CONDITION IN THE CONTROL FILE',2) !2 = USER ERROR
          
          IF(FIXITY(1).EQ.2 .OR. FIXITY(2).EQ.2 .OR. FIXITY(3).EQ.2) THEN !HAS NONZERO BC VALUES, MULTI-STEPS (GC)
              READ(10,*) !COMMENT STEPS GC
              READ(10,*) !FIXITY2_STEPS: READ IN ALREADY
              IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE PRESCRIBED DISPLACEMENT STEPS IN THE CARD *BOUNDARY CONDITION IN THE CONTROL FILE',2) !2 = USER ERROR
              
              READ(10,*) !COMMENT DISPLACEMENT
              DO I = 1,FIXITY2_STEPS                
                READ(10,*,iostat=IERROR) FIXITY2_NONZERO_EBC(I,:)
                IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE PRESCRIBED DISPLACEMENT IN THE CARD *BOUNDARY CONDITION IN THE CONTROL FILE',2) !2 = USER ERROR
              END DO

              READ(10,*) !COMMENT TIME
              DO I = 1,FIXITY2_STEPS                
                READ(10,*,iostat=IERROR) FIXITY2_TIME(I,:)
                IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE PRESCRIBED DISPLACEMENT DURATION IN THE CARD *BOUNDARY CONDITION IN THE CONTROL FILE',2) !2 = USER ERROR
              END DO
          ENDIF
          
          DO I = 1,MODEL_NODE_SET_LENGTH(SET_INDEX)
            
            NODE_ID = MODEL_NODE_SET_LIST(SET_INDEX,I)
            
            CALL ASSIGN_INT_VALUE_TO_NODE( MODEL_EBC(1,NODE_ID),FIXITY(1), &
            'INPUT ERROR: MULTIPLE VALUES OF X-FIXITY ASSIGNED TO SINGLE NODE')
            CALL ASSIGN_DBL_VALUE_TO_NODE( MODEL_NONZERO_EBC(1,NODE_ID),D3_TEMP(1), &
            'INPUT ERROR: MULTIPLE VALUES OF A NON-ZERO X DISP BC CONDITION ASSIGNED TO SINGLE NODE')
            
            CALL ASSIGN_INT_VALUE_TO_NODE( MODEL_EBC(2,NODE_ID),FIXITY(2), &
            'INPUT ERROR: MULTIPLE VALUES OF Y-FIXITY ASSIGNED TO SINGLE NODE')
            CALL ASSIGN_DBL_VALUE_TO_NODE( MODEL_NONZERO_EBC(2,NODE_ID),D3_TEMP(2), &
            'INPUT ERROR: MULTIPLE VALUES OF A NON-ZERO Y DISP BC CONDITION ASSIGNED TO SINGLE NODE')
            
            CALL ASSIGN_INT_VALUE_TO_NODE( MODEL_EBC(3,NODE_ID),FIXITY(3), &
            'INPUT ERROR: MULTIPLE VALUES OF Z-FIXITY ASSIGNED TO SINGLE NODE')
            CALL ASSIGN_DBL_VALUE_TO_NODE( MODEL_NONZERO_EBC(3,NODE_ID),D3_TEMP(3), &
            'INPUT ERROR: MULTIPLE VALUES OF A NON-ZERO Z DISP BC CONDITION ASSIGNED TO SINGLE NODE')            
          
          END DO
        
        
          ELSEIF (CTEMP.eq.'*ORDERED BOUNDARY CONDITION') THEN
          READ_ORDERED_BC = .TRUE.
          
          READ(10,*) !COMMENT
          !
          ! NOW, DETERMINE IF A NUMBER OF NAME WAS GIVEN
          !
          CALL READ_SET_NO(10,MODEL_SET_NAMES,SET_INDEX, NUM_NODESET,MODEL_NODE_SET_ID,'*BOUNDARY CONDITION')
          
          READ(10,*) !COMMENT
          READ(10,*,iostat=IERROR) FIXITY
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE VARIABLE FIXITY IN THE CARD *ORDERED BOUNDARY CONDITION IN THE CONTROL FILE',2) !2 = USER ERROR
          
          IF(FIXITY(1).EQ.2 .OR. FIXITY(2).EQ.2 .OR. FIXITY(3).EQ.2) THEN !HAS NONZERO BC VALUES, MULTI-STEPS (GC)
              READ(10,*) !COMMENT STEPS GC
              READ(10,*) !FIXITY2_STEPS: READ IN ALREADY
              IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE PRESCRIBED DISPLACEMENT STEPS IN THE CARD *BOUNDARY CONDITION IN THE CONTROL FILE',2) !2 = USER ERROR
              
              READ(10,*) !COMMENT DISPLACEMENT
              DO I = 1,FIXITY2_STEPS                
                READ(10,*,iostat=IERROR) FIXITY2_NONZERO_EBC(I,:)
                IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE PRESCRIBED DISPLACEMENT IN THE CARD *BOUNDARY CONDITION IN THE CONTROL FILE',2) !2 = USER ERROR
              END DO

              READ(10,*) !COMMENT TIME
              DO I = 1,FIXITY2_STEPS                
                READ(10,*,iostat=IERROR) FIXITY2_TIME(I,:)
                IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE PRESCRIBED DISPLACEMENT DURATION IN THE CARD *BOUNDARY CONDITION IN THE CONTROL FILE',2) !2 = USER ERROR
              END DO
              
             !GC: CHECK INPUT OF PRESCRIBED DISPLACEMENT
              DO I = 1,FIXITY2_STEPS
                DO J = 1,2
                  IF(FIXITY2_TIME(I,J) .GT. TIME_END) CALL EXIT_PROGRAM('THERE WAS AN ERROR THE PRESCRIBED DISPLACEMENT TIME DURATION. TIME DURATION IS GREAT THAN TIME END',2)
                END DO
                
                IF (FIXITY2_TIME(I,1) .GT. FIXITY2_TIME(I,2)) CALL EXIT_PROGRAM('THERE WAS AN ERROR THE PRESCRIBED DISPLACEMENT TIME DURATION. TIME START IS GREAT THAN TIME END',2)

              END DO             
              
          ENDIF
          
          DO I = 1,MODEL_NODE_SET_LENGTH(SET_INDEX)
            
            NODE_ID = MODEL_NODE_SET_LIST(SET_INDEX,I)
            MODEL_EBC(1:3,NODE_ID) = FIXITY(1:3)
            MODEL_NONZERO_EBC(1:3,NODE_ID) = D3_TEMP(1:3)
            
          END DO
        
          WRITE(CHAR41,'(A56,I12)') 'USING ORDERED BOUNDARY CONDITION, APPLYING SET_INDEX = ', SET_INDEX
          CALL WARN(CHAR41)
          
          
          
          ELSEIF (CTEMP.eq.'*BODY DEFINITION') THEN
          READ_BODY = .TRUE.
          
          READ(10,*) !COMMENT
          CALL READ_SET_NO(10,MODEL_SET_NAMES,SET_INDEX, NUM_NODESET,MODEL_NODE_SET_ID,'*BODY DEFINITION')
          
          READ(10,*) !COMMENT
          READ(10,*,iostat=IERROR) BODY_ID
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE VARIABLE BODY_ID IN THE CARD *BODY DEFINITION IN THE CONTROL FILE',2) !2 = USER ERROR
          IF (BODY_ID.GT.1000) THEN
            CALL EXIT_PROGRAM('INPUT ERROR: PLEASE CHOOSE A BODY ID THAT IS LESS THAN OR EQUAL TO 1000',2)
          END IF
          
          !0329 READ BODYFORCE
          READ(10,*) !COMMENT
          READ(10,*,iostat=IERROR) BODYFORCE
          
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE VARIABLE BODYFORCE (NEW VARIABLE, DBL(3) !) IN THE CARD *BODY DEFINITION IN THE CONTROL FILE',2) !2 = USER ERROR
          
          !KC
          !MODEL_BODYFORCE = BODYFORCE
          DO I = 1,MODEL_NODE_SET_LENGTH(SET_INDEX)
            
            NODE_ID = MODEL_NODE_SET_LIST(SET_INDEX,I)
            
            CALL ASSIGN_INT_VALUE_TO_NODE( MODEL_BODY_ID(NODE_ID),BODY_ID, &
            'INPUT ERROR: MULTIPLE VALUES OF BODY_ID ASSIGNED TO SINGLE NODE')
            
            !0702 KC
            CALL ASSIGN_DBL_VALUE_TO_NODE( MODEL_BODYFORCE(1,NODE_ID),BODYFORCE(1), &
            'INPUT ERROR: MULTIPLE VALUES OF A NON-ZERO VX INITIAL CONDITION ASSIGNED TO SINGLE NODE')
            
            CALL ASSIGN_DBL_VALUE_TO_NODE( MODEL_BODYFORCE(2,NODE_ID),BODYFORCE(2), &
            'INPUT ERROR: MULTIPLE VALUES OF A NON-ZERO VY INITIAL CONDITION ASSIGNED TO SINGLE NODE')
            
            CALL ASSIGN_DBL_VALUE_TO_NODE( MODEL_BODYFORCE(3,NODE_ID),BODYFORCE(3), &
            'INPUT ERROR: MULTIPLE VALUES OF A NON-ZERO VZ INITIAL CONDITION ASSIGNED TO SINGLE NODE')
            
          
          END DO
          
          ELSEIF (CTEMP.eq.'*INITIAL CONDITION') THEN
          READ_IC = .TRUE.

          READ(10,*) !COMMENT
          !
          ! NOW, DETERMINE IF A NUMBER OF NAME WAS GIVEN
          !
          CALL READ_SET_NO(10,MODEL_SET_NAMES,SET_INDEX, NUM_NODESET,MODEL_NODE_SET_ID,'*INITIAL CONDITION')
          !
          READ(10,*) !COMMENT
          READ(10,*,iostat=IERROR) VELOCITIES
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE VARIABLE VELOCITIES IN THE CARD *INITIAL CONDITION IN THE CONTROL FILE',2) !2 = USER ERROR
          
          DO I = 1,MODEL_NODE_SET_LENGTH(SET_INDEX)
            
            NODE_ID = MODEL_NODE_SET_LIST(SET_INDEX,I)
            
            CALL ASSIGN_DBL_VALUE_TO_NODE( MODEL_VINIT(1,NODE_ID),VELOCITIES(1), &
            'INPUT ERROR: MULTIPLE VALUES OF A NON-ZERO VX INITIAL CONDITION ASSIGNED TO SINGLE NODE')
            
            CALL ASSIGN_DBL_VALUE_TO_NODE( MODEL_VINIT(2,NODE_ID),VELOCITIES(2), &
            'INPUT ERROR: MULTIPLE VALUES OF A NON-ZERO VY INITIAL CONDITION ASSIGNED TO SINGLE NODE')
            
            CALL ASSIGN_DBL_VALUE_TO_NODE( MODEL_VINIT(3,NODE_ID),VELOCITIES(3), &
            'INPUT ERROR: MULTIPLE VALUES OF A NON-ZERO VZ INITIAL CONDITION ASSIGNED TO SINGLE NODE')
            
          
          END DO
          
        ELSEIF (CTEMP.eq.'*WINDOW DEFINITION') THEN
        
          READ_WIN = .TRUE.

          READ(10,*) !COMMENT
          !
          ! NOW, DETERMINE IF A NUMBER OF NAME WAS GIVEN
          !
          CALL READ_SET_NO(10,MODEL_SET_NAMES,SET_INDEX, NUM_NODESET,MODEL_NODE_SET_ID,'*WINDOW DEFINITION')
          !
          READ(10,*) !COMMENT
          READ(10,*,iostat=IERROR) NORM_WINDOW
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE VARIABLE NORM_WINDOW IN THE CARD *WINDOW DEFINITION IN THE CONTROL FILE',2) !2 = USER ERROR
          
          DO I = 1,MODEL_NODE_SET_LENGTH(SET_INDEX)
            
            NODE_ID = MODEL_NODE_SET_LIST(SET_INDEX,I)
            
            CALL ASSIGN_DBL_VALUE_TO_NODE( MODEL_NORM_WIN(NODE_ID),NORM_WINDOW, &
            'INPUT ERROR: MULTIPLE VALUES OF A NORMALIZED WINDOW WAS ASSIGNED TO SINGLE NODE')
            
          END DO
          
          
        ELSEIF (CTEMP.eq.'*MATERIALS') THEN
          READ_MAT = .TRUE.
          
          READ(10,*) !COMMENT
          !
          ! NOW, DETERMINE IF A NUMBER OF NAME WAS GIVEN
          !
          CALL READ_SET_NO(10,MODEL_SET_NAMES,SET_INDEX, NUM_NODESET,MODEL_NODE_SET_ID,'*MATERIALS')
          !
          READ(10,*) !COMMENT
          READ(10,*,iostat=IERROR) MAT_TYPE
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE VARIABLE MAT_TYPE IN THE CARD *MATERIALS IN THE CONTROL FILE',2) !2 = USER ERROR
          READ(10,*) !COMMENT
          !
          ! READ IN ALL THE PROPERTIES
          !
          READ(10,*,iostat=IERROR) PROPS(1:5)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *MATERIALS, VALUE ROW 1',2)
          READ(10,*,iostat=IERROR) PROPS(6:10)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *MATERIALS, VALUE ROW 2',2)
          READ(10,*,iostat=IERROR) PROPS(11:15)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *MATERIALS, VALUE ROW 3',2)
          READ(10,*,iostat=IERROR) PROPS(16:20)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *MATERIALS, VALUE ROW 4',2)
          READ(10,*,iostat=IERROR) PROPS(21:25)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *MATERIALS, VALUE ROW 5',2)
          READ(10,*,iostat=IERROR) PROPS(26:30)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *MATERIALS, VALUE ROW 6',2)
          
          
          DO I = 1,MODEL_NODE_SET_LENGTH(SET_INDEX)
            
            NODE_ID = MODEL_NODE_SET_LIST(SET_INDEX,I)
            
            CALL ASSIGN_INT_VALUE_TO_NODE(  MODEL_MAT_TYPE(NODE_ID),MAT_TYPE, &
            'INPUT ERROR: MULTIPLE VALUES OF A MATERIAL ID ASSIGNED TO SINGLE NODE')
            
            CALL GET_NUM_PROP(MAT_TYPE,NUM_PROP)
            !
            ! PREP THE DATA: NOTE, ONLY DRUCKER PRAGER NEEDS MODIFICATION
            !
			  IF (MAT_TYPE.EQ.3) THEN !DRUCKER PRAGER WITH TENSION CUTTOFF AND DAMAGE
              
	              !POISS = PROPS(1)
	              !YOUNG = PROPS(2)
	              !Q_PHI = PROPS(4)
	              !K_PHI = PROPS(5)
	              !Q_PSI = PROPS(6)
	              !T_CUT = PROPS(7)
	              !K_I = PROPS(8) !DAM PAR #1
	              !K_C = PROPS(9) !DAM PAR #2
	              !DAM_MAX = PROPS(10)
                  DP_TYPE = PROPS(4)
				  PHI = PROPS(5)
				  C = PROPS(6)
				  PSI = PROPS(7)
				  !
				  PHI_RAD = PHI*3.141592653589793d0/180.0d0
				  PSI_RAD = PSI*3.141592653589793d0/180.0d0
				  SRT32 = SQRT(3.0d0/2.d0)
                  
                  IF (DP_TYPE.EQ.1) THEN
                !  
                !  Circumscribes the Mohr–Coulomb yield surface
                !  
				  Q_PHI = 2.0d0*SIN(PHI_RAD)/(SRT32*(3.0d0 - SIN(PHI_RAD))) *3.d0
				  K_PHI = 6.0d0*C*COS(PHI_RAD)/(SRT32*(3.0d0 - SIN(PHI_RAD)))
				  !Q_PSI = 2.0d0*SIN(PSI_RAD)/(SRT3*(3.0d0 - SIN(PSI_RAD))) 

                  
                ELSEIF (DP_TYPE.EQ.2) THEN
                !
                !  Middle circumscribes the Mohr–Coulomb yield surface
                !
				  Q_PHI = 2.0d0*SIN(PHI_RAD)/(SRT32*(3.0d0 + SIN(PHI_RAD)))  *3.d0
				  K_PHI = 6.0d0*C*COS(PHI_RAD)/(SRT32*(3.0d0 + SIN(PHI_RAD)))
				  !Q_PSI = 2.0d0*SIN(PSI_RAD)/(SRT3*(3.0d0 + SIN(PSI_RAD)))

                  
                ELSEIF (DP_TYPE.EQ.3) THEN
                !
                !   Inscribes the Mohr–Coulomb yield surface
                !
				  Q_PHI = SIN(PHI_RAD)/SQRT(9.0d0 + 3.0d0*SIN(PHI_RAD)*SIN(PHI_RAD))  /SQRT(0.5d0)   *3.d0
				  K_PHI = 3.0d0*C*COS(PHI_RAD)/SQRT(9.0d0 + 3.0d0*SIN(PHI_RAD)*SIN(PHI_RAD)) /SQRT(0.5d0)
				  !Q_PSI = SIN(PSI_RAD)/SQRT(9.0d0 + 3.0d0*SIN(PSI_RAD)*SIN(PSI_RAD))   
                  
                 ELSEIF (DP_TYPE.EQ.4) THEN
                !  
                !  Directly input A and B values
                !                   
                  Q_PHI = PROPS(5)   *3.d0
                  K_PHI = PROPS(6)
                  Q_PSI = PROPS(7)
                  
                ENDIF
                
                 !DOES NOT ALLOW DILATATION FOR NOW
                 IF(ABS(Q_PSI) .GT. 1.0E-09) THEN
                   CALL WARN('DILATATION EFFECT IS NOT IMPLEMENTED IN DRUCKER-PRAGER MATERIAL MODEL, Q_PSI IS RESET TO ZERO')
                 END IF
                
                 Q_PSI = 0.D0
                
				!
                !  write(*,*) DP_TYPE,Q_PHI,K_PHI
                !  pause
				  PROPS(25) = Q_PHI
				  PROPS(26) = K_PHI
				  PROPS(27) = Q_PSI
                  
                  
				  !PROPS(7) = Q_PHI
				  !PROPS(8) = K_PHI
				  !PROPS(9) = Q_PSI                  
                  
			  END IF
			  
	              EPOISS = PROPS(1)
	              EYOUNG = PROPS(2)
	              
	              LAMDA = EYOUNG*EPOISS/((1.0d0+EPOISS)*(1.0d0-2.0d0*EPOISS))
	              
	              MU = EYOUNG/(2.0d0*(1.0d0+EPOISS))
	              
	              LAMDA_PLUS_2MU = LAMDA + 2.0d0*MU
	              
				  PROPS(28) = LAMDA
				  PROPS(29) = MU
				  PROPS(30) = LAMDA_PLUS_2MU
				  
                  MODEL_PROP(:,NODE_ID)=PROPS
          
          END DO
          
          
          
        END IF !IDENTIFY CARD IF STRUCTURE
        
        
        
	  END DO !READ FILE WITH ARBITRARY DO
	  
	  CLOSE(10)
      
      
      !
      ! CHECK ASSIGNMENTS FOR MODEL FOR CONSISTENCY
      !
      
      IF ((READ_BC).AND.(READ_ORDERED_BC)) THEN
      
        CALL EXIT_PROGRAM('EITHER ORDERED OR NON-ORDERED BOUNDARY CONDITIONS MUST BE USED',2)
        
      ELSEIF ((.NOT.READ_BC).AND.(.NOT.READ_ORDERED_BC)) THEN
        CALL WARN('BOUNDARY CONDITIONS WERE NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
        CALL LOG_APPEND_SPACE('ALL DOFS = FREE')
        WRITE(*,*) 'ALL DOFS = FREE'
      ELSE
      
      
      NUM_MISSING = 0
      DO I=1, NUMP
        DO J=1,3
        EBCTEMP = MODEL_EBC(J,I)
        IF (EBCTEMP.EQ.(-1)) THEN
           NUM_MISSING = NUM_MISSING + 1
           MODEL_EBC(J,I) = 0
        END IF
        END DO
        
      END DO
      
      IF (NUM_MISSING.GT.0) THEN
        WRITE(SNUM_MISSING,'(I8)') NUM_MISSING
        CALL LOG_APPEND('ESSENTIAL BOUNDARY CONDITIONS ARE MISSING FOR PART OF THE MODEL, NUMBER MISSING =')
        CALL LOG_APPEND(SNUM_MISSING)
        CALL WARN('SOME ESSENTIAL BOUNDARY CONDITIONS WERE NOT FOUND IN THE CONTROL FILE, DEFAULT OF FREE WILL BE USED')
        CALL LOG_APPEND_SPACE('NODE DOFS WITH UNASSIGNED ESSENTIAL BOUNDARY CONDITIONS = FREE')
      END IF
      
      END IF
      
      
      !
      ! CHECK ASSIGNMENT OF MODEL_BODY_ID
      !
      IF (.NOT.READ_BODY) THEN
      
        CALL WARN('BODY IDS WERE NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
        CALL LOG_APPEND_SPACE('ALL DOFS = BODY 1')
        WRITE(*,*) ' '
        WRITE(*,*) 'ALL DOFS = BODY 1'
        WRITE(*,*) ' '
        
      ELSE
      
      NUM_MISSING = 0
      DO I=1, NUMP
        BODY_ID = MODEL_BODY_ID(I)
        IF (BODY_ID.EQ.(-1)) THEN
           NUM_MISSING = NUM_MISSING + 1
           MODEL_BODY_ID(I) = 1
        END IF
      END DO
      
      IF (NUM_MISSING.GT.0) THEN
        WRITE(SNUM_MISSING,'(I8)') NUM_MISSING
        CALL LOG_APPEND('BODY ID DEFINITIONS ARE MISSING FOR PART OF THE MODEL, NUMBER MISSING =')
        CALL LOG_APPEND(SNUM_MISSING)
        CALL EXIT_PROGRAM('BODY ID DEFINITIONS ARE IS MISSING FOR PART OF THE MODEL, THEY MUST BE SPECIFIED (OPTIONS:SPECIFY VALUES FOR ALL NODES, OR DONT SPECIFY AT ALL)',2)
      END IF
      
      END IF
      
      
      
      IF (.NOT.READ_IC) THEN
        CALL WARN('INITIAL CONDITIONS WERE NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
        CALL LOG_APPEND_SPACE('ALL DOFS = AT REST')
        WRITE(*,*) ' '
        WRITE(*,*) 'ALL DOFS = AT REST'
        WRITE(*,*) ' '
      END IF
      
      IF (.NOT.READ_WIN) THEN
        CALL WARN('WINDOW DEFINITIONS WERE NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
        CALL LOG_APPEND_SPACE('ALL DOF NORMALIZED WINDOWS = 1.8')
        WRITE(*,*) ' '
        WRITE(*,*) 'ALL DOF NORMALIZED WINDOWS = 1.8'
        WRITE(*,*) ' '
        MODEL_NORM_WIN = 1.8d0
      ELSE
      
      NUM_MISSING = 0
      DO I=1, NUMP
        NORM_WINDOW = MODEL_NORM_WIN(I)
        IF (NORM_WINDOW.LT.(0.00001d0)) THEN
           NUM_MISSING = NUM_MISSING + 1
        END IF
      END DO
      
      IF (NUM_MISSING.GT.0) THEN
        WRITE(SNUM_MISSING,'(I8)') NUM_MISSING
        CALL LOG_APPEND('WINDOW DEFINITIONS ARE MISSING FOR PART OF THE MODEL, NUMBER MISSING =')
        CALL LOG_APPEND(SNUM_MISSING)
        CALL EXIT_PROGRAM('WINDOW DEFINITIONS ARE IS MISSING FOR PART OF THE MODEL, THEY MUST BE SPECIFIED (ALL NODES OR NONE)',2)
      END IF
      
      
      END IF
      
      
      
      
      IF (.NOT.READ_MAT) THEN
        CALL EXIT_PROGRAM('MATERIAL DEFINITIONS WERE NOT FOUND IN THE CONTROL FILE, THEY MUST BE SPECIFIED',2)
      ELSE
      
      NUM_MISSING = 0
      DO I=1, NUMP
        MAT_TYPE = MODEL_MAT_TYPE(I)
        IF (MAT_TYPE.EQ.(-1)) THEN
           NUM_MISSING = NUM_MISSING + 1
        END IF
      END DO
      
      
      IF (NUM_MISSING.GT.0) THEN
        WRITE(SNUM_MISSING,'(I8)') NUM_MISSING
        CALL LOG_APPEND('A MATERIAL DEFINITION IS MISSING FOR PART OF THE MODEL, NUMBER MISSING =')
        CALL LOG_APPEND(SNUM_MISSING)
        CALL EXIT_PROGRAM('A MATERIAL DEFINITION IS MISSING FOR PART OF THE MODEL, THEY MUST BE SPECIFIED',2)
      END IF
      
      END IF
      
      
	  CONTINUE
	  RETURN
	  END SUBROUTINE
      
      
      
          SUBROUTINE READ_SET_NO(FID,MODEL_SET_NAMES,SET_INDEX, NUM_NODESET,MODEL_NODE_SET_ID,CARD)
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! READ EITHER A STRING OR INTEGER FROM CONTROL AND RETURN THE SET INDEX
	  !
          IMPLICIT NONE
          
          INTEGER, INTENT(IN):: FID, NUM_NODESET
          INTEGER, INTENT(IN):: MODEL_NODE_SET_ID(NUM_NODESET)
	      CHARACTER(100), INTENT(IN):: MODEL_SET_NAMES(1000)
	      CHARACTER(*), INTENT(IN):: CARD
          
          INTEGER, INTENT(OUT):: SET_INDEX
          
	      CHARACTER(100):: CTEMP,CTEMP2
	      INTEGER:: IERROR
	      CHARACTER(100):: SET_NAME
          INTEGER:: I, SET_NO
          LOGICAL:: USE_STRING
	      CHARACTER(100):: CARD_WARN
          !
          USE_STRING = .FALSE.
          READ(FID,'(A50)') CTEMP
          READ(CTEMP,'(I1000)',iostat=IERROR) SET_NO
          IF (IERROR.NE.0) THEN 
            !
            ! READ IN WALL NAME
            !
            READ(CTEMP,'(A50)',iostat=IERROR) CTEMP2
            CTEMP = CTEMP2
            USE_STRING = .TRUE.
            IF (IERROR.NE.0) THEN
            WRITE(CARD_WARN,'(A29,A20)') 'A FATAL ERROR WAS DETECTED IN CARD' , CARD
            CALL WRITE_OUT(CARD_WARN)
            CALL EXIT_PROGRAM('INVALID SET NAME',1) 
            END IF
          END IF
          !
          ! DETERMINE SET NUMBER
          !
          IF (USE_STRING) THEN
          
          SET_INDEX = 0
            DO I=1,1000
              IF (MODEL_SET_NAMES(I).EQ.CTEMP) SET_INDEX = I
            END DO
              IF (SET_INDEX.EQ.0) THEN
                !THE SET WAS NOT IN THE LIST!
                WRITE(CARD_WARN,'(A29,A20)') 'A FATAL ERROR WAS DETECTED IN CARD' , CARD
                CALL WRITE_OUT(CARD_WARN)
	            CALL EXIT_PROGRAM('SET NAME NOT FOUND',2)
              END IF
            
            CONTINUE
            
          ELSE
          
              SET_INDEX = 0
              DO I=1,NUM_NODESET
                IF (MODEL_NODE_SET_ID(I).EQ.SET_NO) SET_INDEX = I
              END DO
              IF (SET_INDEX.EQ.0) THEN
                !THE SET WAS NOT IN THE LIST!
                WRITE(CARD_WARN,'(A29,A20)') 'A FATAL ERROR WAS DETECTED IN CARD' , CARD
                CALL WRITE_OUT(CARD_WARN)
	            CALL EXIT_PROGRAM('SET NUMBER NOT FOUND',2)
              END IF
          
          END IF
          
          RETURN
          
          END SUBROUTINE
