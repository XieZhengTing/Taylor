
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

      SUBROUTINE PRE_MODEL
	  
      !
      ! PURPOSE OF THIS LOGIC:
	  !
	  ! READ IN THE *.k FILE (FOR NOW USE LS-DYNA FILES) AND ASSIGN
	  ! THE INFO TO MEGA ARRAYS
	  !
      USE MODEL
      !
	  IMPLICIT NONE
      ! 
      ! GLOBAL IN-OUT
      !
      !N/A
	  !
	  ! LOCAL I/O FLAGS ETC
	  !
	  CHARACTER(100):: CTEMP
	  INTEGER:: IERROR, LARG1
      !
      !MISC COUNTERS ETC
	  !
      
	  INTEGER:: I, J, K
      INTEGER:: CNT_NODESET, NODE_SET_POSITION, NUM_READ
      !
      INTEGER:: I_TEMP, I2_TEMP(2), I3_TEMP(3), I4_TEMP(4), I5_TEMP(5)
      INTEGER:: I6_TEMP(6), I7_TEMP(7), I8_TEMP(8), I9_TEMP(9), I10_TEMP(8), K_TEMP(1000)
      !
      DOUBLE PRECISION:: D_TEMP, D2_TEMP(2), D3_TEMP(3), D4_TEMP(4)
      DOUBLE PRECISION:: D5_TEMP(5), D6_TEMP(6), D7_TEMP(7), D8_TEMP(8)
	  
	  
	  OPEN(10,FILE='Control.dat',STATUS='OLD',ACTION='READ',iostat=IERROR)

	  IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR OPENING THE CONTROL FILE "Control.dat"',2)
	  
	  
	  !
	  ! GET THE NAME OF THE MODEL FROM THE CONTROL FILE
	  !
	  DO 
	  
	    READ(10,'(A100)',iostat=IERROR) CTEMP
	  
	    IF (IERROR.EQ.(-1)) EXIT
	    CTEMP=TRIM(CTEMP)
        
	    IF (CTEMP.eq.'*MODEL') THEN !*MODEL FILE CARD

	      READ (10,*) 
	      READ (10,'(A50)',iostat=IERROR) MODEL_FILE
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *MODEL IN THE CONTROL FILE',2)
	      CTEMP=TRIM(MODEL_FILE)
	    END IF
	  END DO
	  !
	  ! CHECK FOR ERRORS IN OPENING THE FILE
	  !
	  OPEN(10,FILE=MODEL_FILE,STATUS='OLD',ACTION='READ',iostat=IERROR)
	  !
	  IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR OPENING THE MODEL FILE',2)
      !
      CLOSE(10)
      !
      ! COUNT THE NUMBER OF NODES
      !
      CALL COUNT_NODES(MODEL_FILE, MODEL_NUMP)
      !
      ! COUNT THE NUMBER OF ELEMENTS
      !
      CALL COUNT_ELEMENTS(MODEL_FILE, MODEL_NUMEL, MODEL_NUMBLOCK, MODEL_NUMEL_BLOCKS)
      !
      ! COUNT THE NUMBER OF NODE SETS
      !
      CALL COUNT_NODE_SETS(MODEL_FILE, NUM_NODESET,MAX_NINODESET)
      !
      ! COUNT THE NUMBER OF PRESCRIBED STEPS
      !
      CALL COUNT_FIXITY2_STEPS(MODEL_FILE, FIXITY2_STEPS)
      !
      ! ALLOCATE NODAL AND ELEMENT ARRAYS
      !
      ALLOCATE(MODEL_COO(3,MODEL_NUMP))
      ALLOCATE(MODEL_ELCON(8,MODEL_NUMEL))
      ALLOCATE(MODEL_ELBID(MODEL_NUMEL))
      ALLOCATE(MODEL_SM_LEN(6,MODEL_NUMP))
      ALLOCATE(MODEL_SM_AREA(3,MODEL_NUMP))
      ALLOCATE(MODEL_SM_VOL(MODEL_NUMP))
      ALLOCATE(MODEL_WIN(3,MODEL_NUMP))
      ALLOCATE(MODEL_VOL(MODEL_NUMP))
      ALLOCATE(MODEL_NSNI_FAC(3,MODEL_NUMP))
      ALLOCATE(MODEL_NODE_SET_LIST(NUM_NODESET,MAX_NINODESET))
      ALLOCATE(MODEL_NODE_SET_ID(NUM_NODESET))
      ALLOCATE(MODEL_NODE_SET_LENGTH(NUM_NODESET))
      ALLOCATE(MODEL_EBC(3,MODEL_NUMP))
      ALLOCATE(MODEL_NONZERO_EBC(3,MODEL_NUMP))
      ALLOCATE(MODEL_VINIT(3,MODEL_NUMP))
      ALLOCATE(MODEL_MAT_TYPE(MODEL_NUMP))
      ALLOCATE(MODEL_PROP(30,MODEL_NUMP))
      ALLOCATE(MODEL_BODY_ID(MODEL_NUMP))
      ALLOCATE(MODEL_NORM_WIN(MODEL_NUMP))
      ALLOCATE(MODEL_MASS(MODEL_NUMP*3))
	  ALLOCATE(MODEL_X_MOM(MODEL_NUMP),MODEL_Y_MOM(MODEL_NUMP),MODEL_Z_MOM(MODEL_NUMP))
	  
	  !FOR OUTPUT... TRANSFER LOCAL TO GLOBAL MODEL-SIZE ARRAYS... PROBABLY NOT THE BEST WAY TO DO THINGS!
      ALLOCATE(MODEL_ACL(3,MODEL_NUMP))
      ALLOCATE(MODEL_VEL(3,MODEL_NUMP))
      ALLOCATE(MODEL_DSP(3,MODEL_NUMP))
      ALLOCATE(MODEL_DSP_TOT(3,MODEL_NUMP))
      ALLOCATE(MODEL_COO_CURRENT(3,MODEL_NUMP))
      ALLOCATE(MODEL_FINT(3*MODEL_NUMP))
      ALLOCATE(MODEL_NODE_IDS(MODEL_NUMP))
      !0702
      ALLOCATE(MODEL_BODYFORCE(3,MODEL_NUMP))

      ALLOCATE(FIXITY2_TIME(FIXITY2_STEPS,2)) !GC
      ALLOCATE(FIXITY2_NONZERO_EBC(FIXITY2_STEPS,3)) 

      
      DO I=1,MODEL_NUMP
      MODEL_NODE_IDS(I) = I
      END DO
	  
      !
      ! READ IN THE COORDINATES AND ELEMENTS
      !
      CNT_NODESET = 0
      MODEL_NODE_SET_LIST = 0
      MODEL_TYPEL_BLOCKS=0
      K_TEMP=0
      !
	  OPEN(10,FILE=MODEL_FILE,STATUS='OLD',ACTION='READ',iostat=IERROR)
      
	  DO 
	  
100     READ(10,'(A100)',iostat=IERROR) CTEMP
	  
	    IF (IERROR.EQ.(-1)) EXIT
        
	    CTEMP=TRIM(CTEMP)
        
	    IF (CTEMP.eq.'*NODE') THEN
	  
	      READ (10,*) 
	      DO I=1,MODEL_NUMP
             READ(10,'(I8,3F16.0,2I8)',IOSTAT=IERROR) I_TEMP, MODEL_COO(:,I_TEMP), I2_TEMP
	      END DO
          
        ELSEIF (CTEMP.eq.'*ELEMENT_SOLID') THEN
	  
          K = 0
        
	        DO I=1,MODEL_NUMBLOCK
            
              !SKIP A LINE IF NEW BLOCK
              IF (I.GT.1)  READ(10,*)
            
	          DO J=1,MODEL_NUMEL_BLOCKS(I)
              
                K = K + 1
                READ(10,'(10I8)',IOSTAT=IERROR) I_TEMP, I_TEMP, MODEL_ELCON(:,K)
              MODEL_ELBID(K)=I
              END DO
              K_TEMP(MODEL_ELBID(K))=K
            END DO
	  
	    ELSEIF (CTEMP.eq.'*SET_NODE_LIST') THEN
	    
          CNT_NODESET= CNT_NODESET + 1
          NODE_SET_POSITION = 0
          READ(10,*) MODEL_NODE_SET_ID(CNT_NODESET) !SET ID

	        DO 
200           READ(10,'(A100)',iostat=IERROR) CTEMP
	          IF (CTEMP.eq.'*SET_NODE_LIST') THEN
                 NODE_SET_POSITION = 0
                 CNT_NODESET= CNT_NODESET + 1
                 READ(10,*) MODEL_NODE_SET_ID(CNT_NODESET) !SET ID
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
                
                MODEL_NODE_SET_LIST(CNT_NODESET,NODE_SET_POSITION+1:NODE_SET_POSITION+NUM_READ) = I8_TEMP(1:NUM_READ)
                
                MODEL_NODE_SET_LENGTH(CNT_NODESET) = NODE_SET_POSITION+NUM_READ
                
                !MOVE THE "HEAD"
                NODE_SET_POSITION = NODE_SET_POSITION + NUM_READ
                
                
              ELSE
                GOTO 100
              END IF
            
            END DO
            
	    END IF 
	  
	  END DO !READ FILE WITH ARBITRARY DO
      
      CLOSE(10)
      
      !DETERMINE THE ELEMENT TYPE (HEX OR TET) OF EACH ELEMENT
      	        DO I=1,MODEL_NUMBLOCK
                IF (MODEL_ELCON(7,K_TEMP(I)).EQ.MODEL_ELCON(8,K_TEMP(I))) THEN
                MODEL_TYPEL_BLOCKS(I)=1
                ELSE
                MODEL_TYPEL_BLOCKS(I)=2
                END IF
                END DO
      !
      ! NOW PREPAIRE THE OTHER DATA
      !
      CONTINUE
      !
      ! GET THE DILATIONS (WINDOWS)
      !
      CALL GET_ADIAL(MODEL_WIN,MODEL_NUMP,MODEL_NUMEL,MODEL_ELCON,MODEL_COO,MODEL_TYPEL_BLOCKS,MODEL_ELBID)
      !
      ! GET THE VOLUMES
      !
      CALL GET_VOLUMES(MODEL_VOL,MODEL_NUMP,MODEL_NUMEL,MODEL_ELCON,MODEL_COO,MODEL_TYPEL_BLOCKS,MODEL_ELBID)
      !
      ! GET THE SMOOTHING LENGTHS
      

      !
      CALL GET_SMOOTHING_LENGTHS(MODEL_VOL,MODEL_NUMP,MODEL_COO,MODEL_WIN, &
                                MODEL_SM_LEN,MODEL_SM_AREA,MODEL_SM_VOL,MODEL_XDIST_MAX, &
                                MODEL_YDIST_MAX, MODEL_ZDIST_MAX, &
                                MODEL_X_MOM,MODEL_Y_MOM,MODEL_Z_MOM) 
      !
      ! GET THE APPROXIMATE 2ND MOMENTS OF INERTIA (/VOL)
      !
      CALL SECOND_MOMENT(MODEL_VOL,MODEL_NUMP,MODEL_COO,MODEL_WIN, &   !IN
                                MODEL_NSNI_FAC) !OUT
      !
      ! ECHO THE MODEL TO THE LOG FILE
      !
      CALL ECHO_MODEL(MODEL_WIN,MODEL_VOL,MODEL_NUMP,MODEL_NUMEL,MODEL_ELCON,MODEL_COO, &
                            MODEL_ELBID,MODEL_FILE, MODEL_NUMBLOCK, MODEL_NUMEL_BLOCKS, &
                             MODEL_SM_LEN,MODEL_SM_AREA,MODEL_SM_VOL,MODEL_NSNI_FAC)
      
	  RETURN
	  
	  END SUBROUTINE
	  
      
      
      
      
      
      
      
      

      
