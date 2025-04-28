
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
    
    

	  SUBROUTINE ECHO_MODEL(MODEL_WIN,MODEL_VOL,MODEL_NUMP,MODEL_NUMEL,MODEL_ELCON,MODEL_COO, &
                            MODEL_ELBID,MODEL_FILE, MODEL_NUMBLOCK, MODEL_NUMEL_BLOCKS, &
                             MODEL_SM_LEN,MODEL_SM_AREA,MODEL_SM_VOL,MODEL_NSNI_FAC)
      
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! ECHO THE MODEL
	  !
      IMPLICIT NONE
      
	  DOUBLE PRECISION, INTENT(IN):: MODEL_VOL(*)
	  DOUBLE PRECISION, INTENT(IN):: MODEL_WIN(3,*)
	  DOUBLE PRECISION, INTENT(IN):: MODEL_COO(3,*)
      INTEGER, INTENT(IN):: MODEL_ELCON(8,*)
      INTEGER, INTENT(IN):: MODEL_ELBID(*)
      INTEGER, INTENT(IN):: MODEL_NUMP, MODEL_NUMEL
      INTEGER, INTENT(IN):: MODEL_NUMBLOCK
      INTEGER, INTENT(IN):: MODEL_NUMEL_BLOCKS(1000)
      
      
	  DOUBLE PRECISION, INTENT(IN):: MODEL_NSNI_FAC(3,*)
	  DOUBLE PRECISION, INTENT(IN):: MODEL_SM_LEN(6,*)
	  DOUBLE PRECISION, INTENT(IN):: MODEL_SM_AREA(3,*)
	  DOUBLE PRECISION, INTENT(IN):: MODEL_SM_VOL(*)
      
      CHARACTER(50), INTENT(IN):: MODEL_FILE
      
      !LOCAL
      INTEGER:: I, J, K
      
      OPEN(50,FILE='MEGA.log',STATUS='OLD',POSITION='APPEND')
      
      WRITE(50,*)
      WRITE(50,*) 'CONTROL FILE:'
      WRITE(50,'(A)') MODEL_FILE 
      WRITE(50,*)
      
      WRITE(50,*)
      WRITE(50,*) 'NUMBER OF BLOCKS:'
      WRITE(50,'(I10)') MODEL_NUMBLOCK 
      WRITE(50,*)
      
      WRITE(50,*)
      WRITE(50,*) ' ************************ BLOCKS ************************ '
      WRITE(50,*)
      WRITE(50,'(A10,A30)') 'BLOCK ID', 'NUMBER OF ELEMENTS'
      DO I=1,MODEL_NUMBLOCK
      WRITE(50,'(I10,I30)') I, MODEL_NUMEL_BLOCKS(I)
      END DO
      WRITE(50,*)
      
      
      WRITE(50,*)
      WRITE(50,*) '************************ NODES ************************'
      WRITE(50,*)
      WRITE(50,'(A9,20A15)') 'NODE ID','X-COORDINATE','Y-COORDINATE','Z-COORDINATE', 'VOLUME','X-WINDOW','Y-WINDOW','Z-WINDOW', &
      'SMTH VOL','+X SMTH DIST','-X SMTH DIST','+Y SMTH DIST','-Y SMTH DIST','+Z SMTH DIST','-Z SMTH DIST', &
      'X SMTH AREA','Y SMTH AREA','Z SMTH AREA'!,'X NSNI FAC','Y NSNI FAC','Z NSNI FAC'
      
      DO I=1,MODEL_NUMP
        WRITE(50,'(I9,20E15.7)') I, MODEL_COO(:,I), MODEL_VOL(I), MODEL_WIN(:,I), &
         MODEL_SM_VOL(I), MODEL_SM_LEN(:,I), MODEL_SM_AREA(:,I)!, MODEL_NSNI_FAC(:,I)
      END DO
      WRITE(50,*)
      
      WRITE(50,*)
      WRITE(50,*) '************************ ELEMENTS ************************'
      WRITE(50,*)
      
      WRITE(50,'(10A10)') 'ELEM-ID','BLOCK-ID','NODE-1','NODE-2','NODE-3','NODE-4', &
                                                 'NODE-5','NODE-6','NODE-7','NODE-8'
      
          K = 0
        
	        DO I=1,MODEL_NUMBLOCK
            
	          DO J=1,MODEL_NUMEL_BLOCKS(I)
              
                K = K + 1
                WRITE(50,'(10I10)') K, MODEL_ELBID(K), MODEL_ELCON(:,K)
              
              END DO
            END DO
	  
      
      
      
      
      CLOSE(50)
      
      
      RETURN
      
      END SUBROUTINE
      