
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

    SUBROUTINE GET_ADIAL(ADIAL,NUMPTS,NUMELS,IDELM,XYZ,MODEL_TYPEL_BLOCKS,MODEL_ELBID)
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! GENERATE THE WINDOW FUNCTIONS (NORMALIZED TO ONE)
	  !
    IMPLICIT NONE
    !
    !**GLOBAL**
    !
    !NUMBER OF POINTS
    INTEGER:: NUMPTS
    !
    !DILATIONS
    DOUBLE PRECISION:: ADIAL(3,NUMPTS)
    !
    !NUMBER OF ELEMENTS
    INTEGER:: NUMELS
    !
    !
    INTEGER::MODEL_TYPEL_BLOCKS(1000),MODEL_ELBID(NUMELS)
    !ID CONNECTIVITY ARRAY
    INTEGER:: IDELM(8,NUMELS)
    !
    !COORDINATE ARRAY
    DOUBLE PRECISION:: XYZ(3,NUMPTS)
    !
    !**LOCAL**
    !
    INTEGER:: I, J, K
    !
    !TEMPORARY VARIABLES FOR COMPUTING DILATIONS
    DOUBLE PRECISION:: DIST_TEMP
    DOUBLE PRECISION:: XDIST_MAX, YDIST_MAX, ZDIST_MAX 
    !
    !LOCAL IDS OF NODES TO AN ELEMENT
    INTEGER:: I8_TEMP(8)
    !
    !TEMPORARY ELEMENT COORDINATE INFORMATION
    DOUBLE PRECISION:: XYZEL(3,8)
    
    
    !INITIALIZE VALUES
    ADIAL = 0.0d0
    
  
    
    !LOOP OVER ELEMENTS AND CHECK THE MAX X, Y, Z DIMENSIONS
    DO I=1,NUMELS
    
       !ELEMENT IDS
       I8_TEMP=IDELM(:,I)
       
       !ELEMENT COORDINATES
       DO J=1,3
         DO K=1,8
           XYZEL(J,K) = XYZ(J,I8_TEMP(K))
         END DO
       END DO
       
       !FIND THE MAX DIMS OF THIS ELEMENT
       XDIST_MAX = 0.0d0
       YDIST_MAX = 0.0d0
       ZDIST_MAX = 0.0d0
       IF (MODEL_TYPEL_BLOCKS(MODEL_ELBID(I)).EQ.2) THEN
       DO K=1,8
         DO J=1,8
           IF (J.NE.K) THEN
           
             DIST_TEMP = DABS(XYZEL(1,K) - XYZEL(1,J))
             IF (DIST_TEMP.GT.XDIST_MAX) XDIST_MAX = DIST_TEMP
             
             DIST_TEMP = DABS(XYZEL(2,K) - XYZEL(2,J))
             IF (DIST_TEMP.GT.YDIST_MAX) YDIST_MAX = DIST_TEMP
             
             DIST_TEMP = DABS(XYZEL(3,K) - XYZEL(3,J))
             IF (DIST_TEMP.GT.ZDIST_MAX) ZDIST_MAX = DIST_TEMP
             
           END IF
         END DO
       END DO
       
       !CHECK THE CURRENT DILATIONS FOR THE ATTACHED NODES,
       !AND IF IT'S SMALLER THAN THE MAX DIMS OF THE ELEMENT,
       !ASSIGN IT TO THE NODE
       DO K=1,8
       
          IF (XDIST_MAX.GT.ADIAL(1,I8_TEMP(K))) THEN
            ADIAL(1,I8_TEMP(K)) = XDIST_MAX
          END IF
       
          IF (YDIST_MAX.GT.ADIAL(2,I8_TEMP(K))) THEN
            ADIAL(2,I8_TEMP(K)) = YDIST_MAX
          END IF
       
          IF (ZDIST_MAX.GT.ADIAL(3,I8_TEMP(K))) THEN
            ADIAL(3,I8_TEMP(K)) = ZDIST_MAX
          END IF
       END DO
       ELSE
              DO K=1,4
         DO J=1,4
           IF (J.NE.K) THEN
           
             DIST_TEMP = DABS(XYZEL(1,K) - XYZEL(1,J))
             IF (DIST_TEMP.GT.XDIST_MAX) XDIST_MAX = DIST_TEMP
             
             DIST_TEMP = DABS(XYZEL(2,K) - XYZEL(2,J))
             IF (DIST_TEMP.GT.YDIST_MAX) YDIST_MAX = DIST_TEMP
             
             DIST_TEMP = DABS(XYZEL(3,K) - XYZEL(3,J))
             IF (DIST_TEMP.GT.ZDIST_MAX) ZDIST_MAX = DIST_TEMP
             
           END IF
         END DO
       END DO
       
       !CHECK THE CURRENT DILATIONS FOR THE ATTACHED NODES,
       !AND IF IT'S SMALLER THAN THE MAX DIMS OF THE ELEMENT,
       !ASSIGN IT TO THE NODE
       DO K=1,4
       
          IF (XDIST_MAX.GT.ADIAL(1,I8_TEMP(K))) THEN
            ADIAL(1,I8_TEMP(K)) = XDIST_MAX
          END IF
       
          IF (YDIST_MAX.GT.ADIAL(2,I8_TEMP(K))) THEN
            ADIAL(2,I8_TEMP(K)) = YDIST_MAX
          END IF
       
          IF (ZDIST_MAX.GT.ADIAL(3,I8_TEMP(K))) THEN
            ADIAL(3,I8_TEMP(K)) = ZDIST_MAX
          END IF
       END DO
       END IF
      
       
     END DO
    
    
    RETURN
    
    END SUBROUTINE
    