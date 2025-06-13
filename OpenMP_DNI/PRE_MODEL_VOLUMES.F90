
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

    SUBROUTINE GET_VOLUMES(VOLUME,NUMPTS,NUMELS,IDELM,XYZ,MODEL_TYPEL_BLOCKS,MODEL_ELBID)
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! RETRUN THE VOLUMES ASSOCATED WITH NODES, GIVEN CONNECTIVITY
	  !
    IMPLICIT NONE
    !
    !**GLOBAL**
    !
    !NUMBER OF POINTS
    INTEGER:: NUMPTS
    !
    !VOLUMES
    DOUBLE PRECISION:: VOLUME(NUMPTS)
    !
    !NUMBER OF ELEMENTS
    INTEGER:: NUMELS
    !
    !ID CONNECTIVITY ARRAY
    INTEGER:: IDELM(8,NUMELS)
    !
    !COORDINATE ARRAY
    DOUBLE PRECISION:: XYZ(3,NUMPTS)
    !
    !MODEL TYPE =1 TET, =2 HEX
    INTEGER:: MODEL_TYPEL_BLOCKS(1000)
    !
    INTEGER:: MODEL_ELBID(NUMELS)
    !**LOCAL**
    !
    INTEGER:: I, J, K
    !
    !HEX SPLITTING TET ID MAP
    INTEGER::T(12,4)
    !
    !LOCAL IDS OF NODES TO AN ELEMENT
    INTEGER:: I8_TEMP(8)
    !
    !TEMPORARY ELEMENT COORDINATE INFORMATION
    DOUBLE PRECISION:: XYZEL(3,9)
    !
    !TEMPORARY VOLUME
    DOUBLE PRECISION:: DVOL
    
    
    
    !INITIALIZE THE VOLUMES
    VOLUME = 0.0d0
    
    !MAP OF HEX NODES INTO 12 TETS
    CALL GET_TET_IDS(T)
    
    !LOOP OVER ELEMENTS, COMPUTE VOLUMES
    DO I=1,NUMELS
    
       !ELEMENT IDS
       I8_TEMP=IDELM(:,I)
       
       DO J=1,3
         DO K=1,8
           XYZEL(J,K) = XYZ(J,I8_TEMP(K))
         END DO
       END DO
       IF (MODEL_TYPEL_BLOCKS(MODEL_ELBID(I)).EQ.2) THEN
       XYZEL(:,9)=0.0d0
       DO J=1,3
         DO K=1,8
           XYZEL(J,9) = XYZEL(J,9) + XYZEL(J,K)
         END DO
       END DO
       
       !NODE 9 IS THE CENTER
       XYZEL(:,9)=XYZEL(:,9)*0.125d0
       
       !GET THE TOTAL VOLUME FOR THIS HEX
       CALL GET_VOL2(DVOL,XYZEL)
       !ASSIGN 1/8 OF THE VOLUME TO THE NODES
       DO K=1,8
         VOLUME(I8_TEMP(K)) = VOLUME(I8_TEMP(K)) + 0.125d0*DVOL
       END DO
       ELSE
              XYZEL(:,9)=0.0d0
       DO J=1,3
         DO K=1,4
           XYZEL(J,9) = XYZEL(J,9) + XYZEL(J,K)
         END DO
       END DO
       
       !NODE 9 IS THE CENTER
       XYZEL(:,9)=XYZEL(:,9)*0.25d0
       
       !GET THE TOTAL VOLUME FOR THIS HEX
       CALL GET_VOL1(DVOL,XYZEL)
       
       !ASSIGN 1/8 OF THE VOLUME TO THE NODES
       DO K=1,4
         VOLUME(I8_TEMP(K)) = VOLUME(I8_TEMP(K)) + 0.25d0*DVOL
       END DO
       END IF
       
    END DO
    
    RETURN
    
    END SUBROUTINE
    
    