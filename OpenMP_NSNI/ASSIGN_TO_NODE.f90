
    
    
            SUBROUTINE ASSIGN_INT_VALUE_TO_NODE( EXISTING_VAL,TRIAL_VAL,EXIT_STRING)
            
            INTEGER:: EXISTING_VAL,TRIAL_VAL
	        CHARACTER(*):: EXIT_STRING
            
            
            IF (EXISTING_VAL.NE.(-1)) THEN 
              !THIS VALUE HAS BEEN ASIGNED!
              CALL EXIT_PROGRAM(EXIT_STRING,2)
            ELSE
              !THIS VALUE HAS NOT BEEN ASIGNED YET
              EXISTING_VAL = TRIAL_VAL
            END IF
            
            RETURN
            
            END SUBROUTINE
            
            
            SUBROUTINE ASSIGN_DBL_VALUE_TO_NODE( EXISTING_VAL,TRIAL_VAL,EXIT_STRING)
            
            DOUBLE PRECISION:: EXISTING_VAL,TRIAL_VAL
	        CHARACTER(*):: EXIT_STRING
            
            
            IF (DABS(EXISTING_VAL).GT.(1.0D-16)) THEN
              !THIS VALUE HAS BEEN ASIGNED!
              CALL EXIT_PROGRAM(EXIT_STRING,2)
            ELSE
              !THIS VALUE HAS NOT BEEN ASIGNED YET
              EXISTING_VAL = TRIAL_VAL
            END IF
            
            END SUBROUTINE
            