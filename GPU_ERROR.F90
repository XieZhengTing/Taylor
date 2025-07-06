MODULE GPU_ERROR
    IMPLICIT NONE
    
    ! GPU error flag and information
    INTEGER :: GPU_ERROR_FLAG = 0
    INTEGER :: GPU_ERROR_TYPE = 0
    CHARACTER(200) :: GPU_ERROR_REASON = ''
    
    !$ACC DECLARE CREATE(GPU_ERROR_FLAG, GPU_ERROR_TYPE, GPU_ERROR_REASON)
    
CONTAINS
    
    ! GPU-safe error setting routine
    SUBROUTINE SET_GPU_ERROR(REASON, ITYPE)
        !$ACC ROUTINE SEQ
        IMPLICIT NONE
        CHARACTER(*), INTENT(IN) :: REASON
        INTEGER, INTENT(IN) :: ITYPE
        
        ! Set error flag and type
        GPU_ERROR_FLAG = 1
        GPU_ERROR_TYPE = ITYPE
        
        ! Cannot copy string on GPU easily, just set a code
        ! We'll map error codes to messages on CPU side
    END SUBROUTINE SET_GPU_ERROR
    
    ! CPU-side error checking routine
    SUBROUTINE CHECK_GPU_ERROR()
        IMPLICIT NONE
        INTEGER :: LOCAL_FLAG, LOCAL_TYPE
        
        ! Copy error flag from GPU to CPU
        !$ACC UPDATE HOST(GPU_ERROR_FLAG, GPU_ERROR_TYPE)
        
        LOCAL_FLAG = GPU_ERROR_FLAG
        LOCAL_TYPE = GPU_ERROR_TYPE
        
        IF (LOCAL_FLAG .NE. 0) THEN
            ! Map error codes to messages
            SELECT CASE(LOCAL_FLAG)
                CASE(1)
                    CALL EXIT_PROGRAM('PLASTICITY MATERIAL MODEL DID NOT CONVERGE', LOCAL_TYPE)
                CASE(2)
                    CALL EXIT_PROGRAM('INVALID MATERIAL TYPE IN SUBROUTINE CONSTITUTION', LOCAL_TYPE)
                CASE(3)
                    CALL EXIT_PROGRAM('NOT A HYPERELASTIC MATERIAL', LOCAL_TYPE)
                CASE DEFAULT
                    CALL EXIT_PROGRAM('UNKNOWN GPU ERROR', LOCAL_TYPE)
            END SELECT
        END IF
        
        ! Reset error flag
        GPU_ERROR_FLAG = 0
        !$ACC UPDATE DEVICE(GPU_ERROR_FLAG)
        
    END SUBROUTINE CHECK_GPU_ERROR
    
END MODULE GPU_ERROR