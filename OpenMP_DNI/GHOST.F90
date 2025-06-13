
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
    
    
    
    
    
    SUBROUTINE GHOSTER(HPC_SCHEME)
                                
      IMPLICIT NONE
                                
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! ASSIGN INITIAL GHOST NODES
	  !
    INTEGER, INTENT(IN):: HPC_SCHEME
    
	  SELECT CASE (HPC_SCHEME)
	  
	  CASE(1)
      
	      !
	      ! SERIAL: THERE ARE NO GHOSTS
          !
      
	  CASE DEFAULT
	      !
	      !SOMETHING WENT WRONG
	      !
	      CALL EXIT_PROGRAM('INVALID HPC_SCHEME TYPE IN SUBROUTINE ASSEMBLER',0)
	      !
      END SELECT
		
      RETURN
		
      END SUBROUTINE
      
      
    SUBROUTINE GHOST_INIT(HPC_SCHEME, LOCAL_NUMP, GHOST_NUMP, LOCAL_GHOST)
                                
      IMPLICIT NONE
                                
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! ASSIGN INITIAL GHOST NODES
	  !
    INTEGER, INTENT(IN):: HPC_SCHEME
    INTEGER, INTENT(IN):: LOCAL_NUMP
    INTEGER, INTENT(OUT):: GHOST_NUMP
    INTEGER, INTENT(OUT):: LOCAL_GHOST(LOCAL_NUMP)
    !
    ! LOCAL
    !
    INTEGER:: I
    
	  SELECT CASE (HPC_SCHEME)
	  
	  CASE(1)
	      !
	      ! SERIAL
	      !
      GHOST_NUMP = 0
          !
          ! THERE ARE NO GHOSTS
          !
      DO I=1,LOCAL_NUMP
      
         LOCAL_GHOST(I) = 0
         
      END DO
      
      
	  CASE DEFAULT
	      !
	      !SOMETHING WENT WRONG
	      !
	      CALL EXIT_PROGRAM('INVALID HPC_SCHEME TYPE IN SUBROUTINE ASSEMBLER',0)
	      !
      END SELECT
		
      RETURN
		
      END SUBROUTINE
      
      
    SUBROUTINE GHOST_INIT_MAP(HPC_SCHEME, GHOST_NUMP, LOCAL_NUMP, TOTAL_MODEL_MAP)
                                
      IMPLICIT NONE
                                
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! ASSIGN INITIAL GHOST NODES
	  !
    INTEGER, INTENT(IN):: HPC_SCHEME
    INTEGER, INTENT(IN):: LOCAL_NUMP
    INTEGER, INTENT(IN):: GHOST_NUMP
    INTEGER, INTENT(OUT):: TOTAL_MODEL_MAP(LOCAL_NUMP+GHOST_NUMP)
    !
    ! LOCAL
    !
    INTEGER:: I, K
    
    
    
    
	  SELECT CASE (HPC_SCHEME)
	  
	  CASE(1)
	      !
	      ! SERIAL
	      !
          K = LOCAL_NUMP+GHOST_NUMP
          !
          ! THERE ARE NO GHOSTS
          !
      DO I=1,K
      
         TOTAL_MODEL_MAP(I) = I
         
      END DO
      
      
	  CASE DEFAULT
	      !
	      !SOMETHING WENT WRONG
	      !
	      CALL EXIT_PROGRAM('INVALID HPC_SCHEME TYPE IN SUBROUTINE ASSEMBLER',0)
	      !
      END SELECT
		
      RETURN
		
      END SUBROUTINE
      