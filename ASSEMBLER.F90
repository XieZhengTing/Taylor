
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
	  
	  
	  SUBROUTINE ASSEMBLER(NUMP,FINT,HPC_SCHEME)
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! COMPUTE THE INTERNAL FORCE FOR ALL LOCALLY OWNED NODES
	  !
	  
	  IMPLICIT NONE
	  
	  INTEGER, INTENT(IN):: NUMP
	  DOUBLE PRECISION, INTENT(INOUT):: FINT(NUMP*3)
	  INTEGER, INTENT(IN):: HPC_SCHEME
	  
	  SELECT CASE (HPC_SCHEME)
	  
	  CASE(1)
	      !
	      ! SERIAL
	      !
	      !
	      ! DO NOTHING
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
      
      
	  
      subroutine OUTPUT_ASSEMBLER(HPC_SCHEME, LOCAL_NUMP, MODEL_NUMP, TOTAL_MODEL_MAP, &
			                         LOCAL_ACL,LOCAL_VEL,LOCAL_DSP,LOCAL_DSP_TOT,LOCAL_COO_CURRENT,LOCAL_FINT, &
			                         MODEL_ACL,MODEL_VEL,MODEL_DSP,MODEL_DSP_TOT,MODEL_COO_CURRENT,MODEL_FINT)

	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! COMPUTE THE INTERNAL FORCE FOR ALL LOCALLY OWNED NODES
	  !
	  
	  IMPLICIT NONE
	  
	  INTEGER, INTENT(IN):: HPC_SCHEME, LOCAL_NUMP, MODEL_NUMP, TOTAL_MODEL_MAP(MODEL_NUMP)
	  DOUBLE PRECISION, INTENT(IN):: LOCAL_ACL(3,LOCAL_NUMP),LOCAL_VEL(3,LOCAL_NUMP), &
	                                 LOCAL_DSP(3,LOCAL_NUMP),LOCAL_DSP_TOT(3,LOCAL_NUMP), &
	                                 LOCAL_COO_CURRENT(3,LOCAL_NUMP),LOCAL_FINT(3*LOCAL_NUMP)
	  DOUBLE PRECISION, INTENT(OUT)::MODEL_ACL(3,LOCAL_NUMP),MODEL_VEL(3,LOCAL_NUMP), &
	                                 MODEL_DSP(3,LOCAL_NUMP),MODEL_DSP_TOT(3,LOCAL_NUMP), &
	                                 MODEL_COO_CURRENT(3,LOCAL_NUMP),MODEL_FINT(3*LOCAL_NUMP)

	  SELECT CASE (HPC_SCHEME)
	  
	  CASE(1)
	      !
	      ! SERIAL
	      !
	      !ASSIGN VALUES WITH NO MAP
		  MODEL_ACL=LOCAL_ACL
		  MODEL_VEL=LOCAL_VEL
		  MODEL_DSP=LOCAL_DSP
		  MODEL_DSP_TOT=LOCAL_DSP_TOT
		  MODEL_COO_CURRENT=LOCAL_COO_CURRENT
		  MODEL_FINT=LOCAL_FINT
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
	  