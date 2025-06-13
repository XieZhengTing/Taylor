
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
    
      MODULE MODEL
      !
      ! PURPOSE OF THIS LOGIC:
      !
      ! this is a module file that contains declarations for the "model" which is
      ! all the data which is shared globally and seen by all processors. This 
      ! module should only be included in very top-level routines such as MAIN
      !
      IMPLICIT NONE
    

      
	  !DATA
      CHARACTER(50):: MODEL_FILE
	  DOUBLE PRECISION, ALLOCATABLE:: MODEL_VOL(:)
	  DOUBLE PRECISION, ALLOCATABLE:: MODEL_WIN(:,:)
	  DOUBLE PRECISION, ALLOCATABLE:: MODEL_COO(:,:)
      INTEGER, ALLOCATABLE:: MODEL_ELCON(:,:)
      INTEGER, ALLOCATABLE:: MODEL_ELBID(:)
      INTEGER:: MODEL_NUMP
      INTEGER:: MODEL_NUMEL
      INTEGER:: MODEL_NUMBLOCK
      INTEGER:: MODEL_NUMEL_BLOCKS(1000)
      INTEGER:: MODEL_TYPEL_BLOCKS(1000)
      INTEGER:: NUM_NODESET
      INTEGER:: MAX_NINODESET
      INTEGER, ALLOCATABLE:: MODEL_NODE_SET_LIST(:,:)
      INTEGER, ALLOCATABLE:: MODEL_NODE_SET_ID(:)
      INTEGER, ALLOCATABLE:: MODEL_NODE_SET_LENGTH(:)
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_SM_LEN(:,:)
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_SM_AREA(:,:)
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_SM_VOL(:)
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_NSNI_FAC(:,:)
	  CHARACTER(100):: MODEL_SET_NAMES(1000)
      !CONTROL OF MODEL
      INTEGER, ALLOCATABLE:: MODEL_EBC(:,:)
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_NONZERO_EBC(:,:)
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_VINIT(:,:)
      INTEGER, ALLOCATABLE:: MODEL_MAT_TYPE(:)
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_PROP(:,:)
      INTEGER, ALLOCATABLE:: MODEL_BODY_ID(:)
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_NORM_WIN(:) 
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_X_MOM(:),MODEL_Y_MOM(:),MODEL_Z_MOM(:)       
      DOUBLE PRECISION, SAVE:: MODEL_XDIST_MAX, MODEL_YDIST_MAX,MODEL_ZDIST_MAX
      DOUBLE PRECISION, SAVE:: TIME, TIME_END, DLT
      INTEGER,SAVE::FIXITY2_STEPS 
      DOUBLE PRECISION,ALLOCATABLE::FIXITY2_TIME(:,:),FIXITY2_NONZERO_EBC(:,:)
      
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_MASS(:)
      
	  
	  !FOR OUTPUT... TRANSFER LOCAL TO GLOBAL MODEL-SIZE ARRAYS... PROBABLY NOT THE BEST WAY TO DO THINGS!
	  
	  
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_ACL(:,:)
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_VEL(:,:)
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_DSP(:,:)
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_DSP_TOT(:,:)
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_COO_CURRENT(:,:)
      DOUBLE PRECISION, ALLOCATABLE:: MODEL_FINT(:)
      INTEGER, ALLOCATABLE:: MODEL_NODE_IDS(:)
      
      !0329
     ! DOUBLE PRECISION, ALLOCATABLE:: MODEL_BODYFORCE(:)
      !0702 KC
	  DOUBLE PRECISION, ALLOCATABLE:: MODEL_BODYFORCE(:,:)
      !Output exodus format
      CHARACTER(len=50), SAVE::exodusVars   
      INTEGER, SAVE::exodusVarsCount,exodusVarsCountMax,exodusStep
      DIMENSION exodusVars(:) 
      ALLOCATABLE exodusVars  
      CHARACTER(LEN=50)exodusFileName      
      
      END MODULE MODEL
      
