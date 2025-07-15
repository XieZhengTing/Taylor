
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

    PROGRAM MEGA
    !
    USE MODEL
    USE CONTROL
    USE GPU_ERROR
    !
    IMPLICIT NONE
    !
    INTEGER:: I, J, M, ITEMP, I_STEP, STEP_NUM, NUM_EBC


    !
    !
    INTEGER:: GHOST_NUMP, TOTAL_LOCAL_NUMP, TOTAL_LOCAL_SIZE
    !
    ! ARRAYS THE SIZE OF LOCAL NODES
    !
    INTEGER:: LOCAL_NUMP                
    INTEGER , ALLOCATABLE:: MODEL_MAP(:)                    !LOCAL-GLOBAL MODEL MAP
    DOUBLE PRECISION , ALLOCATABLE:: LOCAL_SM_LEN(:,:)      !NODAL SMOOTHING LENGTHS
    DOUBLE PRECISION , ALLOCATABLE:: LOCAL_SM_AREA(:,:)     !NODAL SMOOTHING AREAS
    DOUBLE PRECISION , ALLOCATABLE:: LOCAL_SM_VOL(:)        !NODAL SMOOTHING VOLUMES
    DOUBLE PRECISION , ALLOCATABLE:: LOCAL_WIN(:,:)         !NODAL WINDOWS
    DOUBLE PRECISION , ALLOCATABLE:: LOCAL_VOL(:)           !NODAL VOLUMES
    INTEGER , ALLOCATABLE:: LOCAL_EBC(:,:)                  !NODAL BOUNDARY CONDITION FLAGS
    DOUBLE PRECISION , ALLOCATABLE:: LOCAL_NONZERO_EBC(:,:) !NODAL NON-ZERO BOUNDARY CONDITIONS
    DOUBLE PRECISION , ALLOCATABLE:: LOCAL_VINIT(:,:)       !NODAL INITIAL VELOCITIES
    INTEGER , ALLOCATABLE:: LOCAL_MAT_TYPE(:)               !NODAL MATERIAL TYPE
    DOUBLE PRECISION , ALLOCATABLE:: LOCAL_PROP(:,:)        !NODAL MATERIAL PROPERTIES
    INTEGER , ALLOCATABLE:: LOCAL_BODY_ID(:)                !NODAL BODY IDS
    INTEGER, ALLOCATABLE:: TOTAL_MODEL_MAP(:)               !LOCAL-GLOBAL MODEL MAP
    DOUBLE PRECISION , ALLOCATABLE:: LOCAL_CHAR_DIST(:)     !NODAL CHARACTERISTIC DISTANCE
    DOUBLE PRECISION , ALLOCATABLE:: LOCAL_WAVE_VEL(:)      !NODAL WAVE VELOCITY ESTIMATE
    DOUBLE PRECISION:: LOCAL_XDIST_MAX                      !MAX SMOOTHING DISTANCE IN X DIRECTION
    DOUBLE PRECISION:: LOCAL_YDIST_MAX                      !MAX SMOOTHING DISTANCE IN Y DIRECTION
    DOUBLE PRECISION:: LOCAL_ZDIST_MAX                      !MAX SMOOTHING DISTANCE IN Z DIRECTION
    INTEGER , ALLOCATABLE:: LOCAL_IJKSPC(:,:)               !IJK SPACE OF NODES IN A BIN (OUTPUT)
    !
    DOUBLE PRECISION , ALLOCATABLE:: LOCAL_NSNI_FAC(:,:)    !
    !
    INTEGER, ALLOCATABLE:: LOCAL_GHOST(:)
    !
    ! ARRAYS THE SIZE OF GHOST + LOCAL NODES
    !
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_STATE(:,:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_STRESS(:,:)


    !NSNI
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_DX_STRESS(:,:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_DY_STRESS(:,:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_DZ_STRESS(:,:)

    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_X_MOM(:), LOCAL_Y_MOM(:), LOCAL_Z_MOM(:) !NODAL MOMENTS FOR NSNI


    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_STRAIN(:,:), LOCAL_STRAIN_EQ(:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_MASS(:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_FEXT(:),LOCAL_FEXT_NMO(:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_FINT(:),LOCAL_FINT_NMO(:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_ACL(:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_VEL(:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_DSP(:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_DSP_TOT_PHY(:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_DSP_TOT(:)
 DOUBLE PRECISION, ALLOCATABLE:: GDINC(:)
 DOUBLE PRECISION, ALLOCATABLE:: GVEL(:) 
 DOUBLE PRECISION, ALLOCATABLE:: GACL(:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_COO(:,:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_COO_CURRENT(:,:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_PRFORCE(:,:)

    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_H_STRESS(:,:)!GC
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_S_STRESS(:,:)!GC
    
    DOUBLE PRECISION, ALLOCATABLE:: TOTAL_FORCE(:,:) !KC

    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_ACL_PHY(:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_VEL_PHY(:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_DINC_PHY(:)

    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_DX_STRAIN(:,:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_DY_STRAIN(:,:)
    DOUBLE PRECISION, ALLOCATABLE:: LOCAL_DZ_STRAIN(:,:)

    LOGICAL, ALLOCATABLE:: LOCAL_EBC_NODES(:)
	
	DOUBLE PRECISION:: MPDC !MASS PROPORTIAL DAMPING MATRIX (DIAG)

    LOGICAL:: DO_INTERP


    !
    DOUBLE PRECISION:: TIME_COUNTER,DLOCAL_INT_ENERGY,DLOCAL_KIN_ENERGY,DLOCAL_TOTAL_ENERGY, &
        MAXV, MAXD, DLOCAL_EXT_ENERGY, LOCAL_INT_WORK
    !
    LOGICAL:: LINIT, LFORCE,LINIT_TIME
    !
    INTEGER:: STEPS, TIMER_STEPS
    !
    REAL:: REAL_TIME_0, REAL_TIME_1, REAL_TIME_2, REAL_TIME_TEMP
    !
    DOUBLE PRECISION:: SIM_TIME_1, SIM_TIME_2, SIM_TIME_TEMP
    DOUBLE PRECISION:: SPEED, SIM_TIME_LEFT, REAL_TIME_REMAINING
    LOGICAL::TIMER_FLAG
    CHARACTER*13:: CTIME_ALL

    DOUBLE PRECISION:: TOTAL_REAL_TIME
    CHARACTER*40:: CTOTAL_REAL_TIME

    DOUBLE PRECISION:: LOCAL_DLT

    LOGICAL:: LINIT_TEMP
    !
    !
    !*********************************************************
    !******************** EXECUTABLE CODE ********************
    !*********************************************************
    !
    !
    CALL PRINT_OPENING
    !
    CALL LOG_HEADER
    !
    CALL PRE_OMP
    !
    CALL PRE_MODEL
    !
    CALL LOG_APPEND_SPACE('MODEL FILE READ SUCSESSFULLY')
    !
    CALL PRE_CONTROL
    !
    CALL LOG_APPEND_SPACE('CONTROL FILE READ SUCSESSFULLY')
 ! Debug and sync LFEM_OUTPUT for OpenACC
 PRINT *, 'MAIN: LFEM_OUTPUT after PRE_CONTROL =', LFEM_OUTPUT
 !$ACC UPDATE DEVICE(LFEM_OUTPUT)
   
   ! Size check
   PRINT *, '=== SIZE CHECK ==='
   PRINT *, 'MODEL_NUMP:', MODEL_NUMP
    !
    CALL ASSIGN_PARALLEL(HPC_SCHEME, MODEL_NUMP, LOCAL_NUMP)
    !
    ALLOCATE(MODEL_MAP(LOCAL_NUMP))
    !
    CALL ASSIGN_PARALLEL_MAP(HPC_SCHEME, MODEL_NUMP, LOCAL_NUMP, MODEL_MAP)
    !
    ! ALLOCATE LOCAL MODEL ARRAYS (IF LOCALLY OWNED NODES CHANGE, THEY WILL HAVE TO BE REALLOCATED)
    !
    ALLOCATE(LOCAL_SM_LEN(6,LOCAL_NUMP), LOCAL_SM_AREA(3,LOCAL_NUMP), LOCAL_SM_VOL(LOCAL_NUMP))
    ALLOCATE(LOCAL_WIN(3,LOCAL_NUMP),    LOCAL_VOL(LOCAL_NUMP),       LOCAL_NSNI_FAC(3,LOCAL_NUMP))
    ALLOCATE(LOCAL_VINIT(3,LOCAL_NUMP),  LOCAL_MAT_TYPE(LOCAL_NUMP),  LOCAL_PROP(30,LOCAL_NUMP))
    ALLOCATE(LOCAL_BODY_ID(LOCAL_NUMP),  LOCAL_CHAR_DIST(LOCAL_NUMP), LOCAL_WAVE_VEL(LOCAL_NUMP))
    ALLOCATE(LOCAL_X_MOM(LOCAL_NUMP),    LOCAL_Y_MOM(LOCAL_NUMP),     LOCAL_Z_MOM(LOCAL_NUMP))
	ALLOCATE(LOCAL_IJKSPC(3,LOCAL_NUMP))
    
    ALLOCATE(TOTAL_FORCE(3,2))!KC FOR 2 BODIES
    !
    LOCAL_CHAR_DIST = 0.0d0
    LOCAL_WAVE_VEL = 0.0d0
    !
    ! ASSIGN LOCAL MODEL ARRAYS
    !
    CALL PARALLEL_MODEL(    MODEL_NUMP,    LOCAL_NUMP,     MODEL_MAP,       &
        MODEL_SM_LEN,  MODEL_SM_AREA,  MODEL_SM_VOL,    &
        MODEL_WIN,     MODEL_VOL,      MODEL_NSNI_FAC,  &
        MODEL_VINIT,   MODEL_MAT_TYPE,                  &
        MODEL_PROP,    MODEL_BODY_ID,  LOCAL_SM_LEN,    &
        LOCAL_SM_AREA, LOCAL_SM_VOL,   LOCAL_WIN,       &
        LOCAL_VOL,     LOCAL_NSNI_FAC,                  &
        LOCAL_VINIT,   LOCAL_MAT_TYPE, LOCAL_PROP,      &
        MODEL_XDIST_MAX, MODEL_YDIST_MAX, MODEL_ZDIST_MAX, &
        LOCAL_XDIST_MAX, LOCAL_YDIST_MAX, LOCAL_ZDIST_MAX, &
        MODEL_X_MOM, MODEL_Y_MOM, MODEL_Z_MOM , &
        LOCAL_X_MOM, LOCAL_Y_MOM, LOCAL_Z_MOM , &
        LOCAL_BODY_ID)
    !
    ! DEALLOCATE MODEL ARRAYS (IF LOCALLY OWNED NODES CHANGE, DO NOT DO THIS)
    !
    DEALLOCATE(MODEL_SM_LEN,MODEL_SM_AREA,MODEL_MAP,MODEL_SM_VOL,MODEL_WIN,MODEL_VOL,MODEL_NSNI_FAC,MODEL_MAT_TYPE, &
    MODEL_PROP,MODEL_BODY_ID,MODEL_X_MOM,MODEL_Y_MOM, MODEL_Z_MOM)

    !
    ! GET INITIAL GHOST NODES
    !
    ALLOCATE(LOCAL_GHOST(LOCAL_NUMP))
    !
    CALL GHOST_INIT(HPC_SCHEME, LOCAL_NUMP, GHOST_NUMP, LOCAL_GHOST)
    !
    TOTAL_LOCAL_NUMP = LOCAL_NUMP + GHOST_NUMP
    !
    ! GET MAP FOR THE GHOST INFORMATION (IT'S JUST FOR ASSIGNING INITIAL VELOCITY AND COORDINATES...)
    !
    ALLOCATE(TOTAL_MODEL_MAP(LOCAL_NUMP + GHOST_NUMP))
    !
    CALL GHOST_INIT_MAP(HPC_SCHEME, GHOST_NUMP, LOCAL_NUMP, TOTAL_MODEL_MAP)
    !
    ! GENERATE STATE AND FIELD VARIABLE ARRAYS (FOR GHOSTS, SOME QUANTITIES WILL BE ZERO)
    !

    TOTAL_LOCAL_SIZE = LOCAL_NUMP + GHOST_NUMP + GHOST_BUFFER

   PRINT *, 'TOTAL_LOCAL_SIZE:', TOTAL_LOCAL_SIZE
   PRINT *, 'TOTAL_LOCAL_NUMP:', TOTAL_LOCAL_NUMP
   PRINT *, 'LOCAL_NUMP:', LOCAL_NUMP
   IF (TOTAL_LOCAL_SIZE .NE. MODEL_NUMP) THEN
       PRINT *, 'WARNING: Size mismatch!'
       PRINT *, '  TOTAL_LOCAL_SIZE =', TOTAL_LOCAL_SIZE
       PRINT *, '  MODEL_NUMP =', MODEL_NUMP
   END IF

    ALLOCATE(LOCAL_STATE(20,TOTAL_LOCAL_SIZE),    LOCAL_STRESS(6,TOTAL_LOCAL_SIZE))
    ALLOCATE(LOCAL_DX_STRESS(6,TOTAL_LOCAL_SIZE), LOCAL_DY_STRESS(6,TOTAL_LOCAL_SIZE),   LOCAL_DZ_STRESS(6,TOTAL_LOCAL_SIZE))
    ALLOCATE(LOCAL_DX_STRAIN(6,TOTAL_LOCAL_SIZE), LOCAL_DY_STRAIN(6,TOTAL_LOCAL_SIZE),   LOCAL_DZ_STRAIN(6,TOTAL_LOCAL_SIZE))
    ALLOCATE(LOCAL_STRAIN(6,TOTAL_LOCAL_SIZE),    LOCAL_COO(3,TOTAL_LOCAL_SIZE),         LOCAL_MASS(3*TOTAL_LOCAL_SIZE), LOCAL_STRAIN_EQ(TOTAL_LOCAL_SIZE))
    ALLOCATE(LOCAL_FEXT(3*TOTAL_LOCAL_SIZE),      LOCAL_FEXT_NMO(3*TOTAL_LOCAL_SIZE),    LOCAL_FINT(3*TOTAL_LOCAL_SIZE),LOCAL_FINT_NMO(3*TOTAL_LOCAL_SIZE))
    ALLOCATE(LOCAL_ACL(3*TOTAL_LOCAL_SIZE),       LOCAL_VEL(3*TOTAL_LOCAL_SIZE),         LOCAL_DSP(3*TOTAL_LOCAL_SIZE))
    ALLOCATE(LOCAL_PRFORCE(3,TOTAL_LOCAL_SIZE),   LOCAL_ACL_PHY(3*TOTAL_LOCAL_SIZE),     LOCAL_VEL_PHY(3*TOTAL_LOCAL_SIZE),LOCAL_DINC_PHY(3*TOTAL_LOCAL_SIZE))
    ALLOCATE(LOCAL_DSP_TOT(3*TOTAL_LOCAL_SIZE),   LOCAL_DSP_TOT_PHY(3*TOTAL_LOCAL_SIZE), LOCAL_COO_CURRENT(3,TOTAL_LOCAL_SIZE))
    ALLOCATE(LOCAL_EBC(3,TOTAL_LOCAL_SIZE),       LOCAL_EBC_NODES(TOTAL_LOCAL_SIZE),     LOCAL_NONZERO_EBC(3,TOTAL_LOCAL_SIZE))
    ALLOCATE(LOCAL_H_STRESS(6,TOTAL_LOCAL_SIZE),  LOCAL_S_STRESS(6,TOTAL_LOCAL_SIZE))
 ALLOCATE(GDINC(3*TOTAL_LOCAL_SIZE))
 ALLOCATE(GVEL(3*TOTAL_LOCAL_SIZE))
 ALLOCATE(GACL(3*TOTAL_LOCAL_SIZE))
 
 ! 初始化
 GDINC = 0.0d0
 GVEL = 0.0d0
 GACL = 0.0d0   
    !
    ! GET INITIAL VALUES
    !
    CALL STATE_FEILD_INIT(TOTAL_LOCAL_SIZE,TOTAL_LOCAL_NUMP, LOCAL_NUMP, MODEL_NUMP, MODEL_VINIT, TOTAL_MODEL_MAP, MODEL_COO,MODEL_MASS,MODEL_EBC, &
        MODEL_NONZERO_EBC,LOCAL_STATE,LOCAL_STRESS,LOCAL_STRAIN,LOCAL_H_STRESS,LOCAL_S_STRESS,LOCAL_COO,LOCAL_FEXT,LOCAL_FINT, &
        LOCAL_ACL,LOCAL_VEL,LOCAL_DSP,LOCAL_DSP_TOT,LOCAL_DSP_TOT_PHY,LOCAL_MASS,LOCAL_COO_CURRENT, LOCAL_EBC, LOCAL_NONZERO_EBC,LOCAL_EBC_NODES, &
        LOCAL_PRFORCE, LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, LOCAL_DX_STRAIN, LOCAL_DY_STRAIN, LOCAL_DZ_STRAIN,LOCAL_STRAIN_EQ)

   ! Check LOCAL_EBC after initialization
   PRINT *, '=== LOCAL_EBC CHECK after STATE_INIT ==='

   ! OpenACC特定修正：修正節點集重疊導致的分類錯誤
   #ifdef _OPENACC
   PRINT *, '=== Applying OpenACC node classification fix ==='
   DO I = 1, LOCAL_NUMP
       ! BAR節點 (body_id=1) 不應有約束，應有初速 -373 m/s
       IF (LOCAL_BODY_ID(I) .EQ. 1) THEN
           IF (ANY(LOCAL_EBC(:,I) .NE. 0)) THEN
               PRINT '(A,I5,A)', 'Fixing BAR node ', I, ' (removing constraints)'
           END IF
           ! 清除所有約束
           LOCAL_EBC(:,I) = 0
           LOCAL_EBC_NODES(I) = .FALSE.
           ! 設定正確的初速
           LOCAL_VEL((I-1)*3+1) = 0.0d0      ! X velocity
           LOCAL_VEL((I-1)*3+2) = 0.0d0      ! Y velocity
           LOCAL_VEL((I-1)*3+3) = -373.0d0   ! Z velocity
       END IF
   END DO
   
   ! 同步修正到GPU
   !$ACC UPDATE DEVICE(LOCAL_VEL, LOCAL_EBC, LOCAL_EBC_NODES)
   
   ! 驗證節點5
   IF (LOCAL_NUMP .GE. 5) THEN
       PRINT *, 'Node 5 after OpenACC fix:'
       PRINT '(A,3F10.2)', '  Velocity:', LOCAL_VEL(13:15)
       PRINT '(A,3I3)', '  EBC:', LOCAL_EBC(:,5)
   END IF
   #endif

   NUM_EBC = 0
   DO I = 1, LOCAL_NUMP
       DO J = 1, 3
           IF (LOCAL_EBC(J,I) .NE. 0) THEN
               NUM_EBC = NUM_EBC + 1
               IF (NUM_EBC .LE. 5) THEN
                   PRINT '(A,I4,A,I1,A,I2)', 'Node', I, ' Dir', J, ' EBC=', LOCAL_EBC(J,I)
               END IF
           END IF
       END DO
   END DO
   PRINT *, 'Total EBC entries:', NUM_EBC
   ! Special check for nodes 1-5
   DO I = 1, MIN(5, LOCAL_NUMP)
       IF (ANY(LOCAL_EBC(:,I) .NE. 0)) THEN
           PRINT '(A,I4,A,3I2)', 'Node', I, ' EBC:', LOCAL_EBC(:,I)
       END IF
   END DO
    !
    DEALLOCATE(MODEL_VINIT,MODEL_MASS,MODEL_EBC,MODEL_COO,MODEL_NONZERO_EBC)

    !
    TIME = 0.0d0
    STEP_NUM = 0

    ! OpenACC: Begin GPU data region for main time integration loop
    !
!$ACC DATA COPYIN(                                          &!← 把 DLT 及所有初始化陣列搬到 GPU
!$ACC&     DLT,                                             &
!$ACC&     LOCAL_COO, LOCAL_COO_CURRENT, LOCAL_MASS,       &
!$ACC&     LOCAL_NONZERO_EBC, LOCAL_EBC_NODES)             &
!$ACC& COPY(LOCAL_EBC,                                  & ! 改為 COPY 以便雙向同步
!$ACC&     LOCAL_SM_LEN, LOCAL_SM_AREA, LOCAL_SM_VOL,       &
!$ACC&     LOCAL_WIN, LOCAL_VOL, LOCAL_NSNI_FAC,           &
!$ACC&     LOCAL_MAT_TYPE, LOCAL_PROP, LOCAL_BODY_ID,       &
!$ACC&     LOCAL_CHAR_DIST, LOCAL_WAVE_VEL,                &
!$ACC&     MODEL_BODYFORCE,                                &
!$ACC&     LINIT,                                         &
!$ACC&     LOCAL_X_MOM, LOCAL_Y_MOM, LOCAL_Z_MOM,          &
!$ACC&     LFINITE_STRAIN, LLAGRANGIAN,                    &
!$ACC&     MODEL_ELCON, MODEL_NUMEL)                       &
!$ACC& COPY(                                               &!← 在離開 region 時自動拷回以下更新結果
!$ACC&     LOCAL_STATE, LOCAL_STRESS, LOCAL_STRAIN, LOCAL_STRAIN_EQ, &
!$ACC&     GDINC, GVEL, GACL,                              &
!$ACC&     LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS,     &
!$ACC&     LOCAL_H_STRESS, LOCAL_S_STRESS,                         &
!$ACC&     LOCAL_FINT, LOCAL_FEXT, LOCAL_FINT_NMO, LOCAL_FEXT_NMO, &
!$ACC&     LOCAL_ACL, LOCAL_DSP,                                   &
!$ACC&     LOCAL_ACL_PHY, LOCAL_VEL_PHY, LOCAL_DINC_PHY,           &
!$ACC&     LOCAL_DSP_TOT, LOCAL_DSP_TOT_PHY,                       &
!$ACC&     LOCAL_PRFORCE, TOTAL_FORCE)                      &
!$ACC& COPYIN(LOCAL_VEL)  ! Ensure initial velocity is copied to GPU

   ! Check LOCAL_EBC inside DATA region
   PRINT *, '=== LOCAL_EBC CHECK inside DATA region ==='
   DO I = 1, MIN(5, LOCAL_NUMP)
       PRINT '(A,I4,A,3I2)', 'Node', I, ' EBC:', LOCAL_EBC(:,I)
   END DO
   ! Check specific nodes that should have constraints
   DO I = LOCAL_NUMP-5, LOCAL_NUMP
       IF (I .GT. 0) THEN
           PRINT '(A,I4,A,3I2)', 'Node', I, ' EBC:', LOCAL_EBC(:,I)
       END IF
   END DO

   ! Debug: Check initial velocities before time integration
   PRINT *, '=== INITIAL VELOCITY CHECK ==='
   PRINT *, 'First 5 nodes Z-velocity (should be -373 for Taylor bar):'
   DO I = 1, MIN(5, LOCAL_NUMP)
       PRINT '(A,I4,A,3E15.5)', 'Node', I, ' Vel:', &
           LOCAL_VEL((I-1)*3+1), LOCAL_VEL((I-1)*3+2), LOCAL_VEL((I-1)*3+3)
   END DO
   PRINT *, 'Max absolute velocity:', MAXVAL(ABS(LOCAL_VEL))
   
   ! CRITICAL: Ensure initial velocities are on GPU
   !$ACC UPDATE DEVICE(LOCAL_VEL)

    !
    !******************************************OUTPUT******************************************
    !
    ! MAKE HEADER FOR EXODUS FILE
    !
    !(COMMENT THE FOLLOWING IF EXODUS NOT INSTALLED, DO NOT INCLUDE OUTPUT.F90) #NOEXODUS
    !call output_Header(LOCAL_NUMP,LOCAL_COO) !TEMP
    !
    CALL CLEAN_VTKS
    !
    !CALL OUTPUT_ASSEMBLER(HPC_SCHEME, LOCAL_NUMP, MODEL_NUMP, TOTAL_MODEL_MAP, &
    !	                         LOCAL_ACL,LOCAL_VEL,LOCAL_DSP,LOCAL_DSP_TOT,LOCAL_COO_CURRENT,LOCAL_FINT, &
    !	                         MODEL_ACL,MODEL_VEL,MODEL_DSP,MODEL_DSP_TOT,MODEL_COO_CURRENT,MODEL_FINT)
    !
    OPEN (120,FILE='TIME.OUT',STATUS='REPLACE')
    
    OPEN (121,FILE='ENERGY_INFO.OUT',STATUS='REPLACE')
    
    OPEN (122,FILE='TOTAL_FORCE.OUT',STATUS='REPLACE')
    
    WRITE(121,'(5A20)') 'STEP','INTERNAL ENERGY','KINETIC ENERGY','EXTERNAL ENERGY','TOTAL ENERGY'
    
    WRITE(122,'(10A15)') 'Time','Body_ID','Fx','Fy','Fz','Body_ID','Fx','Fy','Fz'
    exodusStep=0
    !
    IF (HPC_SCHEME.EQ.1) THEN
        !
        IF (UNF_OUTPUT) THEN
            call UNF_OUTPUT_STEP_VTK(exodusStep,LOCAL_NUMP,MODEL_NUMEL, MODEL_ELCON, LOCAL_COO_CURRENT,MODEL_NODE_IDS,LOCAL_DSP_TOT,LOCAL_VEL, &
                LOCAL_ACL, LOCAL_FINT, LOCAL_EBC,  LOCAL_BODY_ID,  LOCAL_MAT_TYPE,  LOCAL_COO,    &
                LOCAL_VINIT,  MODEL_NORM_WIN,  LOCAL_WIN,  LOCAL_VOL,  LOCAL_MASS,  LOCAL_PROP,    &
                LOCAL_CHAR_DIST, LOCAL_WAVE_VEL, LOCAL_PRFORCE, LOCAL_STRESS, LOCAL_STRAIN,LOCAL_STATE, LOCAL_STRAIN_EQ, &
				LOCAL_IJKSPC)
        ELSE
            call OUTPUT_STEP_VTK(exodusStep,LOCAL_NUMP,MODEL_NUMEL, MODEL_ELCON, LOCAL_COO_CURRENT,MODEL_NODE_IDS,LOCAL_DSP_TOT,LOCAL_VEL, &
                LOCAL_ACL, LOCAL_FINT, LOCAL_EBC,  LOCAL_BODY_ID,  LOCAL_MAT_TYPE,  LOCAL_COO,    &
                LOCAL_VINIT,  MODEL_NORM_WIN,  LOCAL_WIN,  LOCAL_VOL,  LOCAL_MASS,  LOCAL_PROP,    &
                LOCAL_CHAR_DIST, LOCAL_WAVE_VEL, LOCAL_PRFORCE, LOCAL_STRESS, LOCAL_STRAIN,LOCAL_STATE, LOCAL_STRAIN_EQ, &
				LOCAL_IJKSPC)
        END IF

    END IF
    !
    exodusStep = exodusStep + 1
    !
    STEPS = 0
    TIMER_STEPS = 0
    !
    !******************************************OUTPUT******************************************
    !
    TIME_COUNTER = 0.0d0
    LINIT = .TRUE.
    !$ACC UPDATE DEVICE(LINIT)
    LINIT_TIME = .TRUE.
    SIM_TIME_1 = 0.0d0
    CALL CPU_TIME(REAL_TIME_0)
    CALL CPU_TIME(REAL_TIME_1)
    !
    ! START TIME INTEGRATION
    !
    CALL WRITE_OUT('STARTING TIME INTEGRATION')
    !

    !WRITE(50,'(2A8,2A15,4A15)') 'Out Step','Step','Time','Delta T','Int Engery', 'Kin Energy', 'Ext Energy', 'Tot Energy', 'Est Time Rem'
    !WRITE(*,'(2A8,2A15,4A15)')  'Out Step','Step','Time','Delta T','Int Engery', 'Kin Energy', 'Ext Energy', 'Tot Energy', 'Est Time Rem'
    !
    WRITE(50,'(2A8,2A15,1A30)') 'Out Step','Step','Time','Delta T', 'Estimated Time Remaining'
    WRITE(*,'(2A8,2A15,1A30)')  'Out Step','Step','Time','Delta T', 'Estimated Time Remaining'
    
	!
    DLOCAL_INT_ENERGY = 0.0d0
    DLOCAL_KIN_ENERGY = 0.0d0
    DLOCAL_EXT_ENERGY = 0.0d0
    LOCAL_INT_WORK = 0.0d0
    TOTAL_FORCE(:,:)=0.0d0
    !
    DO
        !
        !
        !

       ! Verify LOCAL_EBC at start of time loop
       IF (STEPS .EQ. 0) THEN
           PRINT *, '=== LOCAL_EBC at start of time loop ==='
           ! Sync from GPU to ensure we have latest values
           !$ACC UPDATE HOST(LOCAL_EBC)
           NUM_EBC = 0
           DO I = 1, LOCAL_NUMP
               DO J = 1, 3
                   IF (LOCAL_EBC(J,I) .NE. 0) THEN
                       NUM_EBC = NUM_EBC + 1
                       IF (NUM_EBC .LE. 5) THEN
                           PRINT *, 'Node', I, 'Dir', J, 'EBC=', LOCAL_EBC(J,I)
                       END IF
                   END IF
               END DO
           END DO
           PRINT *, 'Total EBC in time loop:', NUM_EBC
           ! Check node 5 specifically
           PRINT *, 'Node 5 EBC in time loop:', LOCAL_EBC(:,5)
       END IF
       

        ! ASSIGN NEW GHOSTS (RIGHT NOW, THIS SUBROUTINE DOES NOTHING)
        !
        CALL GHOSTER(HPC_SCHEME)
        !
        ! PERFORM NECESSARY GHOSTING TASKS (REALLOCATE, ASSIGN, ETC.)
        !
        CONTINUE
        !
        ! CHECK THE TIME_STEP
        !


        IF(LINIT) THEN
            LOCAL_DSP = 0.0d0 !TO GET ZERO FINT, JUST TO GET SHP
            DO_INTERP = .FALSE.
            CALL HANDELER(       LOCAL_WIN,      LOCAL_VOL,         LOCAL_NUMP,     LOCAL_COO,      LOCAL_COO_CURRENT,     &
                LOCAL_SM_LEN,    LOCAL_SM_AREA,  LOCAL_SM_VOL,      LOCAL_NSNI_FAC, LOCAL_GHOST,       &
                LOCAL_PROP,      LOCAL_STATE,    LOCAL_STRESS,      LOCAL_STRAIN,   LOCAL_H_STRESS,LOCAL_S_STRESS,   LOCAL_DSP, LOCAL_DSP_TOT,        &
                LOCAL_FINT,      LOCAL_MAT_TYPE, TOTAL_LOCAL_SIZE,  LINIT,          LOCAL_DLT, &
                LOCAL_FEXT,      LOCAL_CHAR_DIST, LOCAL_WAVE_VEL, LOCAL_EBC_NODES, &
                DO_INTERP,       LOCAL_DINC_PHY, LOCAL_VEL_PHY, LOCAL_VEL, LOCAL_ACL, LOCAL_ACL_PHY, &
                LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
                LOCAL_X_MOM,     LOCAL_Y_MOM,    LOCAL_Z_MOM, &
                LOCAL_XDIST_MAX, LOCAL_YDIST_MAX, LOCAL_ZDIST_MAX, MODEL_BODYFORCE,DLT, LOCAL_INT_WORK, LOCAL_BODY_ID, LOCAL_STRAIN_EQ, &
				LOCAL_IJKSPC)
        ! Check for GPU errors
        CALL CHECK_GPU_ERROR()
            !LINIT = .FALSE.  !SET IT TO FALSE AFTER THE SHP/DSHP COMPUTATION
        ENDIF
        !

        !DO I=1,LOCAL_NUMP

        !
        !REDUCE VALUES TO ALL PROCS HERE
        !

        !
        !  PREDICTOR
        !
        !CALL PREDICTOR(TOTAL_LOCAL_SIZE,LOCAL_ACL,LOCAL_VEL,LOCAL_DSP,DLT)

     IF (LINIT .AND. AUTO_TS) THEN
         LOCAL_DSP = 0.0d0
         !$ACC UPDATE DEVICE(LOCAL_DSP)           !← 把主機的零值傳到 GPU
     ELSE
        !$ACC PARALLEL LOOP PRESENT(LOCAL_DSP, LOCAL_VEL, LOCAL_ACL, LOCAL_DSP_TOT) ASYNC(1)

         DO I = 1, 3*TOTAL_LOCAL_SIZE
             LOCAL_DSP(I)     = DLT * LOCAL_VEL(I) + 0.5d0*DLT**2 * LOCAL_ACL(I)
             LOCAL_DSP_TOT(I) = LOCAL_DSP_TOT(I) + LOCAL_DSP(I)
             LOCAL_VEL(I)     = LOCAL_VEL(I) + 0.5d0*DLT      * LOCAL_ACL(I)
         END DO
        !$ACC END PARALLEL LOOP
        !$ACC WAIT(1)                             !← 確保非同步計算完成
        !$ACC UPDATE HOST(LOCAL_DSP, LOCAL_DSP_TOT, LOCAL_VEL)  !← 把 GPU 計算結果拷回主機
     ! Debug: Check if predictor is working
        IF (STEPS .LE. 2) THEN
            PRINT *, '=== PREDICTOR CHECK (Step', STEPS, ') ==='
            PRINT *, 'Time step DLT:', DLT
            PRINT *, 'First node after predictor:'
            PRINT *, '  Vel:', LOCAL_VEL(1:3)
            PRINT *, '  Accel:', LOCAL_ACL(1:3)
            PRINT *, '  Disp incr:', LOCAL_DSP(1:3)
            PRINT *, '  Total disp:', LOCAL_DSP_TOT(1:3)
        END IF
     END IF

 ! CRITICAL: Copy displacement increments to GDINC for interpolation
 GDINC = LOCAL_DSP
 GVEL = LOCAL_VEL  
 GACL = LOCAL_ACL
 
 ! Ensure these arrays are synchronized to GPU before interpolation
 !$ACC UPDATE DEVICE(GDINC, GVEL, GACL)
 
 ! Debug: Verify GDINC has correct values
 IF (STEPS .LE. 2) THEN
     PRINT *, '=== GDINC CHECK before interpolation ==='
     PRINT *, 'GDINC(1:3):', GDINC(1:3)
     PRINT *, 'Should match LOCAL_DSP:', LOCAL_DSP(1:3)
 END IF
        !
        !COMPUTE THE PREDICTED CURRENT COORD FOR SEMI-LAG SHP CACULATION
        !

        DO_INTERP = .TRUE.
        LINIT_TEMP = .FALSE.
        CALL HANDELER(      LOCAL_WIN,       LOCAL_VOL,      LOCAL_NUMP,        LOCAL_COO,      LOCAL_COO_CURRENT,     &
            LOCAL_SM_LEN,    LOCAL_SM_AREA,  LOCAL_SM_VOL,      LOCAL_NSNI_FAC, LOCAL_GHOST,       &
            LOCAL_PROP,      LOCAL_STATE,    LOCAL_STRESS,      LOCAL_STRAIN,   LOCAL_H_STRESS, LOCAL_S_STRESS,   LOCAL_DSP,  LOCAL_DSP_TOT,      &
            LOCAL_FINT,      LOCAL_MAT_TYPE, TOTAL_LOCAL_SIZE,  LINIT_TEMP,          LOCAL_DLT, &
            LOCAL_FEXT,      LOCAL_CHAR_DIST, LOCAL_WAVE_VEL, LOCAL_EBC_NODES, &
            DO_INTERP,       LOCAL_DINC_PHY, LOCAL_VEL_PHY, LOCAL_VEL, LOCAL_ACL, LOCAL_ACL_PHY, &
            LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
            LOCAL_X_MOM,     LOCAL_Y_MOM,    LOCAL_Z_MOM, &
            LOCAL_XDIST_MAX, LOCAL_YDIST_MAX, LOCAL_ZDIST_MAX,MODEL_BODYFORCE,DLT, LOCAL_INT_WORK, LOCAL_BODY_ID, LOCAL_STRAIN_EQ, &
				LOCAL_IJKSPC)
        ! Check for GPU errors
        CALL CHECK_GPU_ERROR()


        !$ACC UPDATE HOST(LOCAL_DINC_PHY, LOCAL_VEL_PHY, LOCAL_ACL_PHY)
 ! Debug: Check PHY variables
 IF (STEPS .LE. 2) THEN
     PRINT *, '=== PHY VARIABLES CHECK (Step', STEPS, ') ==='
     PRINT *, 'LOCAL_DINC_PHY (first 3):', LOCAL_DINC_PHY(1:3)
     PRINT *, 'Before update - LOCAL_DSP_TOT_PHY:', LOCAL_DSP_TOT_PHY(1:3)
 END IF
        LOCAL_DSP_TOT_PHY = LOCAL_DSP_TOT_PHY + LOCAL_DINC_PHY

 IF (STEPS .LE. 2) THEN
     PRINT *, 'After update - LOCAL_DSP_TOT_PHY:', LOCAL_DSP_TOT_PHY(1:3)
 END IF

        ! Ensure LOCAL_DSP_TOT_PHY is on GPU for coordinate update
        !$ACC UPDATE DEVICE(LOCAL_DSP_TOT_PHY)
        
        !
        ! LIKELY HAVE TO UPDATE GHOSTS HERE FOR GHOST SCHEMES, THESE ARRAYS
        ! SHOULD BE DIFFERENT SIZES THEN #TODO
        !
        !$ACC PARALLEL LOOP PRESENT(LOCAL_COO_CURRENT, LOCAL_COO, LOCAL_DSP_TOT_PHY)
        DO I=1,LOCAL_NUMP
            DO J=1,3
                LOCAL_COO_CURRENT(J,I) = LOCAL_COO(J,I) + LOCAL_DSP_TOT_PHY((I-1)*3+J)
            END DO
        END DO
        !$ACC END PARALLEL LOOP

        ! Sync updated coordinates back to host
        !$ACC UPDATE HOST(LOCAL_COO_CURRENT)

 ! Debug: Verify coordinate update
 IF (MOD(STEPS, 16) .EQ. 0) THEN
     PRINT *, '=== COORDINATE UPDATE CHECK (Step', STEPS, ') ==='
     PRINT *, 'Sample displacement:', LOCAL_DSP_TOT_PHY(1:3)
     PRINT *, 'Original coord:', LOCAL_COO(:,1)
     PRINT *, 'Current coord:', LOCAL_COO_CURRENT(:,1)
 END IF

 ! Update all arrays needed for VTK output
 !$ACC UPDATE HOST(LOCAL_DSP_TOT_PHY)
 !$ACC UPDATE HOST(LOCAL_STRESS, LOCAL_STRAIN, LOCAL_STATE, LOCAL_STRAIN_EQ)
 !$ACC UPDATE HOST(LOCAL_H_STRESS, LOCAL_S_STRESS)
 ! CRITICAL: Update MODEL_ELCON for element connectivity output
 !$ACC UPDATE HOST(MODEL_ELCON)

 ! Debug: Check coordinate values before VTK output
 PRINT *, '=== COORDINATE CHECK BEFORE VTK OUTPUT ==='
 PRINT *, 'First 3 nodes coordinates:'
 DO I = 1, 3
     PRINT '(A,I4,A,3E15.5)', 'Node', I, ':', LOCAL_COO_CURRENT(:,I)
 END DO
 PRINT *, 'Min/Max X:', MINVAL(LOCAL_COO_CURRENT(1,:)), MAXVAL(LOCAL_COO_CURRENT(1,:))
 PRINT *, 'Min/Max Y:', MINVAL(LOCAL_COO_CURRENT(2,:)), MAXVAL(LOCAL_COO_CURRENT(2,:))
 PRINT *, 'Min/Max Z:', MINVAL(LOCAL_COO_CURRENT(3,:)), MAXVAL(LOCAL_COO_CURRENT(3,:))

        !TEST HUGHS-WINDET ROTATION ALGORITHM, LATER SHOULD BE REMOVED
        !CALL ROTATION_TEST(LOCAL_DSP,LOCAL_COO,LOCAL_NUMP,TIME,DLT)
        !
        !
        !WRITE(*,*) SQRT((LOCAL_COO_CURRENT(1,1)-LOCAL_COO_CURRENT(1,305))**2+(LOCAL_COO_CURRENT(2,1)-LOCAL_COO_CURRENT(2,305))**2)
        !
        ! GET THE INTERNAL FORCE AND NEIGHBORS ETC
        !
        DO_INTERP = .FALSE.
        LINIT_TEMP = .FALSE.
        IF(LINIT) THEN
            LOCAL_FINT_NMO = 0.0d0
            LOCAL_FEXT_NMO = 0.0d0
        ELSE
            LOCAL_FINT_NMO = LOCAL_FINT
            LOCAL_FEXT_NMO = LOCAL_FEXT !NOT USED! NEED TO FILL THIS OUT IN HANDLER
        END IF

        CALL HANDELER(      LOCAL_WIN,       LOCAL_VOL,      LOCAL_NUMP,        LOCAL_COO,      LOCAL_COO_CURRENT,     &
            LOCAL_SM_LEN,    LOCAL_SM_AREA,  LOCAL_SM_VOL,      LOCAL_NSNI_FAC, LOCAL_GHOST,       &
            LOCAL_PROP,      LOCAL_STATE,    LOCAL_STRESS,      LOCAL_STRAIN,   LOCAL_H_STRESS,LOCAL_S_STRESS,   LOCAL_DSP, LOCAL_DSP_TOT,        &
            LOCAL_FINT,      LOCAL_MAT_TYPE, TOTAL_LOCAL_SIZE,  LINIT_TEMP,     LOCAL_DLT, &
            LOCAL_FEXT,      LOCAL_CHAR_DIST, LOCAL_WAVE_VEL, LOCAL_EBC_NODES, &
            DO_INTERP,       LOCAL_DINC_PHY, LOCAL_VEL_PHY, LOCAL_VEL, LOCAL_ACL, LOCAL_ACL_PHY, &
            LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
            LOCAL_X_MOM,     LOCAL_Y_MOM,    LOCAL_Z_MOM, &
            LOCAL_XDIST_MAX, LOCAL_YDIST_MAX, LOCAL_ZDIST_MAX,MODEL_BODYFORCE,DLT, LOCAL_INT_WORK, LOCAL_BODY_ID, LOCAL_STRAIN_EQ, &
				LOCAL_IJKSPC)
        ! Check for GPU errors
        CALL CHECK_GPU_ERROR()

        !
        ! GET THE INTERNAL FORCE SUFFICIENT FOR CONTINUEING ONTO TIME INTEGRATION
        !(INTERNAL FORCE FOR LOCALLY OWNED NODES)
        !
        ! NOTE THE BELOW WILL HAVE TO BE MODIFIED FOR HPC, I DIDNT PUT ENOUGH HOOKS IN THE CODE...
        ! #TODO
        !
        CALL ASSEMBLER(LOCAL_NUMP,LOCAL_FINT,HPC_SCHEME)
        !LOCAL_FINT = 0.D0 !TEMP FOR TESTING ROTATION
        IF (AUTO_TS) DLT = LOCAL_DLT
        
        IF (PDSTIME.NE.0.0D0) PDSEARCH=CEILING(PDSTIME/DLT) 

       ! Final check before boundary enforcement
       IF (STEPS .LE. 2) THEN
           PRINT *, '=== LOCAL_EBC before boundary enforcement ==='
           !$ACC UPDATE HOST(LOCAL_EBC)  ! Ensure latest values from GPU
           DO I = 1, MIN(5, LOCAL_NUMP)
               IF (ANY(LOCAL_EBC(:,I) .NE. 0)) THEN
                   PRINT *, 'Node', I, 'EBC:', LOCAL_EBC(:,I)
               END IF
           END DO
       END IF

        DO I=1,LOCAL_NUMP

            ! Note: This loop contains reductions (TOTAL_FORCE)
            ! Will be handled in Step 3 with proper reduction clauses

            DO J=1,3

                M = (I-1)*3+J

                IF (LOCAL_EBC(J,I).EQ.1) THEN

                    LOCAL_VEL(M) = 0.0d0
                    LOCAL_ACL(M) = 0.0d0

                    LOCAL_PRFORCE(J,I) = - LOCAL_FINT(M)

                    LOCAL_FEXT(M) = - LOCAL_FINT(M)

                ELSEIF (LOCAL_EBC(J,I).EQ.2) THEN !NON-ZERO ESSENTIAL BC

                    ! FIND WHICH TIME_DURATION THE CURRENT TIME FALL INTO
                    DO I_STEP = 1, FIXITY2_STEPS
                        IF (FIXITY2_TIME(I_STEP,1) .LT. TIME .AND. TIME .LT. FIXITY2_TIME(I_STEP,2)) THEN
                            STEP_NUM = I_STEP
                            EXIT
                        END IF
                    END DO

                    ! ASSIGN VELOCITY
                    IF (1 .LE. STEP_NUM .AND. STEP_NUM .LE. FIXITY2_STEPS) THEN
                        LOCAL_VEL(M) = FIXITY2_NONZERO_EBC(STEP_NUM,J)/(FIXITY2_TIME(STEP_NUM,2) - FIXITY2_TIME(STEP_NUM,1))
                    ELSE
                        LOCAL_VEL(M) = 0.0d0
                    END IF

                    STEP_NUM = 0 ! REINITIALIZE

                    LOCAL_ACL(M) = 0.0d0

                    LOCAL_PRFORCE(J,I) = - LOCAL_FINT(M)

                    LOCAL_FEXT(M) = - LOCAL_FINT(M)

                ELSE   !FREE

				MPDC = LOCAL_PROP(21,I)*LOCAL_MASS(M) !MASS PROPORTIAL DAMPING MATRIX (DIAG)
				
                    LOCAL_ACL(M) = (LOCAL_FEXT(M) - LOCAL_FINT(M) - MPDC*LOCAL_VEL(M))/ &
					(LOCAL_MASS(M)+0.5d0*DLT*MPDC)
					
                    LOCAL_VEL(M) = LOCAL_VEL(M) + 0.5d0*DLT*LOCAL_ACL(M)

                    LOCAL_PRFORCE(J,I) = 0.0d0

                END IF

               ! Debug: Check forces for first few nodes
               IF ((I .EQ. 5 .OR. I .LE. 3) .AND. STEPS .LE. 2) THEN
                   IF (J .EQ. 1) THEN
                       PRINT '(A,I3,A)', '=== Node ', I, ' forces ==='
                       PRINT *, '  EBC:', LOCAL_EBC(:,I)
                       PRINT *, '  FINT:', LOCAL_FINT((I-1)*3+1:I*3)
                       PRINT *, '  FEXT:', LOCAL_FEXT((I-1)*3+1:I*3)
                       PRINT *, '  Mass:', LOCAL_MASS((I-1)*3+1:I*3)
                       PRINT *, '  Accel:', LOCAL_ACL((I-1)*3+1:I*3)
                   END IF
               END IF

            TOTAL_FORCE(J,LOCAL_BODY_ID(I))=TOTAL_FORCE(J,LOCAL_BODY_ID(I))+LOCAL_PRFORCE(J,I)

            END DO
        END DO
        ! Synchronize CPU updates to GPU for next iteration
        !$ACC UPDATE DEVICE(LOCAL_VEL, LOCAL_ACL)
        WRITE(122,'(E15.5,I8,3(E15.5),I8,3(E15.5))') TIME,LOCAL_BODY_ID(1), TOTAL_FORCE(1,LOCAL_BODY_ID(1)), TOTAL_FORCE(2,LOCAL_BODY_ID(1)),TOTAL_FORCE(3,LOCAL_BODY_ID(1)), LOCAL_BODY_ID(LOCAL_NUMP), TOTAL_FORCE(1,LOCAL_BODY_ID(LOCAL_NUMP)),TOTAL_FORCE(2,LOCAL_BODY_ID(LOCAL_NUMP)),TOTAL_FORCE(3,LOCAL_BODY_ID(LOCAL_NUMP))


        !
        ! GET THE ACCELERATION FROM EQUATION OF MOTION
        !
        !CALL EOM(TOTAL_LOCAL_NUMP,LOCAL_FINT,LOCAL_FEXT,LOCAL_MASS,LOCAL_ACL)
        !
        ! ENFORCE LOCAL_EBCS
        ! LIKELY HAVE TO DO IT FOR ALL GHOSTS AS WELL #TODO
        !
        !CALL BOUNDARY(LOCAL_NUMP, LOCAL_FINT, LOCAL_ACL, LOCAL_VEL, LOCAL_DSP, LOCAL_EBC)
        !
        ! CORRECTOR ALGORITHM
        !
        !CALL CORRECTOR(TOTAL_LOCAL_NUMP,LOCAL_VEL,LOCAL_ACL,DLT)
        !
        ! UPDATE TIME QUANTITIES
        !
        !
        ! GET THE CURRENT REMAINING TIME ESTIMATE
        !
        CALL EST_TIME(NCORES_INPUT,TIMER_STEPS,SIM_TIME_1,SIM_TIME_2, &
            REAL_TIME_1,REAL_TIME_2, &
            SIM_TIME_LEFT,TIME,TIME_END, &
            REAL_TIME_REMAINING,LINIT)
        !
        ! OUTPUT OPERATIONS
        !
        !
        ! CALCULATE ENERGIES
        ! NEED TO CALCULATE THE EXTERNAL ENERGY SO THAT THE ENERGIES ADD UP #TODO
        !
        !
        DLOCAL_INT_ENERGY = 0.0d0
        DLOCAL_KIN_ENERGY = 0.0d0
        DO I=1,TOTAL_LOCAL_NUMP*3

            IF (.NOT.LFINITE_STRAIN) THEN
                DLOCAL_INT_ENERGY = DLOCAL_INT_ENERGY + 0.50d0*(LOCAL_FINT(I))*LOCAL_DSP_TOT(I)
            END IF

            !THE BELOW IS ONE FORMULA IN BELYSCHKO'S BOOK
            !DLOCAL_INT_ENERGY = DLOCAL_INT_ENERGY + 0.50d0*(LOCAL_FINT(I)+LOCAL_FINT_NMO(I))*LOCAL_DSP(I)

            !THE BELOW IS ONE FORMULA IN BELYSCHKO'S BOOK
            !DLOCAL_INT_ENERGY = DLOCAL_INT_ENERGY + 0.50d0*DLT*(LOCAL_FINT(I)+LOCAL_FINT_NMO(I))*LOCAL_VEL(I)

            DLOCAL_KIN_ENERGY = DLOCAL_KIN_ENERGY + 0.50d0*LOCAL_MASS(I)*LOCAL_VEL(I)**2
            DLOCAL_EXT_ENERGY = DLOCAL_EXT_ENERGY + 0.50d0*(LOCAL_FEXT(I)+LOCAL_FEXT_NMO(I))*LOCAL_DSP(I)

        END DO

        IF (LFINITE_STRAIN) THEN
            DLOCAL_INT_ENERGY = DLOCAL_INT_ENERGY + LOCAL_INT_WORK
        END IF

        !
        !BELOW IS THE WORK WITH INTEGRATION IN THE INTERNAL FORCE ROUTINE, DOES NOT WORK ATM, CHECK
        !BELYSCHKO'S BOOK. THIS SHOULD BE THE CORRECT WAY TO DO IT FOR FINITE-STRAIN NONLINEAR PROBLEMS
        !
        !DLOCAL_INT_ENERGY = DLOCAL_INT_ENERGY + LOCAL_INT_WORK
        !
        DLOCAL_TOTAL_ENERGY = DLOCAL_INT_ENERGY + DLOCAL_KIN_ENERGY - DLOCAL_EXT_ENERGY

        WRITE(121,'(I20,4E20.5)') STEPS, DLOCAL_INT_ENERGY ,DLOCAL_KIN_ENERGY,DLOCAL_EXT_ENERGY,DLOCAL_TOTAL_ENERGY
        !
        IF (TIME_COUNTER.GT.TIME_OUTPUT) THEN
            !
            !
            IF (LINIT_TIME) THEN

                CALL EST_TIME(NCORES_INPUT,TIMER_STEPS,SIM_TIME_1,SIM_TIME_2, &
                    REAL_TIME_1,REAL_TIME_2, &
                    SIM_TIME_LEFT,TIME,TIME_END, &
                    REAL_TIME_REMAINING,LINIT_TIME)
                LINIT_TIME = .FALSE.
            END IF
            ! ASSEMBLE LOCAL QUANTATIES INTO GLOBAL FOR OUTPUT, LATER THIS LIKELY
            ! NEEDS TO BE MODIFIED
            !
            !CALL OUTPUT_ASSEMBLER(HPC_SCHEME, LOCAL_NUMP, MODEL_NUMP, TOTAL_MODEL_MAP, &
            !                         LOCAL_ACL,LOCAL_VEL,LOCAL_DSP,LOCAL_DSP_TOT,LOCAL_COO_CURRENT,LOCAL_FINT, &
            !                         MODEL_ACL,MODEL_VEL,MODEL_DSP,MODEL_DSP_TOT,MODEL_COO_CURRENT,MODEL_FINT)
            !
            !
            !******************************************OUTPUT******************************************
            !
            ! CALL THE SUBROUTINE TO OUTPUT TO THE EXODUS FILE
            !
            !  (COMMENT THE FOLLOWING IF EXODUS NOT INSTALLED, DO NOT INCLUDE OUTPUT.F90) #NOEXODUS
            !call output_Step(        LOCAL_NUMP,LOCAL_COO_CURRENT,MODEL_NODE_IDS,LOCAL_DSP_TOT_PHY,LOCAL_VEL_PHY, &
            !                         LOCAL_ACL_PHY, LOCAL_FINT, LOCAL_EBC,  LOCAL_BODY_ID,  LOCAL_MAT_TYPE,  LOCAL_COO,    &
            !                         LOCAL_VINIT,  LOCAL_WIN,  LOCAL_VOL,  LOCAL_MASS,  LOCAL_PROP,    &
            !                         LOCAL_CHAR_DIST, LOCAL_WAVE_VEL)  !TEMP
            !
            ! CALL THE SUBROUTINE TO OUTPUT TO THE VTK FILE
            !
            IF (HPC_SCHEME.EQ.1) THEN
     ! Ensure all data is synchronized before output
     !$ACC UPDATE HOST(LOCAL_DSP_TOT_PHY, LOCAL_VEL_PHY, LOCAL_ACL_PHY)
     !$ACC UPDATE HOST(LOCAL_STRESS, LOCAL_STRAIN, LOCAL_STATE, LOCAL_STRAIN_EQ)
     !$ACC UPDATE HOST(MODEL_ELCON)
                IF (UNF_OUTPUT) THEN
                    call UNF_OUTPUT_STEP_VTK(exodusStep,LOCAL_NUMP,MODEL_NUMEL, MODEL_ELCON, LOCAL_COO_CURRENT,MODEL_NODE_IDS,LOCAL_DSP_TOT_PHY,LOCAL_VEL_PHY, &
                        LOCAL_ACL_PHY, LOCAL_FINT, LOCAL_EBC,  LOCAL_BODY_ID,  LOCAL_MAT_TYPE,  LOCAL_COO,    &
                        LOCAL_VINIT,  MODEL_NORM_WIN,  LOCAL_WIN,  LOCAL_VOL,  LOCAL_MASS,  LOCAL_PROP,    &
                        LOCAL_CHAR_DIST, LOCAL_WAVE_VEL, LOCAL_PRFORCE, LOCAL_STRESS, LOCAL_STRAIN,LOCAL_STATE, LOCAL_STRAIN_EQ, &
				LOCAL_IJKSPC)
                ELSE
                    call OUTPUT_STEP_VTK(exodusStep,LOCAL_NUMP,MODEL_NUMEL, MODEL_ELCON, LOCAL_COO_CURRENT,MODEL_NODE_IDS,LOCAL_DSP_TOT_PHY,LOCAL_VEL, &
                        LOCAL_ACL, LOCAL_FINT, LOCAL_EBC,  LOCAL_BODY_ID,  LOCAL_MAT_TYPE,  LOCAL_COO,    &
                        LOCAL_VINIT,  MODEL_NORM_WIN,  LOCAL_WIN,  LOCAL_VOL,  LOCAL_MASS,  LOCAL_PROP,    &
                        LOCAL_CHAR_DIST, LOCAL_WAVE_VEL, LOCAL_PRFORCE, LOCAL_STRESS, LOCAL_STRAIN,LOCAL_STATE, LOCAL_STRAIN_EQ, &
				LOCAL_IJKSPC)
                END IF

            END IF
            !
            !******************************************OUTPUT******************************************
            !
            TIME_COUNTER = 0.0d0
            !
            ! OUTPUT TO LOG FILE AND SCREEN
            !
            OPEN(50,FILE='MEGA.log',STATUS='OLD',POSITION='APPEND')
            !
            !MAXV = MAXVAL(DABS(LOCAL_VEL))
            !MAXD = MAXVAL(DABS(LOCAL_DSP_TOT))
            !WRITE(50,'(2F15.5,3E20.10)') MAXV, MAXD, &
            !                    DLOCAL_INT_ENERGY, DLOCAL_KIN_ENERGY, DLOCAL_TOTAL_ENERGY
            !WRITE(*,'(2F15.5,3E20.10)') MAXV, MAXD, &
            !                    DLOCAL_INT_ENERGY, DLOCAL_KIN_ENERGY, DLOCAL_TOTAL_ENERGY

            CALL MAKE_CTIME(CTIME_ALL,REAL_TIME_REMAINING,TIMER_FLAG)

            IF (TIMER_FLAG) CALL WARN('TIME REMAINING > 10000 DAYS')

            !WRITE(50,'(2I8,2E15.5,4E15.5,A15)') exodusStep, STEPS, TIME, DLT, &
            !    DLOCAL_INT_ENERGY, DLOCAL_KIN_ENERGY, DLOCAL_EXT_ENERGY, DLOCAL_TOTAL_ENERGY,CTIME_ALL
            !WRITE(*,'(2I8,2E15.5,4E15.5,A15)') exodusStep, STEPS, TIME, DLT, &
            !    DLOCAL_INT_ENERGY, DLOCAL_KIN_ENERGY, DLOCAL_EXT_ENERGY, DLOCAL_TOTAL_ENERGY,CTIME_ALL
            WRITE(50,'(2I8,2E15.5,A30)') exodusStep, STEPS, TIME, DLT, &
                CTIME_ALL
            WRITE(*,'(2I8,2E15.5,A30)') exodusStep, STEPS, TIME, DLT, &
                CTIME_ALL
				
            CLOSE(50)

            !
            exodusStep=exodusStep+1
        END IF
        

        !
        ! WRAP-UP OPERATIONS FOR THIS TIME-STEP
        !
        TIME = TIME + DLT
        TIME_COUNTER = TIME_COUNTER + DLT
        !
        STEPS = STEPS + 1
        TIMER_STEPS = TIMER_STEPS + 1
        !
        LINIT = .FALSE.
        !$ACC UPDATE DEVICE(LINIT)


        WRITE(120,*) TIME


        IF (TIME.GT.TIME_END) EXIT
        !
    END DO

    !
    ! End OpenACC data region - copy results back to host
    !
    !$ACC END DATA
    !

	!GC
	DEALLOCATE(MODEL_ELCON)
	
    CLOSE(120)
    CLOSE(121)
    CLOSE(122) !KC FOR TOTATL FORCES

    CALL CPU_TIME(REAL_TIME_1)

    TOTAL_REAL_TIME = REAL_TIME_1 - REAL_TIME_0

    CALL MAKE_CTIME(CTIME_ALL,TOTAL_REAL_TIME,TIMER_FLAG)


    CALL WRITE_OUT('MEGA COMPLETE')
    WRITE(CTOTAL_REAL_TIME, '(A27,A13)') 'TOTAL TIME FOR SIMULATION: ', CTIME_ALL
    CALL WRITE_OUT(CTOTAL_REAL_TIME)
    !

    CONTINUE

    END PROGRAM
