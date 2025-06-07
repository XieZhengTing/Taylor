
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
    !
    IMPLICIT NONE
    !
    INTEGER:: I, J, M, ITEMP, I_STEP, STEP_NUM
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
    INTEGER :: ierr_handeler ! Error flag from HANDELER
    DOUBLE PRECISION :: FEXT_M_AT_ENTRY ! Declare FEXT_M_AT_ENTRY in the main program's declaration block



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
    !
    ! GET INITIAL VALUES
    !
    CALL STATE_FEILD_INIT(TOTAL_LOCAL_SIZE,TOTAL_LOCAL_NUMP, LOCAL_NUMP, MODEL_NUMP, MODEL_VINIT, TOTAL_MODEL_MAP, MODEL_COO,MODEL_MASS,MODEL_EBC, &
        MODEL_NONZERO_EBC,LOCAL_STATE,LOCAL_STRESS,LOCAL_STRAIN,LOCAL_H_STRESS,LOCAL_S_STRESS,LOCAL_COO,LOCAL_FEXT,LOCAL_FINT, &
        LOCAL_ACL,LOCAL_VEL,LOCAL_DSP,LOCAL_DSP_TOT,LOCAL_DSP_TOT_PHY,LOCAL_MASS,LOCAL_COO_CURRENT, LOCAL_EBC, LOCAL_NONZERO_EBC,LOCAL_EBC_NODES, &
        LOCAL_PRFORCE, LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, LOCAL_DX_STRAIN, LOCAL_DY_STRAIN, LOCAL_DZ_STRAIN,LOCAL_STRAIN_EQ)

    !
    DEALLOCATE(MODEL_VINIT,MODEL_MASS,MODEL_EBC,MODEL_COO,MODEL_NONZERO_EBC)

    !
    TIME = 0.0d0
    STEP_NUM = 0
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
    LINIT_TIME = .TRUE.
    SIM_TIME_1 = 0.0d0
    CALL CPU_TIME(REAL_TIME_0)
    CALL CPU_TIME(REAL_TIME_1)
    !
    ! START TIME INTEGRATION
    !
    CALL WRITE_OUT('STARTING TIME INTEGRATION')
    !
    WRITE(*,*)
    WRITE(*,*) 'DEBUG: MAIN - Array size checks before ACC DATA region:'
    WRITE(*,*) '  TOTAL_LOCAL_SIZE = ', TOTAL_LOCAL_SIZE
    WRITE(*,*) '  Expected FEXT/FINT size = ', 3*TOTAL_LOCAL_SIZE
    
    IF (ALLOCATED(LOCAL_FEXT)) THEN
        WRITE(*,*) '  Actual LOCAL_FEXT size = ', SIZE(LOCAL_FEXT)
        IF (SIZE(LOCAL_FEXT) < 3*TOTAL_LOCAL_SIZE) THEN ! Use < for insufficient size
            WRITE(*,*) 'FATAL ERROR: MAIN - LOCAL_FEXT size mismatch! Allocated: ', SIZE(LOCAL_FEXT), ' Required: ', 3*TOTAL_LOCAL_SIZE
            CALL EXIT_PROGRAM('LOCAL_FEXT size error', -33)
        END IF
    ELSE
        WRITE(*,*) 'FATAL ERROR: MAIN - LOCAL_FEXT not allocated before ACC DATA region!'
        CALL EXIT_PROGRAM('LOCAL_FEXT not allocated', -34)
    END IF
    
    IF (ALLOCATED(LOCAL_FINT)) THEN
        WRITE(*,*) '  Actual LOCAL_FINT size = ', SIZE(LOCAL_FINT)
        IF (SIZE(LOCAL_FINT) < 3*TOTAL_LOCAL_SIZE) THEN ! Use < for insufficient size
            WRITE(*,*) 'FATAL ERROR: MAIN - LOCAL_FINT size mismatch! Allocated: ', SIZE(LOCAL_FINT), ' Required: ', 3*TOTAL_LOCAL_SIZE
            CALL EXIT_PROGRAM('LOCAL_FINT size error', -35)
        END IF
    ELSE
        WRITE(*,*) 'FATAL ERROR: MAIN - LOCAL_FINT not allocated before ACC DATA region!'
        CALL EXIT_PROGRAM('LOCAL_FINT not allocated', -36)
    END IF
    WRITE(*,*)
    ! Ensure critical arrays are allocated before entering ACC DATA region
    IF (.NOT. ALLOCATED(LOCAL_COO)) THEN
        WRITE(*,*) 'FATAL ERROR: MAIN - LOCAL_COO not allocated before ACC DATA region.'
        CALL EXIT_PROGRAM('Array LOCAL_COO not allocated', -30)
    END IF
    IF (.NOT. ALLOCATED(LOCAL_FINT)) THEN
        WRITE(*,*) 'FATAL ERROR: MAIN - LOCAL_FINT not allocated before ACC DATA region.'
        CALL EXIT_PROGRAM('Array LOCAL_FINT not allocated', -31)
    END IF
    IF (.NOT. ALLOCATED(LOCAL_FEXT)) THEN
        WRITE(*,*) 'FATAL ERROR: MAIN - LOCAL_FEXT not allocated before ACC DATA region.'
        CALL EXIT_PROGRAM('Array LOCAL_FEXT not allocated', -32)
    END IF
    ! Add checks for other critical allocatable arrays as needed...

    ! Initialize NMO arrays on host before entering data region
    LOCAL_FINT_NMO = 0.0D0
    LOCAL_FEXT_NMO = 0.0D0
    LOCAL_FINT = 0.0D0
    LOCAL_FEXT = 0.0D0
    LOCAL_DLT = 0.0d0 ! Corresponds to DLT_FINT from HANDELER/CONSTRUCT_FINT
    LOCAL_CHAR_DIST = 0.0d0
    LOCAL_WAVE_VEL = 0.0d0
    DLOCAL_INT_ENERGY = 0.0d0
    DLOCAL_KIN_ENERGY = 0.0d0
    DLOCAL_EXT_ENERGY = 0.0d0
    LOCAL_INT_WORK = 0.0d0
    TOTAL_FORCE = 0.0d0

!$ACC DATA &
      !$ACC COPYIN(LOCAL_COO, MODEL_BODYFORCE, MODEL_ELCON, MODEL_NODE_IDS, MODEL_NORM_WIN, &
      !$ACC        LOCAL_SM_LEN, LOCAL_SM_AREA, LOCAL_SM_VOL, LOCAL_WIN, LOCAL_VOL, &
      !$ACC        LOCAL_NSNI_FAC, LOCAL_VINIT, LOCAL_MAT_TYPE, LOCAL_PROP, LOCAL_BODY_ID, &
      !$ACC        LOCAL_X_MOM, LOCAL_Y_MOM, LOCAL_Z_MOM, LOCAL_IJKSPC, &
      !$ACC        LOCAL_XDIST_MAX, LOCAL_YDIST_MAX, LOCAL_ZDIST_MAX, &
      !$ACC        LOCAL_EBC, LOCAL_NONZERO_EBC, LOCAL_GHOST) &
      !$ACC COPY(LOCAL_STATE, LOCAL_STRESS, LOCAL_STRAIN, LOCAL_H_STRESS, LOCAL_S_STRESS, &
      !$ACC      LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
      !$ACC      LOCAL_ACL, LOCAL_VEL, LOCAL_DSP, LOCAL_DSP_TOT, LOCAL_DSP_TOT_PHY, &
      !$ACC      LOCAL_COO_CURRENT, LOCAL_PRFORCE, LOCAL_MASS, LOCAL_STRAIN_EQ, &
      !$ACC      LOCAL_EBC_NODES, LOCAL_CHAR_DIST, LOCAL_WAVE_VEL, &
      !$ACC      LOCAL_FINT, LOCAL_FEXT, &
      !$ACC      TIME, DLT, LINIT, LINIT_TIME, DO_INTERP, TIME_COUNTER, STEPS, TIMER_STEPS, exodusStep) &
      !$ACC CREATE(LOCAL_FINT_NMO, LOCAL_FEXT_NMO, LOCAL_ACL_PHY, LOCAL_VEL_PHY, LOCAL_DINC_PHY) &
      !$ACC COPYOUT(LOCAL_DLT, DLOCAL_INT_ENERGY, DLOCAL_KIN_ENERGY, DLOCAL_EXT_ENERGY, &
      !$ACC          DLOCAL_TOTAL_ENERGY, TOTAL_FORCE)
    !WRITE(50,'(2A8,2A15,4A15)') 'Out Step','Step','Time','Delta T','Int Engery', 'Kin Energy', 'Ext Energy', 'Tot Energy', 'Est Time Rem'
    !WRITE(*,'(2A8,2A15,4A15)')  'Out Step','Step','Time','Delta T','Int Engery', 'Kin Energy', 'Ext Energy', 'Tot Energy', 'Est Time Rem'
    !
    WRITE(50,'(2A8,2A15,1A30)') 'Out Step','Step','Time','Delta T', 'Estimated Time Remaining'
    WRITE(*,'(2A8,2A15,1A30)')  'Out Step','Step','Time','Delta T', 'Estimated Time Remaining'
    
	!

    !
    DO
        !
        !$ACC UPDATE DEVICE(TIME, DLT, LINIT, LINIT_TIME, DO_INTERP, TIME_COUNTER, STEPS, TIMER_STEPS, exodusStep) ! Ensure scalars are up-to-date on device

        !
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
!            !$ACC UPDATE DEVICE(LINIT)
            LOCAL_DSP = 0.0d0 ! TO GET ZERO FINT, JUST TO GET SHP
            ! Explicitly initialize FINT and FEXT on the device since they are now CREATE'd
            !$ACC PARALLEL LOOP DEFAULT(PRESENT)
            DO I = 1, 3*TOTAL_LOCAL_SIZE
                LOCAL_FINT(I) = 0.0D0
                LOCAL_FEXT(I) = 0.0D0
            END DO
            !$ACC END PARALLEL LOOP
            DO_INTERP = .FALSE.
            ! LOCAL_DSP is already on device and modified. DO_INTERP is scalar, will be passed by value.


            ! Detailed debug information and parameter checks
            WRITE(*,*)
            WRITE(*,*) 'DEBUG: MAIN - Before initial HANDELER call (LINIT=.TRUE.):'
            WRITE(*,*) '  LINIT              = ', LINIT
            WRITE(*,*) '  LOCAL_NUMP         = ', LOCAL_NUMP
            WRITE(*,*) '  TOTAL_LOCAL_SIZE   = ', TOTAL_LOCAL_SIZE
            IF (ALLOCATED(LOCAL_FINT)) THEN
                WRITE(*,*) '  ALLOCATED(LOCAL_FINT) = .TRUE., SIZE = ', SIZE(LOCAL_FINT)
            ELSE
                WRITE(*,*) '  ALLOCATED(LOCAL_FINT) = .FALSE.'
            END IF
            IF (ALLOCATED(LOCAL_FEXT)) THEN
                WRITE(*,*) '  ALLOCATED(LOCAL_FEXT) = .TRUE., SIZE = ', SIZE(LOCAL_FEXT)
            ELSE
                WRITE(*,*) '  ALLOCATED(LOCAL_FEXT) = .FALSE.'
            END IF
            WRITE(*,*)
            ! Critical check before calling HANDELER
            IF (LOCAL_NUMP <= 0) THEN
                WRITE(*,*) 'FATAL ERROR: MAIN - LOCAL_NUMP <= 0 before initial HANDELER call. LOCAL_NUMP = ', LOCAL_NUMP
                CALL EXIT_PROGRAM('Invalid LOCAL_NUMP in MAIN before HANDELER', -20)
            END IF
            ! Check FEXT array validity
            IF (.NOT. ALLOCATED(LOCAL_FEXT)) THEN
                WRITE(*,*) 'FATAL ERROR: MAIN - LOCAL_FEXT not allocated before HANDELER call.'
                CALL EXIT_PROGRAM('LOCAL_FEXT allocation error in MAIN', -21)
            END IF
            IF (ALLOCATED(LOCAL_FEXT) .AND. (SIZE(LOCAL_FEXT) < 3*TOTAL_LOCAL_SIZE)) THEN ! Assuming FEXT should be 3*TOTAL_LOCAL_SIZE
                WRITE(*,*) 'FATAL ERROR: MAIN - LOCAL_FEXT size insufficient.'
                WRITE(*,*) '  Required at least: ', 3*TOTAL_LOCAL_SIZE, ' Available: ', SIZE(LOCAL_FEXT)
                CALL EXIT_PROGRAM('LOCAL_FEXT size error in MAIN', -22)
            END IF
            CALL HANDELER(       LOCAL_WIN,      LOCAL_VOL,         LOCAL_NUMP,     LOCAL_COO,      LOCAL_COO_CURRENT,     &
                LOCAL_SM_LEN,    LOCAL_SM_AREA,  LOCAL_SM_VOL,      LOCAL_NSNI_FAC, LOCAL_GHOST,       &
                LOCAL_PROP,      LOCAL_STATE,    LOCAL_STRESS,      LOCAL_STRAIN,   LOCAL_H_STRESS,LOCAL_S_STRESS,   LOCAL_DSP, LOCAL_DSP_TOT,        &
                LOCAL_FINT,      LOCAL_MAT_TYPE, TOTAL_LOCAL_SIZE,  LINIT,          LOCAL_DLT, &
                LOCAL_FEXT,      LOCAL_CHAR_DIST, LOCAL_WAVE_VEL, LOCAL_EBC_NODES, &
                DO_INTERP,       LOCAL_DINC_PHY, LOCAL_VEL_PHY, LOCAL_VEL, LOCAL_ACL, LOCAL_ACL_PHY, &
                LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
                LOCAL_X_MOM,     LOCAL_Y_MOM,    LOCAL_Z_MOM,  &
                LOCAL_XDIST_MAX, LOCAL_YDIST_MAX, LOCAL_ZDIST_MAX,MODEL_BODYFORCE,DLT, LOCAL_INT_WORK, LOCAL_BODY_ID, LOCAL_STRAIN_EQ, &
				LOCAL_IJKSPC, ierr_handeler)

            ! After host execution of HANDELER (LINIT=.TRUE.), update device with newly computed/allocated SAVE arrays

            ! LOCAL_DLT is now COPYOUT, will be updated at END DATA or via explicit UPDATE HOST if needed sooner.


           IF (ierr_handeler .NE. 0) THEN
                WRITE(*,*) 'FATAL ERROR: MAIN - Error returned from HANDELER during LINIT. ierr_handeler = ', ierr_handeler

                CALL EXIT_PROGRAM('Error during HANDELER initialization', ierr_handeler)
            END IF
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
!        !$ACC PARALLEL LOOP DEFAULT(PRESENT) PRIVATE(I) COLLAPSE(1) &
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) PRIVATE(M) &
        !$ACC IF( (.NOT.LINIT) .OR. (.NOT.AUTO_TS) )
!        IF (LINIT.AND.(AUTO_TS)) THEN
!            !WE DONT HAVE A TIME STEP ESTIMATE YET
!            LOCAL_DSP = 0.0d0
!        ELSE
!            !PREDICT THE DISPLACEMENT INCREMENT AND VELOCITY FROM THE PREVIOUS ACCELERATION
!            LOCAL_DSP = DLT * LOCAL_VEL + DLT**2 * 0.5d0 * LOCAL_ACL
!            LOCAL_DSP_TOT = LOCAL_DSP_TOT + LOCAL_DSP
!            LOCAL_VEL = LOCAL_VEL + DLT * 0.5d0 * LOCAL_ACL
!        END IF
!         !$ACC END PARALLEL LOOP IF( (.NOT.LINIT) .OR. (.NOT.AUTO_TS) )
        DO M=1, 3*TOTAL_LOCAL_SIZE ! Loop over all components of vector arrays
            IF (LINIT.AND.(AUTO_TS)) THEN
                !WE DONT HAVE A TIME STEP ESTIMATE YET
                LOCAL_DSP(M) = 0.0d0
            ELSE
                !PREDICT THE DISPLACEMENT INCREMENT AND VELOCITY FROM THE PREVIOUS ACCELERATION
                LOCAL_DSP(M) = DLT * LOCAL_VEL(M) + DLT**2 * 0.5d0 * LOCAL_ACL(M)
                LOCAL_DSP_TOT(M) = LOCAL_DSP_TOT(M) + LOCAL_DSP(M)
                LOCAL_VEL(M) = LOCAL_VEL(M) + DLT * 0.5d0 * LOCAL_ACL(M)
            END IF
        END DO
        !$ACC END PARALLEL LOOP ! Removed the IF clause from END PARALLEL LOOP


        !
        !COMPUTE THE PREDICTED CURRENT COORD FOR SEMI-LAG SHP CACULATION
        !
        DO_INTERP = .TRUE.
        LINIT_TEMP = .FALSE.
        ! HANDELER's DO_INTERP part will run its own !$ACC PARALLEL LOOP
        ! Ensure inputs to this HANDELER call are correct on host/device
        ! These are already on device via COPY or modified in previous steps. Scalars passed by value.

        CALL HANDELER(      LOCAL_WIN,       LOCAL_VOL,      LOCAL_NUMP,        LOCAL_COO,      LOCAL_COO_CURRENT,     &
            LOCAL_SM_LEN,    LOCAL_SM_AREA,  LOCAL_SM_VOL,      LOCAL_NSNI_FAC, LOCAL_GHOST,       &
            LOCAL_PROP,      LOCAL_STATE,    LOCAL_STRESS,      LOCAL_STRAIN,   LOCAL_H_STRESS, LOCAL_S_STRESS,   LOCAL_DSP,  LOCAL_DSP_TOT,      &
            LOCAL_FINT,      LOCAL_MAT_TYPE, TOTAL_LOCAL_SIZE,  LINIT_TEMP,          LOCAL_DLT, &
            LOCAL_FEXT,      LOCAL_CHAR_DIST, LOCAL_WAVE_VEL, LOCAL_EBC_NODES, &
            DO_INTERP,       LOCAL_DINC_PHY, LOCAL_VEL_PHY, LOCAL_VEL, LOCAL_ACL, LOCAL_ACL_PHY, &
            LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
            LOCAL_X_MOM,     LOCAL_Y_MOM,    LOCAL_Z_MOM, &
            LOCAL_XDIST_MAX, LOCAL_YDIST_MAX, LOCAL_ZDIST_MAX,MODEL_BODYFORCE,DLT, LOCAL_INT_WORK, LOCAL_BODY_ID, LOCAL_STRAIN_EQ, &
			LOCAL_IJKSPC, ierr_handeler)

        ! LOCAL_DINC_PHY is COPY, will be updated at END DATA. If needed sooner: !$ACC UPDATE HOST(LOCAL_DINC_PHY)

        !$ACC PARALLEL LOOP DEFAULT(PRESENT)
        DO I = 1, 3*TOTAL_LOCAL_SIZE
            LOCAL_DSP_TOT_PHY(I) = LOCAL_DSP_TOT_PHY(I) + LOCAL_DINC_PHY(I)
        END DO
        !$ACC END PARALLEL LOOP
        !
        ! LIKELY HAVE TO UPDATE GHOSTS HERE FOR GHOST SCHEMES, THESE ARRAYS
        ! SHOULD BE DIFFERENT SIZES THEN #TODO
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) PRIVATE(J) COLLAPSE(1)
        DO I=1,TOTAL_LOCAL_SIZE ! Should this be LOCAL_NUMP or TOTAL_LOCAL_SIZE? Original was LOCAL_NUMP

            DO J=1,3
                LOCAL_COO_CURRENT(J,I) = LOCAL_COO(J,I) + LOCAL_DSP_TOT_PHY((I-1)*3+J)
            END DO
        END DO
        !$ACC END PARALLEL LOOP
        ! LOCAL_COO_CURRENT is COPY, changes on device are kept.


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
!        !$ACC UPDATE DEVICE(DO_INTERP, LINIT_TEMP)
        IF(LINIT) THEN
            LOCAL_FINT_NMO = 0.0d0
            LOCAL_FEXT_NMO = 0.0d0
        ELSE
            LOCAL_FINT_NMO = LOCAL_FINT
            LOCAL_FEXT_NMO = LOCAL_FEXT !NOT USED! NEED TO FILL THIS OUT IN HANDLER
        END IF        
        ! Arguments are passed. LOCAL_FINT_NMO, LOCAL_FEXT_NMO are CREATE, modified on device.

        CALL HANDELER(      LOCAL_WIN,       LOCAL_VOL,      LOCAL_NUMP,        LOCAL_COO,      LOCAL_COO_CURRENT,     &
            LOCAL_SM_LEN,    LOCAL_SM_AREA,  LOCAL_SM_VOL,      LOCAL_NSNI_FAC, LOCAL_GHOST,       &
            LOCAL_PROP,      LOCAL_STATE,    LOCAL_STRESS,      LOCAL_STRAIN,   LOCAL_H_STRESS,LOCAL_S_STRESS,   LOCAL_DSP, LOCAL_DSP_TOT,        &
            LOCAL_FINT,      LOCAL_MAT_TYPE, TOTAL_LOCAL_SIZE,  LINIT_TEMP,     LOCAL_DLT, &
            LOCAL_FEXT,      LOCAL_CHAR_DIST, LOCAL_WAVE_VEL, LOCAL_EBC_NODES, &
            DO_INTERP,       LOCAL_DINC_PHY, LOCAL_VEL_PHY, LOCAL_VEL, LOCAL_ACL, LOCAL_ACL_PHY, &
            LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
            LOCAL_X_MOM,     LOCAL_Y_MOM,    LOCAL_Z_MOM,  &
            LOCAL_XDIST_MAX, LOCAL_YDIST_MAX, LOCAL_ZDIST_MAX,MODEL_BODYFORCE,DLT, LOCAL_INT_WORK, LOCAL_BODY_ID, LOCAL_STRAIN_EQ, &
			LOCAL_IJKSPC, ierr_handeler)

        ! After HANDELER, LOCAL_FINT, LOCAL_FEXT, LOCAL_DLT, LOCAL_STATE, LOCAL_STRESS etc. are updated on device.
        ! If DLT needs to be used on host immediately: !$ACC UPDATE HOST(LOCAL_DLT)



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
        ! DLT and PDSEARCH are scalars, passed by value to kernels if needed.
        ! If they are part of COPY, their device versions are updated.


!$ACC PARALLEL LOOP DEFAULT(PRESENT) &
        !$ACC PRIVATE(J,M,I_STEP,STEP_NUM,MPDC) &
        !$ACC GANG WORKER VECTOR_LENGTH(32)
        DO I=1,LOCAL_NUMP
            !$ACC LOOP SEQ PRIVATE(M,I_STEP,STEP_NUM,MPDC)
            DO J=1,3
                M = (I-1)*3+J

                IF (LOCAL_EBC(J,I).EQ.1) THEN
                    ! 固定邊界條件
                    LOCAL_VEL(M) = 0.0d0
                    LOCAL_ACL(M) = 0.0d0
                    LOCAL_PRFORCE(J,I) = - LOCAL_FINT(M)
                    LOCAL_FEXT(M) = - LOCAL_FINT(M)

                ELSEIF (LOCAL_EBC(J,I).EQ.2) THEN 
                    ! 非零位移邊界條件
                    STEP_NUM = 0
                    ! 找出當前時間所在的時間區間
                    DO I_STEP = 1, FIXITY2_STEPS
                        IF (FIXITY2_TIME(I_STEP,1) .LT. TIME .AND. TIME .LT. FIXITY2_TIME(I_STEP,2)) THEN
                            STEP_NUM = I_STEP
                            EXIT
                        END IF
                    END DO

                    ! 指定速度
                    IF (1 .LE. STEP_NUM .AND. STEP_NUM .LE. FIXITY2_STEPS) THEN
                        LOCAL_VEL(M) = FIXITY2_NONZERO_EBC(STEP_NUM,J)/(FIXITY2_TIME(STEP_NUM,2) - FIXITY2_TIME(STEP_NUM,1))
                    ELSE
                        LOCAL_VEL(M) = 0.0d0
                    END IF

                    STEP_NUM = 0 ! 重新初始化

                    LOCAL_ACL(M) = 0.0d0
                    LOCAL_PRFORCE(J,I) = - LOCAL_FINT(M)
                    LOCAL_FEXT(M) = - LOCAL_FINT(M)

                ELSE   
                    ! 自由節點
                    MPDC = LOCAL_PROP(21,I)*LOCAL_MASS(M) ! 質量比例阻尼
                    
                    LOCAL_ACL(M) = (LOCAL_FEXT(M) - LOCAL_FINT(M) - MPDC*LOCAL_VEL(M))/ &
                                   (LOCAL_MASS(M)+0.5d0*DLT*MPDC)
                                   
                    LOCAL_VEL(M) = LOCAL_VEL(M) + 0.5d0*DLT*LOCAL_ACL(M)
                    
                    LOCAL_PRFORCE(J,I) = 0.0d0
                END IF
                
                ! 累加總力（所有邊界條件類型）
                !$ACC ATOMIC UPDATE
                TOTAL_FORCE(J,LOCAL_BODY_ID(I)) = TOTAL_FORCE(J,LOCAL_BODY_ID(I)) + LOCAL_PRFORCE(J,I)

            END DO ! J loop
        END DO ! I loop
        !$ACC END PARALLEL LOOP

        ! Original single loop structure (commented out for reference after applying the split)
        ! !$ACC PARALLEL LOOP DEFAULT(PRESENT) &
        ! !$ACC PRIVATE(J,M,I_STEP,STEP_NUM,MPDC, FEXT_M_AT_ENTRY) REDUCTION(+:TOTAL_FORCE)
        ! DO I=1,LOCAL_NUMP
        !     !$ACC LOOP SEQ ! Explicitly mark inner J loop as sequential for this I
        !     DO J=1,3
        !         M = (I-1)*3+J
        !         FEXT_M_AT_ENTRY = LOCAL_FEXT(M) ! Store the value of LOCAL_FEXT(M) at the beginning of this (I,J) iteration
        !         IF (LOCAL_EBC(J,I).EQ.1) THEN
        !             ...
        !         ELSEIF (LOCAL_EBC(J,I).EQ.2) THEN !NON-ZERO ESSENTIAL BC
        !             ...
        !         ELSE   !FREE
        !             LOCAL_FEXT(M) = FEXT_M_AT_ENTRY ! Explicitly assign LOCAL_FEXT(M) to its entry value in this branch
        !             ...
        !         END IF
        !     TOTAL_FORCE(J,LOCAL_BODY_ID(I))=TOTAL_FORCE(J,LOCAL_BODY_ID(I))+LOCAL_PRFORCE(J,I)
        !     END DO
        ! END DO
        ! !$ACC END PARALLEL LOOP
        ! 再處理 FREE 類型的節點的自由度 (LOCAL_EBC(J,I) == 0)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) PRIVATE(J,M,MPDC)
        DO I=1,LOCAL_NUMP
            !$ACC LOOP SEQ
            DO J=1,3
                M = (I-1)*3+J
                IF (LOCAL_EBC(J,I).EQ.0) THEN ! FREE node degree of freedom

				MPDC = LOCAL_PROP(21,I)*LOCAL_MASS(M) !MASS PROPORTIAL DAMPING MATRIX (DIAG)
				
                    LOCAL_ACL(M) = (LOCAL_FEXT(M) - LOCAL_FINT(M) - MPDC*LOCAL_VEL(M))/ &
					(LOCAL_MASS(M)+0.5d0*DLT*MPDC)
					
                    LOCAL_VEL(M) = LOCAL_VEL(M) + 0.5d0*DLT*LOCAL_ACL(M)

                    LOCAL_PRFORCE(J,I) = 0.0d0
                    !$ACC ATOMIC UPDATE
                    TOTAL_FORCE(J,LOCAL_BODY_ID(I)) = TOTAL_FORCE(J,LOCAL_BODY_ID(I)) + LOCAL_PRFORCE(J,I)

                END IF
                            
!            TOTAL_FORCE(J,LOCAL_BODY_ID(I))=TOTAL_FORCE(J,LOCAL_BODY_ID(I))+LOCAL_PRFORCE(J,I)

            END DO
        END DO

        !$ACC END PARALLEL LOOP


        ! TIME and TOTAL_FORCE are COPYOUT or COPY, will be updated at END DATA. If needed sooner:
    ! 只在需要輸出時才同步
    IF (MOD(STEPS, 100) == 0) THEN  ! 每 100 步輸出一次
        !$ACC UPDATE HOST(TOTAL_FORCE)    

        WRITE(122,'(E15.5,I8,3(E15.5),I8,3(E15.5))') TIME,LOCAL_BODY_ID(1), TOTAL_FORCE(1,LOCAL_BODY_ID(1)), TOTAL_FORCE(2,LOCAL_BODY_ID(1)),TOTAL_FORCE(3,LOCAL_BODY_ID(1)), LOCAL_BODY_ID(LOCAL_NUMP), TOTAL_FORCE(1,LOCAL_BODY_ID(LOCAL_NUMP)),TOTAL_FORCE(2,LOCAL_BODY_ID(LOCAL_NUMP)),TOTAL_FORCE(3,LOCAL_BODY_ID(LOCAL_NUMP))
    END IF

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

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) PRIVATE(I) &
        !$ACC REDUCTION(+:DLOCAL_INT_ENERGY, DLOCAL_KIN_ENERGY, DLOCAL_EXT_ENERGY)

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
        !$ACC END PARALLEL LOOP
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

        ! Energies and STEPS, TIME_COUNTER are COPYOUT or COPY.
        ! If needed for immediate host logic: !$ACC UPDATE HOST(STEPS, DLOCAL_INT_ENERGY, ...)

        !$ACC UPDATE HOST(DLOCAL_INT_ENERGY, DLOCAL_KIN_ENERGY, DLOCAL_EXT_ENERGY, DLOCAL_TOTAL_ENERGY)

        WRITE(121,'(I20,4E20.5)') STEPS, DLOCAL_INT_ENERGY ,DLOCAL_KIN_ENERGY,DLOCAL_EXT_ENERGY,DLOCAL_TOTAL_ENERGY
        !
!        !$ACC UPDATE HOST(TIME_COUNTER)
        IF (TIME_COUNTER >= TIME_OUTPUT) THEN ! Use >= for safer comparison with floating point

            !
            IF (LINIT_TIME) THEN

                CALL EST_TIME(NCORES_INPUT,TIMER_STEPS,SIM_TIME_1,SIM_TIME_2, &
                    REAL_TIME_1,REAL_TIME_2, &
                    SIM_TIME_LEFT,TIME,TIME_END, &
                    REAL_TIME_REMAINING,LINIT_TIME)
                LINIT_TIME = .FALSE.
                ! LINIT_TIME is COPY, device version updated.

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
                   ! 在呼叫 VTK 輸出副程式前，將裝置上的資料更新回主機
                   !$ACC UPDATE HOST(exodusStep, LOCAL_COO_CURRENT, LOCAL_DSP_TOT_PHY, LOCAL_VEL_PHY, LOCAL_ACL_PHY, & ! REVIEW THIS LIST
                   !$ACC&                  LOCAL_VEL, LOCAL_ACL, LOCAL_FINT, LOCAL_MASS, LOCAL_CHAR_DIST, LOCAL_WAVE_VEL, & ! REVIEW THIS LIST
                   !$ACC&                  LOCAL_PRFORCE, LOCAL_STRESS, LOCAL_STRAIN, LOCAL_STATE, LOCAL_STRAIN_EQ, & ! REVIEW THIS LIST
                   !$ACC&                  LOCAL_EBC, LOCAL_BODY_ID, LOCAL_MAT_TYPE, LOCAL_COO, LOCAL_VINIT, MODEL_NORM_WIN, & ! REVIEW THIS LIST
                   !$ACC&                  LOCAL_WIN, LOCAL_VOL, LOCAL_PROP, LOCAL_IJKSPC) ! REVIEW THIS LIST

                    ! All these variables are part of the main DATA region with COPY or COPYOUT.
                    ! If their most up-to-date values are needed for UNF_OUTPUT_STEP_VTK on host *before*
                    ! the main DATA region ends, then explicit UPDATE HOST would be needed here.
                    ! Otherwise, rely on END DATA synchronization.

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
            WRITE(50,'(2I8,2E15.5,A30)') exodusStep, STEPS, TIME, DLT, & ! These are host values

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

        LINIT = .FALSE.
        ! All these (TIME, TIME_COUNTER, STEPS, TIMER_STEPS, LINIT, exodusStep) are COPY, device versions updated.


        WRITE(120,*) TIME


        IF (TIME.GT.TIME_END) EXIT
        !
    END DO
    !$ACC END DATA

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
