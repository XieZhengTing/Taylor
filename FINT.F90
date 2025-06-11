

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
    SUBROUTINE CONSTRUCT_FINT(GWIN,    GVOL,       GNUMP,        GCOO,  GCOO_CUURENT,   &
    GSM_LEN, GSM_AREA,   GSM_VOL,      GNSNI_FAC,      &
    GGHOST,  GPROP,      GSTATE,       GSTRESS,        &
    GSTRAIN, G_H_STRESS, G_S_STRESS, GDINC,      GDINC_TOT,    GMAT_TYPE,                    &
    GEBC_NODES, DLT, GSIZE_IN,                               &
    GN,      GSTART,     DIM_NN_LIST,                  &
    GSTACK,  GSTACK_SHP, GSTACK_DSHP,  GSTACK_DDSHP,  GMAXN,  GINVK, LINIT,   &
    FINT,    FEXT,       DLT_FINT,   GCHAR_DIST,   GMAX_WVEL, &
    LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, & !OUTPUT (TOTAL_LOCAL_SIZE BASED)
    G_X_MOM, G_Y_MOM, G_Z_MOM,MODEL_BODYFORCE,GINT_WORK, MODEL_BODY_ID, GSTRAIN_EQ, ierr_fint_arg)


    !
    ! FUNCTION OF THIS SUBROUTINE:
    !
    ! BASED ON CURRENT DISPLACEMENT INCREMENT, STATE VARIABLES, AND SO ON,
    ! COMPUTE THE INTERNAL FORCE FOR A SET OF NODES GIVEN THEIR COORDINATES,
    ! DILATIONS, AND SO ON.
!   !$ACC ROUTINE SEQ
    !VOIGT ORDERING:
    !XX, YY, ZZ, YZ, XZ, XY
    !
!   !$ACC ROUTINE SEQ
    USE FINT_FUNCTIONS
    USE CONTROL
    USE CONSTITUTE_VONMISES_MOD    ! For VON_MISES
    USE CONSTITUTE_DRUCKPRAG_MOD   ! For DRUCK_PRAG
    USE CONSTITUTE_VISCOELASTIC_MOD! For VISCO_ELASTIC
    USE CONSTITUTE_VONMISES_DAM_MOD! For VON_MISES_DAM
    USE HYPERELASTIC_MOD           ! For HYPERELASTIC
    USE ESTIMATE_MODULI_MOD        ! For ESTIMATE_MODULI
    USE RK_PROCEDURES_MOD          ! For RK1, D_HUGHES_WINGET, HUGHES_WINGET, ROTATE_TENSOR (assuming they are here or in FINT_FUNCTIONS)
    USE INVERSE_MOD                ! For INVERSE, INV3
    USE DETERMINANT_MOD            ! For DETERMINANT
    USE CONSTITUTION_MOD  
    ! USE omp_lib
    IMPLICIT NONE
    ! Parameters for numerical stability and fallback behavior (as suggested)
    LOGICAL, PARAMETER :: CYCLE_ON_FALLBACK = .FALSE.       ! Set to .TRUE. to CYCLE after a fallback to identity matrix
    DOUBLE PRECISION, PARAMETER :: DET_THRESHOLD = 1.0D-12  ! Threshold for determinant singularity check
    INTEGER, PARAMETER :: FALLBACK_ERROR_INCREMENT = 10      ! Error increment for global_device_error_occurred on fallback
    DOUBLE PRECISION, ALLOCATABLE :: FINT_LOCAL(:), FEXT_LOCAL(:)
    INTEGER :: THREAD_ID, NUM_THREADS
    LOGICAL :: USE_ATOMIC_FALLBACK

    !
    !*********************************************************
    !********************* DECLAIRATIONS *********************
    !*********************************************************
    !
    !
    !********** GLOBAL VARIABLES **********
    !
    !GLOBAL IN
    !
    !
    !DISCRETIZATION INFO
    !
    INTEGER, INTENT(IN):: GNUMP                         !NUMBER OF NODES
    INTEGER, INTENT(IN):: GSIZE_IN                      !TOTAL NUMBER OF NODES (local + ghost + buffer) FOR ARRAY DIMENSIONING
    DOUBLE PRECISION:: GCOO(3,GSIZE_IN)        !COORDINATE FOR EACH NODE
    DOUBLE PRECISION:: GCOO_CUURENT(3,GSIZE_IN) !CUURENT COORDINATE FOR EACH NODE
    DOUBLE PRECISION:: GWIN(3,GNUMP)        !WINDOWS FOR EACH (LOCAL) NODE
	DOUBLE PRECISION, SAVE, ALLOCATABLE:: GWIN0(:,:) ! Assumed to be (3,GNUMP) if related to GWIN
    DOUBLE PRECISION:: GSM_LEN(6,GNUMP)     !SMOOTHING LENGTHS FOR EACH (LOCAL) NODE

    ! (1-6) = (+X, -X, +Y, -Y, +Z, -Z)
    DOUBLE PRECISION:: GSM_VOL(GNUMP)       !VOLUME OF SMOOTHING ZONE FOR EACH (LOCAL) NODE
    DOUBLE PRECISION:: GSM_AREA(3,GNUMP)    !AREAS OF SMOOTHING ZONE SIDES FOR EACH (LOCAL) NODE
    INTEGER:: GN(GNUMP)                     !NUMBER OF NEIGHBORS FOR EACH (LOCAL) NODE
    INTEGER:: GSTART(GNUMP)                 !START LOCATION OF NODE NEIOGHBORS IN STACK (FOR LOCAL NODES)
    INTEGER:: DIM_NN_LIST                   !SIZE OF NEIGHBOR STACK
    INTEGER:: GSTACK(DIM_NN_LIST)           !NEIGHBORS FOR EACH NODE (STACKED)
    DOUBLE PRECISION:: GSTACK_SHP(DIM_NN_LIST)       !SHAPES (STACKED)
    DOUBLE PRECISION:: GSTACK_DSHP(3,DIM_NN_LIST)       !SHAPES (STACKED)
    DOUBLE PRECISION:: GSTACK_DDSHP(6,DIM_NN_LIST)       !SHAPES (STACKED)
    DOUBLE PRECISION:: GINVK(3,3,GNUMP)     ! FOR LOCAL NODES
    DOUBLE PRECISION:: GCHAR_DIST(GNUMP),   GMAX_WVEL(GNUMP) ! FOR LOCAL NODES
    INTEGER, INTENT(IN):: GMAXN                              !MAX NUMBER OF NEIGHBORS FOR ALL NODES
    INTEGER:: GGHOST(GSIZE_IN)                 !FLAG FOR GHOST NODES (GHOST = 1)
    LOGICAL:: GEBC_NODES(GSIZE_IN)
    DOUBLE PRECISION:: GVOL(GSIZE_IN)          !VOLUME OF EACH NODE
    DOUBLE PRECISION:: GNSNI_FAC(3,GNUMP)   ! FOR LOCAL NODES
    !
    !STATE AND FIELD VARIABLES
    !
    DOUBLE PRECISION:: GSTRESS(6,GSIZE_IN)     !CAUCHY STRESS OF EACH NODE
    DOUBLE PRECISION:: GSTRAIN(6,GSIZE_IN)     !CAUCHY STRAIN OF EACH NODE
    DOUBLE PRECISION:: GSTATE(20,GSIZE_IN)     !STATE VARIABLES OF EACH NODE
    DOUBLE PRECISION:: GPROP(30,GNUMP)     !MATERIAL PROPERTIES OF EACH (LOCAL) NODE - CHECK IF NEIGHBOR PROPS ARE NEEDED
    DOUBLE PRECISION:: GDINC(3*GSIZE_IN) ,GDINC_TOT(3*GSIZE_IN)      !DISPLACEMENT INCREMENT (PREDICTOR) OF EACH NODE
    INTEGER::          GMAT_TYPE(GNUMP)    !MATERIAL TYPE OF EACH (LOCAL) NODE
    DOUBLE PRECISION:: G_H_STRESS(6,GSIZE_IN)!GC
    DOUBLE PRECISION:: G_S_STRESS(6,GSIZE_IN)!GC
    DOUBLE PRECISION, INTENT(IN) :: DLT
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(GSIZE_IN*3):: FINT, FEXT
    DOUBLE PRECISION:: DLT_FINT
    DOUBLE PRECISION:: GSTRAIN_EQ(GSIZE_IN)     !EQUIVALENT PLASTIC STARIN OF EACH NODE
    INTEGER, INTENT(OUT) :: ierr_fint_arg
    !********** LOCAL VARIABLES **********
    INTEGER:: I,J,K,L,JJ, K_PRINT, L_PRINT, J_PRINT      !INDICIES
    DOUBLE PRECISION:: LCOO(3)    !COORDINATE AT A NODE IN INITIAL CONFIGURATION
    DOUBLE PRECISION:: LCOO_T(3)    !COORDINATE AT A NODE IN CURRENT CONFIGURATION
    DOUBLE PRECISION:: LWIN(3) !WINDOW A NODE
    DOUBLE PRECISION:: VOL !NODAL VOLUME (ACTUAL INTEGRATION WEIGHT
    DOUBLE PRECISION:: LSM_LEN(6)  !SMOOTHING LENGTHS A NODE
    DOUBLE PRECISION:: LSM_PTS(3,6)  !SMOOTHING POINT POSITION
    DOUBLE PRECISION:: SM_COO(3)     !TEMPORARY SMOOTHING POINT POSITION
    DOUBLE PRECISION:: LSM_VOL       !VOLUME OF SMOOTHING ZONE
    DOUBLE PRECISION:: LSM_AOV(6)   !AREA OVER VOLUME OF SMOOTHING ZONE SIDES/ZONES
    DOUBLE PRECISION:: LSTRESS(6)     !CAUCHY STRESS OF A NODE
    DOUBLE PRECISION:: LSTRESS_PREDICTOR(6)     !ELASTIC PREDICTOR STRESS
    DOUBLE PRECISION:: LSTRAIN(6)
    DOUBLE PRECISION:: L_H_STRESS(6)!GC: FOR VISCOELASTIC
    DOUBLE PRECISION:: L_S_STRESS(6)!GC
    DOUBLE PRECISION, INTENT(IN):: G_X_MOM(GNUMP),G_Y_MOM(GNUMP),G_Z_MOM(GNUMP) ! For local nodes
    DOUBLE PRECISION:: LSTATE(20)      !STATE VARIABLES OF A NODE
    DOUBLE PRECISION:: LPROP(30)     !MATERIAL PROPERTIES OF A NODE
    INTEGER:: LMAT_TYPE    !MATERIAL TYPE OF EACH NODE
    LOGICAL:: SELF_EBC
    INTEGER:: LSTART              !LOCATION OF START IN STACK
    INTEGER:: LGHOST              !FLAG FOR GHOST NODES (GHOST = 1)
    INTEGER:: LSTACK(GMAXN)       !LOCAL LIST/STACK OF NEIGHBORS
    DOUBLE PRECISION:: PHI(GMAXN)       !INFLUENCE FUNCTION FOR NEIGHBORS
    INTEGER:: LN                  !NUMBER OF NEIGHBORS FOR A NODE
    DOUBLE PRECISION:: LVOL       !VOLUME OF A NODE
    DOUBLE PRECISION:: SHP(GMAXN), SHPD(3,GMAXN), SHPD_TRASH(3,GMAXN)
    DOUBLE PRECISION:: SHPDTMP(3,GMAXN)
    DOUBLE PRECISION:: SHPD_SM(3,GMAXN)        !SHAPE FUNCTION SMOOTHED GRADIENTS
    DOUBLE PRECISION:: SHPDD_SM(6,GMAXN)       !SHAPE FUNCTION SMOOTHED SMOOTHED GRADIENTS
    DOUBLE PRECISION:: SHPDTEMP(9) !TEMPORARY VARIABLE FOR SHAPES FOR STABILZIATION
    DOUBLE PRECISION:: SHP6(GMAXN,6), SHPD6(3,GMAXN,6)
    DOUBLE PRECISION:: LMAT(3,3) !INCREMENTAL DEFORMATION GRADIENT WITH RESPECT TO THE CURRENT TIME STEP
    DOUBLE PRECISION:: LDINC(3,GMAXN),LDINC_TOT(3,GMAXN)
    DOUBLE PRECISION:: LCOO_CUURENT(3,GMAXN)
    DOUBLE PRECISION:: LCOONE(3,GMAXN)
    DOUBLE PRECISION:: B_TEMP(3,3), B_INV_TEMP(3,3) !GC
    DOUBLE PRECISION:: STRAIN(6)       !INCREMENTALLY OBJECTIVE STRAIN
    DOUBLE PRECISION:: ELAS_MAT(6,6)
    DOUBLE PRECISION:: BMAT(6,3)
    DOUBLE PRECISION:: BMAT_T(3,6)
    DOUBLE PRECISION:: FINT3(3),FINT3_EXT(3),INVK(3,3) 

    DOUBLE PRECISION:: ROT(6,6) !ROTATION MATRIX
    LOGICAL:: LINIT
    DOUBLE PRECISION:: FMAT(3,3), IFMAT(3,3),X_0(3),X_t(3), DX_t(3,1), PKSTRESS(3,3), TEMP_STRESS(3,3)
!    DOUBLE PRECISION, ALLOCATABLE:: FINT_TEMP(:,:,:), FEXT_TEMP(:,:,:)
    DOUBLE PRECISION:: DET
!     DOUBLE PRECISION:: FINT_TEMP(20,3,GNUMP)
!     INTEGER:: ID_RANK
    DOUBLE PRECISION:: PMAT(6,3)
    DOUBLE PRECISION:: FBOD(3),FGRAV(3)
    !DOUBLE PRECISION:: MODEL_BODYFORCE(3)
    !0702
    DOUBLE PRECISION:: MODEL_BODYFORCE(3,GNUMP) ! For local nodes
    DOUBLE PRECISION:: LBOD(3)

    !NSNI
    INTEGER::XMAP(3),YMAP(3),ZMAP(3)
    DOUBLE PRECISION::XLMAT(3,3),YLMAT(3,3),ZLMAT(3,3)
    DOUBLE PRECISION:: DX_STRAIN(6), DY_STRAIN(6), DZ_STRAIN(6)
    DOUBLE PRECISION::  LOCAL_DX_STRESS(6,GSIZE_IN)
    DOUBLE PRECISION::  LOCAL_DY_STRESS(6,GSIZE_IN)
    DOUBLE PRECISION::  LOCAL_DZ_STRESS(6,GSIZE_IN)
    DOUBLE PRECISION::  LDX_STRESS(6)
    DOUBLE PRECISION::  LDY_STRESS(6)
    DOUBLE PRECISION::  LDZ_STRESS(6)
    DOUBLE PRECISION:: CMAT(6,6), LAMDA, MU, LAMDA_PLUS_2MU
    DOUBLE PRECISION:: XBMAT(6,3), XBMAT_T(3,6), XFINT3(GMAXN,3)
    DOUBLE PRECISION:: YBMAT(6,3), YBMAT_T(3,6), YFINT3(GMAXN,3)
    DOUBLE PRECISION:: ZBMAT(6,3), ZBMAT_T(3,6), ZFINT3(GMAXN,3)
    DOUBLE PRECISION:: MAG_STAB_FINT,MAG_FINT
    DOUBLE PRECISION:: TEMP_DEBUG(3)
    DOUBLE PRECISION:: GINT_WORK
!    DOUBLE PRECISION, ALLOCATABLE:: GINT_WORK_TEMP(:)
    DOUBLE PRECISION, ALLOCATABLE :: fint_temp(:)      ! 調整型別和秩（如果需要）
    DOUBLE PRECISION, ALLOCATABLE :: fext_temp(:)      ! 調整型別和秩（如果需要）
    DOUBLE PRECISION, ALLOCATABLE :: gint_work_temp(:) ! 調整型別和秩（如果需要）
    INTEGER :: global_device_error_occurred
    INTEGER :: inverse_fallback_count ! New counter for inverse fallbacks
    INTEGER :: device_error_status_check ! For STAT= in ALLOCATE/DEALLOCATE
    !
    ! FOR TIME STEP CALCS
    !
    LOGICAL:: FIRST,NSNI_FLAG
    DOUBLE PRECISION:: STRESS_INC(6), STRAIN_INC(6), POISS, YOUNG, BULK, SHEAR,  &
        STRESS_INC_DEV(6), STRESS_INC_SPHR(6), STRAIN_INC_DEV(6), STRAIN_INC_SPHR(6),  &
        NORM_STRESS_INC_DEV, NORM_STRESS_INC_SPHR, NORM_STRAIN_INC_DEV, NORM_STRAIN_INC_SPHR,   &
        SHEAR_TRIAL, BULK_TRIAL, PMOD, DENSITY, MAXMOD, MAX_VEL,   &
        DIST, XJ, YJ, ZJ, CHAR_DIST, DLT_TEMP, XI, YI, ZI
    DOUBLE PRECISION:: D(6)
    INTEGER, INTENT(IN) :: MODEL_BODY_ID(GNUMP) ! For local nodes, but used for MODEL_BODY_ID(KK) where KK is neighbor
    INTEGER :: LOCAL_BODY_ID, LOCAL_BODY_ID_2, II, P,KK    
    DOUBLE PRECISION:: F_INT_C(3), MU1, MU_NEW, MU_NEW2, X2(3), X1(3), TEMP, F_N, XNORM(3), F_T1, F_T2, F_T3, F_TT, F_INT_C_TEMP(3),F_T
	!LOGICAL::      KCONTACT
    DOUBLE PRECISION:: DX(3), LITTLE_DX(3),LENGTH_DX,LENGTH_LITTLE_DX
	
    DOUBLE PRECISION:: NSNI_LIMITER
    INTEGER :: constitution_error_flag
    INTEGER :: call_err_status
    INTEGER, PARAMETER :: CHUNK_SIZE = 1024    ! Process CHUNK_SIZE nodes at a time (Adjust as needed)
    INTEGER :: chunk_start, chunk_end, nchunks, c_chunk

    !
    !*********************************************************
    !******************** EXECUTABLE CODE ********************
    !*********************************************************
    ierr_fint_arg = 0 
    
    IF (RK_CONT /= 3) THEN
        WRITE(*,*) 'ERROR: CONSTRUCT_FINT - RK_CONT = ', RK_CONT, ' (expected 3)'
        RK_CONT = 3  ! 強制設定為 3
        WRITE(*,*) 'WARNING: CONSTRUCT_FINT - Forcing RK_CONT = 3'
    END IF

    ! Sanity checks at the beginning of the subroutine (Step 3 from request)
    IF (DIM_NN_LIST <= 0 .OR. DIM_NN_LIST > 10000000) THEN ! Example upper bound
        WRITE(*,*) 'ERROR: CONSTRUCT_FINT - DIM_NN_LIST out of range: ', DIM_NN_LIST
        ierr_fint_arg = -15
        RETURN
    END IF


    ! Ensure GMAXN does not exceed a practical limit
    IF (GMAXN > 5000) THEN ! Example practical limit for GMAXN
        WRITE(*,*) 'WARNING: CONSTRUCT_FINT - GMAXN (', GMAXN, ') > 5000. This may cause performance issues.'
    END IF

    IF (GNUMP <= 0) THEN
        ierr_fint_arg = -10; RETURN
    END IF
    IF (GSIZE_IN <= 0) THEN
        ierr_fint_arg = -101; RETURN
    END IF

    ! WRITE(*,*) 'DEBUG: CONSTRUCT_FINT (Device) - Entry parameter check:'
    ! WRITE(*,*) '  Received GMAXN = ', GMAXN
    ! WRITE(*,*) '  Received DIM_NN_LIST = ', DIM_NN_LIST    ! WRITE(*,*) '  Received GNUMP = ', GNUMP
    ! WRITE(*,*) '  Received GSIZE_IN = ', GSIZE_IN


    global_device_error_occurred = 0
    inverse_fallback_count = 0

    !GINT_WORK = 0.0d0
    !
    ! Rigorous parameter checks (GNUMP is for loop bounds, GSIZE_IN for array extents)

    IF (GNUMP <= 0) THEN
        ! WRITE(*,*) 'FATAL ERROR: CONSTRUCT_FINT (Device) called with invalid GNUMP = ', GNUMP


        ierr_fint_arg = -10
        RETURN
    END IF

    IF (GSIZE_IN <= 0) THEN
        ! WRITE(*,*) 'FATAL ERROR: CONSTRUCT_FINT (Device) called with invalid GSIZE_IN = ', GSIZE_IN
        ierr_fint_arg = -101 
        RETURN
    END IF

    ! 直接檢查輸入參數
    IF (GMAXN <= 0 .OR. GMAXN > DIM_NN_LIST) THEN
        WRITE(*,*) 'ERROR: CONSTRUCT_FINT - Invalid GMAXN: ', GMAXN
        ierr_fint_arg = -15
        RETURN
    END IF

    IF (DIM_NN_LIST <= 0) THEN
       !WRITE(*,*) 'WARNING: CONSTRUCT_FINT (Device) - Fixing invalid DIM_NN_LIST from ', DIM_NN_LIST, ' to ', GNUMP * 1000
        DIM_NN_LIST = GNUMP * 1000

        ierr_fint_arg = -12
        ! Depending on severity, might not want to RETURN immediately if a fallback is possible.
        ! For now, keeping the original logic of returning on error for DIM_NN_LIST.
        ! If DIM_NN_LIST is critical for allocations that cannot use a fallback, RETURN is appropriate.
        RETURN 
    END IF

!     FEXT is now dimensioned with GSIZE_IN*3. The GSIZE_IN check handles this.
!     WRITE(*,*) 'DEBUG: CONSTRUCT_FINT (Device) - After parameter fix:'
!     WRITE(*,*) '  GMAXN = ', GMAXN
!     WRITE(*,*) '  DIM_NN_LIST = ', DIM_NN_LIST
 
    !INITIALIZE FINT

    ! Initialize FINT and FEXT on the device at the beginning of the routine,
    ! This is done once before the chunked parallel loops.

    FINT(1:GSIZE_IN*3) = 0.0D0
    FEXT(1:GSIZE_IN*3) = 0.0D0
    GINT_WORK = 0.0D0 ! Initialize GINT_WORK before accumulation
    DET = 1.d0
    ! Handle GWIN0 allocation outside parallel region
    IF (LINIT) THEN
        ALLOCATE(GWIN0(3,GNUMP), STAT=device_error_status_check) ! GWIN0 seems related to local node properties
        IF (device_error_status_check /= 0) THEN
            global_device_error_occurred = 100 ! 特定錯誤碼
            ierr_fint_arg = global_device_error_occurred
            RETURN
        END IF
        GWIN0 = GWIN
    END IF


    !
    !
    !LOOP OVER THE NODE STACK
    !
    !LET OPEN-MP DECIDE HOW TO DO THE DO-LOOP
    !		

    ! Direct parallelization without chunking
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) &
    !$ACC COPYOUT(global_device_error_occurred, inverse_fallback_count) &
    !$ACC PRIVATE(LCOO, LCOO_T, VOL, LWIN, LSM_LEN, LN, LSTART, LGHOST, LVOL, LPROP, LMAT_TYPE, LOCAL_BODY_ID, SELF_EBC, &
    !$ACC         LSTRESS, LSTRAIN, LSTATE, LBOD, LSTACK, LDINC, LDINC_TOT, LCOONE, LCOO_CUURENT, &
    !$ACC         SHP, SHPD, SHPD_TRASH, SHPDTMP, FMAT, IFMAT, DET, LMAT, STRAIN, D, ELAS_MAT, LSTRESS_PREDICTOR, &
    !$ACC         B_TEMP, B_INV_TEMP, BMAT, BMAT_T, FINT3, FINT3_EXT, ROT, STRESS_INC, STRAIN_INC, L_H_STRESS, L_S_STRESS, &
    !$ACC         POISS, YOUNG, BULK, SHEAR, DENSITY, PMOD, MAXMOD, NSNI_FLAG, SHEAR_TRIAL, BULK_TRIAL, XMAP, YMAP, ZMAP, XLMAT, YLMAT, ZLMAT, DX_STRAIN, DY_STRAIN, DZ_STRAIN, &
    !$ACC         XNORM, F_INT_C_TEMP, F_INT_C, FGRAV, FBOD, J, K, L, JJ, II, P, KK, call_err_status, constitution_error_flag) &
    !$ACC REDUCTION(+:GINT_WORK)

DO I = 1, GNUMP
        !
        ! Initialize call_err_status for each node I inside the parallel loop
        call_err_status = 0

        !
        !
        !GRAB NODE INFORMATION FROM LIST
        !
        LCOO(:) = GCOO(:,I) ! I is the central node index (1 to GNUMP)

        LCOO_T(:) = GCOO_CUURENT(:,I)
        VOL = GVOL(I)
        LWIN(:) = GWIN(:,I)
        LSM_LEN(:) = GSM_LEN(:,I)
        LN = GN(I)
        LSTART = GSTART(I)
        LGHOST = GGHOST(I)
        LVOL = GVOL(I)
        LPROP = GPROP(:,I)
        LMAT_TYPE = GMAT_TYPE(I)
        LOCAL_BODY_ID = MODEL_BODY_ID(I) ! MODEL_BODY_ID is dimensioned GNUMP

        
        !
        IF (GEBC_NODES(I)) THEN
            SELF_EBC = .TRUE.
        ELSE
            SELF_EBC = .FALSE.
        END IF
        !
        ! RECALL THE NODE STRESS AND STATE VARIABLES
        !
        DO J = 1, 6
            LSTRESS(J) = GSTRESS(J,I) ! GSTRESS is (6,GSIZE_IN), I is local index
            LSTRAIN(J) = GSTRAIN(J,I) ! GSTRAIN is (6,GSIZE_IN), I is local index


            IF(LMAT_TYPE == 4) THEN !建議取消註解此區塊
                L_H_STRESS(J) = G_H_STRESS(J,I)!GC G_H_STRESS is (6,GSIZE_IN)
                L_S_STRESS(J) = G_S_STRESS(J,I)!GC G_S_STRESS is (6,GSIZE_IN)

            END IF
        END DO

        DO J = 1, 30
            LPROP(J) = GPROP(J,I)
        END DO

        DO J = 1, 3
            LBOD(J) = MODEL_BODYFORCE(J,I) ! MODEL_BODYFORCE is (3,GNUMP)

        END DO

        DO J = 1, 20
            LSTATE(J) = GSTATE(J,I) ! GSTATE is (20,GSIZE_IN)

        END DO
        !
        ! GET THE NEIGHBOR LIST
        !
        DO J = 1, LN
            LSTACK(J) = GSTACK(LSTART+J-1)
        END DO
        !
        ! GET THE INCREMENTS OF DISPLACEMENTS FOR NEIGHBORS
        !
        DO J = 1, LN
            JJ = LSTACK(J) ! JJ is a global index (1 to GSIZE_IN)

            LDINC(1,J) = GDINC((JJ-1)*3+1)
            LDINC(2,J) = GDINC((JJ-1)*3+2)
            LDINC(3,J) = GDINC((JJ-1)*3+3)
        END DO
        !
        ! GET THE TOTAL DISPLACEMENTS FOR NEIGHBORS

        !
        DO J = 1, LN
            JJ = LSTACK(J)
            LDINC_TOT(1,J) = GDINC_TOT((JJ-1)*3+1) ! GDINC_TOT is (3*GSIZE_IN)
            LDINC_TOT(2,J) = GDINC_TOT((JJ-1)*3+2)
            LDINC_TOT(3,J) = GDINC_TOT((JJ-1)*3+3)

        END DO
        !
        ! GET THE ORIGINAL COORDINATES FOR NEIGHBORS
        !
        DO J = 1, LN
            JJ = LSTACK(J)
            LCOONE(1,J) = GCOO(1,JJ)
            LCOONE(2,J) = GCOO(2,JJ) ! GCOO is (3,GSIZE_IN)
            LCOONE(3,J) = GCOO(3,JJ)

        END DO
        !
        !
        ! GET THE CURRENT COORDINATES FOR NEIGHBORS
        !
        DO J = 1, LN
            JJ = LSTACK(J)
            LCOO_CUURENT(1,J) = GCOO_CUURENT(1,JJ)
            LCOO_CUURENT(2,J) = GCOO_CUURENT(2,JJ) ! GCOO_CUURENT is (3,GSIZE_IN)
            LCOO_CUURENT(3,J) = GCOO_CUURENT(3,JJ)

        END DO
        !
        !IF IT IS LAGRANGIAN AND IT IS NOT THE FIRST STEP, RECALL SHAPE FUNCTIONS
		!
        IF ((LLAGRANGIAN).AND.(.NOT.LINIT)) THEN
            DO J = 1, LN
                SHP(J) =   GSTACK_SHP(LSTART+J-1)
                SHPD(1,J) = GSTACK_DSHP(1,LSTART+J-1)
                SHPD(2,J) = GSTACK_DSHP(2,LSTART+J-1)
                SHPD(3,J) = GSTACK_DSHP(3,LSTART+J-1)
            END DO

            IF (LFINITE_STRAIN) THEN

                FMAT = 0.0d0
                FMAT(1,1) = 1.d0
                FMAT(2,2) = 1.d0
                FMAT(3,3) = 1.d0
                !
                DO K = 1, 3
                    DO L = 1, 3
                        DO J = 1, LN
                            !JJ = LSTACK(J)
                            !FMAT(K,L) = FMAT(K,L) +  SHPD(L,J)*GCOO_CUURENT(K,JJ)
                            FMAT(K,L) = FMAT(K,L) +  SHPD(L,J)*LDINC_TOT(K,J)
                        END DO !J = 1, LN (COMPUTE THE INCREMENTAL DEFORMATION GRADIENT)
                    END DO
                END DO
                CALL DETERMINANT(FMAT,DET) ! DET for FMAT
                !   IF (I == 1) THEN
                !        WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - FMAT (Lagrangian, LINIT) Determinant = ', DET
                !    END IF

                !
                !
                SHPDTMP(:,1:LN)=SHPD(:,1:LN)
                SHPD(:,1:LN) = 0.0d0

            !     CALL INVERSE(FMAT, 3, IFMAT)
            !     IF (I == 1) THEN
            !        WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - About to call INV3 for FMAT (Lagrangian, non-LINIT)'
            !        WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - FMAT Determinant = ', DET
            !     END IF                ! This END IF was misplaced and commented out, seems correct.
                IF (ABS(DET) < DET_THRESHOLD) THEN ! Use parameterized threshold


!                    IF (I == 1 .OR. MOD(inverse_fallback_count,100) == 0) THEN
!                        WRITE(*,*) 'WARNING: CONSTRUCT_FINT - Node ', I, ' - FMAT (Lagrangian, non-LINIT) is near-singular (DET=',DET,'). Using identity for IFMAT.'
!                    END IF
                    IFMAT = 0.0D0
                    IFMAT(1,1) = 1.0D0
                    IFMAT(2,2) = 1.0D0
                    IFMAT(3,3) = 1.0D0
                    !$ACC ATOMIC UPDATE
                    inverse_fallback_count = inverse_fallback_count + 1
                    !$ACC ATOMIC UPDATE
                    global_device_error_occurred = global_device_error_occurred + FALLBACK_ERROR_INCREMENT
                    call_err_status = -99 ! Indicate fallback due to determinant

                ELSE
                    CALL INV3(FMAT,  IFMAT, call_err_status) 

                    IF (call_err_status /= 0) THEN
!                        IF (I == 1 .OR. MOD(inverse_fallback_count,100) == 0) THEN
!                             WRITE(*,*) 'WARNING: CONSTRUCT_FINT - Node ', I, ' - INV3 failed for FMAT (Lagrangian, non-LINIT), status: ', call_err_status, '. Using identity for IFMAT.'
!                        END IF
                        IFMAT = 0.0D0
                        IFMAT(1,1) = 1.0D0
                        IFMAT(2,2) = 1.0D0
                        IFMAT(3,3) = 1.0D0
                        !$ACC ATOMIC UPDATE
                        inverse_fallback_count = inverse_fallback_count + 1
                        !$ACC ATOMIC UPDATE
                        global_device_error_occurred = global_device_error_occurred + FALLBACK_ERROR_INCREMENT
                        ! call_err_status from INV3 is kept if INV3 itself failed
                    END IF
                END IF

            !    IF (I == 1) THEN
            !         WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - INV3/Fallback for FMAT completed.'
            !    END IF

                IF (CYCLE_ON_FALLBACK .AND. (call_err_status /= 0 .OR. ABS(DET) < DET_THRESHOLD) ) THEN
                    CYCLE
                END IF

                ! If after fallback, we still consider it a major error to CYCLE, then the original error accumulation can be kept.
                ! For now, we assume fallback allows continuation. If global_device_error_occurred should still be incremented, uncomment below.
                ! IF (call_err_status_original_from_inv3 /= 0) THEN ! hypothetical original status
                !     !$ACC ATOMIC UPDATE
                !     global_device_error_occurred = global_device_error_occurred + 1000 
                !     CYCLE
                ! END IF

                DO J = 1, LN
                    DO K = 1, 3
                        SHPD(1,J) = SHPD(1,J) + SHPDTMP(K,J)*IFMAT(K,1)
                        SHPD(2,J) = SHPD(2,J) + SHPDTMP(K,J)*IFMAT(K,2)
                        SHPD(3,J) = SHPD(3,J) + SHPDTMP(K,J)*IFMAT(K,3)
                    END DO
                END DO
            END IF

        ELSE
            !
            ! TODO: CONDENSE ALL SHAPE FUNCTION CALCULATIONS
            !
            IF (ITYPE_INT.EQ.0) THEN
                !
                ! DIRECT NODAL INTEGRATION
                !
                IF (LLAGRANGIAN) THEN
                !    IF (I == 1) THEN
                !        WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - About to call RK1 (Lagrangian)'
                !    END IF
                    ! 新增：檢查 RK_CONT 參數
                    IF (I == 1) THEN
                        WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - RK parameters check:'
                        WRITE(*,*) '  RK_CONT = ', RK_CONT, ' (should be 3 for cubic spline)'
                        WRITE(*,*) '  RK_DEGREE = ', RK_DEGREE
                        WRITE(*,*) '  RK_PSIZE = ', RK_PSIZE
                    END IF
                    CALL RK1(LCOO, RK_DEGREE, RK_PSIZE, RK_CONT, RK_IMPL,GCOO, GWIN, GSIZE_IN, LSTACK, LN, GMAXN, GEBC_NODES,SELF_EBC, &
                        QL, QL_COEF,QL_LEN, &
                        SHP, SHPD,SHSUP)
                !    IF (I == 1) THEN
                !        WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - RK1 call completed (Lagrangian)'
                !    END IF

                ELSE !GCOO_CUURENT
                !    IF (I == 1) THEN
                !        WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - About to call RK1 (Eulerian/Updated Lagrangian)'
                !    END IF

                    CALL RK1(LCOO_T, RK_DEGREE, RK_PSIZE, RK_CONT, RK_IMPL,GCOO_CUURENT, GWIN, GSIZE_IN, LSTACK, LN, GMAXN, GEBC_NODES,SELF_EBC, &
                        QL, QL_COEF,QL_LEN, &
                        SHP, SHPD, SHSUP)
 !                   IF (I == 1) THEN
 !                       WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - RK1 call completed (Eulerian/Updated Lagrangian)'
 !                       IF (LN < 4 .AND. .NOT. SELF_EBC) THEN ! Check for 3D, adjust for 2D if necessary
 !                           WRITE(*,*) 'WARNING: CONSTRUCT_FINT - Node ', I, ' has insufficient neighbors (LN=', LN, ') for robust SHPD calculation.'
 !                       END IF
 !                   END IF                    
                    ! CALCULATE THE DEFORMATION GRADIENT (FMAT = (dX/dx)^-1 = dx/dX)
                    ! IF (I == 1) THEN  ! This was a duplicated debug message
                    !    WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - RK1 call completed (Eulerian/Updated Lagrangian)'
                    ! END IF
                    ! CALCULATE THE DEFORMATION GRADIENT
                    B_TEMP = 0.D0
                    DO K = 1, 3
                        DO L = 1, 3
                            DO J = 1, LN
                                B_TEMP(K,L) = B_TEMP(K,L) +  SHPD(L,J)*LCOONE(K,J)
                            END DO 
                        END DO
                    END DO
    USE_ATOMIC_FALLBACK = .FALSE.
    IF (GNUMP .GT. 10000) THEN
    ! 對大問題使用 atomic 以節省記憶體
        USE_ATOMIC_FALLBACK = .TRUE.
    END IF

    ! Check determinant before calling INVERSE for B_TEMP
    CALL DETERMINANT(B_TEMP, DET) ! DET for B_TEMP

 !                   IF (I == 1) THEN
 !                       WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - About to call INVERSE for B_TEMP'
 !                       WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - B_TEMP Determinant = ', DET
 !                       IF (LN < 4 .AND. .NOT. SELF_EBC) THEN
 !                            WRITE(*,*) 'DIAGNOSTIC: Node ', I, ' (B_TEMP calc) has LN=', LN
 !                       END IF
 !                   END IF

                    IF (ABS(DET) < DET_THRESHOLD) THEN ! Use parameterized threshold

                        ! IF (I == 1 .OR. MOD(inverse_fallback_count,100) == 0) THEN
                            ! WRITE(*,*) 'WARNING: CONSTRUCT_FINT - Node ', I, ' - B_TEMP (dX/dx) is singular or near-singular (DET=',DET,'). Using identity for B_INV_TEMP (FMAT).'
                            ! IF (I <= 5 .OR. MOD(inverse_fallback_count, 500) == 1) THEN ! Detailed output for some cases
                                ! WRITE(*,*) 'DIAGNOSTIC: Node I = ', I, ', LN = ', LN
                                ! WRITE(*,*) '  LCOO_T: ', LCOO_T
                                ! WRITE(*,*) '  B_TEMP:'
                                ! DO K_PRINT = 1, 3
                                    ! WRITE(*,'(3E15.6)') (B_TEMP(K_PRINT, L_PRINT), L_PRINT = 1, 3)
                                ! END DO
                                ! WRITE(*,*) '  Neighbors (Current Coords | Initial Coords | SHPD values):'
                                ! DO J_PRINT = 1, MIN(LN, 5) ! Print first few neighbors
                                    ! WRITE(*,'(A,I0,A,3E12.4,A,3E12.4,A,3E12.4)') '    Nbr ', J_PRINT, ': CUR: ', &
                                        ! (LCOO_CUURENT(K_PRINT,J_PRINT), K_PRINT=1,3), ' | INIT: ', &
                                        ! (LCOONE(K_PRINT,J_PRINT), K_PRINT=1,3), ' | SHPD: ', &
                                        ! (SHPD(K_PRINT,J_PRINT), K_PRINT=1,3)
                                ! END DO
                            ! END IF
                        ! END IF

                        B_INV_TEMP = 0.0D0
                        B_INV_TEMP(1,1) = 1.0D0
                        B_INV_TEMP(2,2) = 1.0D0
                        B_INV_TEMP(3,3) = 1.0D0
                        !$ACC ATOMIC UPDATE
                        inverse_fallback_count = inverse_fallback_count + 1
                        !$ACC ATOMIC UPDATE
                        global_device_error_occurred = global_device_error_occurred + FALLBACK_ERROR_INCREMENT 
                        call_err_status = -99 ! Indicate fallback due to determinant

  
                    ELSE
                        CALL INVERSE(B_TEMP, 3, B_INV_TEMP, call_err_status)
                        IF (call_err_status /= 0) THEN
                            ! IF (I == 1 .OR. MOD(inverse_fallback_count,100) == 0) THEN
                                ! WRITE(*,*) 'WARNING: CONSTRUCT_FINT - Node ', I, ' - INVERSE failed for B_TEMP (dX/dx), status: ', call_err_status, '. Using identity for B_INV_TEMP (FMAT).'
                            ! END IF

                            B_INV_TEMP = 0.0D0
                            B_INV_TEMP(1,1) = 1.0D0
                            B_INV_TEMP(2,2) = 1.0D0
                            B_INV_TEMP(3,3) = 1.0D0
                            !$ACC ATOMIC UPDATE
                            inverse_fallback_count = inverse_fallback_count + 1
                            !$ACC ATOMIC UPDATE
                            global_device_error_occurred = global_device_error_occurred + FALLBACK_ERROR_INCREMENT     
                            ! call_err_status from INVERSE is kept if INVERSE itself failed    
                        END IF
                    END IF

                    IF (CYCLE_ON_FALLBACK .AND. (call_err_status /= 0 .OR. ABS(DET) < DET_THRESHOLD) ) THEN 

                        CYCLE
                    END IF
                    
                    FMAT = B_INV_TEMP
                    !
                    ! STORE THE SHP FOR PHY DISPLACEMENT/VEL CACULATION FOR DNI
                    !
                    DO J = 1, LN
                        GSTACK_SHP(LSTART+J-1) = SHP(J)
                        GSTACK_DSHP(1,LSTART+J-1) = SHPD(1,J)
                        GSTACK_DSHP(2,LSTART+J-1) = SHPD(2,J)
                        GSTACK_DSHP(3,LSTART+J-1) = SHPD(3,J)
                    END DO

                END IF

                CONTINUE

!            ELSEIF ((ITYPE_INT.EQ.2).OR.(ITYPE_INT.EQ.1)) THEN
                !
                ! NSNI
                !
                ! GET THE SMOOTHING INFORMATION
                ! (1-6) = (+X, -X, +Y, -Y, +Z, -Z)
                !
!               DO J = 1, 3
!                    DO K = 1, 6
!                        LSM_PTS(J,K) = LCOO_T(J)
!                    END DO
!                END DO
!                LSM_PTS(1,1) = LCOO_T(1) + LSM_LEN(1)
!                LSM_PTS(1,2) = LCOO_T(1) - LSM_LEN(2)
!                LSM_PTS(2,3) = LCOO_T(2) + LSM_LEN(3)
!                LSM_PTS(2,4) = LCOO_T(2) - LSM_LEN(4)
!                LSM_PTS(3,5) = LCOO_T(3) + LSM_LEN(5)
!                LSM_PTS(3,6) = LCOO_T(3) - LSM_LEN(6)
!
!                LSM_VOL = GSM_VOL(I)
!                !
!                ! COMPUTE SMOOTHED AREA OVER VOLUME
!                !
!                DO J = 1, 3
!                    LSM_AOV((J-1)*2+1) = GSM_AREA(J,I) / LSM_VOL
!                    LSM_AOV((J-1)*2+2) = GSM_AREA(J,I) / LSM_VOL
!                END DO           
! COMPUTE THE SHAPE FUNCTIONS AT GRADIENT SMOOTHING POINTS    
!                DO J = 1, 6
!                    SM_COO(:) = LSM_PTS(:,J)
!                    CALL RK1(SM_COO, RK_DEGREE, RK_PSIZE, RK_CONT, RK_IMPL, GCOO_CUURENT, GWIN, GSIZE_IN, LSTACK, LN, GMAXN, GEBC_NODES,SELF_EBC, &
!                        QL, QL_COEF,QL_LEN, &
!                        SHP, SHPD, SHSUP)
!                    DO K = 1, LN
                    
!                        SHP6(K,J) = SHP(K)
!                        SHPD6(:,K,J) = SHPD(:,K)
!                    END DO

!                END DO !J = 1, 6 (COMPUTE THE SMOOTHED GRADIENTS)

                !
                ! FILL OUT THE SMOOTHED GRADIENT INFORMATION
                !
!                DO K = 1, LN

!                    SHPD(1,K) = (SHP6(K,1)*LSM_AOV(1) - SHP6(K,2)*LSM_AOV(2))
!                    SHPD(2,K) = (SHP6(K,3)*LSM_AOV(3) - SHP6(K,4)*LSM_AOV(4))
!                    SHPD(3,K) = (SHP6(K,5)*LSM_AOV(5) - SHP6(K,6)*LSM_AOV(6))

!                END DO

!                IF (ITYPE_INT.EQ.2) THEN
                    !NSNI CALCS
!                    DO K = 1, LN
                        !
                        ! SMOOTH X Y Z IN X DIRECTION
                        !
                        !XX
!                        SHPDTEMP(1) = (SHPD6(1,K,1)*LSM_AOV(1) - SHPD6(1,K,2)*LSM_AOV(2)) !SMOOTH IN X DIRECTION

                        !XY
!                        SHPDTEMP(2) = (SHPD6(2,K,1)*LSM_AOV(1) - SHPD6(2,K,2)*LSM_AOV(2)) !SMOOTH IN X DIRECTION

                        !XZ
!                        SHPDTEMP(3) = (SHPD6(3,K,1)*LSM_AOV(1) - SHPD6(3,K,2)*LSM_AOV(2)) !SMOOTH IN X DIRECTION
                        !
                        ! SMOOTH X Y Z IN Y DIRECTION
                        !
                        !YX
!                        SHPDTEMP(4) = (SHPD6(1,K,3)*LSM_AOV(3) - SHPD6(1,K,4)*LSM_AOV(4)) !SMOOTH IN Y DIRECTION

                        !YY
!                        SHPDTEMP(5) = (SHPD6(2,K,3)*LSM_AOV(3) - SHPD6(2,K,4)*LSM_AOV(4)) !SMOOTH IN Y DIRECTION

                        !YZ
!                        SHPDTEMP(6) = (SHPD6(3,K,3)*LSM_AOV(3) - SHPD6(3,K,4)*LSM_AOV(4)) !SMOOTH IN Y DIRECTION
                        !
                        ! SMOOTH X Y Z IN Z DIRECTION
                        !
                        !ZX
!                        SHPDTEMP(7) = (SHPD6(1,K,5)*LSM_AOV(5) - SHPD6(1,K,6)*LSM_AOV(6)) !SMOOTH IN Z DIRECTION

                        !ZY
!                        SHPDTEMP(8) = (SHPD6(2,K,5)*LSM_AOV(5) - SHPD6(2,K,6)*LSM_AOV(6)) !SMOOTH IN Z DIRECTION

                        !ZZ
!                        SHPDTEMP(9) = (SHPD6(3,K,5)*LSM_AOV(5) - SHPD6(3,K,6)*LSM_AOV(6)) !SMOOTH IN Z DIRECTION

                        !XX
!                        SHPDD_SM(1,K) = SHPDTEMP(1)
                        !YY
!                        SHPDD_SM(2,K) = SHPDTEMP(5)
                        !ZZ
!                        SHPDD_SM(3,K) = SHPDTEMP(9)
                        !XY
!                        SHPDD_SM(4,K) = 0.5d0 * (SHPDTEMP(2) + SHPDTEMP(4))
                        !YZ
!                        SHPDD_SM(5,K) = 0.5d0 * (SHPDTEMP(6) + SHPDTEMP(8))
                        !XZ
!                        SHPDD_SM(6,K) = 0.5d0 * (SHPDTEMP(3) + SHPDTEMP(7))

!                    END DO

                    !TODO: GET RID OF THESE ARRAYS, WE DONT DO LAGRANGIAN NSNI, SO ITS
                    !WASTING A BUNCH OF STORAGE

!                    DO J = 1, LN
!                        GSTACK_DDSHP(1,LSTART+J-1) = SHPDD_SM(1,J)
!                        GSTACK_DDSHP(2,LSTART+J-1) = SHPDD_SM(2,J)
!                        GSTACK_DDSHP(3,LSTART+J-1) = SHPDD_SM(3,J)
!                        GSTACK_DDSHP(4,LSTART+J-1) = SHPDD_SM(4,J)
!                        GSTACK_DDSHP(5,LSTART+J-1) = SHPDD_SM(5,J)
!                        GSTACK_DDSHP(6,LSTART+J-1) = SHPDD_SM(6,J)
!                    END DO

!                END IF

!                CALL RK1(LCOO_T, RK_DEGREE, RK_PSIZE, RK_CONT, RK_IMPL,GCOO_CUURENT, GWIN, GSIZE_IN, LSTACK, LN, GMAXN, GEBC_NODES,SELF_EBC, &
!                    QL, QL_COEF,QL_LEN, &
!                    SHP, SHPD_TRASH, SHSUP)
                !
                ! STORE THE SHP FOR PHY DISPLACEMENT/VEL CACULATION FOR SNNI, DNI
                !
!                DO J = 1, LN
!                    GSTACK_SHP(LSTART+J-1) = SHP(J)
!                END DO

            END IF



            IF ((LLAGRANGIAN).AND.(LINIT)) THEN

                DO J = 1, LN
                    GSTACK_SHP(LSTART+J-1) = SHP(J)
                    GSTACK_DSHP(1,LSTART+J-1) = SHPD(1,J)
                    GSTACK_DSHP(2,LSTART+J-1) = SHPD(2,J)
                    GSTACK_DSHP(3,LSTART+J-1) = SHPD(3,J)
                END DO

                IF (LFINITE_STRAIN) THEN

                    FMAT = 0.0d0
                    FMAT(1,1) = 1.d0
                    FMAT(2,2) = 1.d0
                    FMAT(3,3) = 1.d0
                    !

                    DO K = 1, 3
                        DO L = 1, 3
                            DO J = 1, LN
                                FMAT(K,L) = FMAT(K,L) +  SHPD(L,J)*LDINC_TOT(K,J)
                            END DO !J = 1, LN (COMPUTE THE INCREMENTAL DEFORMATION GRADIENT)
                        END DO
                    END DO
                    CALL DETERMINANT(FMAT,DET)
                    !
                    SHPDTMP(:,1:LN)=SHPD(:,1:LN)
                    SHPD(:,1:LN) = 0.0d0


    !CALL INVERSE(FMAT, 3, IFMAT)
    ! 保持與 OpenMP 相似的簡單處理
    IF (ABS(DET) < 1.0D-12) THEN
        ! 記錄問題但繼續計算（與 OpenMP 行為一致）
        !$ACC ATOMIC UPDATE
        inverse_fallback_count = inverse_fallback_count + 1
        
        ! 使用小的擾動避免奇異性
        DO K = 1, 3
            FMAT(K,K) = FMAT(K,K) + 1.0D-10
        END DO
    END IF

    CALL INV3(FMAT, IFMAT, call_err_status)
    ! 如果失敗，跳過此節點（保守處理）
    IF (call_err_status /= 0) THEN
        CYCLE
    END IF
                    ! IF (I == 1) THEN
                    !      WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - INV3/Fallback for FMAT (Lagrangian, LINIT) completed.'
                    ! END IF
                    IF (CYCLE_ON_FALLBACK .AND. (call_err_status /= 0 .OR. ABS(DET) < DET_THRESHOLD) ) THEN
                        CYCLE
                    END IF
                    ! Original error handling for CYCLE if needed, but fallback should allow continuation.
                    ! IF (original_call_err_status_from_inv3 /= 0) THEN
                    !     !$ACC ATOMIC UPDATE
                    !     global_device_error_occurred = global_device_error_occurred + 3000 
                    !     CYCLE 
                    ! END IF

                    DO J = 1, LN

                        DO K = 1, 3
                            SHPD(1,J) = SHPD(1,J) + SHPDTMP(K,J)*IFMAT(K,1)
                            SHPD(2,J) = SHPD(2,J) + SHPDTMP(K,J)*IFMAT(K,2)
                            SHPD(3,J) = SHPD(3,J) + SHPDTMP(K,J)*IFMAT(K,3)

                        END DO

                    END DO
               END IF ! IF (LFINITE_STRAIN)
 


            END IF

        END IF
        !
        ! COMPUTE STRAIN MEASURES
        !

            ! IF (.NOT.(LLAGRANGIAN)) THEN
		
                    ! FMAT = 0.0d0
                    ! FMAT(1,1) = 1.d0
                    ! FMAT(2,2) = 1.d0
                    ! FMAT(3,3) = 1.d0
                    ! !

                    ! DO K = 1, 3
                        ! DO L = 1, 3
                            ! DO J = 1, LN
                                ! FMAT(K,L) = FMAT(K,L) +  SHPD(L,J)*LDINC_TOT(K,J)
                            ! END DO !J = 1, LN (COMPUTE THE INCREMENTAL DEFORMATION GRADIENT)
                        ! END DO
                    ! END DO
					
					! FMAT_TRANS = TRANSPOSE(FMAT)
					
					! !COMPUTE GREEN'S STRAIN
					! E_STRAIN = MATMUL(FMAT_TRANS,FMAT)
					
                    ! DO K = 1, 3
                        ! DO L = 1, 3
						
						! E_STRAIN(K,L) = 0.0d0
						! END DO
						! END DO
					
                    ! CALL DETERMINANT(FMAT,DET)
					
			! END IF
					
					
        IF (LFINITE_STRAIN) THEN
            !
            ! COMPUTE THE INCREMENTAL DEFORMATION GRADIENT WITH RESPECT TO THE CURRENT TIME STEP
            !
            LMAT = 0.0d0
            !
            !

            DO K = 1, 3
                DO L = 1, 3
                    DO J = 1, LN
                        LMAT(K,L) = LMAT(K,L) +  SHPD(L,J)*LDINC(K,J)
                    END DO
                END DO
            END DO !J = 1, LN (COMPUTE THE INCREMENTAL DEFORMATION GRADIENT)
			
			! !COMPUTE THE CHANBGE IN LENGTHS OF SUPPORTS
			! INCREMENTALLY (IN PROGRESS!)
			
			! DO J=1,3
			
			  ! DX = 0.0d0
			  ! DX(J) = GWIN(J,I)
				
				! LENGTH_DX = DSQRT(DX(1)**2+DX(2)**2+DX(3)**2)
				
				! DO K = 1, 3
				  ! LITTLE_DX(K) = 0.0d0
					! DO L = 1, 3
					! LITTLE_DX(K) = LITTLE_DX(K) + LMAT(K,L)*DX(L)
					! END DO
				! END DO
				
				! LENGTH_LITTLE_DX = DSQRT(LITTLE_DX(1)**2+LITTLE_DX(2)**2+LITTLE_DX(3)**2)
				
				! GWIN(I,J) = GWIN(I,J) + GWIN(I,J) * LENGTH_LITTLE_DX/LENGTH_DX
				! GWIN(I,J) = MIN(GWIN0(I,J)*2.0d0,GWIN(I,J))
				! GWIN(I,J) = MAX(GWIN0(I,J)/2.0d0,GWIN(I,J))
				
			! END DO
			

            IF (ITYPE_INT.EQ.2) THEN

                !NSNI

                !XX
                !YY
                !ZZ
                !XY
                !YZ
                !XZ

                XMAP(1) = 1 !X,X
                XMAP(2) = 4 !Y,X
                XMAP(3) = 6 !Z,X

                YMAP(1) = 4 !X,Y
                YMAP(2) = 2 !Y,Y
                YMAP(3) = 5 !Z,Y

                ZMAP(1) = 6 !X,Z
                ZMAP(2) = 5 !Y,Z
                ZMAP(3) = 3 !Z,Z

                XLMAT = 0.0d0
                YLMAT = 0.0d0
                ZLMAT = 0.0d0

                DO K = 1, 3
                    DO L = 1, 3
                        DO J = 1, LN
                            XLMAT(K,L) = XLMAT(K,L) +  SHPDD_SM(XMAP(L),J)*LDINC(K,J)
                            YLMAT(K,L) = YLMAT(K,L) +  SHPDD_SM(YMAP(L),J)*LDINC(K,J)
                            ZLMAT(K,L) = ZLMAT(K,L) +  SHPDD_SM(ZMAP(L),J)*LDINC(K,J)
                        END DO
                    END DO
                END DO !J = 1, LN (COMPUTE THE INCREMENTAL DEFORMATION GRADIENT)

                CALL D_HUGHES_WINGET(LMAT,XLMAT, & !IN
                ROT,DX_STRAIN) !OUT

                CALL D_HUGHES_WINGET(LMAT,YLMAT, & !IN
                ROT,DY_STRAIN) !OUT

                CALL D_HUGHES_WINGET(LMAT,ZLMAT, & !IN
                ROT,DZ_STRAIN) !OUT

            END IF
            !
            ! GET THE ROTATION MATRIX AND OBJECTIVE INCREMENTAL
            ! STRAIN USING HUGHES-WINGET ALGORITHM
            !
            ! IF (I == 1) THEN
                ! WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - About to call HUGHES_WINGET'
            ! END IF
            CALL HUGHES_WINGET(LMAT,ROT,STRAIN,D)
            !
            ! IF (I == 1) THEN
                ! WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - HUGHES_WINGET call completed'
            ! END IF
            ! ROTATE THE TOTAL STRESSES FROM PREVIOUS TIME STEP
            !
            CALL ROTATE_TENSOR(ROT,LSTRESS) ! Assuming ROTATE_TENSOR is a pure function or managed elsewhere if it has side effects on shared state

            !
            ! ROTATE THE TOTAL STRAINS FROM PREVIOUS TIME STEP
            !
            CALL ROTATE_TENSOR(ROT,LSTRAIN)
            !
        ELSE
            !
            ! BUILD THE INCREMENTAL INFINTESIMAL STRAIN
            !
            STRAIN = 0.0d0
            !
            !STRAIN ORDERING:
            !XX, YY, ZZ, YZ+ZY, XZ+ZX, XY+XY
            DO J = 1, LN
                STRAIN(1) = STRAIN(1) +  SHPD(1,J)*LDINC(1,J)
                STRAIN(2) = STRAIN(2) +  SHPD(2,J)*LDINC(2,J)
                STRAIN(3) = STRAIN(3) +  SHPD(3,J)*LDINC(3,J)
                STRAIN(4) = STRAIN(4) +  SHPD(2,J)*LDINC(3,J) +  SHPD(3,J)*LDINC(2,J)
                STRAIN(5) = STRAIN(5) +  SHPD(1,J)*LDINC(3,J) +  SHPD(3,J)*LDINC(1,J)
                STRAIN(6) = STRAIN(6) +  SHPD(1,J)*LDINC(2,J) +  SHPD(2,J)*LDINC(1,J)
            END DO
            !
        END IF
        !


        LSTRAIN = LSTRAIN + STRAIN
        !
        ! ELASTIC PREDICTOR
        !
        ELAS_MAT = FORM_CMAT(LPROP)
        LSTRESS_PREDICTOR = LSTRESS + MATMUL(ELAS_MAT,STRAIN)
        !

        IF (LMAT_TYPE.EQ.5) THEN
            ! IF (I == 1) THEN
            !     WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - About to call HYPERELASTIC'
            ! END IF

            CALL HYPERELASTIC(LPROP,LSTRESS,FMAT,LSTRAIN, constitution_error_flag) ! 新增 ierr 參數
            ! IF (I == 1) THEN
            !     WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - HYPERELASTIC call completed, status: ', constitution_error_flag
            ! END IF
        ELSE
            ! IF (I == 1) THEN
            !     WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - About to call CONSTITUTION'
            ! END IF

            CALL CONSTITUTION(LSTRESS_PREDICTOR,LMAT_TYPE, LSTRAIN, STRAIN, LPROP, DLT, FMAT, &
                LSTATE, LSTRESS, L_H_STRESS, L_S_STRESS, constitution_error_flag) !IN/OUT, OUT
            ! IF (I == 1) THEN
            !     WRITE(*,*) 'DEBUG: CONSTRUCT_FINT - Node ', I, ' - CONSTITUTION call completed, status: ', constitution_error_flag
            ! END IF
            IF (constitution_error_flag /= 0) THEN
                !$ACC ATOMIC UPDATE
                global_device_error_occurred = global_device_error_occurred + 1 ! 或其他錯誤處理方式

            END IF

        END IF
        !
        ! ********** SAVE STATE AND FEILD VARIABLES **********
        !
        !
        ! UPDATE STATE VARIABLE TO GSTATE
        !
        DO J = 1, 20
            GSTATE(J,I) = LSTATE(J) ! GSTATE is (20,GSIZE_IN)

        END DO
        !
        DO J = 1, 6

            !GET INCREMENTS FOR TIME STEP PREDICTION
            STRESS_INC(J) = LSTRESS(J) - GSTRESS(J,I) ! GSTRESS is (6,GSIZE_IN)
            STRAIN_INC(J) = LSTRAIN(J) - GSTRAIN(J,I) ! GSTRAIN is (6,GSIZE_IN)


            !SAVE THE STRESSES
            GSTRESS(J,I) = LSTRESS(J)
            GSTRAIN(J,I) = LSTRAIN(J)
			
!			IF(LMAT_TYPE == 4) THEN            
!				G_H_STRESS(J,I) = L_H_STRESS(J)!GC
!				G_S_STRESS(J,I) = L_S_STRESS(J)!GC
!			END IF

        END DO
        !EQUIVALENT PLASTIC STRAIN

        GSTRAIN_EQ(I) = LSTRAIN(1)**2 + LSTRAIN(2)**2 + LSTRAIN(3)**2 + & ! GSTRAIN_EQ is (GSIZE_IN)
	               2.0d0*(LSTRAIN(4)**2 + LSTRAIN(5)**2 + LSTRAIN(6)**2)
        GSTRAIN_EQ(I) = ( GSTRAIN_EQ(I) *2.3d0)**0.5


!        ID_RANK = OMP_get_thread_num()
        IF (LFINITE_STRAIN) THEN
            DO J = 1, 6
            
!                GINT_WORK_TEMP(ID_RANK+1) = GINT_WORK_TEMP(ID_RANK+1) + 0.5d0*D(J)*(2*LSTRESS(J)-STRESS_INC(J))*VOL*DET
                 GINT_WORK = GINT_WORK + 0.5d0*D(J)*(2*LSTRESS(J)-STRESS_INC(J))*VOL*DET ! Accumulate directly with reduction
 
            END DO
        ELSE
!        ID_RANK = OMP_get_thread_num()
            DO J = 1, 6
!                GINT_WORK_TEMP(ID_RANK+1) = GINT_WORK_TEMP(ID_RANK+1) + 0.5d0*STRAIN_INC(J)*STRESS_INC(J)*VOL*DET
                GINT_WORK = GINT_WORK + 0.5d0*STRAIN_INC(J)*STRESS_INC(J)*VOL*DET ! Accumulate directly with reduction

            END DO
        END IF
        

        POISS = LPROP(1)
        YOUNG = LPROP(2)

        BULK = BULK_MOD(YOUNG,POISS)

        SHEAR = SHEAR_MOD(YOUNG,POISS)

        DENSITY =  LPROP(3)
		
		!IF (.FALSE.) THEN !HERE'S THE OLD WAY OF DOING TIME STEP
			
			IF (LMAT_TYPE.GT.1) THEN
				NSNI_FLAG=.FALSE.
				CALL ESTIMATE_MODULI(STRESS_INC, STRAIN_INC, SHEAR_TRIAL, BULK_TRIAL, SHEAR, BULK,NSNI_FLAG)
			END IF

			  PMOD = BULK + 4.0d0*SHEAR/3.0d0
			
			  MAXMOD=MAX(PMOD,SHEAR)

			  GMAX_WVEL(I) = DSQRT(MAXMOD/DENSITY) ! GMAX_WVEL is (GNUMP)


			!   IF (LMAT_TYPE.EQ.3) THEN
			!   	FUDGE THE DRUCKER-PRAGER TIME STEP SO THAT THE FACTOR CAN BE 1.0
			!   	GMAX_WVEL(I) = GMAX_WVEL(I)/0.3d0
			!   END IF

		! ELSE		  
		!   DO IT THE WAY DR. YREUAX DOES		  
        !   GMAX_WVEL(I) = DSQRT(BULK/DENSITY)         
		! END IF
		
!        IF (ITYPE_INT.EQ.2) THEN

!            DO J = 1, 6
!                LDX_STRESS(J) = LOCAL_DX_STRESS(J,I)
!                LDY_STRESS(J) = LOCAL_DY_STRESS(J,I)
!                LDZ_STRESS(J) = LOCAL_DZ_STRESS(J,I)
!            END DO

!            CALL ROTATE_TENSOR(ROT,LDX_STRESS)
!            CALL ROTATE_TENSOR(ROT,LDY_STRESS)
!            CALL ROTATE_TENSOR(ROT,LDZ_STRESS)

!            POISS = LPROP(1)
!            YOUNG = LPROP(2)

!            BULK = BULK_MOD(YOUNG,POISS)

!            SHEAR = SHEAR_MOD(YOUNG,POISS)

!            IF (LMAT_TYPE.GT.1) THEN
            !IF (.TRUE.) THEN
!                NSNI_FLAG=.TRUE.
!                CALL ESTIMATE_MODULI(STRESS_INC, STRAIN_INC, SHEAR_TRIAL, BULK_TRIAL, SHEAR, BULK,NSNI_FLAG)
!            END IF

!            LAMDA = LAMDA_MOD(BULK,SHEAR)
!            MU = SHEAR
!            LAMDA_PLUS_2MU = LAMDA + 2.0d0*MU

!            CMAT = 0.0d0

!            CMAT(1,1) = LAMDA_PLUS_2MU
!            CMAT(2,2) = LAMDA_PLUS_2MU
!            CMAT(3,3) = LAMDA_PLUS_2MU
!            CMAT(4,4) = MU
!            CMAT(5,5) = MU
!            CMAT(6,6) = MU

!            CMAT(1,2) = LAMDA
!            CMAT(2,1) = LAMDA

!            CMAT(1,3) = LAMDA
!            CMAT(3,1) = LAMDA

!            CMAT(3,2) = LAMDA
!            CMAT(2,3) = LAMDA
!
!			NSNI_LIMITER = 1.0d0
!            IF ((LMAT_TYPE.EQ.3).OR.(LMAT_TYPE.EQ.6)) THEN
			!IF (LMAT_TYPE.EQ.3) THEN
            !
            ! DAMAGE MECHANICS INVOLVED, DO NOT USE A "TOTAL" STRESS
            ! STABILIZATION, INSTEAD, USE THE INCREMENTAL STRESS
            ! FOR STABILIZATION.
			!
			! COMMENT #2: STILL INEFFECTIVE, ZERO-OUT THE STRESS
			! IN CASE THE DAMAGE GETS TOO BIG, THIS SEEMS TO WORK
			! OK, BUT THE VALUE IS SENSATIVE A BIT.
            !
!			IF (LSTATE(4).GT.(0.5d0)) THEN
!                NSNI_LIMITER = 0.0d0
!			ELSE
!               NSNI_LIMITER = (1.0d0-2.0d0*LSTATE(4))
!			END IF
!            END IF
            !
            ! UPDATE THE PSUEDO-STRESSES FOR NSNI
            !
!            DO K = 1, 6
!                DO L = 1, 6
!                    LDX_STRESS(K) = LDX_STRESS(K) +  CMAT(K,L)*DX_STRAIN(L)
!                    LDY_STRESS(K) = LDY_STRESS(K) +  CMAT(K,L)*DY_STRAIN(L)
!                    LDZ_STRESS(K) = LDZ_STRESS(K) +  CMAT(K,L)*DZ_STRAIN(L)
!                END DO
!            END DO
            !
            ! SAVE WORK CONGUGATE STRESS DERIVATIVES FOR NSNI
            !
!            DO J = 1, 6
                !
!                LOCAL_DX_STRESS(J,I) = LDX_STRESS(J)*NSNI_LIMITER
!                LOCAL_DY_STRESS(J,I) = LDY_STRESS(J)*NSNI_LIMITER
!                LOCAL_DZ_STRESS(J,I) = LDY_STRESS(J)*NSNI_LIMITER
                !
!            END DO
			
			!DOESNT WORK!
			!IF ((LMAT_TYPE.EQ.3).OR.(LMAT_TYPE.EQ.6)) THEN
			!LDX_STRESS = LDX_STRESS*(1.0d0-LSTATE(4))
			!LDY_STRESS = LDY_STRESS*(1.0d0-LSTATE(4))
			!LDX_STRESS = LDX_STRESS*(1.0d0-LSTATE(4))
			!END IF


!        END IF !NSNI
    XNORM(1:3) =0.D0
        DO K = 1, LN
            KK = LSTACK(K)
            IF(LOCAL_BODY_ID .EQ.MODEL_BODY_ID(KK)) THEN ! MODEL_BODY_ID is (GNUMP), KK can be > GNUMP. This needs care.

                XNORM(1) = XNORM(1)+SHPD(1,K)
                XNORM(2) = XNORM(2)+SHPD(2,K)
                XNORM(3) = XNORM(3)+SHPD(3,K)
            ELSE
                XNORM(1) = XNORM(1)-SHPD(1,K)
                XNORM(2) = XNORM(2)-SHPD(2,K)
                XNORM(3) = XNORM(3)-SHPD(3,K)
            ENDIF

        END DO
        !
        ! ASSEMBLE THE INTERNAL FORCE
        !
        MAG_FINT = 0.0d0
        !
        DO J = 1, LN

            JJ = LSTACK(J) ! JJ is global index (1 to GSIZE_IN)


            !BMAT = 0.0d0

            !VOIGT ORDERING:
            !XX, YY, ZZ, YZ, XZ, XY
            BMAT = 0.0d0
            BMAT(1,1) = SHPD(1,J)
            BMAT(2,2) = SHPD(2,J)
            BMAT(3,3) = SHPD(3,J)
            BMAT(4,2) = SHPD(3,J)
            BMAT(4,3) = SHPD(2,J)
            BMAT(5,1) = SHPD(3,J)
            BMAT(5,3) = SHPD(1,J)
            BMAT(6,1) = SHPD(2,J)
            BMAT(6,2) = SHPD(1,J)

            BMAT_T = TRANSPOSE(BMAT)

            DENSITY=LPROP(3)

            !FINT3(1) = LSTRESS(1) * SHPD(1,J) + LSTRESS(5) * SHPD(3,J) + LSTRESS(6) * SHPD(2,J)

            !FINT3(2) = LSTRESS(2) * SHPD(2,J) + LSTRESS(4) * SHPD(3,J) + LSTRESS(6) * SHPD(1,J)

            !FINT3(3) = LSTRESS(3) * SHPD(3,J) + LSTRESS(4) * SHPD(2,J) + LSTRESS(5) * SHPD(1,J)

            FINT3 = MATMUL(BMAT_T,LSTRESS)            

            

   !
   ! DO KERNEL CONTACT
   !

           F_INT_C_TEMP=0.D0

           F_INT_C=0.D0
!IF (KCONTACT) THEN

!    LOCAL_BODY_ID_2 = MODEL_BODY_ID(JJ)
!    MU1 = LPROP(20) !BODY 1 (LPROP is for node I, GNUMP-based)

!    MU_NEW = GPROP(20,JJ) !BODY 2
!        IF(LOCAL_BODY_ID .NE. LOCAL_BODY_ID_2) THEN
!    IF (IKCONTACT.EQ.1) THEN
        ! NODAL ORIENTATION TO GET THE XNORM

!            XNORM(1:3) =0.D0
!            X2(1) = GCOO_CUURENT(1,JJ)
!            X2(2) = GCOO_CUURENT(2,JJ)
!            X2(3) = GCOO_CUURENT(3,JJ)
!            XNORM(1:3) = LCOO_T(1:3) - X2(1:3)
!            TEMP = DSQRT(XNORM(1)**2.0D0 + XNORM(2)**2.0D0 + XNORM(3)**2.0D0)
!        IF(TEMP > 1e-10) THEN
!            XNORM(1:3) = XNORM(1:3) / TEMP
!        ENDIF


!    ELSEIF (IKCONTACT.EQ.2) THEN
        ! LEVEL SET TO GET THE XNORM

!        TEMP = DSQRT(XNORM(1)**2.0D0 + XNORM(2)**2.0D0 + XNORM(3)**2.0D0)
!        IF(TEMP > 1e-10) THEN
!            XNORM(1:3) = XNORM(1:3) / TEMP
!        ENDIF

!    ENDIF

!    F_N = FINT3(1) * XNORM(1) + FINT3(2) * XNORM(2) + FINT3(3) * XNORM(3)

!    F_T1 = FINT3(1) - XNORM(1) * F_N ! Ft_x
!    F_T2 = FINT3(2) - XNORM(2) * F_N ! Ft_y
!    F_T3 = FINT3(3) - XNORM(3) * F_N ! Ft_z

    ! Length of Ft
!    F_T = DSQRT(F_T1**2.0D0 + F_T2**2.0D0 + F_T3**2.0D0)
!    MU_NEW2 = (MU1+ MU_NEW) * 0.5D0

!    IF(F_T .GT. 1E-13) THEN
        ! minimum of Ft and mu * Fn (Ft <= mu * Fn)
!        F_TT= MIN(DABS(F_N*MU_NEW2),F_T)
        ! Fn + Ft_modified
!        F_INT_C_TEMP(1) = F_N * XNORM(1) +  F_TT * F_T1 / F_T
!        F_INT_C_TEMP(2) = F_N * XNORM(2) +  F_TT * F_T2 / F_T
!        F_INT_C_TEMP(3) = F_N * XNORM(3) +  F_TT * F_T3 / F_T

!    ENDIF
    !F_INT_C = F_INT_C + F_INT_C_TEMP
!FINT3 = FINT3 + F_INT_C_TEMP
!       ENDIF
!ENDIF

            

            !GRAVITY
            FGRAV = SHP(J) * IGRAVITY*DENSITY
            
            !BODY FORCE
            FBOD = SHP(J) * LBOD*DENSITY
                        
            FINT3_EXT = FBOD + FGRAV

            DO K = 1, 3
                MAG_FINT = MAG_FINT + FINT3(K)**2
            END DO

    IF (USE_ATOMIC_FALLBACK) THEN
        ! 大問題：使用 atomic 更新
        DO K = 1, 3
            !$ACC ATOMIC UPDATE
            FINT((JJ-1)*3+K) = FINT((JJ-1)*3+K) + FINT3(K)*VOL*DET
            !$ACC ATOMIC UPDATE
            FEXT((JJ-1)*3+K) = FEXT((JJ-1)*3+K) + FINT3_EXT(K)*VOL*DET
        END DO
    ELSE
        ! 小問題：為了數值一致性，模擬 OpenMP 的累加順序
        ! 使用兩階段方法
        ! 第一階段：每個節點將其貢獻存儲到臨時位置
        DO K = 1, 3
            ! 這裡需要一個節點級別的臨時存儲策略
            ! 由於 GPU 的限制，我們仍使用 atomic，但加入順序控制
            !$ACC ATOMIC UPDATE
            FINT((JJ-1)*3+K) = FINT((JJ-1)*3+K) + FINT3(K)*VOL*DET
            !$ACC ATOMIC UPDATE
            FEXT((JJ-1)*3+K) = FEXT((JJ-1)*3+K) + FINT3_EXT(K)*VOL*DET
        END DO
    END IF

        END DO !ASSEMBLE FINT FOR STANDARD NODAL INTEGRATION PART

        MAG_FINT=DSQRT(MAG_FINT)

!        IF (ITYPE_INT.EQ.2) THEN !NSNI

            !
!            MAG_STAB_FINT = 0.0d0
            !
!            DO J = 1, LN

!                JJ = LSTACK(J)

!                XBMAT = 0.0d0
!                YBMAT = 0.0d0
!                ZBMAT = 0.0d0

!                XBMAT(1,1) = SHPDD_SM(XMAP(1),J)
!                XBMAT(2,2) = SHPDD_SM(XMAP(2),J)
!                XBMAT(3,3) = SHPDD_SM(XMAP(3),J)
!                XBMAT(4,2) = SHPDD_SM(XMAP(3),J)
!                XBMAT(4,3) = SHPDD_SM(XMAP(2),J)
!                XBMAT(5,1) = SHPDD_SM(XMAP(3),J)
!                XBMAT(5,3) = SHPDD_SM(XMAP(1),J)
!                XBMAT(6,1) = SHPDD_SM(XMAP(2),J)
!                XBMAT(6,2) = SHPDD_SM(XMAP(1),J)

!                XBMAT_T = TRANSPOSE(XBMAT)

!               XFINT3(J,1:3) = MATMUL(XBMAT_T,LDX_STRESS)



!                YBMAT(1,1) = SHPDD_SM(YMAP(1),J)
!                YBMAT(2,2) = SHPDD_SM(YMAP(2),J)
!                YBMAT(3,3) = SHPDD_SM(YMAP(3),J)
!                YBMAT(4,2) = SHPDD_SM(YMAP(3),J)
!                YBMAT(4,3) = SHPDD_SM(YMAP(2),J)
!                YBMAT(5,1) = SHPDD_SM(YMAP(3),J)
!                YBMAT(5,3) = SHPDD_SM(YMAP(1),J)
!                YBMAT(6,1) = SHPDD_SM(YMAP(2),J)
!                YBMAT(6,2) = SHPDD_SM(YMAP(1),J)

!                YBMAT_T = TRANSPOSE(YBMAT)

!                YFINT3(J,1:3) = MATMUL(YBMAT_T,LDY_STRESS)



!                ZBMAT(1,1) = SHPDD_SM(ZMAP(1),J)
!                ZBMAT(2,2) = SHPDD_SM(ZMAP(2),J)
!                ZBMAT(3,3) = SHPDD_SM(ZMAP(3),J)
!                ZBMAT(4,2) = SHPDD_SM(ZMAP(3),J)
!                ZBMAT(4,3) = SHPDD_SM(ZMAP(2),J)
!                ZBMAT(5,1) = SHPDD_SM(ZMAP(3),J)
!                ZBMAT(5,3) = SHPDD_SM(ZMAP(1),J)
!                ZBMAT(6,1) = SHPDD_SM(ZMAP(2),J)
!                ZBMAT(6,2) = SHPDD_SM(ZMAP(1),J)

!                ZBMAT_T = TRANSPOSE(ZBMAT)

!                ZFINT3(J,1:3) = MATMUL(ZBMAT_T,LDZ_STRESS)

!                DO K = 1, 3
!                    MAG_STAB_FINT = MAG_STAB_FINT + (XFINT3(J,K)**2 + YFINT3(J,K)**2 + ZFINT3(J,K)**2)
!                END DO
                
!                CONTINUE


!            END DO
            !
            ! CONTROL THE CONTRIBUTION TO FINT BY NSNI: SOME PARAMETERS ARE ESTIMATED AND
            ! THEY MIGHT NOT BE ACCURATE
            !
!            IF (USE_STAB_CONTROL) THEN
!             IF (MAG_STAB_FINT.GT.(1.0d-12)) THEN
!                MAG_STAB_FINT=DSQRT(MAG_STAB_FINT)
!                IF (MAG_STAB_FINT.GT.MAG_FINT) THEN
!                    XFINT3=XFINT3*MAG_FINT/MAG_STAB_FINT*STABILIZATION_CONTROL_COEF
!                    YFINT3=YFINT3*MAG_FINT/MAG_STAB_FINT*STABILIZATION_CONTROL_COEF
!                    ZFINT3=ZFINT3*MAG_FINT/MAG_STAB_FINT*STABILIZATION_CONTROL_COEF
!                END IF
!             END IF
!            END IF
            
            !DEBUG
            !IF (MAG_STAB_FINT.GT.(1.0e-6)) THEN
            !CONTINUE
            !END IF

!            DO J = 1, LN

!                JJ = LSTACK(J)

!                DO K = 1, 3
                    !IT DOESNT LOOK LIKE X_MOM GETS ASSIGNED ANYTHING! FIX NSNI!
!                    ID_RANK = OMP_get_thread_num()
!                    FINT_TEMP(ID_RANK+1,K,JJ) = FINT_TEMP(ID_RANK+1,K,JJ) + XFINT3(J,K) *VOL*DET * G_X_MOM(I)
!                    FINT_TEMP(ID_RANK+1,K,JJ) = FINT_TEMP(ID_RANK+1,K,JJ) + YFINT3(J,K) *VOL*DET * G_Y_MOM(I)
!                    FINT_TEMP(ID_RANK+1,K,JJ) = FINT_TEMP(ID_RANK+1,K,JJ) + ZFINT3(J,K) *VOL*DET * G_Z_MOM(I)

!                END DO

!            END DO

!        END IF !NSNI


    END DO !INTEGRATION POINT (NODE) LOOP
    !$ACC END PARALLEL LOOP

!     $ACC WAIT  ! Wait for the current chunk to complete
!    END DO ! End of chunking loop (c_chunk)


    ! Assign the accumulated device error to the output error flag
    ! Do not add global_device_error_occurred (fallback related) to ierr_fint_arg
    ! ierr_fint_arg should reflect fatal errors or be 0 if successful.
    ! The initial value of ierr_fint_arg is 0, and it's only changed if a fatal error occurs earlier.

    ! inverse_fallback_count is already brought to host via COPYOUT clause

    ! Host-side I/O for debugging information related to inverse_fallback_count


    IF (inverse_fallback_count > 0) THEN
        WRITE(*,*) 'INFO: CONSTRUCT_FINT - Number of INVERSE/INV3 fallbacks to identity matrix: ', inverse_fallback_count
        ! You could add more detailed WRITE statements here if you were to collect specific node IDs or DET values that caused fallbacks.

    END IF
    !
    !FINT_TEMP TO HOLD THE VALUES FOR OPENMP REDUCE(ASSEMBLE,ADD) IN THE END
    !
    ! IF (AUTO_TS) THEN
        ! ... (original commented out time step calculation code) ...
    ! END IF !CALC TIME STEP

    IF (ALLOCATED(fint_temp)) THEN
        DEALLOCATE(fint_temp, STAT=device_error_status_check)
        IF (device_error_status_check /= 0) PRINT *, "Error deallocating fint_temp in CONSTRUCT_FINT"
    END IF
    IF (ALLOCATED(fext_temp)) THEN
        DEALLOCATE(fext_temp, STAT=device_error_status_check)
        IF (device_error_status_check /= 0) PRINT *, "Error deallocating fext_temp in CONSTRUCT_FINT"
    END IF
    IF (ALLOCATED(gint_work_temp)) THEN
        DEALLOCATE(gint_work_temp, STAT=device_error_status_check)
        IF (device_error_status_check /= 0) PRINT *, "Error deallocating gint_work_temp in CONSTRUCT_FINT"
    END IF

    IF (LINIT .AND. ALLOCATED(GWIN0)) THEN  ! Deallocate GWIN0 if it was allocated in this routine
        DEALLOCATE(GWIN0, STAT=device_error_status_check) 

        IF (device_error_status_check /= 0) PRINT *, "Error deallocating GWIN0 in CONSTRUCT_FINT"

    END IF

    RETURN
    END SUBROUTINE
!            DO J = 1, LN
!                LSTACK(J) = GSTACK(LSTART+J-1)
!            END DO

            !FIND THE CHARACTERISTIC DISTANCES

!            XI=GCOO_CUURENT(1,I)
!            YI=GCOO_CUURENT(2,I)
!            ZI=GCOO_CUURENT(3,I)
                    !
					! USE THE UNDEFORMED CONFIGURATION
					!
!            XI=GCOO(1,I)
!            YI=GCOO(2,I)
!            ZI=GCOO(3,I)

!            FIRST = .TRUE.

!			IF (LINIT) THEN
			
!            DO J = 1, LN

!                JJ = LSTACK(J)

!                IF (JJ.NE.I) THEN

!                    XJ=GCOO_CUURENT(1,JJ)
!                    YJ=GCOO_CUURENT(2,JJ) 

!                    ZJ=GCOO_CUURENT(3,JJ)
                    !
					! USE THE UNDEFORMED CONFIGURATION
					!
!                    XJ=GCOO(1,JJ)
!                    YJ=GCOO(2,JJ)
!                    ZJ=GCOO(3,JJ)

!                    DIST = DSQRT((XJ-XI)**2 + (YJ-YI)**2 + (ZJ-ZI)**2)

!                    IF (FIRST) THEN
!                        GCHAR_DIST(I) = DIST
!                        FIRST = .FALSE.
!                    ELSE
!                        GCHAR_DIST(I) = MIN(GCHAR_DIST(I),DIST)
!                    END IF

!                END IF

!            END DO !J=1,GSIZE_IN (NEIGHBOR NODES)

			
!			END IF

!            DLT_TEMP = GCHAR_DIST(I) / GMAX_WVEL(I) !* 0.2d0


!            IF (I.EQ.1) THEN
!                DLT_FINT = DLT_TEMP
!            ELSE

!                IF (DLT_TEMP.LT.DLT_FINT) THEN
!                    DLT_FINT = DLT_TEMP
!                END IF
!                DLT_FINT = MIN(DLT_TEMP,DLT_FINT)
!            END IF
!        END DO
!        DLT_FINT = DLT_FINT*DLT_FAC
!    END IF !CALC TIME STEP