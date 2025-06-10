
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
    
      SUBROUTINE HANDELER(GWIN,    GVOL,      GNUMP,     GCOO,      GCOO_CUURENT,      &   !FROM MAIN
                          GSM_LEN, GSM_AREA,  GSM_VOL,   GNSNI_FAC, GGHOST,       &   !FROM MAIN
                          GPROP,   GSTATE,    GSTRESS,   GSTRAIN,   G_H_STRESS, G_S_STRESS,   GDINC,     GDINC_TOT,     &   !FROM MAIN
                          GFINT,   GMAT_TYPE, GSIZE,     LINIT,     DLT_FINT, &
                          GFEXT,      GCHAR_DIST,   GMAX_WVEL, GEBC_NODES, &
						  DO_INTERP, GDINC_PHY, GVEL_PHY, GVEL, GACL, GACL_PHY, &
						  LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
                          G_X_MOM,     G_Y_MOM,    G_Z_MOM, &  
                          GXDIST_MAX, GYDIST_MAX, GZDIST_MAX, MODEL_BODYFORCE,DLT, GINT_WORK, MODEL_BODY_ID, GSTRAIN_EQ, &
						  GIJKSPC, H_IERR)
                          
      
      USE CONTROL
      !USE THE CONTROL FILE SINCE WE WILL NEED IT FOR NODE SEARCH ALGORITHMS
      
      IMPLICIT NONE
      
      INTEGER:: GSIZE
      
	  !DATA
      INTEGER:: GNUMP
	  DOUBLE PRECISION:: GVOL(GNUMP)
	  DOUBLE PRECISION:: GWIN(3,GNUMP)
	  DOUBLE PRECISION:: GCOO(3,GSIZE)
	  DOUBLE PRECISION:: GCOO_CUURENT(3,GSIZE)
      DOUBLE PRECISION:: GSM_LEN(6,GNUMP)
      DOUBLE PRECISION:: GSM_AREA(3,GNUMP)
      DOUBLE PRECISION:: GSM_VOL(GNUMP)
      DOUBLE PRECISION:: GNSNI_FAC(3,GNUMP)
      DOUBLE PRECISION:: DLT_FINT
      
	  DOUBLE PRECISION:: GCHAR_DIST(GNUMP),   GMAX_WVEL(GNUMP)
      INTEGER:: GIJKSPC(3,GNUMP)
	  
      INTEGER:: GGHOST(GNUMP)
      DOUBLE PRECISION:: GPROP(30,GNUMP)
      DOUBLE PRECISION:: GSTATE(20,GSIZE)
      DOUBLE PRECISION:: GSTRESS(6,GSIZE)
      DOUBLE PRECISION:: GSTRAIN(6,GSIZE)
      DOUBLE PRECISION:: GDINC(3*GSIZE),GDINC_TOT(3*GSIZE)
      DOUBLE PRECISION:: GDINC_PHY(3*GSIZE)
      
      DOUBLE PRECISION:: GACL_PHY(3*GSIZE)
      DOUBLE PRECISION:: GACL(3*GSIZE)
      
      DOUBLE PRECISION:: GVEL_PHY(3*GSIZE)
      DOUBLE PRECISION:: GVEL(3*GSIZE)
      
      DOUBLE PRECISION:: GFINT(3*GSIZE), GFEXT(3*GSIZE)
	  LOGICAL:: GEBC_NODES(GSIZE)
      INTEGER:: GMAT_TYPE(GNUMP)
      !
    DOUBLE PRECISION:: G_H_STRESS(6,GNUMP)!GC
    DOUBLE PRECISION:: G_S_STRESS(6,GNUMP)!GC
    DOUBLE PRECISION:: DLT
      
	  INTEGER,ALLOCATABLE, SAVE:: GN(:)                     !NUMBER OF NEIGHBORS FOR EACH NODE
     !$ACC DECLARE LINK(GN)
	  INTEGER,ALLOCATABLE, SAVE:: GSTART(:)                 !START LOCATION OF NODE NEIOGHBORS IN STACK
     !$ACC DECLARE LINK(GSTART)
	  INTEGER, SAVE:: DIM_NN_LIST                           !SIZE OF NEIGHBOR STACK
     !$ACC DECLARE COPY(DIM_NN_LIST)
	  INTEGER,ALLOCATABLE, SAVE:: GSTACK(:)                 !NEIGHBORS FOR EACH NODE (STACKED)
     !$ACC DECLARE LINK(GSTACK)
	  DOUBLE PRECISION,ALLOCATABLE, SAVE:: GSTACK_SHP(:)                 !NEIGHBORS FOR EACH NODE (STACKED)
     !$ACC DECLARE LINK(GSTACK_SHP)
	  DOUBLE PRECISION,ALLOCATABLE, SAVE:: GSTACK_DSHP(:,:)                 !NEIGHBORS FOR EACH NODE (STACKED)
     !$ACC DECLARE LINK(GSTACK_DSHP)
	  DOUBLE PRECISION,ALLOCATABLE, SAVE:: GSTACK_DDSHP(:,:)                 !NEIGHBORS FOR EACH NODE (STACKED)
     !$ACC DECLARE LINK(GSTACK_DDSHP)
	  DOUBLE PRECISION,ALLOCATABLE, SAVE:: GINVK(:,:,:)      
     !$ACC DECLARE LINK(GINVK)
	  INTEGER, SAVE:: GMAXN                                 !MAX NUMBER OF NEIGHBORS FOR ALL NODES
     !$ACC DECLARE COPY(GMAXN)
	  INTEGER, SAVE:: SEARCHCOUNT
           !$ACC DECLARE COPY(SEARCHCOUNT)
      !BINNING INFORMATION
      
      DOUBLE PRECISION,SAVE:: XBIN,YBIN,ZBIN,XMAX,YMAX,ZMAX,XMIN,YMIN,ZMIN,XBIN_TEMP,YBIN_TEMP,ZBIN_TEMP
        !$ACC DECLARE COPY(XBIN,YBIN,ZBIN,XMAX,YMAX,ZMAX,XMIN,YMIN,ZMIN)
 

      DOUBLE PRECISION,SAVE::DOM_XLENGTH,DOM_YLENGTH,DOM_ZLENGTH
     !$ACC DECLARE COPY(DOM_XLENGTH,DOM_YLENGTH,DOM_ZLENGTH)


      INTEGER,SAVE::NBINSX,NBINSY,NBINSZ,NBINS,MAX_NEIGH
     !$ACC DECLARE COPY(NBINSX,NBINSY,NBINSZ,NBINS,MAX_NEIGH)

      INTEGER,SAVE,ALLOCATABLE::ISPACE(:),JSPACE(:),KSPACE(:),NODES_IN_BIN(:),NODELIST_IN_BIN(:,:)
     !$ACC DECLARE LINK(ISPACE,JSPACE,KSPACE,NODES_IN_BIN,NODELIST_IN_BIN)

      
      ! Scalar SAVE variables that are modified and used in logic or passed to kernels

      DOUBLE PRECISION, SAVE:: CNT_SEARCH
        !$ACC DECLARE COPY(CNT_SEARCH)
      
      
	  LOGICAL:: LINIT !TRUE = FIRST TIME STEP
	  
	  LOGICAL:: DO_INTERP
      	  
	  INTEGER:: I, J, K, M, LSTART, MM
	  DOUBLE PRECISION:: SHPT
	  
      !0329
	  DOUBLE PRECISION:: MODEL_BODYFORCE(3,GNUMP)
	  INTEGER, INTENT(IN):: MODEL_BODY_ID(GNUMP)
       DOUBLE PRECISION:: GSTRAIN_EQ(GSIZE)
	  !NSNI
      DOUBLE PRECISION::  LOCAL_DX_STRESS(6,GSIZE)
      DOUBLE PRECISION::  LOCAL_DY_STRESS(6,GSIZE)
      DOUBLE PRECISION::  LOCAL_DZ_STRESS(6,GSIZE)
	  
      DOUBLE PRECISION, INTENT(IN)::GXDIST_MAX, GYDIST_MAX, GZDIST_MAX
      DOUBLE PRECISION, INTENT(IN)::G_X_MOM(GNUMP),G_Y_MOM(GNUMP),G_Z_MOM(GNUMP)
	  
	    INTEGER:: BINS_SAFETY_FAC
      DOUBLE PRECISION::  GINT_WORK      

      INTEGER :: ierr_fint
      INTEGER, INTENT(OUT) :: H_IERR
      INTEGER :: temp_gmaxn_device_check ! Declare the temporary variable
      INTEGER(KIND=8) :: temp_dim_long
      LOGICAL :: NEED_SEARCH         ! MOVED: Declaration for search logic
      LOGICAL :: NEED_BIN_UPDATE     ! MOVED: Declaration for bin update logic
      
IF (DO_INTERP) THEN
    IF(.NOT. PERIDYNAMICS) THEN  !RKPM
        ! 先確保 GSTACK_SHP 已經計算
        IF (LINIT) THEN
            WRITE(*,*) 'WARNING: DO_INTERP called but shape functions may not be ready during LINIT'
            ! 設定預設值避免錯誤
            !$ACC PARALLEL LOOP DEFAULT(PRESENT)
            DO I = 1, 3*GNUMP
                GDINC_PHY(I) = 0.0d0
                GVEL_PHY(I) = 0.0d0
                GACL_PHY(I) = 0.0d0
            END DO
            !$ACC END PARALLEL LOOP
        ELSE
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) &
            !$ACC PRIVATE(LSTART, K, M, J, SHPT, MM)
            DO I=1,GNUMP
                LSTART = GSTART(I)
                DO K = 1, 3
                    M=(I-1)*3+K
                    GDINC_PHY(M) = 0.0d0
                    GVEL_PHY(M) = 0.0d0
                    GACL_PHY(M) = 0.0d0
                    IF (GN(I) > 0) THEN  ! 檢查是否有鄰居
                        DO J = 1, GN(I)
                            SHPT=GSTACK_SHP(LSTART+J-1)
                            MM = (GSTACK(LSTART+J-1) - 1)*3 + K
                            GDINC_PHY(M) = GDINC_PHY(M) +  SHPT*GDINC(MM)
                            GVEL_PHY(M) = GVEL_PHY(M) +  SHPT*GVEL(MM)
                            GACL_PHY(M) = GACL_PHY(M) +  SHPT*GACL(MM)
                        END DO
                    END IF
                END DO
            END DO
            !$ACC END PARALLEL LOOP
        END IF
    ELSE  !PERIDYNAMICS
            GDINC_PHY = GDINC
            GVEL_PHY = GVEL
            GACL_PHY = GACL
        
        END IF
				
				
      ELSE
      
      !CALL THE NODE SEARCH ROUTINE HERE IF NEEDED
      IF ((LINIT).OR.(.NOT.LLAGRANGIAN)) THEN
      H_IERR = 0
      IF (GNUMP <= 0) THEN
          WRITE(*,*) 'ERROR: HANDELER - GNUMP <= 0 before calculating DIM_NN_LIST. GNUMP = ', GNUMP
          H_IERR = -40 
          RETURN
      END IF    
      WRITE(*,*) 'DEBUG: HANDELER - GNUMP = ', GNUMP, ' (before DIM_NN_LIST calculation)'
      
      ! 使用與 OpenMP 版本相同的計算方式
      DIM_NN_LIST = GNUMP * 1000
      GMAXN = GNUMP
      WRITE(*,*) 'DEBUG: HANDELER - Host DIM_NN_LIST set to = ', DIM_NN_LIST
      ! Immediately update DIM_NN_LIST on the device
      !$ACC UPDATE DEVICE(DIM_NN_LIST)

      ! GMAXN is typically determined by search routines.
      ! If GMAXN has a default or pre-search value, it can be set here.
      ! For now, GMAXN will be set by SOFT_SEARCH.
      ! WRITE(*,*) 'DEBUG: HANDELER - Initial GMAXN (before search) = ', GMAXN  ! Debug print if GMAXN has an initial value

      ! Update GMAXN on the device AFTER it has been set on the host
      ! If GMAXN and DIM_NN_LIST are used in device kernels within HANDELER itself,
      ! an !$ACC UPDATE DEVICE(GMAXN, DIM_NN_LIST) would be needed here.



      IF (LINIT) THEN
          ALLOCATE(GN(GNUMP)) 
          ALLOCATE(GSTART(GNUMP)) 
          
          ! DIM_NN_LIST 必須與 OpenMP 版本完全相同
          DIM_NN_LIST=GNUMP*1000  ! 與 OpenMP 版本一致
          GMAXN=GNUMP             ! 與 OpenMP 版本一致
          
          ALLOCATE(GSTACK(DIM_NN_LIST))                
          ALLOCATE(GSTACK_SHP(DIM_NN_LIST))                    
          ALLOCATE(GSTACK_DSHP(3,DIM_NN_LIST))  
          ALLOCATE(GSTACK_DDSHP(6,DIM_NN_LIST)) 
          ALLOCATE(GINVK(3,3,GNUMP))
          
          ! 初始化陣列為 0（與 OpenMP 版本的隱含行為一致）
          GN = 0
          GSTART = 0
          GSTACK = 0
          GSTACK_SHP = 0.0D0
          GSTACK_DSHP = 0.0D0
          GSTACK_DDSHP = 0.0D0
          GINVK = 0.0D0
          
          ! 將初始化的資料複製到 GPU
          !$ACC ENTER DATA COPYIN(GN, GSTART, GSTACK, GSTACK_SHP, GSTACK_DSHP, GSTACK_DDSHP, GINVK)
          
          CNT_SEARCH = 0.0d0
          SEARCHCOUNT=PDSEARCH
          
          ! 更新標量到 GPU
          !$ACC UPDATE DEVICE(CNT_SEARCH, SEARCHCOUNT, DIM_NN_LIST, GMAXN)
      END IF
      
      !
      ! HARD_SEARCH
      !
      ! WORKS BY OPENING THE FOLLOWING COMMENT, ESPECIALLY FOR DEBUGGING PURPOSE
      !CALL HARD_SEARCH(GNUMP,GCOO_CUURENT,GN,GSTART,DIM_NN_LIST,GSTACK,GMAXN,GWIN,GXDIST_MAX, GYDIST_MAX,GZDIST_MAX)

	  CNT_SEARCH = CNT_SEARCH + DLT
      ! If CNT_SEARCH is used in device kernels within HANDELER for logic:
      ! !$ACC UPDATE DEVICE(CNT_SEARCH) 

	  
      !IF ((CNT_SEARCH.GE.TIME_SEARCH).OR.(LINIT)) THEN
	  
	  IF (.TRUE.) THEN !THE TIMED BIN SEARCH (INTERMIETENT) DOESN'T WORK!
	  

	  CNT_SEARCH = 0.0d0
      ! 宣告搜尋控制變數

      NEED_SEARCH = .FALSE.
      NEED_BIN_UPDATE = .FALSE.

      !!!
      !!! BUCKET_SEARCH
      !!!      
         
        
         
         IF (LINIT) THEN
         !
		 ! GET THE BIN SIZE, CONSIDERING THE MAX WINDOW SIZE, 
		 ! AND MAX SMOOTHING SIZE, AND BOUNDING BOX FOR THE DOMAIN.
		 ! SINCE WE DONT REGENERATE THE BINS, JUST DO THIS ONCE!
		 !
         XBIN = 0.0d0
         YBIN = 0.0d0
         ZBIN = 0.0d0
         
         I=1
         XMAX = GCOO_CUURENT(1,I)
         YMAX = GCOO_CUURENT(2,I)
         ZMAX = GCOO_CUURENT(3,I)
         
         XMIN = GCOO_CUURENT(1,I)
         YMIN = GCOO_CUURENT(2,I)
         ZMIN = GCOO_CUURENT(3,I)
         
         DO I=1,GNUMP
         
           IF (GWIN(1,I).GT.XBIN) XBIN=GWIN(1,I) 
           IF (GWIN(2,I).GT.YBIN) YBIN=GWIN(2,I)
           IF (GWIN(3,I).GT.ZBIN) ZBIN=GWIN(3,I)
           
           IF (GWIN(1,I)+GXDIST_MAX.GT.XBIN) XBIN=GWIN(1,I)+GXDIST_MAX
           IF (GWIN(2,I)+GYDIST_MAX.GT.YBIN) YBIN=GWIN(2,I)+GYDIST_MAX
           IF (GWIN(3,I)+GZDIST_MAX.GT.ZBIN) ZBIN=GWIN(3,I)+GZDIST_MAX           
           
           
           IF (GCOO_CUURENT(1,I).GT.XMAX) XMAX=GCOO_CUURENT(1,I)
           IF (GCOO_CUURENT(2,I).GT.YMAX) YMAX=GCOO_CUURENT(2,I)
           IF (GCOO_CUURENT(3,I).GT.ZMAX) ZMAX=GCOO_CUURENT(3,I)
           
           IF (GCOO_CUURENT(1,I).LT.XMIN) XMIN=GCOO_CUURENT(1,I)
           IF (GCOO_CUURENT(2,I).LT.YMIN) YMIN=GCOO_CUURENT(2,I)
           IF (GCOO_CUURENT(3,I).LT.ZMIN) ZMIN=GCOO_CUURENT(3,I)
           
         END DO
       !  SHOULD BE RELATED TO SMOOTHING FACTOR
         XBIN = XBIN !*1.1d0
         YBIN = YBIN !*1.1d0
         ZBIN = ZBIN !*1.1d0
         
        ! ADD SOME SAFE LENGTH TO THE EDGES
         XMIN = XMIN - XBIN * 10.0D0
         YMIN = YMIN - YBIN * 10.0D0
         ZMIN = ZMIN - ZBIN * 10.0D0
         
         XMAX = XMAX + XBIN * 10.0D0
         YMAX = YMAX + YBIN * 10.0D0
         ZMAX = ZMAX + ZBIN * 10.0D0        
         
         
         DOM_XLENGTH = XMAX - XMIN
         DOM_YLENGTH = YMAX - YMIN
         DOM_ZLENGTH = ZMAX - ZMIN
        
         NBINSX = CEILING(DOM_XLENGTH/XBIN) !+1
         NBINSY = CEILING(DOM_YLENGTH/YBIN) !+1
         NBINSZ = CEILING(DOM_ZLENGTH/ZBIN) !+1
         
         NBINS = NBINSX*NBINSY*NBINSZ
         
         MAX_NEIGH= 1000 !MAX NODES IN A BIN, SET TO LARGE VALUE FOR NOW
		 BINS_SAFETY_FAC = 5 !NOT USED

           ALLOCATE(ISPACE(GNUMP),JSPACE(GNUMP),KSPACE(GNUMP))
           ALLOCATE(NODES_IN_BIN(NBINSX*NBINSY*NBINSZ))
           ALLOCATE(NODELIST_IN_BIN(NBINSX*NBINSY*NBINSZ,MAX_NEIGH))
           !$ACC ENTER DATA CREATE(ISPACE,JSPACE,KSPACE,NODES_IN_BIN,NODELIST_IN_BIN)
           !$ACC UPDATE DEVICE(XBIN,YBIN,ZBIN,XMAX,YMAX,ZMAX,XMIN,YMIN,ZMIN)
           !$ACC UPDATE DEVICE(DOM_XLENGTH,DOM_YLENGTH,DOM_ZLENGTH)
           !$ACC UPDATE DEVICE(NBINSX,NBINSY,NBINSZ,NBINS,MAX_NEIGH)


           ! If these binning arrays/scalars are used in device kernels within HANDELER:
           ! !$ACC UPDATE DEVICE(ISPACE,JSPACE,KSPACE,NODES_IN_BIN,NODELIST_IN_BIN)
           ! !$ACC UPDATE DEVICE(XBIN,YBIN,ZBIN,XMAX,YMAX,ZMAX,XMIN,YMIN,ZMIN)
           ! !$ACC UPDATE DEVICE(DOM_XLENGTH,DOM_YLENGTH,DOM_ZLENGTH)
           ! !$ACC UPDATE DEVICE(NBINSX,NBINSY,NBINSZ,NBINS,MAX_NEIGH)



         END IF !(LINIT) 
         ! 決定是否需要進行鄰居搜尋
         ! 三種情況需要搜尋：
         ! 1. 初始化時（LINIT = .TRUE.）
         ! 2. 非拉格朗日方法且達到搜尋間隔
         ! 3. 強制更新（未來擴充用）
         
         IF (LINIT) THEN
             ! === 情況 1：初始化必須搜尋 ===
             NEED_SEARCH = .TRUE.
             NEED_BIN_UPDATE = .TRUE.  ! 初始化時也要更新 bin
             WRITE(*,*) 'DEBUG: HANDELER - Initial search required (LINIT=.TRUE.)'
             
         ELSEIF (.NOT.LLAGRANGIAN) THEN
             ! === 情況 2：歐拉或 ALE 方法 ===
             ! 檢查是否達到搜尋間隔
             IF (SEARCHCOUNT.GE.PDSEARCH) THEN
                 NEED_SEARCH = .TRUE.
                 NEED_BIN_UPDATE = .TRUE.  ! 達到間隔時更新 bin
                 SEARCHCOUNT = 1
                 WRITE(*,*) 'DEBUG: HANDELER - Periodic search triggered, SEARCHCOUNT reset to 1'
             ELSE
                 ! 還沒到搜尋時間，計數器加一
                 SEARCHCOUNT = SEARCHCOUNT + 1
                 NEED_SEARCH = .FALSE.
             END IF
             
         ELSE
             ! === 情況 3：拉格朗日方法（初始化後不需要搜尋）===
             NEED_SEARCH = .FALSE.
             WRITE(*,*) 'DEBUG: HANDELER - Lagrangian method, no search needed'
         END IF
         ! ==== 修改區塊結束 ====
         
         IF (NEED_SEARCH) THEN
             WRITE(*,*) 'DEBUG: HANDELER - Executing search...'
             
             ! --- Step 1: 完整同步所有需要的資料到 CPU ---
             !$ACC UPDATE HOST(GCOO_CUURENT, GWIN, GSM_LEN)
             
             ! --- Step 2: 更新空間分割（與 OpenMP 版本完全相同的邏輯）---
             IF (NEED_BIN_UPDATE) THEN
                 CALL GET_BINS(GNUMP,GCOO_CUURENT,NODES_IN_BIN,MAX_NEIGH,NODELIST_IN_BIN, &
                               NBINS,NBINSX,NBINSY,NBINSZ,ISPACE,JSPACE,KSPACE,XMIN,YMIN,ZMIN,XBIN,YBIN,ZBIN)
                 
                 DO I=1,GNUMP
                     GIJKSPC(1,I) = ISPACE(I)
                     GIJKSPC(2,I) = JSPACE(I)
                     GIJKSPC(3,I) = KSPACE(I)
                 END DO
             END IF
            
             ! --- Step 3: 執行鄰居搜尋（與 OpenMP 版本相同）---
             CALL SOFT_SEARCH(GNUMP,GCOO_CUURENT,GN,GSTART,DIM_NN_LIST,GSTACK,GMAXN,GWIN, &
                              NODES_IN_BIN,MAX_NEIGH,NODELIST_IN_BIN, &
                              NBINS,NBINSX,NBINSY,NBINSZ,ISPACE,JSPACE,KSPACE, &
                              GXDIST_MAX, GYDIST_MAX, GZDIST_MAX,GSM_LEN)
             
! --- Step 4: 完整同步所有修改過的資料回 GPU ---
! 關鍵：確保所有形狀函數相關資料都同步
!$ACC UPDATE DEVICE(GN, GSTART, GSTACK, GMAXN, DIM_NN_LIST)
!$ACC UPDATE DEVICE(GSTACK_SHP, GSTACK_DSHP, GSTACK_DDSHP)  ! 新增：同步形狀函數
!$ACC UPDATE DEVICE(GINVK)  ! 新增：同步逆矩陣
!$ACC UPDATE DEVICE(NODES_IN_BIN, NODELIST_IN_BIN)
!$ACC UPDATE DEVICE(ISPACE, JSPACE, KSPACE, GIJKSPC)
!$ACC UPDATE DEVICE(GCOO_CUURENT, GWIN, GSM_LEN)

! 等待所有更新完成
!$ACC WAIT



      WRITE(*,*) 'DEBUG: HANDELER - After SOFT_SEARCH, Host GMAXN = ', GMAXN

      !$ACC UPDATE DEVICE(GMAXN) ! Ensure GMAXN is updated on device if SOFT_SEARCH changed it on host

      WRITE(*,*) 'DEBUG: HANDELER - Final Host GMAXN before CONSTRUCT_FINT call (after SOFT_SEARCH) = ', GMAXN
      ! GMAXN is updated on the device after SOFT_SEARCH modifies its host value.

      ! Also update GN, GSTART, GSTACK as they are modified by SOFT_SEARCH on host.

      ! It's good practice to update all relevant variables that might have changed on host before a device kernel call.
      ! If GMAXN, GN, GSTART, GSTACK are used in device kernels within HANDELER after SOFT_SEARCH:
      ! !$ACC UPDATE DEVICE(GMAXN, GN, GSTART, GSTACK) 
      ! Otherwise, their values will be passed as arguments to CONSTRUCT_FINT.



      WRITE(*,*) 'DEBUG: HANDELER - After UPDATE DEVICE(GMAXN), Host GMAXN (should be unchanged) = ', GMAXN
      
      ! Verify by copying back from device
      ! This debug step might be misleading if GMAXN is not explicitly on device
      ! or not in a data region. For now, let's keep it for user's debug flow.
      ! If GMAXN is not on device, this UPDATE HOST might fetch an uninitialized value or cause an error.
      ! Assuming GMAXN might be on device due to being an argument to a previous (hypothetical) kernel or if it's part of an implicit data environment.
      ! The user's plan removes the problematic UPDATE HOST here, which is good.




      temp_gmaxn_device_check = GMAXN
      WRITE(*,*) 'DEBUG: HANDELER - GMAXN value on DEVICE (read back to host via GMAXN) = ', temp_gmaxn_device_check
      ! Restore host GMAXN to device if it was overwritten by UPDATE HOST for checking (original logic)
      ! This step ensures the device has the host's intended value if the UPDATE HOST was just for checking.
      ! !$ACC UPDATE DEVICE(GMAXN) IF (PRESENT(GMAXN))



                            
        !DEALLOCATE(ISPACE,JSPACE,KSPACE,NODES_IN_BIN,NODELIST_IN_BIN)       
          
        !!!CONTINUE
      END IF !CNT_SEARCH=MAX_CNT_SEARCH
END IF	  
      END IF
      ierr_fint = 0    


      ! 確保所有必要的資料都在 GPU 上
      !$ACC UPDATE DEVICE(GCOO, GCOO_CUURENT, GWIN, GVOL)
      !$ACC UPDATE DEVICE(GN, GSTART, GSTACK)
      !$ACC WAIT

      !WRITE(*,*) 'DEBUG: HANDELER - Before calling CONSTRUCT_FINT/PD:' 

      WRITE(*,*) '  LINIT = ', LINIT
      !$ACC UPDATE HOST(GMAXN, DIM_NN_LIST) ! Get current device values for check
      ! This UPDATE HOST for GMAXN and DIM_NN_LIST might be problematic if they are not consistently on device
      ! or if their device values are not what we expect here.
      ! It's better to rely on the host values that were just set/calculated.
      ! The GMAXN here is the one potentially modified by SOFT_SEARCH.


      WRITE(*,*) '  GNUMP = ', GNUMP
      WRITE(*,*) '  GMAXN = ', GMAXN
      WRITE(*,*) '  DIM_NN_LIST = ', DIM_NN_LIST
      ! 確保關鍵變量有效 (as per user request).
      ! Check GMAXN (which will be GMAXN_IN for CONSTRUCT_FINT)
      IF (GMAXN > 2000 .OR. GMAXN <= 0) THEN ! Added GMAXN <= 0 check
          WRITE(*,*) 'WARNING: HANDELER - GMAXN (value: ', GMAXN, ') out of range (1-2000), resetting to 2000 (or GNUMP if smaller).'
          GMAXN = MIN(2000, MAX(1, GNUMP)) ! Ensure GMAXN is at least 1 and not more than GNUMP if GNUMP < 2000
      END IF

      IF (GMAXN <= 0) THEN
          WRITE(*,*) 'WARNING: HANDELER - GMAXN invalid before CONSTRUCT_FINT call, resetting to GNUMP'
          GMAXN = GNUMP
          !$ACC UPDATE DEVICE(GMAXN)

          ! If GMAXN is used in other device kernels within HANDELER before CONSTRUCT_FINT, update it:
          ! !$ACC UPDATE DEVICE(GMAXN) IF (PRESENT(GMAXN))
      END IF

      IF (DIM_NN_LIST <= 0) THEN
          WRITE(*,*) 'WARNING: HANDELER - DIM_NN_LIST invalid before CONSTRUCT_FINT call, resetting'
          DIM_NN_LIST = GNUMP * 1000
          !$ACC UPDATE DEVICE(DIM_NN_LIST)

          ! If DIM_NN_LIST is used in other device kernels within HANDELER before CONSTRUCT_FINT, update it:
          ! !$ACC UPDATE DEVICE(DIM_NN_LIST) IF (PRESENT(DIM_NN_LIST))
      END IF

      WRITE(*,*) 'DEBUG: HANDELER - Final values before CONSTRUCT_FINT call (after potential reset):'
      WRITE(*,*) '  GMAXN = ', GMAXN
      WRITE(*,*) '  DIM_NN_LIST = ', DIM_NN_LIST

      WRITE(*,*) '  PERIDYNAMICS = ', PERIDYNAMICS
      IF (ALLOCATED(GN)) THEN
            WRITE(*,*) '  ALLOCATED(GN) = .TRUE., SIZE(GN,1) = ', SIZE(GN,1)
      ELSE
            WRITE(*,*) '  ALLOCATED(GN) = .FALSE.'
      END IF
      IF (ALLOCATED(GSTACK)) THEN
            WRITE(*,*) '  ALLOCATED(GSTACK) = .TRUE., SIZE(GSTACK,1) = ', SIZE(GSTACK,1)
      ELSE
            WRITE(*,*) '  ALLOCATED(GSTACK) = .FALSE.'
      END IF
      IF (GNUMP <= 0) THEN
          WRITE(*,*) 'ERROR: HANDELER - GNUMP <= 0 before calling CONSTRUCT_FINT/PD. GNUMP = ', GNUMP
          H_IERR = -41 ! Specific error code
          RETURN
      END IF
IF (.NOT. ALLOCATED(GN)) THEN
    H_IERR = -42
    WRITE(*,*) 'ERROR: HANDELER - GN not allocated.'
    RETURN
END IF

IF (.NOT. ALLOCATED(GSTACK)) THEN
    H_IERR = -43
    WRITE(*,*) 'ERROR: HANDELER - GSTACK not allocated.'
    RETURN
END IF

      ! 在調用 CONSTRUCT_FINT 之前，確保 FEXT 陣列在設備上是完整的
      ! Ensure GFEXT and GFINT are clean on the device before CONSTRUCT_FINT computes them.

      ! 檢查 GFEXT 陣列的有效性 (GFEXT 是 CONSTRUCT_FINT 的虛擬參數，對應 MAIN 中的 LOCAL_FEXT)
      ! 這些檢查主要應在 MAIN 中進行，因為 HANDELER 中的 GFEXT 是虛擬參數，其分配狀態由實參決定。
      ! 但如果 GSIZE 是 HANDELER 內部計算的，並且與 GFEXT 的預期大小相關，則此處檢查 GSIZE。
      IF (GSIZE <= 0) THEN
          WRITE(*,*) 'ERROR: HANDELER - GSIZE <= 0, cannot determine valid size for GFEXT. GSIZE = ', GSIZE
          H_IERR = -52 
          RETURN
      END IF
      !$ACC PARALLEL LOOP DEFAULT(PRESENT)
      DO I = 1, 3*GSIZE
          GFINT(I) = 0.0D0
          GFEXT(I) = 0.0D0
      END DO
      !$ACC END PARALLEL LOOP
      
      ! 驗證一致性的檢查點
      IF (LINIT) THEN
          ! 在 SOFT_SEARCH 後檢查結果
          !$ACC UPDATE HOST(GN(1:MIN(10,GNUMP)), GMAXN)
      IF (LINIT) THEN
          !$ACC UPDATE HOST(GN, GSTART, GSTACK, GSTACK_SHP, GSTACK_DSHP)
          
          WRITE(*,*) 'DEBUG: First 10 GN values: ', GN(1:MIN(10,GNUMP))
          WRITE(*,*) 'DEBUG: GMAXN after search: ', GMAXN
          
          ! 檢查是否有節點沒有鄰居
          DO I = 1, MIN(10, GNUMP)
              IF (GN(I) <= 0) THEN
                  WRITE(*,*) 'WARNING: Node ', I, ' has no neighbors!'
              END IF
          END DO
          
          IF (GNUMP > 0 .AND. GN(1) > 0) THEN
              WRITE(*,*) 'DEBUG: Node 1 neighbors: ', GSTACK(GSTART(1):GSTART(1)+GN(1)-1)
              
              ! 檢查 shape function 是否已計算
              WRITE(*,*) 'DEBUG: First 10 SHP values for node 1: ', &
                        GSTACK_SHP(GSTART(1):MIN(GSTART(1)+9, GSTART(1)+GN(1)-1))
          END IF
          
          ! 確保所有資料都同步到 GPU
          !$ACC UPDATE DEVICE(GN, GSTART, GSTACK, GSTACK_SHP, GSTACK_DSHP, GSTACK_DDSHP)
      END IF
      END IF

      IF (PERIDYNAMICS) THEN

      !GET THE INTERNAL FORCE
      CALL CONSTRUCT_FINT_PD(GWIN,    GVOL,        GNUMP,     GCOO, GCOO_CUURENT, &   !FROM MAIN
                            GSM_LEN, GSM_AREA,    GSM_VOL,   GNSNI_FAC,          &   !FROM MAIN
                            GGHOST,  GPROP,       GSTATE,    GSTRESS,            &   !FROM MAIN
                            GSTRAIN, G_H_STRESS, G_S_STRESS, GDINC,  GDINC_TOT,     GMAT_TYPE,                 &   !FROM MAIN
                            GEBC_NODES,                                         &   !FROM MAIN
                            GN, GSTART,  DIM_NN_LIST, GSTACK,                        &   !FROM HANDLER
                          GSTACK_SHP,  GSTACK_DSHP,  GSTACK_DDSHP,  GMAXN, GINVK, LINIT,          &   !FROM HANDLER ! GMAXN, DIM_NN_LIST
                            GFINT, DLT_FINT,   GCHAR_DIST,   GMAX_WVEL, &   !delet GFEXT first (0627)
                            LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
                            G_X_MOM, G_Y_MOM, G_Z_MOM, MODEL_BODYFORCE, GINT_WORK,DLT, ierr_fint) !OUTPUT
      IF (ierr_fint /= 0) THEN
          WRITE(*,*) 'ERROR: HANDELER - Error returned from CONSTRUCT_FINT_PD. ierr_fint = ', ierr_fint
          H_IERR = ierr_fint
          RETURN
      END IF      
      H_IERR = ierr_fint                            

	  ELSE      
      ! 在呼叫 CONSTRUCT_FINT 前，先把 H_IERR 初始化
      H_IERR = 0


      ! Comprehensive !$ACC DATA region for CONSTRUCT_FINT call
      ! Consolidate all COPYIN and COPYOUT clauses into a single DATA directive
      ! Wenn HANDELER von MAIN.F90's Datenregion aufgerufen wird, ist diese Datenregion redundant.
      ! Die Daten werden bereits auf dem Gerät verwaltet.
      ! !$ACC DATA COPYIN(GWIN, GVOL, GNUMP, GCOO, GCOO_CUURENT, &
      ! !$ACC             GSM_LEN, GSM_AREA, GSM_VOL, GNSNI_FAC, &
      ! !$ACC             GGHOST, GPROP, GSTATE, GSTRESS, GSTRAIN, G_H_STRESS, G_S_STRESS, &
      ! !$ACC             GDINC, GDINC_TOT, GMAT_TYPE, GEBC_NODES, DLT, GSIZE, &
      ! !$ACC             GN, GSTART, DIM_NN_LIST, GSTACK, GSTACK_SHP, GSTACK_DSHP, GSTACK_DDSHP, &
      ! !$ACC             GMAXN, GINVK, LINIT, &
      ! !$ACC             G_X_MOM, G_Y_MOM, G_Z_MOM, MODEL_BODYFORCE, MODEL_BODY_ID) &
      ! !$ACC COPYOUT(GFINT, GFEXT, GINT_WORK, GSTRAIN_EQ, DLT_FINT, GCHAR_DIST, GMAX_WVEL, &
      ! !$ACC           LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
      ! !$ACC           ierr_fint)

      ! Note: GVOL is GNUMP-sized in HANDELER. CONSTRUCT_FINT uses GVOL(I) for I=1,GNUMP.
      ! GSTRAIN_EQ, DLT_FINT, GCHAR_DIST, GMAX_WVEL, LOCAL_*_STRESS are outputs from CONSTRUCT_FINT.


      !GET THE INTERNAL FORCE
      CALL CONSTRUCT_FINT(GWIN,    GVOL,        GNUMP,     GCOO, GCOO_CUURENT, &
                          GSM_LEN, GSM_AREA,    GSM_VOL,   GNSNI_FAC,          &
                          GGHOST,  GPROP,       GSTATE,    GSTRESS,            &
                          GSTRAIN, G_H_STRESS, G_S_STRESS, GDINC,  GDINC_TOT,     GMAT_TYPE,                 &
                          GEBC_NODES, DLT, GSIZE,                                & ! DLT, GSIZE (passed as GSIZE_IN)
                          GN, GSTART,  DIM_NN_LIST, GSTACK,                        &
                          GSTACK_SHP,  GSTACK_DSHP,  GSTACK_DDSHP,  GMAXN, GINVK, LINIT,          &
                          GFINT, GFEXT, DLT_FINT,   GCHAR_DIST,   GMAX_WVEL, &
                          LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
                          G_X_MOM, G_Y_MOM, G_Z_MOM, MODEL_BODYFORCE, GINT_WORK, MODEL_BODY_ID, GSTRAIN_EQ, &
                          ierr_fint)
          ! !$ACC END DATA




      H_IERR = ierr_fint ! Assign the error code from CONSTRUCT_FINT to H_IERR
      ! Now check H_IERR (which holds the result from ierr_fint)
      IF (H_IERR /= 0) THEN
          WRITE(*,*) 'ERROR: HANDELER - CONSTRUCT_FINT 返回錯誤，H_IERR = ', H_IERR
          RETURN
      END IF

	  END IF
	  
      END IF

    !UPDATE THE PHYSICAL DISPLACEMENT VALUES (?)
      
      RETURN
      
      END SUBROUTINE
                    