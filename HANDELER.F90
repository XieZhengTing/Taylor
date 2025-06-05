
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
      
      IF (DO_INTERP) THEN
!        !$ACC PARALLEL LOOP DEFAULT(PRESENT) &
!        !$ACC PRIVATE(LSTART, K, M, J, SHPT, MM)
		IF(.NOT. PERIDYNAMICS) THEN  !RKPM  
            !$ACC PARALLEL LOOP DEFAULT(PRESENT) &
            !$ACC PRIVATE(LSTART, K, M, J, SHPT, MM)
            DO I=1,GNUMP
                LSTART = GSTART(I)
                DO K = 1, 3
                    M=(I-1)*3+K
                    GDINC_PHY(M) = 0.0d0
                    GVEL_PHY(M) = 0.0d0
                    GACL_PHY(M) = 0.0d0
                    DO J = 1, GN(I)
                        SHPT=GSTACK_SHP(LSTART+J-1)
                        
                        MM = (GSTACK(LSTART+J-1) - 1)*3 + K
                        
                        GDINC_PHY(M) = GDINC_PHY(M) +  SHPT*GDINC(MM)
                        GVEL_PHY(M) = GVEL_PHY(M) +  SHPT*GVEL(MM)
                        GACL_PHY(M) = GACL_PHY(M) +  SHPT*GACL(MM)
                    END DO
                END DO

            END DO
            !$ACC END PARALLEL LOOP
        ELSE  !PERIDYNAMICS
            GDINC_PHY = GDINC
            GVEL_PHY = GVEL
            GACL_PHY = GACL
        
        ENDIF
				
				
      ELSE
      
      !CALL THE NODE SEARCH ROUTINE HERE IF NEEDED
      IF ((LINIT).OR.(.NOT.LLAGRANGIAN)) THEN
      H_IERR = 0 ! Initialize H_IERR
      IF (GNUMP <= 0) THEN
          WRITE(*,*) 'ERROR: HANDELER - GNUMP <= 0 before calculating DIM_NN_LIST. GNUMP = ', GNUMP
          H_IERR = -40 

          RETURN
      END IF    
      WRITE(*,*) 'DEBUG: HANDELER - GNUMP = ', GNUMP, ' (before DIM_NN_LIST calculation)'
     


      ! GMAXN is determined by search routines later.
      ! The GMAXN_IN passed to CONSTRUCT_FINT will be this HANDELER's GMAXN.
      ! We will check GMAXN after search routines and before CONSTRUCT_FINT call.
      ! For now, ensure GMAXN (if it has a preliminary value, e.g. from a previous step or LINIT)
      ! is reasonable if it were to be used before search.
      ! However, GMAXN is typically an *output* of search routines.
      ! The example suggests checking GMAXN_IN, which is the GMAXN from HANDELER passed to CONSTRUCT_FINT.
  
     
      ! Calculate and limit DIM_NN_LIST *before* SOFT_SEARCH or any allocation depending on it.
      temp_dim_long = INT(GNUMP, KIND=8) * 1000_8   ! Default: GNUMP*1000
      IF (temp_dim_long > 500000_8) THEN           ! Cap at 500,000
          temp_dim_long = 500000_8
      END IF
      ! Ensure at least GNUMP*10
      temp_dim_long = MAX(temp_dim_long, INT(GNUMP, KIND=8)*10_8)
      DIM_NN_LIST = INT(temp_dim_long, KIND=4) ! Assign back to INTEGER (default kind)
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
          ! DIM_NN_LIST is set above, ensure it's valid before allocating GSTACK
          IF (DIM_NN_LIST <= 0) THEN 
              DIM_NN_LIST = GNUMP * 1000 ! Fallback
          END IF

    
          ALLOCATE(GSTACK(DIM_NN_LIST))                
          ALLOCATE(GSTACK_SHP(DIM_NN_LIST))                    
          ALLOCATE(GSTACK_DSHP(3,DIM_NN_LIST))  
          ALLOCATE(GSTACK_DDSHP(6,DIM_NN_LIST)) 
          ALLOCATE(GINVK(3,3,GNUMP))
          !$ACC ENTER DATA CREATE(GN, GSTART, GSTACK, GSTACK_SHP, GSTACK_DSHP, GSTACK_DDSHP, GINVK)

          ! For newly allocated arrays that are DECLARE CREATE'd, their device counterparts are also allocated.
          ! If they need initial values from host, an UPDATE DEVICE would be needed after host initialization.
          ! Here, they are just allocated.
          ! Without DECLARE CREATE, their data will be copied to device when passed as arguments to ACC ROUTINEs.

          CNT_SEARCH = 0.0d0

          SEARCHCOUNT=PDSEARCH
          ! If CNT_SEARCH or SEARCHCOUNT are used in device kernels within HANDELER,
          ! an !$ACC UPDATE DEVICE(CNT_SEARCH, SEARCHCOUNT) would be needed.


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
         
		 
		 

         IF (SEARCHCOUNT.LT.PDSEARCH) THEN              
             SEARCHCOUNT=SEARCHCOUNT+1

         ELSE             
         CALL GET_BINS(GNUMP,GCOO_CUURENT,NODES_IN_BIN,MAX_NEIGH,NODELIST_IN_BIN, &
                       NBINS,NBINSX,NBINSY,NBINSZ,ISPACE,JSPACE,KSPACE,XMIN,YMIN,ZMIN,XBIN,YBIN,ZBIN)
         ! If ISPACE, JSPACE, etc. are used in device kernels within HANDELER after GET_BINS:
         ! !$ACC UPDATE DEVICE(ISPACE, JSPACE, KSPACE, NODES_IN_BIN, NODELIST_IN_BIN) 


         DO I=1,GNUMP
         GIJKSPC(1,I) = ISPACE(I)
         GIJKSPC(2,I) = JSPACE(I)
         GIJKSPC(3,I) = KSPACE(I)
		 END DO
		 
         SEARCHCOUNT=1


         ENDIF
                
         CALL SOFT_SEARCH(GNUMP,GCOO_CUURENT,GN,GSTART,DIM_NN_LIST,GSTACK,GMAXN,GWIN, &
                            NODES_IN_BIN,MAX_NEIGH,NODELIST_IN_BIN, &
                            NBINS,NBINSX,NBINSY,NBINSZ,ISPACE,JSPACE,KSPACE, &
                            GXDIST_MAX, GYDIST_MAX, GZDIST_MAX,GSM_LEN)
      !$ACC UPDATE DEVICE(GMAXN, GN, GSTART, GSTACK)


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
      ierr_fint = 0    


      WRITE(*,*) 'DEBUG: HANDELER - Before calling CONSTRUCT_FINT/PD:'
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
      !$ACC DATA COPYIN(GWIN, GVOL, GNUMP, GCOO, GCOO_CUURENT, & ! Line 569 or near

      !$ACC             GSM_LEN, GSM_AREA, GSM_VOL, GNSNI_FAC, &
      !$ACC             GGHOST, GPROP, GSTATE, GSTRESS, GSTRAIN, G_H_STRESS, G_S_STRESS, &
      !$ACC             GDINC, GDINC_TOT, GMAT_TYPE, GEBC_NODES, DLT, GSIZE, &
      !$ACC             GN, GSTART, DIM_NN_LIST, GSTACK, GSTACK_SHP, GSTACK_DSHP, GSTACK_DDSHP, &
      !$ACC             GMAXN, GINVK, LINIT, &
      !$ACC             G_X_MOM, G_Y_MOM, G_Z_MOM, MODEL_BODYFORCE, MODEL_BODY_ID) & ! End of COPYIN list
      !$ACC COPYOUT(GFINT, GFEXT, GINT_WORK, GSTRAIN_EQ, DLT_FINT, GCHAR_DIST, GMAX_WVEL, &
      !$ACC           LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
      !$ACC           ierr_fint)
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
      !$ACC END DATA



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
                    