
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
						  GIJKSPC)
                          
      
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
	  INTEGER,ALLOCATABLE, SAVE:: GSTART(:)                 !START LOCATION OF NODE NEIOGHBORS IN STACK
	  INTEGER, SAVE:: DIM_NN_LIST                           !SIZE OF NEIGHBOR STACK
	  INTEGER,ALLOCATABLE, SAVE:: GSTACK(:)                 !NEIGHBORS FOR EACH NODE (STACKED)
	  DOUBLE PRECISION,ALLOCATABLE, SAVE:: GSTACK_SHP(:)                 !NEIGHBORS FOR EACH NODE (STACKED)
	  DOUBLE PRECISION,ALLOCATABLE, SAVE:: GSTACK_DSHP(:,:)                 !NEIGHBORS FOR EACH NODE (STACKED)
	  DOUBLE PRECISION,ALLOCATABLE, SAVE:: GSTACK_DDSHP(:,:)                 !NEIGHBORS FOR EACH NODE (STACKED)
	  DOUBLE PRECISION,ALLOCATABLE, SAVE:: GINVK(:,:,:)      
	  INTEGER, SAVE:: GMAXN                                 !MAX NUMBER OF NEIGHBORS FOR ALL NODES
	  INTEGER, SAVE:: SEARCHCOUNT
      
      !BINNING INFORMATION
      
      DOUBLE PRECISION,SAVE:: XBIN,YBIN,ZBIN,XMAX,YMAX,ZMAX,XMIN,YMIN,ZMIN,XBIN_TEMP,YBIN_TEMP,ZBIN_TEMP
      DOUBLE PRECISION,SAVE::DOM_XLENGTH,DOM_YLENGTH,DOM_ZLENGTH
      INTEGER,SAVE::NBINSX,NBINSY,NBINSZ,NBINS,MAX_NEIGH
      INTEGER,SAVE,ALLOCATABLE::ISPACE(:),JSPACE(:),KSPACE(:),NODES_IN_BIN(:),NODELIST_IN_BIN(:,:)
         
      
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
      DOUBLE PRECISION, SAVE:: CNT_SEARCH
      
      IF (DO_INTERP) THEN
		IF(.NOT. PERIDYNAMICS) THEN  !RKPM

       ! CRITICAL: Ensure neighbor data is current on GPU
       ! This is necessary because initial CREATE may leave data uninitialized
        ! Only update neighbor structure, NOT shape functions (they come from GPU)
        ! ONLY update neighbor connectivity, never shape functions
        !$ACC UPDATE DEVICE(GN, GSTART, GSTACK)
        ! DO NOT include GSTACK_SHP, GSTACK_DSHP, GSTACK_DDSHP here!
        
        ! Check if shape functions are valid before interpolation
        IF (GSTACK_SHP(GSTART(1)) .EQ. 0.0d0) THEN
            PRINT *, 'WARNING: Shape functions are zero before interpolation!'
            PRINT *, 'This may indicate GPU sync issue'
        END IF       
       ! Debug: Verify neighbor data before interpolation
       PRINT *, '=== NEIGHBOR DATA CHECK ==='
       PRINT *, 'First 3 nodes neighbor count (GN):', GN(1:3)       
       IF (GN(1) > 0) THEN
           PRINT *, 'Node 1 first neighbor:', GSTACK(GSTART(1))
           PRINT *, 'Node 1 first shape function:', GSTACK_SHP(GSTART(1))
       END IF

          !$ACC PARALLEL LOOP PRESENT(GDINC_PHY, GVEL_PHY, GACL_PHY, &
          !$ACC&                      GSTACK_SHP, GSTACK_DSHP, GSTACK_DDSHP, GSTACK, GSTART, GN, &
          !$ACC&                      GDINC, GVEL, GACL) &
        !$ACC&              PRIVATE(LSTART, M, MM, SHPT)
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
        !$ACC UPDATE HOST(GDINC_PHY, GVEL_PHY, GACL_PHY)
        ELSE  !PERIDYNAMICS
            !$ACC PARALLEL LOOP PRESENT(GDINC_PHY, GVEL_PHY, GACL_PHY, &
            !$ACC&                      GDINC, GVEL, GACL)
            DO I = 1, 3*GNUMP
                GDINC_PHY(I) = GDINC(I)
                GVEL_PHY(I) = GVEL(I)
                GACL_PHY(I) = GACL(I)
            END DO
            !$ACC END PARALLEL LOOP
            !$ACC UPDATE HOST(GDINC_PHY, GVEL_PHY, GACL_PHY)
        
        ENDIF
				
				
      ELSE
      
      !CALL THE NODE SEARCH ROUTINE HERE IF NEEDED
      IF ((LINIT).OR.(.NOT.LLAGRANGIAN)) THEN
      
        DIM_NN_LIST=GNUMP*1000
        GMAXN=GNUMP

        !$ACC ENTER DATA COPYIN(DIM_NN_LIST, GMAXN)
      
      IF (LINIT) THEN
      ALLOCATE(GN(GNUMP)) 
	  ALLOCATE(GSTART(GNUMP))         
	  ALLOCATE(GSTACK(DIM_NN_LIST))                
	  ALLOCATE(GSTACK_SHP(DIM_NN_LIST))                    
	  ALLOCATE(GSTACK_DSHP(3,DIM_NN_LIST))  
	  ALLOCATE(GSTACK_DDSHP(6,DIM_NN_LIST)) 
      AlLOCATE(GINVK(3,3,GNUMP))

      ! OpenACC: Create GPU data for neighbor lists
      !
   ! Initialize to zero to avoid accessing random values
   GSTACK_SHP = 0.0d0
   GSTACK_DSHP = 0.0d0
   GSTACK_DDSHP = 0.0d0
      ! Use COPYIN to ensure GPU gets initialized values
      !$ACC ENTER DATA COPYIN(GN, GSTART, GSTACK, GSTACK_SHP, GSTACK_DSHP, GSTACK_DDSHP, GINVK)
	  CNT_SEARCH = 0.0d0
      SEARCHCOUNT=PDSEARCH
      END IF
      
      !
      ! HARD_SEARCH
      !
      ! WORKS BY OPENING THE FOLLOWING COMMENT, ESPECIALLY FOR DEBUGGING PURPOSE
      !CALL HARD_SEARCH(GNUMP,GCOO_CUURENT,GN,GSTART,DIM_NN_LIST,GSTACK,GMAXN,GWIN,GXDIST_MAX, GYDIST_MAX,GZDIST_MAX)

	  CNT_SEARCH = CNT_SEARCH + DLT
	  
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
           !$ACC ENTER DATA CREATE(ISPACE, JSPACE, KSPACE, NODES_IN_BIN, NODELIST_IN_BIN)		   
         END IF !(LINIT) 
         
		 
		 
         IF (SEARCHCOUNT.LT.PDSEARCH) THEN              
             SEARCHCOUNT=SEARCHCOUNT+1
         ELSE             
         CALL GET_BINS(GNUMP,GCOO_CUURENT,NODES_IN_BIN,MAX_NEIGH,NODELIST_IN_BIN, &
                       NBINS,NBINSX,NBINSY,NBINSZ,ISPACE,JSPACE,KSPACE,XMIN,YMIN,ZMIN,XBIN,YBIN,ZBIN)
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

         ! Update neighbor lists on GPU after search
         !
         !$ACC UPDATE DEVICE(GN, GSTART, GSTACK, GSTACK_SHP, GSTACK_DSHP, GSTACK_DDSHP)
         !$ACC UPDATE DEVICE(DIM_NN_LIST, GMAXN)                       
                            
        !DEALLOCATE(ISPACE,JSPACE,KSPACE,NODES_IN_BIN,NODELIST_IN_BIN)       
          
        !!!CONTINUE
      END IF !CNT_SEARCH=MAX_CNT_SEARCH
	  
      END IF
      
      IF (PERIDYNAMICS) THEN
      !GET THE INTERNAL FORCE
      CALL CONSTRUCT_FINT_PD(GWIN,    GVOL,        GNUMP,     GCOO, GCOO_CUURENT, &   !FROM MAIN
                            GSM_LEN, GSM_AREA,    GSM_VOL,   GNSNI_FAC,          &   !FROM MAIN
                            GGHOST,  GPROP,       GSTATE,    GSTRESS,            &   !FROM MAIN
                            GSTRAIN, G_H_STRESS, G_S_STRESS, GDINC,  GDINC_TOT,     GMAT_TYPE,                 &   !FROM MAIN
                            GEBC_NODES,                                         &   !FROM MAIN
                            GN, GSTART,  DIM_NN_LIST, GSTACK,                        &   !FROM HANDLER
                            GSTACK_SHP,  GSTACK_DSHP,  GSTACK_DDSHP,  GMAXN, GINVK, LINIT,          &   !FROM HANDLER
                            GFINT, DLT_FINT,   GCHAR_DIST,   GMAX_WVEL, &   !delet GFEXT first (0627)
                            LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
                            G_X_MOM, G_Y_MOM, G_Z_MOM, MODEL_BODYFORCE, GINT_WORK,DLT)              !OUTPUT
	  ELSE
     ! Diagnostic: Verify GDINC_TOT before GPU computation
     IF (LINIT .OR. ITYPE_INT .EQ. 0) THEN
         PRINT *, '=== HANDELER: GDINC_TOT CHECK ==='
         PRINT *, 'First 3 values:', GDINC_TOT(1:3)
         PRINT *, 'Max value:', MAXVAL(ABS(GDINC_TOT))
     END IF
      !GET THE INTERNAL FORCE
      CALL CONSTRUCT_FINT(GWIN,    GVOL,        GNUMP,     GCOO, GCOO_CUURENT, &   !FROM MAIN
                          GSM_LEN, GSM_AREA,    GSM_VOL,   GNSNI_FAC,          &   !FROM MAIN
                          GGHOST,  GPROP,       GSTATE,    GSTRESS,            &   !FROM MAIN
                          GSTRAIN, G_H_STRESS, G_S_STRESS, GDINC,  GDINC_TOT,     GMAT_TYPE,                 &   !FROM MAIN
                          GEBC_NODES,                                         &   !FROM MAIN
                          GN, GSTART,  DIM_NN_LIST, GSTACK,                        &   !FROM HANDLER
                          GSTACK_SHP,  GSTACK_DSHP,  GSTACK_DDSHP,  GMAXN, GINVK, LINIT,          &   !FROM HANDLER
                          GFINT, GFEXT, DLT_FINT,   GCHAR_DIST,   GMAX_WVEL, &
                          LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
                          G_X_MOM, G_Y_MOM, G_Z_MOM, MODEL_BODYFORCE, GINT_WORK, MODEL_BODY_ID, GSTRAIN_EQ, DLT)              !OUTPUT
	  END IF
	  
      ! Shape functions are already synced in CONSTRUCT_FINT
      ! Do NOT update device here as it would overwrite GPU-computed values
      
      ! Debug: Verify shape functions after FINT computation
      IF (LINIT) THEN
          PRINT *, '=== POST-FINT SHAPE FUNCTION CHECK ==='
          PRINT *, 'GSTACK_SHP(first 5):', GSTACK_SHP(1:MIN(5, DIM_NN_LIST))
          IF (MAXVAL(ABS(GSTACK_SHP(1:MIN(100, DIM_NN_LIST)))) .EQ. 0.0d0) THEN
              PRINT *, 'ERROR: Shape functions are still zero after FINT!'
          END IF
      END IF
      
      END IF
      
      !UPDATE THE PHYSICAL DISPLACEMENT VALUES (?)
      
      RETURN
      
      END SUBROUTINE
                    