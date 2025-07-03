

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
      SUBROUTINE CONSTRUCT_FINT_PD(GWIN,    GVOL,       GNUMP,        GCOO,  GCOO_CUURENT,   &   !FROM MAIN
                                GSM_LEN, GSM_AREA,   GSM_VOL,      GNSNI_FAC,      &   !FROM MAIN
                                GGHOST,  GPROP,      GSTATE,       GSTRESS,        &   !FROM MAIN
                                GSTRAIN, G_H_STRESS, G_S_STRESS, GDINC,      GDINC_TOT,    GMAT_TYPE,                    &   !FROM MAIN
                                GEBC_NODES,                                               &   !FROM MAIN
                                GN,      GSTART,     DIM_NN_LIST,                  &   !FROM HANDLER
                                GSTACK,  GSTACK_SHP, GSTACK_DSHP,  GSTACK_DDSHP,  GMAXN,  GINVK, LINIT,   &   !FROM HANDLER
                                FINT,    DLT_FINT,   GCHAR_DIST,   GMAX_WVEL, &
                                LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, & !OUTPUT
                                G_X_MOM, G_Y_MOM, G_Z_MOM,MODEL_BODYFORCE,GINT_WORK,DLT)           
      
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! BASED ON CURRENT DISPLACEMENT INCREMENT, STATE VARIABLES, AND SO ON, 
	  ! COMPUTE THE INTERNAL FORCE FOR A SET OF NODES GIVEN THEIR COORDINATES,
	  ! DILATIONS, AND SO ON.
	  !
	  !VOIGT ORDERING:
      !XX, YY, ZZ, YZ, XZ, XY
	  !
      USE FINT_FUNCTIONS
      USE CONTROL
      USE omp_lib
      !
	  IMPLICIT NONE
	  !
	  !
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
	  INTEGER:: GNUMP                         !NUMBER OF NODES
	  DOUBLE PRECISION:: GCOO(3,GNUMP)        !COORDINATE FOR EACH NODE
	  DOUBLE PRECISION:: GCOO_CUURENT(3,GNUMP) !CUURENT COORDINATE FOR EACH NODE
	  DOUBLE PRECISION:: GWIN(3,GNUMP)        !WINDOWS FOR EACH NODE
	  DOUBLE PRECISION:: GSM_LEN(6,GNUMP)     !SMOOTHING LENGTHS FOR EACH NODE
	    ! (1-6) = (+X, -X, +Y, -Y, +Z, -Z)
	  DOUBLE PRECISION:: GSM_VOL(GNUMP)       !VOLUME OF SMOOTHING ZONE
	  DOUBLE PRECISION:: GSM_AREA(3,GNUMP)    !AREAS OF SMOOTHING ZONE SIDES
	  INTEGER:: GN(GNUMP)                     !NUMBER OF NEIGHBORS FOR EACH NODE
	  INTEGER:: GSTART(GNUMP)                 !START LOCATION OF NODE NEIOGHBORS IN STACK
	  INTEGER:: DIM_NN_LIST                   !SIZE OF NEIGHBOR STACK
	  INTEGER:: GSTACK(DIM_NN_LIST)           !NEIGHBORS FOR EACH NODE (STACKED)
	  DOUBLE PRECISION:: GSTACK_SHP(DIM_NN_LIST)       !SHAPES (STACKED)
	  DOUBLE PRECISION:: GSTACK_DSHP(3,DIM_NN_LIST)       !SHAPES (STACKED)
	  DOUBLE PRECISION:: GSTACK_DDSHP(6,DIM_NN_LIST)       !SHAPES (STACKED)
	  DOUBLE PRECISION:: GINVK(3,3,GNUMP) 
      
	  DOUBLE PRECISION:: GCHAR_DIST(GNUMP),   GMAX_WVEL(GNUMP)
	  
	  INTEGER:: GMAXN                             !MAX NUMBER OF NEIGHBORS FOR ALL NODES
	  INTEGER:: GGHOST(GNUMP)                 !FLAG FOR GHOST NODES (GHOST = 1)
      
      LOGICAL:: GEBC_NODES(GNUMP)
      
	  DOUBLE PRECISION:: GVOL(GNUMP)          !VOLUME OF EACH NODE
	  DOUBLE PRECISION:: GNSNI_FAC(3,GNUMP)
      
      DOUBLE PRECISION::DLT !GC
      
	  !
	  !STATE AND FIELD VARIABLES
	  !
	  DOUBLE PRECISION:: GSTRESS(6,GNUMP)     !CAUCHY STRESS OF EACH NODE
	  DOUBLE PRECISION:: GSTRAIN(6,GNUMP)     !CAUCHY STRESS OF EACH NODE
	  DOUBLE PRECISION:: GSTATE(20,GNUMP)     !STATE VARIABLES OF EACH NODE
	  DOUBLE PRECISION:: GPROP(30,GNUMP)     !MATERIAL PROPERTIES OF EACH NODE
	  DOUBLE PRECISION:: GDINC(3*GNUMP) ,GDINC_TOT(3*GNUMP)      !DISPLACEMENT INCREMENT (PREDICTOR) OF EACH NODE
	  INTEGER::          GMAT_TYPE(GNUMP)    !MATERIAL TYPE OF EACH NODE
	  !
	  DOUBLE PRECISION:: G_H_STRESS(6,GNUMP)!GC
	  DOUBLE PRECISION:: G_S_STRESS(6,GNUMP)!GC

	  !GLOBAL OUT
	  DOUBLE PRECISION:: FINT(GNUMP*3)
	  DOUBLE PRECISION:: DLT_FINT
	  !
	  !********** LOCAL VARIABLES **********
	  INTEGER:: I,J,K,L,JJ               !INDICIES
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


      DOUBLE PRECISION:: G_X_MOM(GNUMP),G_Y_MOM(GNUMP),G_Z_MOM(GNUMP)
      
      
	  DOUBLE PRECISION:: LSTATE(20)      !STATE VARIABLES OF A NODE
	  DOUBLE PRECISION:: LPROP(30)     !MATERIAL PROPERTIES OF A NODE
	  INTEGER:: LMAT_TYPE    !MATERIAL TYPE OF EACH NODE
      LOGICAL:: SELF_EBC
	  !
	  INTEGER:: LSTART              !LOCATION OF START IN STACK
	  INTEGER:: LGHOST              !FLAG FOR GHOST NODES (GHOST = 1)
	  INTEGER:: LSTACK(GMAXN)       !LOCAL LIST/STACK OF NEIGHBORS
	  INTEGER:: LN                  !NUMBER OF NEIGHBORS FOR A NODE
	  DOUBLE PRECISION:: LVOL       !VOLUME OF A NODE
	  DOUBLE PRECISION:: SHP(GMAXN), SHPD(3,GMAXN), SHPD_TRASH(3,GMAXN)       !SHAPE FUNCTIONS AND GRADIENTS
	  DOUBLE PRECISION:: SHPDTMP(3,GMAXN)       !TEMPORARY SHAPE FUNCTIONS AND GRADIENTS
	  
	  DOUBLE PRECISION:: SHPD_SM(3,GMAXN)        !SHAPE FUNCTION SMOOTHED GRADIENTS
	  DOUBLE PRECISION:: SHPDD_SM(6,GMAXN)       !SHAPE FUNCTION SMOOTHED SMOOTHED GRADIENTS
	  DOUBLE PRECISION:: SHPDTEMP(9) !TEMPORARY VARIABLE FOR SHAPES FOR STABILZIATION
	  DOUBLE PRECISION:: SHP6(GMAXN,6), SHPD6(3,GMAXN,6) !SHAPE FUNCTION AND SM. GRAD. AT SMOOTHING POINTS
      DOUBLE PRECISION:: LMAT(3,3) !INCREMENTAL DEFORMATION GRADIENT WITH RESPECT TO THE CURRENT TIME STEP      
	  DOUBLE PRECISION:: LDINC(3,GMAXN),LDINC_TOT(3,GMAXN)       !DISPLACEMENT INCREMENT (PREDICTOR) OF A NODESTRAIN
	  DOUBLE PRECISION:: LCOO_CUURENT(3,GMAXN)  !CURRENT COORDINATES OF THE NEIGBORS
	  DOUBLE PRECISION:: STRAIN(6)       !INCREMENTALLY OBJECTIVE STRAIN
	  DOUBLE PRECISION:: ELAS_MAT(6,6)
	  DOUBLE PRECISION:: BMAT(6,3)
	  DOUBLE PRECISION:: BMAT_T(3,6)
	  DOUBLE PRECISION:: FINT3(3),FINT3_J(3),INVK(3,3),INVK_J(3,3)
	  DOUBLE PRECISION:: ROT(6,6) !ROTATION MATRIX
	  LOGICAL:: LINIT
	  DOUBLE PRECISION:: FMAT(3,3), IFMAT(3,3),X_0(3),X_t(3), DX_t(3,1), PKSTRESS(3,3), TEMP_STRESS(3,3), DX_t_J(3,1)
	  DOUBLE PRECISION, ALLOCATABLE:: FINT_TEMP(:,:,:)
	  DOUBLE PRECISION:: DET     
	  LOGICAL :: INV_OK 
      !DOUBLE PRECISION:: FINT_TEMP(20,3,GNUMP)
      INTEGER:: ID_RANK
      
      !0329
      DOUBLE PRECISION:: PMAT(6,3)
	  DOUBLE PRECISION:: FBOD(3),FGRAV(3) ! TO DO
      DOUBLE PRECISION:: MODEL_BODYFORCE(3,GNUMP)
      
      !NSNI
      INTEGER::XMAP(3),YMAP(3),ZMAP(3)
      DOUBLE PRECISION::XLMAT(3,3),YLMAT(3,3),ZLMAT(3,3)
      DOUBLE PRECISION:: DX_STRAIN(6), DY_STRAIN(6), DZ_STRAIN(6)
      DOUBLE PRECISION::  LOCAL_DX_STRESS(6,GNUMP)
      DOUBLE PRECISION::  LOCAL_DY_STRESS(6,GNUMP)
      DOUBLE PRECISION::  LOCAL_DZ_STRESS(6,GNUMP)
      
      DOUBLE PRECISION::  LDX_STRESS(6)
      DOUBLE PRECISION::  LDY_STRESS(6)
      DOUBLE PRECISION::  LDZ_STRESS(6)
      DOUBLE PRECISION:: CMAT(6,6), LAMDA, MU, LAMDA_PLUS_2MU
      DOUBLE PRECISION:: XBMAT(6,3), XBMAT_T(3,6), XFINT3(3)
      DOUBLE PRECISION:: YBMAT(6,3), YBMAT_T(3,6), YFINT3(3)
      DOUBLE PRECISION:: ZBMAT(6,3), ZBMAT_T(3,6), ZFINT3(3)
      
      DOUBLE PRECISION:: TEMP_DEBUG(3)
	  DOUBLE PRECISION:: GINT_WORK !NOT USED
      DOUBLE PRECISION:: D(6)
      
      LOGICAL:: NSNI_FLAG
      !
      ! FOR TIME STEP CALCS
	  !
	  LOGICAL:: FIRST
	  DOUBLE PRECISION:: STRESS_INC(6), STRAIN_INC(6), POISS, YOUNG, BULK, SHEAR,  &
	                     STRESS_INC_DEV(6), STRESS_INC_SPHR(6), STRAIN_INC_DEV(6), STRAIN_INC_SPHR(6),  &
	                     NORM_STRESS_INC_DEV, NORM_STRESS_INC_SPHR, NORM_STRAIN_INC_DEV, NORM_STRAIN_INC_SPHR,   &
	                     SHEAR_TRIAL, BULK_TRIAL, PMOD, DENSITY, MAXMOD, MAX_VEL,   &
	                     DIST, XJ, YJ, ZJ, CHAR_DIST, DLT_TEMP, XI, YI, ZI 
						 
	  LOGICAL:: HODIV
		
	  GINT_WORK = 0.0d0
	  !
	  !
	  !
	  !
	  !*********************************************************
	  !******************** EXECUTABLE CODE ********************
	  !*********************************************************
	  !
	  HODIV=.FALSE.
	  !
	  !INITIALIZE FINT
	  !
	  !
	  !
	  !LOOP OVER THE NODE STACK
	  !
	  !LET OPEN-MP DECIDE HOW TO DO THE DO-LOOP
	  !
	  FINT = 0.0d0
      DET = 1.d0
	   !REDUCTION(+:FINT)
	  !
      ALLOCATE(FINT_TEMP(NCORES_INPUT,3,GNUMP))
      FINT_TEMP = 0.D0
      
	  
	  
	  
	  
	  
	  
            IF (LINIT) THEN
			
			!WE NEED TO GRAM ALL THE K VALUES BEFORE WE DO ANYTHING
	  
	  
      !$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED( GNUMP, GCOO, GCOO_CUURENT, GWIN, GSM_LEN, GSM_VOL, GSM_AREA, GN, GSTART, &
      !$OMP                                       DIM_NN_LIST, GSTACK, GSTACK_SHP, GSTACK_DSHP, GSTACK_DDSHP, GINVK, & 
      !$OMP                                       GCHAR_DIST,GMAX_WVEL, GMAXN, GGHOST, GEBC_NODES, GVOL, GNSNI_FAC, &
      !$OMP                                       GSTRESS, LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
      !$OMP                                       GSTRAIN, &
      !$OMP                                       GSTATE, GPROP, GDINC,GDINC_TOT, GMAT_TYPE, FINT, DLT_FINT, FINT_TEMP)       
      
      ID_RANK = OMP_get_thread_num()  !OMPJOE

      !$OMP DO   
	  DO I = 1, GNUMP
	    !
	    !
	    !GRAB NODE INFORMATION FROM LIST
		!
		LCOO(:) = GCOO(:,I)
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
		!
        IF (GEBC_NODES(I)) THEN
          SELF_EBC = .TRUE.
        ELSE
          SELF_EBC = .FALSE.
        END IF
		!
		! RECALL THE NODE STRESS AND STATE VARIABLES
		!
		! GET THE NEIGHBOR LIST
		!
		DO J = 1, LN
		  LSTACK(J) = GSTACK(LSTART+J-1)
		END DO
		
        
                    CALL UDFM_SHAPE_TENSOR(LCOO, RK_DEGREE, RK_PSIZE, RK_CONT, RK_IMPL,GCOO, GVOL, GWIN, GNUMP, LSTACK, LN, GMAXN, GEBC_NODES,SELF_EBC, &
                       QL, QL_COEF,QL_LEN, &
                       SHP, INVK) 
					   
            
		        DO J = 1, LN
		          GSTACK_SHP(LSTART+J-1) = SHP(J)  !STORES THE INFLUENCE FUNCTION
		        END DO
                
                GINVK(:,:,I) =  INVK(:,:)
				
	  END DO !INTEGRATION POINT (NODE) LOOP
      !$OMP END DO
      !$OMP END PARALLEL	 
	  
	  
	  END IF
	  
	  
	  
      !$OMP PARALLEL DEFAULT(FIRSTPRIVATE) SHARED( GNUMP, GCOO, GCOO_CUURENT, GWIN, GSM_LEN, GSM_VOL, GSM_AREA, GN, GSTART, &
      !$OMP                                       DIM_NN_LIST, GSTACK, GSTACK_SHP, GSTACK_DSHP, GSTACK_DDSHP, GINVK, & 
      !$OMP                                       GCHAR_DIST,GMAX_WVEL, GMAXN, GGHOST, GEBC_NODES, GVOL, GNSNI_FAC, &
      !$OMP                                       GSTRESS, LOCAL_DX_STRESS, LOCAL_DY_STRESS, LOCAL_DZ_STRESS, &
      !$OMP                                       GSTRAIN, &
      !$OMP                                       GSTATE, GPROP, GDINC,GDINC_TOT, GMAT_TYPE, FINT, DLT_FINT, FINT_TEMP)       
      
      ID_RANK = OMP_get_thread_num()  !OMPJOE

      !$OMP DO   
	  DO I = 1, GNUMP
	    !
	    !
	    !GRAB NODE INFORMATION FROM LIST
		!
		LCOO(:) = GCOO(:,I)
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
		  LSTRESS(J) = GSTRESS(J,I)
		  LSTRAIN(J) = GSTRAIN(J,I)
		  L_H_STRESS(J) = G_H_STRESS(J,I)!GC
          L_S_STRESS(J) = G_S_STRESS(J,I)!GC
		END DO
		
        DO J = 1, 30
          LPROP(J) = GPROP(J,I)
        END DO
        
		DO J = 1, 20
		  LSTATE(J) = GSTATE(J,I)
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
		  JJ = LSTACK(J)
		  LDINC(1,J) = GDINC((JJ-1)*3+1)
		  LDINC(2,J) = GDINC((JJ-1)*3+2)
		  LDINC(3,J) = GDINC((JJ-1)*3+3)
		END DO
		!
		! GET THE GENERALIZED DISPLACEMENTS FOR NEIGHBORS
		!
		DO J = 1, LN
		  JJ = LSTACK(J)
		  LDINC_TOT(1,J) = GDINC_TOT((JJ-1)*3+1)
		  LDINC_TOT(2,J) = GDINC_TOT((JJ-1)*3+2)
		  LDINC_TOT(3,J) = GDINC_TOT((JJ-1)*3+3)
		END DO        
		!
		! GET THE CURRENT COORDINATES FOR NEIGHBORS
		!
		DO J = 1, LN
		  JJ = LSTACK(J)
		  LCOO_CUURENT(1,J) = GCOO_CUURENT(1,JJ)
		  LCOO_CUURENT(2,J) = GCOO_CUURENT(2,JJ)
		  LCOO_CUURENT(3,J) = GCOO_CUURENT(3,JJ)
		END DO        
        
        IF ((LLAGRANGIAN).AND.(.NOT.LINIT)) THEN
		  !IF IT IS LAGRANGIAN AND IT IS NOT THE FIRST STEP, RECALL SHAPE FUNCTIONS
    	  !IF (.FALSE.) THEN

            DO J = 1, LN
		      SHP(J) =   GSTACK_SHP(LSTART+J-1)
		    END DO
              INVK(:,:) =   GINVK(:,:,I)
			  
		   
            IF (LFINITE_STRAIN) THEN
            
		        FMAT = 0.0d0
                FMAT(1,1) = 1.d0
                FMAT(2,2) = 1.d0
                FMAT(3,3) = 1.d0
		        !
                X_0(1:3) = GCOO(1:3,I)
                X_t(1:3) = GCOO_CUURENT(1:3,I)
                CALL DFM_SHAPE_TENSOR(X_0,X_t, RK_DEGREE, RK_PSIZE, RK_CONT, GCOO, GVOL, GWIN, GNUMP, LSTACK, LN, GMAXN, LCOO_CUURENT, &
                   SHP, FMAT)

                   !F = S*invK
                   FMAT = MATMUL(FMAT, INVK)
                   
                CALL DETERMINANT(FMAT,DET)
                CALL INVERSE(FMAT, 3, IFMAT, INV_OK)
                IF (.NOT. INV_OK) CALL WARN('PROBLEM INVERTING MATRIX')
                
		    END IF
    		
    		
		ELSE
		    !
            ! TODO: CONDENSE ALL SHAPE FUNCTION CALCULATIONS
            !
            !
            ! DIRECT NODAL INTEGRATION
            !
            IF (LLAGRANGIAN) THEN

                    CALL UDFM_SHAPE_TENSOR(LCOO, RK_DEGREE, RK_PSIZE, RK_CONT, RK_IMPL,GCOO, GVOL, GWIN, GNUMP, LSTACK, LN, GMAXN, GEBC_NODES,SELF_EBC, &
                       QL, QL_COEF,QL_LEN, &
                       SHP, INVK) 
                
                       
                       
            ELSE !GCOO_CUURENT
			
			! not implemented
			
            END IF
			
			!#TODO: GET RID OF REDUNDANT CALCS

            
            IF ((LLAGRANGIAN).AND.(LINIT)) THEN
            
                
		        DO J = 1, LN
		          GSTACK_SHP(LSTART+J-1) = SHP(J)  !STORES THE INFLUENCE FUNCTION
		        END DO
                
                GINVK(:,:,I) =  INVK(:,:)
                IFMAT = 0.d0
                IFMAT(1,1) = 1.d0
                IFMAT(2,2) = 1.d0
                IFMAT(3,3) = 1.d0
                FMAT = IFMAT
                DET = 1.d0
            
            END IF
        
        END IF
        !
        ! COMPUTE STRAIN MEASURES
        !
		    !
		    ! COMPUTE THE INCREMENTAL DEFORMATION GRADIENT WITH RESPECT TO THE CURRENT TIME STEP
            !
		    LMAT = 0.0d0
		    !
            !
                X_0(1:3) = GCOO(1:3,I)
                X_t(1:3) = GDINC((I-1)*3+1 : (I-1)*3+3)
                
                CALL DFM_SHAPE_TENSOR(X_0,X_t, RK_DEGREE, RK_PSIZE, RK_CONT, GCOO, GVOL, GWIN, GNUMP, LSTACK, LN, GMAXN, LDINC, &
                   SHP, LMAT)
                   
               LMAT = MATMUL(LMAT,INVK)    
         
                !du/dx = du/dX*invF
                IF(.NOT.LINIT) THEN   
                LMAT = MATMUL(LMAT, IFMAT)
                ENDIF
		!
		! GET THE ROTATION MATRIX AND OBJECTIVE INCREMENTAL
		! STRAIN USING HUGHES-WINGET ALGORITHM
		!
		CALL HUGHES_WINGET(LMAT,ROT,STRAIN,D)
		!
		! ROTATE THE TOTAL STRESSES FROM PREVIOUS TIME STEP
		!
		CALL ROTATE_TENSOR(ROT,LSTRESS)
		!
		! ROTATE THE TOTAL STRAINS FROM PREVIOUS TIME STEP
		!
		CALL ROTATE_TENSOR(ROT,LSTRAIN)
        !
		LSTRAIN = LSTRAIN + STRAIN
		!
		! ELASTIC PREDICTOR
		!
		ELAS_MAT = FORM_CMAT(LPROP)
		LSTRESS_PREDICTOR = LSTRESS + MATMUL(ELAS_MAT,STRAIN)
		!

		CALL CONSTITUTION(LSTRESS_PREDICTOR,LMAT_TYPE, LSTRAIN, STRAIN, LPROP, DLT, FMAT, & !IN
		                  LSTATE, LSTRESS, L_S_STRESS, L_H_STRESS) !IN/OUT, OUT
		!
		! ********** SAVE STATE AND FEILD VARIABLES ********** 
		!
		!
        ! UPDATE STATE VARIABLE TO GSTATE
        !
        DO J = 1, 20
		  GSTATE(J,I) = LSTATE(J) 
		END DO    
        !
		DO J = 1, 6
		
		  !GET INCREMENTS FOR TIME STEP PREDICTION
		  STRESS_INC(J) = LSTRESS(J) - GSTRESS(J,I)
		  STRAIN_INC(J) = LSTRAIN(J) - GSTRESS(J,I)
		  
		  !SAVE THE STRESSES
		  GSTRESS(J,I) = LSTRESS(J)
		  GSTRAIN(J,I) = LSTRAIN(J)

		  G_H_STRESS(J,I) = L_H_STRESS(J)!GC
          G_S_STRESS(J,I) = L_S_STRESS(J)!GC
		  
		END DO
		
		
		
	           POISS = LPROP(1)
	           YOUNG = LPROP(2)
	                     
	           BULK = BULK_MOD(YOUNG,POISS)
	           
	           SHEAR = SHEAR_MOD(YOUNG,POISS)
	           
	           IF (LMAT_TYPE.GT.1) THEN
			     NSNI_FLAG=.FALSE.
	             CALL ESTIMATE_MODULI(STRESS_INC, STRAIN_INC, SHEAR_TRIAL, BULK_TRIAL, SHEAR, BULK , NSNI_FLAG)
	           END IF
	           
	           PMOD = BULK + 4.0d0*SHEAR/3.0d0
	           DENSITY =  LPROP(3)
	           MAXMOD=MAX(PMOD,SHEAR)
	           GMAX_WVEL(I) = DSQRT(MAXMOD/DENSITY)
			   
	           IF (LMAT_TYPE.EQ.3) THEN
			     !FUDGE THE DRUCKER-PRAGER TIME STEP SO THAT THE FACTOR CAN BE 1.0
	             GMAX_WVEL(I) = GMAX_WVEL(I)/0.15d0
			   END IF
			   
		
            
		!
		! ASSEMBLE THE INTERNAL FORCE
		!
        DO J = 1, LN
		
		  JJ = LSTACK(J)
		  
		   
               
			   IF (HODIV) THEN
			   
			   
           !
           !LSTRESS: CAUCHY STRESS --> FIRST PK STRESS
           !
           PKSTRESS = VTENSOR_2_TENSOR(LSTRESS)    
		   
           PKSTRESS = DET* MATMUL(PKSTRESS,TRANSPOSE(IFMAT))
           
		   !it should be the divergence of the nominal stress,
		   !which is the first PK transposed. But this doesnt
		   !work, so the theory and Joes implementation
		   !should be checked, Im not quite sure how
		   !Joe implemented the total Lagrangian formulation
		   !and its relation to the new higher order
		   !divergence operation (e.g., there is volume
		   !twice in the internal force)
		   ! -MH
		   !
			   !DIV(P^T)
			   !PKSTRESS=TRANSPOSE(PKSTRESS)
			   
           DX_t(:,1) = GCOO(:,I) - GCOO(:,JJ)            
           
           TEMP_STRESS = MATMUL(PKSTRESS,INVK) 
		   
           DX_t = MATMUL(TEMP_STRESS,DX_t)
           
           FINT3(1:3) = DX_t(1:3,1)  
           
		   
              INVK_J(:,:) =   GINVK(:,:,JJ)
			  
              TEMP_STRESS = MATMUL(PKSTRESS,INVK_J) 
		   
              DX_t(:,1) = GCOO(:,I) - GCOO(:,JJ)  
		   
              DX_t_J = MATMUL(TEMP_STRESS,DX_t)
           
              FINT3_J(1:3) = DX_t_J(1:3,1)  
			  
			  IF (1.eq.0) THEN
			  WRITE(*,*)
			  WRITE(*,*) 'INVK_J=', INVK_J
			  WRITE(*,*)
			  WRITE(*,*) 'DX_t_J=', DX_t_J
			  WRITE(*,*)
			  WRITE(*,*) 'FINT3_J=', FINT3_J
			  WRITE(*,*)
			  END IF
           
		       DO K = 1, 3
                 ID_RANK = OMP_get_thread_num()  !OMPJOE
             
                    !
                    !ASSEMBLE TO JJ
                    !
                    FINT_TEMP(ID_RANK+1,K,JJ) = FINT_TEMP(ID_RANK+1,K,JJ) - FINT3_J(K)*GVOL(I)*SHP(J)  * GVOL(JJ)
                    !
                    !ASSEMBLE TO I ITSELF
                    !       
                    FINT_TEMP(ID_RANK+1,K,I) = FINT_TEMP(ID_RANK+1,K,I) + FINT3(K)*GVOL(JJ)*SHP(J)    *GVOL(I)      
             
             
		       END DO     
			   
			   
			   
			   ELSE
			   
           !
           !LSTRESS: CAUCHY STRESS --> FIRST PK STRESS
           !
           PKSTRESS = VTENSOR_2_TENSOR(LSTRESS)    
		   
           PKSTRESS = DET* MATMUL(PKSTRESS,TRANSPOSE(IFMAT))
           
           DX_t(:,1) = GCOO(:,I) - GCOO(:,JJ)            
           
           TEMP_STRESS = MATMUL(PKSTRESS,INVK) 
		   
           DX_t = MATMUL(TEMP_STRESS,DX_t)
           
           FINT3(1:3) = DX_t(1:3,1)  
           
		   
		       DO K = 1, 3
                 ID_RANK = OMP_get_thread_num()  !OMPJOE
             
                    !
                    !ASSEMBLE TO JJ
                    !
                    FINT_TEMP(ID_RANK+1,K,JJ) = FINT_TEMP(ID_RANK+1,K,JJ) - FINT3(K)*GVOL(I)*SHP(J)  * GVOL(JJ)
                    !
                    !ASSEMBLE TO I ITSELF
                    !       
                    FINT_TEMP(ID_RANK+1,K,I) = FINT_TEMP(ID_RANK+1,K,I) + FINT3(K)*GVOL(JJ)*SHP(J)    *GVOL(I)      
             
             
		       END DO     
			   
			   END IF
             
		   !END IF
          
          
           
		
		END DO
        
        CONTINUE
	  
             
	  END DO !INTEGRATION POINT (NODE) LOOP
      !$OMP END DO
      !$OMP END PARALLEL	  
      !
      !FINT_TEMP TO HOLD THE VALUES FOR OPENMP REDUCE(ASSEMBLE,ADD) IN THE END
      !
      !$OMP PARALLEL PRIVATE(I,K,ID_RANK) SHARED(FINT_TEMP,NCORES_INPUT)
      !$OMP DO      
      DO I = 1, GNUMP
        DO K = 1, 3
		    FINT((I-1)*3+K) =  SUM(FINT_TEMP(1:NCORES_INPUT,K,I)) 
        END DO 
      
      END DO         
      !$OMP END DO
      !$OMP END PARALLEL 

	  
	  
	  

      IF (AUTO_TS) THEN
      !DO TIME STEP CALCS
	           
      DO I = 1, GNUMP
	  
		LN = GN(I)
		LSTART = GSTART(I)
	    !
		! GET THE NEIGHBOR LIST
		!
		DO J = 1, LN
		  LSTACK(J) = GSTACK(LSTART+J-1)
		END DO
		
	           !FIND THE CHARACTERISTIC DISTANCES
                           
	           XI=GCOO_CUURENT(1,I)
	           YI=GCOO_CUURENT(2,I)
	           ZI=GCOO_CUURENT(3,I)
	           
	           FIRST = .TRUE.
                     
		       DO J = 1, LN
		       
		         JJ = LSTACK(J)
		         
		         IF (JJ.NE.I) THEN
		         
		             XJ=GCOO_CUURENT(1,JJ)
		             YJ=GCOO_CUURENT(2,JJ)
		             ZJ=GCOO_CUURENT(3,JJ)
                          
		             DIST = DSQRT((XJ-XI)**2 + (YJ-YI)**2 + (ZJ-ZI)**2)
                                   
		             IF (FIRST) THEN
		               GCHAR_DIST(I) = DIST
		               FIRST = .FALSE.
		             ELSE
		               GCHAR_DIST(I) = MIN(GCHAR_DIST(I),DIST)
		             END IF
		         
		         END IF
         
		       END DO !J=1,GNUMP (NEIGHBOR NODES)
                     
		       DLT_TEMP = GCHAR_DIST(I) / GMAX_WVEL(I)
                 
                     
		       IF (I.EQ.1) THEN
		         DLT_FINT = DLT_TEMP
		       ELSE
		       
		         IF (DLT_TEMP.LT.DLT_FINT) THEN
		           DLT_FINT = DLT_TEMP
		         END IF
		         DLT_FINT = MIN(DLT_TEMP,DLT_FINT)
		       END IF
		       
        END DO
			   
        DLT_FINT = DLT_FINT*DLT_FAC
			   
			   
	  
      END IF !CALC TIME STEP
      
	  
	  DEALLOCATE(FINT_TEMP)
	  RETURN
	  END SUBROUTINE
	  
	  
	  
	  
					  