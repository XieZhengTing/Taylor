
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

      SUBROUTINE PRE_CONTROL
	  
	  USE CONTROL
	  USE MODEL
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! READ IN THE CONTROL VARIABLES IN THE *.MC FILE
	  !
	  IMPLICIT NONE
      !
      !
      ! GLOBAL IN/OUT
      !
      ! N/A
	  !
	  ! LOCAL I/O FLAGS ETC
	  !
	  CHARACTER(50):: CTEMP
	  INTEGER:: IERROR
	  !LOGICALS TO DETECT IF CARDS WERE READ
      LOGICAL:: READ_QL_PARAMETERS, READ_QL_LEN, READ_TIME_PARAMETERS, READ_GHOSTING_PARAMETERS, &
                READ_INTEGRATION_PARAMETERS, READ_RK_PARAMETERS, READ_HPC_PARAMETERS, &
                READ_OUTPUT_PARAMETERS,READ_FORM_PARAMETERS,READ_OUPUT_VARIABLES, &
				READ_NODE_SEARCH_PARAMETERS,READ_PERIODIC_SEARCH,COMPOSITE
      !TEMPORARY INTEGERS TO ASSIGN LOGICALS MEGA USES
      INTEGER:: I, J, IQL, IFINITE_STRAIN, ILAGRANGIAN,IFEM_OUTPUT
	  INTEGER:: IAUTO_TS, IPERIDYNAMICS, SHSUPT, ICOMPOSITE
      !TEMPORARY DOUBLES TO ASSIGN VALUSE MEGA USES
	  DOUBLE PRECISION:: DLT_VAL
      !...
      INTEGER:: IMWID, IMWID2
      INTEGER, ALLOCATABLE:: MIM(:)
      DOUBLE PRECISION, ALLOCATABLE:: XIM(:), YIM(:)
          
      CHARACTER(55):: QLCHAR
      DOUBLE PRECISION:: QL_LEN_REDUCE_AUTO
      
      ! GRAVITY
      LOGICAL:: READ_GRAVITY
      
	  LOGICAL::READ_FEM_OUTPUT,READ_STAB_CONTROL
      INTEGER:: IUSE_STAB_CONTROL
      
      INTEGER, ALLOCATABLE:: MAT_ID_COMPOSITE(:)
      DOUBLE PRECISION:: DIST, DIST_MIN
      INTEGER:: MAT_ID
      DOUBLE PRECISION:: CPROPS(2,30)
      DOUBLE PRECISION::PHI,C,PSI,PHI_RAD,PSI_RAD,SRT32,Q_PHI,K_PHI,Q_PSI
      INTEGER:: DP_TYPE
      DOUBLE PRECISION:: EPOISS, EYOUNG, LAMDA, MU, LAMDA_PLUS_2MU
      DOUBLE PRECISION:: CMAT_TYPE(2)
      
      
      QL_LEN_REDUCE_AUTO = 0.1d0
      
      !0329
      !DOUBLE PRECISION:: BODYFORCE(3)
	  
      !INTEGER:: IGRAVITY
	  
	  OPEN(10,FILE='Control.dat',STATUS='OLD',ACTION='READ',iostat=IERROR)
	  
	  IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR OPENING THE CONTROL FILE "Control.dat"',2)
	  
	  READ_RK_PARAMETERS = .FALSE.
	  READ_QL_PARAMETERS = .FALSE.
	  READ_HPC_PARAMETERS = .FALSE.
	  READ_TIME_PARAMETERS = .FALSE.
	  READ_NODE_SEARCH_PARAMETERS = .FALSE.
	  READ_GHOSTING_PARAMETERS = .FALSE.
      READ_INTEGRATION_PARAMETERS = .FALSE.
      READ_OUTPUT_PARAMETERS = .FALSE.
      READ_FORM_PARAMETERS = .FALSE.
      READ_OUPUT_VARIABLES = .FALSE.
	  READ_PERIODIC_SEARCH=.FALSE.
	  READ_FEM_OUTPUT = .FALSE.	  
      READ_STAB_CONTROL = .FALSE.
      PERIDYNAMICS =  .FALSE.
      READ_GRAVITY =  .FALSE.
      KCONTACT =  .FALSE.
      COMPOSITE=.FALSE.
      
      node_id_OUTPUT = .FALSE.
      displacement_OUTPUT = .FALSE.
      velocity_OUTPUT = .FALSE.
      acceleration_OUTPUT = .FALSE.
      fint_OUTPUT = .FALSE.
      fixity_OUTPUT = .FALSE.
      bodyid_OUTPUT = .FALSE.
      material_type_OUTPUT = .FALSE.
      initial_coordinate_OUTPUT = .FALSE.
      init_velocity_OUTPUT = .FALSE.
      normalized_window_OUTPUT = .FALSE.
      nodal_volume_OUTPUT = .FALSE.
      nodal_mass_OUTPUT = .FALSE.
      Poissons_ratio_OUTPUT = .FALSE.
      Youngs_modulus_OUTPUT = .FALSE.
      density_OUTPUT = .FALSE.
      physical_window_OUTPUT=.FALSE.
      BASIC_OUTPUT=.FALSE.
      ALL_OUTPUT=.FALSE.
      UNF_OUTPUT = .FALSE.
      characteristic_distance_OUTPUT=.FALSE.
      max_wave_vel_OUTPUT=.FALSE.
      pre_disp_force_OUTPUT=.FALSE.
      stress_OUTPUT=.FALSE.
      strain_OUTPUT=.FALSE.
	  eps_OUTPUT=.FALSE.   
      eps_shear_OUTPUT=.FALSE.       
      eps_vol_OUTPUT=.FALSE.   
      dam_OUTPUT=.FALSE.  
      IJKspace_OUTPUT=.FALSE.
      SHSUP=.FALSE.
      PDSEARCH=1
    
	  DO while(IERROR.eq.0) 
	  
	    READ(10,'(A50)',iostat=IERROR) CTEMP
	  
	    IF (IERROR.EQ.(-1)) EXIT
	    CTEMP=TRIM(CTEMP)
        
	    IF (CTEMP.eq.'*OUTPUT VARIABLES') THEN !*MODEL FILE CARD
	    
	    READ_OUPUT_VARIABLES = .TRUE.
      
	      DO WHILE(.TRUE.)
    	  
	        READ(10,'(A50)',iostat=IERROR) CTEMP
    	  
	        IF (IERROR.EQ.(-1)) GOTO 120
	        CTEMP=TRIM(CTEMP)
	        

            IF (CTEMP.eq.'node_id') node_id_OUTPUT = .TRUE.
            IF (CTEMP.eq.'displacement') displacement_OUTPUT = .TRUE.
            IF (CTEMP.eq.'velocity') velocity_OUTPUT = .TRUE.
            IF (CTEMP.eq.'acceleration') acceleration_OUTPUT = .TRUE.
            IF (CTEMP.eq.'fint') fint_OUTPUT = .TRUE.
            IF (CTEMP.eq.'fixity') fixity_OUTPUT = .TRUE.
            IF (CTEMP.eq.'bodyid') bodyid_OUTPUT = .TRUE.
            IF (CTEMP.eq.'material_type') material_type_OUTPUT = .TRUE.
            IF (CTEMP.eq.'initial_coordinate') initial_coordinate_OUTPUT = .TRUE.
            IF (CTEMP.eq.'init_velocity') init_velocity_OUTPUT = .TRUE.
            IF (CTEMP.eq.'normalized_window') normalized_window_OUTPUT = .TRUE.
            IF (CTEMP.eq.'nodal_volume') nodal_volume_OUTPUT = .TRUE.
            IF (CTEMP.eq.'nodal_mass') nodal_mass_OUTPUT = .TRUE.
            IF (CTEMP.eq.'Poissons_ratio') Poissons_ratio_OUTPUT = .TRUE.
            IF (CTEMP.eq.'Youngs_modulus') Youngs_modulus_OUTPUT = .TRUE.
            IF (CTEMP.eq.'density') density_OUTPUT = .TRUE.
            IF (CTEMP.eq.'physical_window') physical_window_OUTPUT = .TRUE.
            IF (CTEMP.eq.'basic') BASIC_OUTPUT = .TRUE.
            IF (CTEMP.eq.'all') ALL_OUTPUT = .TRUE.
            IF (CTEMP.eq.'unformatted') UNF_OUTPUT = .TRUE.
            IF (CTEMP.eq.'characteristic_distance') characteristic_distance_OUTPUT = .TRUE.
            IF (CTEMP.eq.'max_wave_vel') max_wave_vel_OUTPUT = .TRUE.
            IF (CTEMP.eq.'pre_disp_force') pre_disp_force_OUTPUT = .TRUE.
            IF (CTEMP.eq.'stress') stress_OUTPUT = .TRUE.            
            IF (CTEMP.eq.'strain') strain_OUTPUT = .TRUE.
            IF (CTEMP.eq.'eps') eps_OUTPUT = .TRUE.     
            IF (CTEMP.eq.'eps_shear') eps_shear_OUTPUT = .TRUE.     
            IF (CTEMP.eq.'eps_vol') eps_vol_OUTPUT = .TRUE.     
            IF (CTEMP.eq.'dam') dam_OUTPUT = .TRUE.  
            IF (CTEMP.eq.'IJKspace') IJKspace_OUTPUT = .TRUE.
            
    	  END DO
    	  
120       CONTINUE
    	  
    	  IF (ALL_OUTPUT) THEN
         
              node_id_OUTPUT = .TRUE.
              displacement_OUTPUT = .TRUE.
              velocity_OUTPUT = .TRUE.
              acceleration_OUTPUT = .TRUE.
              fint_OUTPUT = .TRUE.
              fixity_OUTPUT = .TRUE.
              bodyid_OUTPUT = .TRUE.
              material_type_OUTPUT = .TRUE.
              initial_coordinate_OUTPUT = .TRUE.
              init_velocity_OUTPUT = .TRUE.
              normalized_window_OUTPUT = .TRUE.
              nodal_volume_OUTPUT = .TRUE.
              nodal_mass_OUTPUT = .TRUE.
              Poissons_ratio_OUTPUT = .TRUE.
              Youngs_modulus_OUTPUT = .TRUE.
              density_OUTPUT = .TRUE.
              physical_window_OUTPUT= .TRUE.
              characteristic_distance_OUTPUT=.TRUE.
              max_wave_vel_OUTPUT=.TRUE.
              pre_disp_force_OUTPUT=.TRUE.
              stress_OUTPUT=.TRUE.
              strain_OUTPUT=.TRUE.
              eps_OUTPUT = .TRUE.     
              eps_shear_OUTPUT = .TRUE.     
              eps_vol_OUTPUT = .TRUE.     
              dam_OUTPUT = .TRUE.
              IJKspace_OUTPUT=.TRUE.
      
            END IF
            
            IF (BASIC_OUTPUT) THEN
              
			  !append basic variables to whatever else the user wants
              node_id_OUTPUT = .TRUE.
              displacement_OUTPUT = .TRUE.
              velocity_OUTPUT = .TRUE.
              bodyid_OUTPUT = .TRUE.
              eps_OUTPUT = .TRUE.     
              
    	  ENDIF
    	  
    	  ENDIF
    	  
      END DO
      
      CLOSE(10)
	  
	  
	  OPEN(10,FILE='Control.dat',STATUS='OLD',ACTION='READ')
	  
	  DO 
	  
	    READ(10,'(A50)',iostat=IERROR) CTEMP
	  
	    IF (IERROR.EQ.(-1)) EXIT
	    CTEMP=TRIM(CTEMP)
        
	    IF (CTEMP.eq.'*PERIDYNAMICS') THEN !*MODEL FILE CARD
          READ (10,*)  ! PERIDYNAMICS SWITCH
	      READ (10,*,iostat=IERROR) IPERIDYNAMICS
          IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *PERIDYNAMICS IN THE CONTROL FILE',2) !2 = USER ERROR
          IF (IPERIDYNAMICS.EQ.1) THEN
              PERIDYNAMICS =  .TRUE.
          ENDIF
        ENDIF
        
	    IF (CTEMP.eq.'*COMPOSITE') THEN !*MODEL FILE CARD
          READ (10,*)  ! COMPOSITE SWITCH
	      READ (10,*,iostat=IERROR) ICOMPOSITE
          IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *COMPOSITE IN THE CONTROL FILE',2) !2 = USER ERROR
          IF (ICOMPOSITE.EQ.1) THEN
              COMPOSITE =  .TRUE.
          ENDIF
          READ (10,*)  ! MATERIAL 1 TYPE
          !
          READ (10,*)  CMAT_TYPE(1)
          !
          READ (10,*)  ! MATERIAL 1 PROPS
          !
          READ(10,*,iostat=IERROR) CPROPS(1,1:5)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *COMPOSITE (MATERIAL 1), VALUE ROW 1',2)
          READ(10,*,iostat=IERROR) CPROPS(1,6:10)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *COMPOSITE (MATERIAL 1), VALUE ROW 2',2)
          READ(10,*,iostat=IERROR) CPROPS(1,11:15)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *COMPOSITE (MATERIAL 1), VALUE ROW 3',2)
          READ(10,*,iostat=IERROR) CPROPS(1,16:20)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *COMPOSITE (MATERIAL 1), VALUE ROW 4',2)
          READ(10,*,iostat=IERROR) CPROPS(1,21:25)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *COMPOSITE (MATERIAL 1), VALUE ROW 5',2)
          READ(10,*,iostat=IERROR) CPROPS(1,26:30)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *COMPOSITE (MATERIAL 1), VALUE ROW 6',2)
          
          READ (10,*)  ! MATERIAL 1 TYPE
          !
          READ (10,*)  CMAT_TYPE(2)
          !
          READ (10,*)  ! MATERIAL 1 PROPS
          !
          !
          READ(10,*,iostat=IERROR) CPROPS(2,1:5)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *COMPOSITE (MATERIAL 2), VALUE ROW 1',2)
          READ(10,*,iostat=IERROR) CPROPS(2,6:10)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *COMPOSITE (MATERIAL 2), VALUE ROW 2',2)
          READ(10,*,iostat=IERROR) CPROPS(2,11:15)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *COMPOSITE (MATERIAL 2), VALUE ROW 3',2)
          READ(10,*,iostat=IERROR) CPROPS(2,16:20)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *COMPOSITE (MATERIAL 2), VALUE ROW 4',2)
          READ(10,*,iostat=IERROR) CPROPS(2,21:25)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *COMPOSITE (MATERIAL 2), VALUE ROW 5',2)
          READ(10,*,iostat=IERROR) CPROPS(2,26:30)
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('ERROR READING CARD *COMPOSITE (MATERIAL 2), VALUE ROW 6',2)
          
          !PREP THE DATA
          
          DO I=1,2
          
          
            IF (CMAT_TYPE(I).EQ.3) THEN 
            ! DRUCKER PRAGER (THE ONLY ONE THAT NEEDS SPECIAL PREP!
                  DP_TYPE = CPROPS(I,4)
				  PHI = CPROPS(I,5)
				  C = CPROPS(I,6)
				  PSI = CPROPS(I,7)
				  !
				  PHI_RAD = PHI*3.141592653589793d0/180.0d0
				  PSI_RAD = PSI*3.141592653589793d0/180.0d0
				  SRT32 = SQRT(3.0d0/2.d0)
                  
                  IF (DP_TYPE.EQ.1) THEN
                  
				  Q_PHI = 2.0d0*SIN(PHI_RAD)/(SRT32*(3.0d0 - SIN(PHI_RAD))) *3.d0
				  K_PHI = 6.0d0*C*COS(PHI_RAD)/(SRT32*(3.0d0 - SIN(PHI_RAD)))
                  
                ELSEIF (DP_TYPE.EQ.2) THEN
				  Q_PHI = 2.0d0*SIN(PHI_RAD)/(SRT32*(3.0d0 + SIN(PHI_RAD)))  *3.d0
				  K_PHI = 6.0d0*C*COS(PHI_RAD)/(SRT32*(3.0d0 + SIN(PHI_RAD)))
                ELSEIF (DP_TYPE.EQ.3) THEN
				  Q_PHI = SIN(PHI_RAD)/SQRT(9.0d0 + 3.0d0*SIN(PHI_RAD)*SIN(PHI_RAD))  /SQRT(0.5d0)   *3.d0
				  K_PHI = 3.0d0*C*COS(PHI_RAD)/SQRT(9.0d0 + 3.0d0*SIN(PHI_RAD)*SIN(PHI_RAD)) /SQRT(0.5d0)
                 ELSEIF (DP_TYPE.EQ.4) THEN         
                  Q_PHI = CPROPS(I,5)   *3.d0
                  K_PHI = CPROPS(I,6)
                  Q_PSI = CPROPS(I,7)
                  
                  ENDIF
                
                 Q_PSI = 0.D0
                
				  CPROPS(I,25) = Q_PHI
				  CPROPS(I,26) = K_PHI
				  CPROPS(I,27) = Q_PSI
                  
			  END IF !DRUCKER PRAGER
              
              !PREP THE REST OF THE DATA
	              EPOISS = CPROPS(I,1)
	              EYOUNG = CPROPS(I,2)
	              
	              LAMDA = EYOUNG*EPOISS/((1.0d0+EPOISS)*(1.0d0-2.0d0*EPOISS))
	              
	              MU = EYOUNG/(2.0d0*(1.0d0+EPOISS))
	              
	              LAMDA_PLUS_2MU = LAMDA + 2.0d0*MU
	              
				  CPROPS(I,28) = LAMDA
				  CPROPS(I,29) = MU
				  CPROPS(I,30) = LAMDA_PLUS_2MU
				  
          END DO
          
          
          
          
        ENDIF
        
        
        
        IF (CTEMP.eq.'*SUPPORT') THEN !*MODEL FILE CARD
          READ (10,*)  ! SPHERICAL SUPPORT SWITCH
	      READ (10,*,iostat=IERROR) SHSUPT
          IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *SUPPORT IN THE CONTROL FILE',2) !2 = USER ERROR
          IF (SHSUPT.EQ.1) THEN
              SHSUP =  .TRUE.
          ENDIF
        ENDIF
        
        IF (CTEMP.eq.'*PERIODIC SEARCH') THEN !*MODEL FILE CARD
            READ (10,*)  ! SPHERICAL SUPPORT SWITCH
            READ (10,*,iostat=IERROR) PDSTIME
            IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *PERIODIC SEARCH IN THE CONTROL FILE',2) !2 = USER ERROR      
            READ_PERIODIC_SEARCH=.TRUE.
        ENDIF
        
	    IF (CTEMP.eq.'*RK PARAMETERS') THEN !*MODEL FILE CARD
	    
	      READ (10,*)  ! DEGREE     CONTINUETY    IMPLICIT      DILATION
	      READ (10,*,iostat=IERROR) RK_DEGREE, RK_CONT, RK_IMPL
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *RK PARAMETERS IN THE CONTROL FILE',2) !2 = USER ERROR
          READ_RK_PARAMETERS = .TRUE.
          
          !INPUT CHECK
          IF ((RK_DEGREE.LT.0)) THEN
	        CALL EXIT_PROGRAM('ERRROR, CHOOSE RK_DEGREE GREATER THAN 0 IN CARD *RK PARAMETERS IN THE CONTROL FILE',2)
          END IF
          
          IF (RK_DEGREE.GT.2) THEN
	        CALL EXIT_PROGRAM('ERRROR, CHOOSE RK_DEGREE = (0,1,2) IN CARD *RK PARAMETERS IN THE CONTROL FILE',-1) !-1 = NOT IMPLEMENTED
          END IF
          
          IF (RK_CONT.LT.(-1)) THEN
	        CALL EXIT_PROGRAM('ERRROR, CHOOSE RK_CONT > -1 IN CARD *RK PARAMETERS IN THE CONTROL FILE',2)
          END IF
          
          IF ((RK_IMPL.LT.0).OR.(RK_IMPL.GT.1)) THEN
	        CALL EXIT_PROGRAM('ERRROR, CHOOSE RK_IMPL = (0,1) IN CARD *RK PARAMETERS IN THE CONTROL FILE',2)
          END IF
          
	    END IF 
        
          
         IF (CTEMP.eq.'*STABILIZATION CONTROL') THEN
              READ_STAB_CONTROL = .TRUE.
              READ(10,*) !COMMENT
              READ(10,*) IUSE_STAB_CONTROL
              IF (IUSE_STAB_CONTROL.EQ.1) THEN
              USE_STAB_CONTROL=.TRUE.
              ELSE
              USE_STAB_CONTROL=.FALSE.
              END IF
              READ(10,*) !COMMENT
              READ(10,*) STABILIZATION_CONTROL_COEF
          END IF


	    IF (CTEMP.eq.'*FORMULATION PARAMETERS') THEN !*MODEL FILE CARD
	    
	      READ (10,*)  ! LFINITE_STRAIN  LLAGRANGIAN
	      READ (10,*,iostat=IERROR) IFINITE_STRAIN, ILAGRANGIAN
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *FORMULATION PARAMETERS IN THE CONTROL FILE',2)
          READ_FORM_PARAMETERS = .TRUE.
          
	      !TODO FOR ALL CONTROL CARDS: CHECK THE VALUES OF INPUT FOR ALL CONTROL CARDS #TODO
	      
          IF (IFINITE_STRAIN.EQ.1) THEN
             LFINITE_STRAIN = .TRUE.
          ELSEIF (IFINITE_STRAIN.EQ.0) THEN
             LFINITE_STRAIN = .FALSE.
          ELSE
	         CALL EXIT_PROGRAM('ERRROR, CHOOSE IFINITE_STRAIN=(0,1) IN CARD *FORMULATION PARAMETERS IN THE CONTROL FILE',2)
          END IF
          
          IF (ILAGRANGIAN.EQ.1) THEN
             LLAGRANGIAN = .TRUE.
          ELSEIF (ILAGRANGIAN.EQ.0) THEN
             LLAGRANGIAN = .FALSE.
          ELSE
	         CALL EXIT_PROGRAM('ERRROR, CHOOSE ILAGRANGIAN=(0,1) IN CARD *FORMULATION PARAMETERS IN THE CONTROL FILE',2)
          END IF
          
	    END IF 
	    
	    !GC
		IF (CTEMP.eq.'*FEM OUTPUT') THEN !*MODEL FILE CARD
			READ(10,*) !# LFEM_OUTPUT
			READ(10,*,iostat=IERROR) IFEM_OUTPUT
			IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *FEMMESH PARAMETER IN THE CONTROL FILE',2)
			READ_FEM_OUTPUT = .TRUE.
			
			IF (IFEM_OUTPUT .EQ. 1) THEN
			   LFEM_OUTPUT = .TRUE.
			ELSE IF (IFEM_OUTPUT .EQ. 0) THEN
			   LFEM_OUTPUT = .FALSE.
			ELSE
			   CALL EXIT_PROGRAM('ERRROR, CHOOSE IFEM_OUTPUT=(0,1) IN CARD *FEM OUTPUT IN THE CONTROL FILE',2)
			END IF
		END IF
        
	    !KC 0717
        IF (CTEMP.eq.'*KERNEL_CONTACT PARAMETERS') THEN !*MODEL FILE CARD
          READ (10,*)  ! KCONTACT (1=nodal orientation)
	      READ (10,*,iostat=IERROR) IKCONTACT
          IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *KERNEL_CONTACT PARAMETERS IN THE CONTROL FILE',2) 
          KCONTACT = .TRUE.
          	IF (IKCONTACT .EQ. 0) THEN
			   KCONTACT = .FALSE.
			ELSE IF (IKCONTACT .EQ. 1) THEN
			   KCONTACT = .TRUE.
            ELSE IF (IKCONTACT .EQ. 2) THEN
               KCONTACT = .TRUE.
			ELSE
			   CALL EXIT_PROGRAM('ERRROR, CHOOSE IIKCONTACT=(0,1) IN CARD *FEM OUTPUT IN THE CONTROL FILE',2)
			END IF
        ENDIF

        
	    IF (CTEMP.eq.'*OUTPUT PARAMETERS') THEN !*MODEL FILE CARD
	    
	      READ (10,*)  ! DEGREE     CONTINUETY    IMPLICIT      DILATION
	      READ (10,*,iostat=IERROR) TIME_OUTPUT
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *OUTPUT PARAMETERS IN THE CONTROL FILE',2)
          READ_OUTPUT_PARAMETERS = .TRUE.
	    END IF 
        
	    IF (CTEMP.eq.'*HPC PARAMETERS') THEN !*MODEL FILE CARD
	    
	      READ (10,*)  ! TYPE
	      READ (10,*,iostat=IERROR) HPC_SCHEME
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *HPC PARAMETERS IN THE CONTROL FILE',2)
          READ_HPC_PARAMETERS = .TRUE.
	    END IF 
	    
	    IF (CTEMP.eq.'*GHOSTING PARAMETERS') THEN !*MODEL FILE CARD
	    
	      READ (10,*)  ! GHOST_BUFFER
	      READ (10,*,iostat=IERROR) GHOST_BUFFER
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *HPC PARAMETERS IN THE CONTROL FILE',2)
          READ_GHOSTING_PARAMETERS = .TRUE.
	    END IF 
        
	    IF (CTEMP.eq.'*TIME PARAMETERS') THEN !*MODEL FILE CARD
	    
	      READ (10,*)  ! TIME_END
	      READ (10,*,iostat=IERROR) TIME_END
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING TIME_END IN THE CARD *TIME PARAMETERS IN THE CONTROL FILE',2)
	      
	      READ (10,*)  ! IAUTO_TS
	      READ (10,*,iostat=IERROR)  IAUTO_TS
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING IAUTO_TS IN THE CARD *TIME PARAMETERS IN THE CONTROL FILE',2)
	      
	      READ (10,*)  !TIME STEP OR FACTOR FOR AUTO
	      READ (10,*,iostat=IERROR)  DLT_VAL
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING DLT_VAL IN THE CARD *TIME PARAMETERS IN THE CONTROL FILE',2)
	      
	      
	      IF (IAUTO_TS.EQ.1) THEN
	        AUTO_TS = .TRUE.
	        DLT_FAC= DLT_VAL
	      ELSEIF (IAUTO_TS.EQ.0) THEN
	        AUTO_TS = .FALSE.
	        DLT = DLT_VAL
	        !DLT_FAC= 0.9d0
	        !DLT=DLT_FAC*DLT
	      ELSE
	        CALL EXIT_PROGRAM('ERRROR, CHOOSE IAUTO_TS=(0,1) IN CARD "*TIME PARAMETERS" IN THE CONTROL FILE',2)
	      END IF
	      
	      !TODO FOR ALL CONTROL CARDS: GIVE THE LINE IN THE ERROR FOR WHICH MEGA COULD NOT READ #TODO
	      
	      
	      
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *TIME PARAMETERS IN THE CONTROL FILE',2)
          READ_TIME_PARAMETERS = .TRUE.
	    END IF 
		
		
	    IF (CTEMP.eq.'*NODE SEARCH PARAMETERS') THEN !*MODEL FILE CARD
	    
	      READ (10,*)  ! TIME_SEARCH
	      READ (10,*,iostat=IERROR) TIME_SEARCH
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING TIME_SEARCH IN THE CARD *NODE SEARCH PARAMETERS IN THE CONTROL FILE',2)
          READ_NODE_SEARCH_PARAMETERS = .TRUE.
	      
	    END IF
	      
		  
	    IF (CTEMP.eq.'*INTEGRATION PARAMETERS') THEN !*MODEL FILE CARD
	    
	      READ (10,*)  ! TIME_END, DLT
	      READ (10,*,iostat=IERROR) ITYPE_INT
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *INTEGRATION PARAMETERS IN THE CONTROL FILE',2)
          READ_INTEGRATION_PARAMETERS = .TRUE.
	    END IF 
        
     	IF (CTEMP.eq.'*GRAVITY PARAMETERS') THEN !*MODEL FILE CARD
	    
	      READ (10,*)  ! VALUES (X, Y, Z)
	      READ (10,*,iostat=IERROR) IGRAVITY
	      IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *GRAVITY PARAMETERS IN THE CONTROL FILE',2)
          READ_GRAVITY = .TRUE.
	    END IF    
        
	    IF (CTEMP.eq.'*QL PARAMETERS') THEN !*MODEL FILE CARD
	      READ_QL_LEN = .TRUE.
	      READ (10,*)  ! USE, COEFFICIENT, LENGTH
	      READ(10,'(A50)',iostat=IERROR) CTEMP
	      READ (CTEMP,*,iostat=IERROR) IQL, QL_COEF,QL_LEN
	      IF (IERROR.NE.0) THEN
              ! TRY TO READ WITHOUT LENGTH
              READ (CTEMP,*,iostat=IERROR) IQL, QL_COEF
              READ_QL_LEN = .FALSE.
	          IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR READING THE CARD *QL PARAMETERS IN THE CONTROL FILE',2)
          END IF
          READ_QL_PARAMETERS = .TRUE.
          IF (IQL.EQ.1) THEN
          QL = .TRUE.
          END IF
          
	    END IF 
	    
	    
	  END DO !READ FILE WITH ARBITRARY DO
	  
	  CLOSE(10)
      
      
      IF (.NOT.READ_TIME_PARAMETERS) THEN
        CALL EXIT_PROGRAM('TIME PARAMETERS NOT READ IN INPUT FILE (CARD *TIME PARAMETERS), FATAL ERROR',2)
      END IF
      
	  
      IF (.NOT.READ_NODE_SEARCH_PARAMETERS) THEN
        CALL WARN('NODE SEARCH PARAMETERS NOT SPECIFIED, WILL UPDATE SEARCH BUCKETS EVERY TIME STEP')
		TIME_SEARCH = -1.0d0
      END IF
	  
      IF (.NOT.READ_OUTPUT_PARAMETERS) THEN
        CALL STRONG_WARN('OUTPUT PARAMETERS NOT SET, DEFAULTING ARBITRARILY TO 100 OUTPUTS')
	    TIME_OUTPUT = TIME_END*0.01d0
      END IF
      
      IF (.NOT.READ_GHOSTING_PARAMETERS) THEN
        CALL WARN('GHOSTING PARAMETERS WERE NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
        !0704
        !GHOST_BUFFER = 1000
        GHOST_BUFFER = 0
        
        !CALL LOG_APPEND('GHOST_BUFFER = 1000')
        CALL LOG_APPEND('GHOST_BUFFER = 0')
        WRITE(*,*) ' '
        WRITE(*,*) 'GHOST_BUFFER = 0'
        WRITE(*,*) ' '
        
      END IF   
      
      IF (.NOT.PERIDYNAMICS) THEN
        CALL WARN('PERIDYNAMICS PARAMETERS WERE NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
        !0704
        !GHOST_BUFFER = 1000
        PERIDYNAMICS = .FALSE.
        
        !CALL LOG_APPEND('GHOST_BUFFER = 1000')
        CALL LOG_APPEND('PERIDYNAMICS = .FALSE.')
        WRITE(*,*) ' '
        WRITE(*,*) 'PERIDYNAMICS = .FALSE.'
        WRITE(*,*) ' '
        
      END IF   
      
      
      IF (.NOT.READ_OUPUT_VARIABLES) THEN
      
        CALL WARN('THE CARD "*OUTPUT VARIABLES" WAS NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')

            node_id_OUTPUT = .TRUE.
            displacement_OUTPUT = .TRUE.
            bodyid_OUTPUT = .TRUE.
            material_type_OUTPUT = .TRUE.
            nodal_volume_OUTPUT = .TRUE.
            density_OUTPUT = .TRUE.

        CALL LOG_APPEND('node_id_OUTPUT = .TRUE.')
        CALL LOG_APPEND('displacement_OUTPUT = .TRUE.')
        CALL LOG_APPEND('bodyid_OUTPUT = .TRUE.')
        CALL LOG_APPEND('material_type_OUTPUT = .TRUE.')
        CALL LOG_APPEND('nodal_volume_OUTPUT = .TRUE.')
        CALL LOG_APPEND('density_OUTPUT = .TRUE.')
        
        WRITE(*,*) 'node_id_OUTPUT = .TRUE.'
        WRITE(*,*) 'displacement_OUTPUT = .TRUE.'
        WRITE(*,*) 'bodyid_OUTPUT = .TRUE.'
        WRITE(*,*) 'material_type_OUTPUT = .TRUE.'
        WRITE(*,*) 'nodal_volume_OUTPUT = .TRUE.'   
        WRITE(*,*) 'density_OUTPUT = .TRUE.'
        
      END IF
      
      
      IF (.NOT.READ_HPC_PARAMETERS) THEN
      
        CALL WARN('THE CARD "*HPC PARAMETERS" WERE NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
        HPC_SCHEME = 1
        !SERIAL
        CALL LOG_APPEND('HPC_SCHEME = 1 (serial)')
        
        WRITE(*,*) ' '
        WRITE(*,*) 'HPC_SCHEME = 1 (serial)'
        WRITE(*,*) ' '
      END IF
      
      !0704
      IF (.NOT.READ_INTEGRATION_PARAMETERS) THEN
      
        CALL WARN('THE CARD "*INTEGRATION PARAMETERS" WERE NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
        ITYPE_INT = 2
        !SERIAL
        CALL LOG_APPEND('ITYPE_INT = 2 (NSNI)')
        
        WRITE(*,*) ' '
        WRITE(*,*) 'ITYPE_INT = 2 (NSNI)'
        WRITE(*,*) ' '
      END IF
      
      IF (.NOT.READ_STAB_CONTROL) THEN
      
        CALL WARN('STABILIZATION CONTROL WAS NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
        CALL LOG_APPEND_SPACE('STABILZIATION CONTROL = ON, COEF = 1.0')
        WRITE(*,*) ' '
        WRITE(*,*) 'STABILZIATION CONTROL = OFF, COEF = N/A'
        WRITE(*,*) ' '
        USE_STAB_CONTROL=.FALSE.
        STABILIZATION_CONTROL_COEF = 1.0d0
              
      END IF
      
      
      
      IF (.NOT.READ_GRAVITY) THEN
      
        CALL WARN('THE CARD "*GRAVITY PARAMETERS" WERE NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
        IGRAVITY = 0.d0
        !SERIAL
        CALL LOG_APPEND('GRAVITY PARAMETERS = 0 0 0 (NO GRAVITY)')
        
        WRITE(*,*) ' '
        WRITE(*,*) 'GRAVITY PARAMETERS = 0 0 0 (NO GRAVITY)'
        WRITE(*,*) ' '
      END IF
      
	  IF (.NOT. READ_FEM_OUTPUT) THEN
	    CALL WARN('THE CARD "*FEM OUTPUT" WERE NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
		LFEM_OUTPUT = .FALSE.
		CALL LOG_APPEND('LFEM_OUTPUT = .FALSE. (NO FEM MESH OUTPUT)')
		
		WRITE(*,*) ' '
        WRITE(*,*) 'LFEM_OUTPUT = .FALSE. (NO FEM MESH OUTPUT)'
        WRITE(*,*) ' '
      END IF
      
      IF (.NOT. READ_PERIODIC_SEARCH) THEN
          CALL WARN('THE CARD "*PERIODIC SEARCH" WERE NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
          PDSEARCH=1
          PDSTIME=0.0D0
          CALL LOG_APPEND('PERIODIC SEARCH = .FALSE. (SEARCH EVERY STEP)')
          WRITE(*,*) ' '
          WRITE(*,*) 'PERIODIC SEARCH = .FALSE. (SEARCH EVERY STEP'
          WRITE(*,*) ' '
      END IF
      
      IF (.NOT.READ_FORM_PARAMETERS) THEN
      
        CALL WARN('THE CARD "*FORMULATION PARAMETERS" WAS NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
        !
        ! FOR NOW, DEFAULT TO LAGRANGIAN AND INFINTESIMAL STRAIN. IN THE FUTURE, DEFAULT TO THE OPOSITE #TODO
        !
       ! IF (.FALSE.) THEN
        
          LFINITE_STRAIN = .TRUE.
          LLAGRANGIAN = .FALSE.
             
          CALL LOG_APPEND('LFINITE_STRAIN = .TRUE. (FINITE STRAIN FORMULATION)')
          CALL LOG_APPEND('LLAGRANGIAN = .FALSE. (SEMI-LAGRANGIAN FORMULATION)')
        
          WRITE(*,*) 'LFINITE_STRAIN = .TRUE. (FINITE STRAIN FORMULATION)'
          WRITE(*,*) 'LLAGRANGIAN = .FALSE. (SEMI-LAGRANGIAN FORMULATION)'
          
       ! ELSE
        
         ! LFINITE_STRAIN = .FALSE.
          !LLAGRANGIAN = .TRUE.
             
        !  CALL LOG_APPEND('LFINITE_STRAIN = .FALSE. (FINITE STRAIN FORMULATION)')
         ! CALL LOG_APPEND('LLAGRANGIAN = .TRUE. (SEMI-LAGRANGIAN FORMULATION)')
        
         ! WRITE(*,*) 'LFINITE_STRAIN = .FALSE. (FINITE STRAIN FORMULATION)'
         ! WRITE(*,*) 'LLAGRANGIAN = .TRUE. (SEMI-LAGRANGIAN FORMULATION)'
          
       ! END IF
      END IF
      
      IF (.NOT. KCONTACT) THEN
        CALL WARN('KERNEL CONTACT PARAMETERS WERE NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
        KCONTACT = .FALSE.
        IKCONTACT = 0
        CALL LOG_APPEND('KCONTACT = .NATURAL CONTACT.')
        WRITE(*,*) ' '
        WRITE(*,*) 'KCONTACT = .NATURAL CONTACT.'
        WRITE(*,*) ' '
        
      END IF 
      
      IF (.NOT.READ_RK_PARAMETERS) THEN
      
        CALL WARN('RK PARAMETERS WERE NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
        
        RK_DEGREE = 1
        RK_CONT = 3
        RK_IMPL = 0
        
        CALL LOG_APPEND('RK_DEGREE = 1')
        CALL LOG_APPEND('RK_CONT = 3')
        CALL LOG_APPEND('RK_IMPL = 0')
        
        WRITE(*,*) 'RK_DEGREE = 1'
        WRITE(*,*) 'RK_CONT = 3'
        WRITE(*,*) 'RK_IMPL = 0'
        
      END IF
      
      
      IF (RK_DEGREE.EQ.0) RK_PSIZE = 1
      IF (RK_DEGREE.EQ.1) RK_PSIZE = 4
      IF (RK_DEGREE.EQ.2) RK_PSIZE = 10
      
      IF (RK_DEGREE.GT.2) CALL EXIT_PROGRAM('(RK_DEGREE.GT.2) NOT IMPLEMENTED YET',0)
      
      IF ((.NOT.READ_QL_LEN).AND.READ_RK_PARAMETERS) THEN
      
        IF (QL) THEN 
        
        !dont bother to spit out warnings if QL is not used
        CALL WARN('QL LENGTH WAS NOT FOUND IN THE CONTROL FILE, WILL BE AUTOMATICALLY GENERATED')
        
        QL_LEN = MODEL_WIN(1,1) 
        
        DO I=1, MODEL_NUMP
        IF (MODEL_WIN(1,I).LT.QL_LEN) QL_LEN = MODEL_WIN(1,I) 
        IF (MODEL_WIN(2,I).LT.QL_LEN) QL_LEN = MODEL_WIN(2,I) 
        IF (MODEL_WIN(3,I).LT.QL_LEN) QL_LEN = MODEL_WIN(3,I) 
        END DO
        
        QL_LEN = QL_LEN * QL_LEN_REDUCE_AUTO
        
        WRITE(QLCHAR,'(A9,E25.15,A20)') 'QL_LEN = ', QL_LEN, '(AUTOMATIC)'
        
        CALL LOG_APPEND(QLCHAR)
        
        WRITE(*,*) QLCHAR
        
        END IF
        
      END IF
      
      IF (.NOT.READ_QL_PARAMETERS) THEN
      
        CALL WARN('QL PARAMETERS WERE NOT FOUND IN THE CONTROL FILE, DEFAULT VALUES WILL BE USED')
        
        QL = .TRUE.
        QL_COEF = 0.001d0
        QL_LEN = MODEL_WIN(1,1)
        
        DO I=1, MODEL_NUMP
          IF (MODEL_WIN(1,I).LT.QL_LEN) QL_LEN = MODEL_WIN(1,I)
          IF (MODEL_WIN(2,I).LT.QL_LEN) QL_LEN = MODEL_WIN(2,I)
          IF (MODEL_WIN(3,I).LT.QL_LEN) QL_LEN = MODEL_WIN(3,I)
        END DO
        
        QL_LEN = QL_LEN * QL_LEN_REDUCE_AUTO
        
        WRITE(QLCHAR,'(A9,E25.15,A20)') 'QL_LEN = ', QL_LEN, '(AUTOMATIC)'
        
        
        CALL LOG_APPEND('QL = TRUE')
        CALL LOG_APPEND('QL_COEF = 0.001')
        CALL LOG_APPEND(QLCHAR)
        
        WRITE(*,*) 'QL = TRUE'
        WRITE(*,*) 'QL_COEF = 0.001'
        WRITE(*,*) QLCHAR
        
      END IF
      
	  !
      ! READ IN THE CONTROL PARAMETERS AND ASSIGN VALUES TO MODEL
      !
      
      
      
      CALL READ_CONTROL(MODEL_NUMP,MODEL_NODE_SET_LIST,MODEL_NODE_SET_LENGTH, &
                               MODEL_NODE_SET_ID,NUM_NODESET,MAX_NINODESET, FIXITY2_STEPS,TIME_END,&
                               MODEL_EBC,MODEL_NONZERO_EBC,FIXITY2_TIME,FIXITY2_NONZERO_EBC,MODEL_VINIT,MODEL_MAT_TYPE,MODEL_PROP,MODEL_BODY_ID, &
                               MODEL_NORM_WIN,MODEL_SET_NAMES,MODEL_BODYFORCE)
                               
      DO I=1, MODEL_NUMP
         MODEL_WIN(:,I) = MODEL_WIN(:,I) * MODEL_NORM_WIN(I)
      END DO
	  
      
      
      
      IF (COMPOSITE) THEN
      
	  OPEN(20,FILE='Composite.dat',STATUS='OLD',ACTION='READ',iostat=IERROR)
	  
	  IF (IERROR.NE.0) CALL EXIT_PROGRAM('THERE WAS AN ERROR OPENING THE FILE "Composite.dat"',2)
	  
      
          READ (20,*)
          READ (20,*)
          READ (20,*)
          READ (20,*) !HEADER
          READ (20,*)
          READ (20,*)
          READ (20,*)
          READ (20,*)
          READ (20,*) IMWID !IMAGE WIDTH
          READ (20,*)
          
          
          IMWID2=IMWID*IMWID
          
          ALLOCATE(XIM(IMWID2),YIM(IMWID2),MIM(IMWID2))
          
          DO I=1,IMWID2
          
          READ (20,*) XIM(I), YIM(I), MIM(I)
         
          END DO
          
          CLOSE(20)
          CONTINUE
          
          
          ALLOCATE(MAT_ID_COMPOSITE(MODEL_NUMP))
          
          !LOOP OVER ALL NODES
          DO I=1, MODEL_NUMP
          
          IF (I.EQ.100) THEN
          CONTINUE
          END IF
          
          DIST_MIN = DSQRT(   (XIM(1)-MODEL_COO(1,I))**2 +   (YIM(1)-MODEL_COO(2,I))**2 )
          MAT_ID_COMPOSITE(I) = MIM(1)
          
          !FIND THE MATERIAL ID BY LOOPING OVER THE TABLE...
          
              DO J=1,IMWID2
          
                DIST = DSQRT(   (XIM(J)-MODEL_COO(1,I))**2 +   (YIM(J)-MODEL_COO(2,I))**2 )
            
                IF (DIST.LT.DIST_MIN) THEN
                    DIST_MIN = DIST
                    MAT_ID_COMPOSITE(I) = MIM(J)
                END IF
              END DO
          
          IF (I.EQ.100) THEN
          CONTINUE
          END IF
          
          
          
          END DO
          
          
          !JUST NEED TO MODIFY: 
          
          !MODEL_MAT_TYPE = 3
          
          DO I=1,MODEL_NUMP
          
            IF (MAT_ID_COMPOSITE(I).EQ.1) THEN
          
            MODEL_MAT_TYPE(I) = CMAT_TYPE(1)
            
            MODEL_PROP(:,I) = CPROPS(1,:)
            
            ELSEIF (MAT_ID_COMPOSITE(I).EQ.2) THEN
          
            MODEL_PROP(:,I) = CPROPS(2,:)
            
            MODEL_MAT_TYPE(I) = CMAT_TYPE(2)
            
            END IF
            
            
          END DO
          
          
      
      END IF
      
      
      
      DO I=1, MODEL_NUMP
          
          MODEL_MASS((I-1)*3+1) = MODEL_VOL(I) * MODEL_PROP(3,I)
          MODEL_MASS((I-1)*3+2) = MODEL_VOL(I) * MODEL_PROP(3,I)
          MODEL_MASS((I-1)*3+3) = MODEL_VOL(I) * MODEL_PROP(3,I)
      
      END DO
      
      
      
	  CALL ECHO_CONTROL
      
      !DEALLOCATE(MODEL_ELCON)
      DEALLOCATE(MODEL_ELBID)
      DEALLOCATE(MODEL_NODE_SET_LIST)
      DEALLOCATE(MODEL_NODE_SET_ID)
      DEALLOCATE(MODEL_NODE_SET_LENGTH)
      

	  CONTINUE
	  
	  RETURN
	  
	  END SUBROUTINE
	  
	  
	  
