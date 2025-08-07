
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


      SUBROUTINE CLEAN_VTKS
	  !
	  ! FUNCTION OF THIS LOGIC:
	  !
	  ! REMOVE OLD VTKS
	  !
      
      USE MODEL
      !
	  IMPLICIT NONE
      ! 
      CHARACTER*16:: CFILE
      !
      !LOCAL VARIABLES
      INTEGER:: I, IFILE_STATUS
      !
      CALL WRITE_OUT('CLEANING VTK FILES...')
      ! 
      DO I=1, 9999
      
          if (I.lt.10) THEN
            WRITE(CFILE,'(A11,I1,A4)') 'MEGAOUT_000',I,'.vtk'
          elseif (I.lt.100) THEN
            WRITE(CFILE,'(A10,I2,A4)') 'MEGAOUT_00',I,'.vtk'
          elseif (I.lt.1000) THEN
            WRITE(CFILE,'(A9,I3,A4)') 'MEGAOUT_0',I,'.vtk'
          elseif (I.lt.10000) THEN
            WRITE(CFILE,'(A8,I4,A4)') 'MEGAOUT_',I,'.vtk'
          end if
          !
          open(unit=50, iostat=IFILE_STATUS, file=CFILE, status='old')
          if (IFILE_STATUS.EQ.0) close(50, status='delete')
      
      END DO
      
      
      RETURN
      END SUBROUTINE
      
      


      SUBROUTINE OUTPUT_STEP_VTK(exodusStep,MODEL_NUMP,MODEL_NUMEL, MODEL_ELCON,MODEL_COO_CURRENT,MODEL_NODE_IDS,MODEL_DSP_TOT,MODEL_VEL,MODEL_ACL,MODEL_FINT, &
      MODEL_EBC,  MODEL_BODY_ID,  MODEL_MAT_TYPE,  MODEL_COO,  MODEL_VINIT,  &
      MODEL_NORM_WIN,  MODEL_WIN,  MODEL_VOL,  MODEL_MASS,  MODEL_PROP,    &
      LOCAL_CHAR_DIST, LOCAL_WAVE_VEL, LOCAL_PRFORCE, LOCAL_STRESS, LOCAL_STRAIN,LOCAL_STATE, LOCAL_STRAIN_EQ, &
				LOCAL_IJKSPC)
	  !
	  ! FUNCTION OF THIS LOGIC:
	  !
	  ! WRITE TO VTK FILE
	  !
      
      USE CONTROL
      USE MODEL, ONLY: MODEL_ELCON_MOD => MODEL_ELCON
      !
	  IMPLICIT NONE
      ! 
      CHARACTER*16:: CFILE
      !
      ! INPUT VARIABLES
      !
      INTEGER:: exodusStep, MODEL_NUMP, MODEL_NUMEL, MODEL_ELCON(8,MODEL_NUMEL)
      !
      DOUBLE PRECISION::MODEL_COO_CURRENT(3,MODEL_NUMP), MODEL_COO(3,MODEL_NUMP), MODEL_VINIT (3,MODEL_NUMP),MODEL_NORM_WIN(MODEL_NUMP), &
                        MODEL_WIN(3,MODEL_NUMP)
      !
      DOUBLE PRECISION::MODEL_DSP_TOT(3*MODEL_NUMP),MODEL_VEL(3*MODEL_NUMP),MODEL_ACL(3*MODEL_NUMP),MODEL_FINT(3*MODEL_NUMP), &
                        MODEL_MASS(3*MODEL_NUMP)
      !
      DOUBLE PRECISION::MODEL_PROP(30,MODEL_NUMP)
      !
      DOUBLE PRECISION::MODEL_VOL(MODEL_NUMP)
      DOUBLE PRECISION:: LOCAL_PRFORCE(3,MODEL_NUMP),LOCAL_STRESS(6,MODEL_NUMP),LOCAL_STRAIN(6,MODEL_NUMP), LOCAL_STRAIN_EQ(MODEL_NUMP)
      !
      INTEGER::MODEL_NODE_IDS(MODEL_NUMP),MODEL_EBC(3,MODEL_NUMP),  MODEL_BODY_ID(MODEL_NUMP),  MODEL_MAT_TYPE(MODEL_NUMP)
      
      DOUBLE PRECISION::LOCAL_CHAR_DIST(MODEL_NUMP)
      DOUBLE PRECISION::LOCAL_WAVE_VEL(MODEL_NUMP)
      DOUBLE PRECISION::LOCAL_STATE(20,MODEL_NUMP)
      INTEGER::LOCAL_IJKSPC(3,MODEL_NUMP)
       
      !
      !LOCAL VARIABLES
      INTEGER:: I
      DOUBLE PRECISION:: LPROP(30)
      INTEGER:: IMAT_TYPE
      LOGICAL:: LOTSOFO
      ! 
 ! Debug output
 PRINT *, '=== VTK OUTPUT DEBUG ==='
 PRINT *, 'LFEM_OUTPUT =', LFEM_OUTPUT
 PRINT *, 'MODEL_NUMEL =', MODEL_NUMEL
 IF (ALLOCATED(MODEL_ELCON_MOD)) THEN
     PRINT *, 'MODEL_ELCON is allocated'
     ! Print sample connectivity data
     IF (MODEL_NUMEL .GT. 0) THEN
         PRINT *, 'First element connectivity:', MODEL_ELCON(:,1)
         PRINT *, 'Min node index:', MINVAL(MODEL_ELCON)
         PRINT *, 'Max node index:', MAXVAL(MODEL_ELCON)
         PRINT *, 'Expected range: 1 to', MODEL_NUMP
     END IF
 ELSE
     PRINT *, 'ERROR: MODEL_ELCON is NOT allocated!'
 END IF
      !exodusStep = CURRENT OUTPUT STEP
      
      
      if (exodusStep.lt.10) THEN
        WRITE(CFILE,'(A11,I1,A4)') 'MEGAOUT_000',exodusStep,'.vtk'
      elseif (exodusStep.lt.100) THEN
        WRITE(CFILE,'(A10,I2,A4)') 'MEGAOUT_00',exodusStep,'.vtk'
      elseif (exodusStep.lt.1000) THEN
        WRITE(CFILE,'(A9,I3,A4)') 'MEGAOUT_0',exodusStep,'.vtk'
      elseif (exodusStep.lt.10000) THEN
        WRITE(CFILE,'(A8,I4,A4)') 'MEGAOUT_',exodusStep,'.vtk'
      else
        CALL EXIT_PROGRAM('THERE WAS AN ERROR DURRING OUTPUT, OUT OF FILE INDICIES',-1)
      end if
      !
      OPEN(50,FILE = CFILE, STATUS='REPLACE')
      !
	  ! MAKE HEADER FOR VTK FILE 
	  !
      CALL VTK_OUTPUT_HEADER(50)
	  !
      !NOW OUTPUT THE POINTS
      !
 ! Debug: Check coordinates being written
 PRINT *, '=== WRITING COORDINATES TO VTK ==='
 PRINT *, 'First 3 nodes being written:'
 DO I = 1, MIN(3, MODEL_NUMP)
     PRINT '(A,I4,A,3E15.5)', 'Node', I, ':', MODEL_COO_CURRENT(:,I)
 END DO

      DO I=1,MODEL_NUMP
        WRITE(50,'(3E13.5)') MODEL_COO_CURRENT(1,I), MODEL_COO_CURRENT(2,I), MODEL_COO_CURRENT(3,I)
      END DO
      !
      WRITE(50,*) ' '
      IF (LFEM_OUTPUT) THEN
		WRITE(50,'(A6, I10, I20)') 'CELLS ', MODEL_NUMP + MODEL_NUMEL, MODEL_NUMP*2 + MODEL_NUMEL*9
	  ELSE
	    WRITE(50,'(A6, I10, I20)') 'CELLS ', MODEL_NUMP , MODEL_NUMP*2 
	  END IF
      
      DO I=1,MODEL_NUMP
        WRITE(50,'(2I8)') 1, I-1
      END DO
	  
	  IF (LFEM_OUTPUT) THEN
	    DO I=1,MODEL_NUMEL ! NUMBER OF ELEMENTS FOR THE BLOCK WITH THE ELE_INDEX
		  WRITE(50,'(9I8)') 8, MODEL_ELCON(:,I)-1
	    END DO
	  END IF
      
      WRITE(50,*) ' '
      IF (LFEM_OUTPUT) THEN
		WRITE(50,'(A10, I8)') 'CELL_TYPES ', MODEL_NUMP + MODEL_NUMEL
	  ELSE
	    WRITE(50,'(A10, I8)') 'CELL_TYPES ', MODEL_NUMP
	  END IF
      
      DO I=1,MODEL_NUMP
        WRITE(50,'(I8)') 1
      END DO
	  
	  IF (LFEM_OUTPUT) THEN
	    DO I=1,MODEL_NUMEL
	      WRITE(50,'(I8)') 12
	    END DO
	  END IF
      
      WRITE(50,*) ' '
      WRITE(50,*) ' '
      WRITE(50,'(A11, I10)') 'POINT_DATA ', MODEL_NUMP
      
            
      !-------------------------DATA OUTPUT--------------------------



      !INTEGERS SIZE (NUMNP)

      IF(material_type_OUTPUT) CALL VTK_WRITE_INT('material_type', MODEL_MAT_TYPE, MODEL_NUMP, 50)
            
      IF(node_id_OUTPUT) CALL VTK_WRITE_INT('node_id', MODEL_NODE_IDS, MODEL_NUMP, 50)
      
      IF(bodyid_OUTPUT) CALL VTK_WRITE_INT('bodyid', MODEL_BODY_ID, MODEL_NUMP, 50)
      
      IF(material_type_OUTPUT) CALL VTK_WRITE_INT('material_type', MODEL_MAT_TYPE, MODEL_NUMP, 50)
      
	  
      
      
      !INTEGERS SIZE (3*NUMNP)

      IF(fint_OUTPUT) CALL VTK_WRITE_INT_THREE('fixity_x', 'fixity_y', 'fixity_z', MODEL_EBC, MODEL_NUMP, 50) 
     
      IF(IJKspace_OUTPUT) CALL VTK_WRITE_INT_THREE('Ispace', 'Jspace','Kspace',LOCAL_IJKSPC, MODEL_NUMP, 50)
      
      !DOUBLES SIZE (3*NUMNP)
      
      IF(displacement_OUTPUT) CALL VTK_WRITE_DBL_THREE('displacement_x', 'displacement_y', 'displacement_z', MODEL_DSP_TOT, MODEL_NUMP, 50)  
      
      IF(velocity_OUTPUT) CALL VTK_WRITE_DBL_THREE('velocity_x', 'velocity_y', 'velocity_z', MODEL_VEL, MODEL_NUMP, 50)  
      
      IF(acceleration_OUTPUT) CALL VTK_WRITE_DBL_THREE('acceleration_x', 'acceleration_y', 'acceleration_z', MODEL_ACL, MODEL_NUMP, 50)  
      
      IF(fint_OUTPUT) CALL VTK_WRITE_DBL_THREE('fint_x', 'fint_y', 'fint_z', MODEL_FINT, MODEL_NUMP, 50)  

      
      !DOUBLES SIZE (NUMNP)

      IF(characteristic_distance_OUTPUT) CALL VTK_WRITE_DBL('characteristic_distance', LOCAL_CHAR_DIST, MODEL_NUMP, 50)
      IF(max_wave_vel_OUTPUT) CALL VTK_WRITE_DBL('max_wave_vel', LOCAL_WAVE_VEL, MODEL_NUMP, 50)

	   !
	   !OUTPUT STATE VARIABLES:
	   !
      IF(eps_shear_OUTPUT) CALL VTK_WRITE_DBL('eps_shear', LOCAL_STATE(1,:), MODEL_NUMP, 50)
      IF(eps_vol_OUTPUT) CALL VTK_WRITE_DBL('eps_vol',     LOCAL_STATE(2,:), MODEL_NUMP, 50)
      IF(eps_OUTPUT) CALL VTK_WRITE_DBL('eps',             LOCAL_STATE(3,:), MODEL_NUMP, 50)
      IF(dam_OUTPUT) CALL VTK_WRITE_DBL('dam',             LOCAL_STATE(4,:), MODEL_NUMP, 50)


      IF(initial_coordinate_OUTPUT) THEN
         CALL VTK_WRITE_DBL('initial_coordinate_x', MODEL_COO(1,:), MODEL_NUMP, 50)
         CALL VTK_WRITE_DBL('initial_coordinate_y', MODEL_COO(2,:), MODEL_NUMP, 50)
         CALL VTK_WRITE_DBL('initial_coordinate_z', MODEL_COO(3,:), MODEL_NUMP, 50)
      END IF
      
      IF(init_velocity_OUTPUT) THEN
         CALL VTK_WRITE_DBL('init_velocity_x', MODEL_VINIT(1,:), MODEL_NUMP, 50)
         CALL VTK_WRITE_DBL('init_velocity_y', MODEL_VINIT(2,:), MODEL_NUMP, 50)
         CALL VTK_WRITE_DBL('init_velocity_z', MODEL_VINIT(3,:), MODEL_NUMP, 50)
      END IF
      
      IF(normalized_window_OUTPUT) CALL VTK_WRITE_DBL('normalized_window', MODEL_NORM_WIN, MODEL_NUMP, 50)
      
      IF (physical_window_OUTPUT) THEN
         CALL VTK_WRITE_DBL('physical_window_x', MODEL_WIN(1,:), MODEL_NUMP, 50)
         CALL VTK_WRITE_DBL('physical_window_y', MODEL_WIN(2,:), MODEL_NUMP, 50)
         CALL VTK_WRITE_DBL('physical_window_z', MODEL_WIN(3,:), MODEL_NUMP, 50)
      END IF
      
      IF(nodal_volume_OUTPUT) CALL VTK_WRITE_DBL('nodal_volume', MODEL_VOL, MODEL_NUMP, 50)
            
            
      IF(pre_disp_force_OUTPUT) THEN
         CALL VTK_WRITE_DBL('pre_disp_force_x', LOCAL_PRFORCE(1,:), MODEL_NUMP, 50)
         CALL VTK_WRITE_DBL('pre_disp_force_y', LOCAL_PRFORCE(2,:), MODEL_NUMP, 50)
         CALL VTK_WRITE_DBL('pre_disp_force_z', LOCAL_PRFORCE(3,:), MODEL_NUMP, 50)
      END IF 
      
      !STRESS OUTPUT
      IF(stress_OUTPUT) THEN
         CALL VTK_WRITE_DBL('stress_xx', LOCAL_STRESS(1,:), MODEL_NUMP, 50)
         CALL VTK_WRITE_DBL('stress_yy', LOCAL_STRESS(2,:), MODEL_NUMP, 50)
         CALL VTK_WRITE_DBL('stress_zz', LOCAL_STRESS(3,:), MODEL_NUMP, 50)
         CALL VTK_WRITE_DBL('stress_yz', LOCAL_STRESS(4,:), MODEL_NUMP, 50)         
         CALL VTK_WRITE_DBL('stress_xz', LOCAL_STRESS(5,:), MODEL_NUMP, 50)         
         CALL VTK_WRITE_DBL('stress_xy', LOCAL_STRESS(6,:), MODEL_NUMP, 50)         
      ENDIF
      
      !STRAIN OUTPUT
      IF(strain_OUTPUT) THEN
         CALL VTK_WRITE_DBL('strain_xx', LOCAL_STRAIN(1,:), MODEL_NUMP, 50)
         CALL VTK_WRITE_DBL('strain_yy', LOCAL_STRAIN(2,:), MODEL_NUMP, 50)
         CALL VTK_WRITE_DBL('strain_zz', LOCAL_STRAIN(3,:), MODEL_NUMP, 50)
         CALL VTK_WRITE_DBL('strain_yz', LOCAL_STRAIN(4,:), MODEL_NUMP, 50)         
         CALL VTK_WRITE_DBL('strain_xz', LOCAL_STRAIN(5,:), MODEL_NUMP, 50)         
         CALL VTK_WRITE_DBL('strain_xy', LOCAL_STRAIN(6,:), MODEL_NUMP, 50)   
         CALL VTK_WRITE_DBL('strain_eq', LOCAL_STRAIN_EQ(:), MODEL_NUMP, 50)   
      ENDIF
      
      
      !DO NODAL MASSES
      
      IF(nodal_mass_OUTPUT) THEN
      
          CALL START_DATA(50)
          WRITE(50,'(A8, A30, A10)')    'SCALARS ',   'nodal_mass',     'double 1'
          CALL MID_DATA(50)
          !
          DO I=1,MODEL_NUMP
            WRITE(50,'(E20.5)') MODEL_MASS((I-1)*3+1)
          END DO
          
      END IF
      
      
      
      !MATERIAL CONSTANTS
      
      IF(Poissons_ratio_OUTPUT) THEN
          CALL START_DATA(50)
          WRITE(50,'(A8, A30, A10)')    'SCALARS ',   'Poissons_ratio',     'double 1'
          CALL MID_DATA(50)
          
          DO I=1,MODEL_NUMP
            WRITE(50,'(E20.5)') MODEL_PROP(1,I)
          END DO
      END IF
      
      
      IF(Youngs_modulus_OUTPUT) THEN
            
          CALL START_DATA(50)
          WRITE(50,'(A8, A30, A10)')    'SCALARS ',   'Youngs_modulus',     'double 1'
          CALL MID_DATA(50)
          
          DO I=1,MODEL_NUMP
            WRITE(50,'(E20.5)') MODEL_PROP(2,I)
          END DO
      END IF
      
      IF(density_OUTPUT) THEN
          CALL START_DATA(50)
          WRITE(50,'(A8, A30, A10)')    'SCALARS ',   'density',     'double 1'
          CALL MID_DATA(50)
          
          DO I=1,MODEL_NUMP
            WRITE(50,'(E20.5)') MODEL_PROP(3,I)
          END DO
      END IF
      
      
      CLOSE(50)
      
      RETURN
      
      
      END SUBROUTINE
      
      
      
	  
      SUBROUTINE VTK_OUTPUT_HEADER(VTKIOUNIT)
	  !
	  ! FUNCTION OF THIS LOGIC:
	  !
	  ! OPEN AND MAKE A HEADER FOR THE VTK FILE
	  !
      USE MODEL
      !
	  IMPLICIT NONE
      !
      !MODULAR INPUT:
      !
      INTEGER, INTENT(IN):: VTKIOUNIT
      !
      WRITE(VTKIOUNIT,'(A26)') '# vtk DataFile Version 2.0 '
      WRITE(VTKIOUNIT,'(A30)') '(MEGA VTK FILE)'
      WRITE(VTKIOUNIT,'(A5)') 'ASCII'
      WRITE(VTKIOUNIT,'(A25)') 'DATASET UNSTRUCTURED_GRID'
      WRITE(VTKIOUNIT,'(A7,I10,A10)') 'POINTS ', MODEL_NUMP, 'double'
      
      RETURN
      
      END SUBROUTINE
      
      
      
      SUBROUTINE VTK_WRITE_DBL(VARIABLE_NAME, VARIABLE, MODEL_NUMP, VTKIOUNIT)
	  !
	  ! FUNCTION OF THIS LOGIC:
	  !
	  ! WRITE A DOUBLE TO THE VTK FILE
	  !
	  IMPLICIT NONE
      !
      !MODULAR INPUT:
      !
	  CHARACTER(*):: VARIABLE_NAME
      INTEGER, INTENT(IN):: MODEL_NUMP, VTKIOUNIT
      DOUBLE PRECISION, INTENT(IN):: VARIABLE(MODEL_NUMP)
      !
      !LOCAL VARIABLES
      !
      INTEGER:: I
      
      CALL START_DATA(VTKIOUNIT)
      
      WRITE(VTKIOUNIT,'(A8, A30, A10)')    'SCALARS ',   VARIABLE_NAME,     'double 1'

      CALL MID_DATA(VTKIOUNIT)
      !
      DO I=1,MODEL_NUMP
        WRITE(50,'(E20.5)') VARIABLE(I)
      END DO
      
      
      RETURN
      
      END SUBROUTINE
      
      
      
      SUBROUTINE VTK_WRITE_DBL_THREE(VARIABLE_NAME_X,VARIABLE_NAME_Y,VARIABLE_NAME_Z, VARIABLE, MODEL_NUMP, VTKIOUNIT)
	  !
	  ! FUNCTION OF THIS LOGIC:
	  !
	  ! WRITE A DOUBLE TO THE VTK FILE WITH A STACKED FORMAT FOR DOFS
	  !
	  IMPLICIT NONE
      !
      !MODULAR INPUT:
      !
	  CHARACTER(*):: VARIABLE_NAME_X,VARIABLE_NAME_Y,VARIABLE_NAME_Z
      INTEGER, INTENT(IN):: MODEL_NUMP, VTKIOUNIT
      DOUBLE PRECISION, INTENT(IN):: VARIABLE(3*MODEL_NUMP)
      !
      !LOCAL
      !
      INTEGER:: I
      
      
          CALL START_DATA(VTKIOUNIT)
          WRITE(VTKIOUNIT,'(A8, A30, A10)')    'SCALARS ',   VARIABLE_NAME_X,     'double 1'
          CALL MID_DATA(VTKIOUNIT)
          DO I=1,MODEL_NUMP
            WRITE(VTKIOUNIT,'(E20.5)') VARIABLE((I-1)*3+1)
          END DO
          
          CALL START_DATA(VTKIOUNIT)
          WRITE(VTKIOUNIT,'(A8, A30, A10)')    'SCALARS ',   VARIABLE_NAME_Y,     'double 1'
          CALL MID_DATA(VTKIOUNIT)
          DO I=1,MODEL_NUMP
            WRITE(VTKIOUNIT,'(E20.5)') VARIABLE((I-1)*3+2)
          END DO
          
          CALL START_DATA(VTKIOUNIT)
          WRITE(VTKIOUNIT,'(A8, A30, A10)')    'SCALARS ',  VARIABLE_NAME_Z,     'double 1'
          CALL MID_DATA(50)
          DO I=1,MODEL_NUMP
            WRITE(VTKIOUNIT,'(E20.5)') VARIABLE((I-1)*3+3)
          END DO
      
      
      RETURN
      
      END SUBROUTINE
      

      
      SUBROUTINE VTK_WRITE_INT_THREE(VARIABLE_NAME_X,VARIABLE_NAME_Y,VARIABLE_NAME_Z, VARIABLE, MODEL_NUMP, VTKIOUNIT)
	  !
	  ! FUNCTION OF THIS LOGIC:
	  !
	  ! WRITE A DOUBLE TO THE VTK FILE WITH A STACKED FORMAT FOR DOFS
	  !
	  IMPLICIT NONE
      !
      !MODULAR INPUT:
      !
	  CHARACTER(*):: VARIABLE_NAME_X,VARIABLE_NAME_Y,VARIABLE_NAME_Z
      INTEGER, INTENT(IN):: MODEL_NUMP, VTKIOUNIT
      INTEGER, INTENT(IN):: VARIABLE(3*MODEL_NUMP)
      !
      !LOCAL
      !
      INTEGER:: I
      
          CALL START_DATA(VTKIOUNIT)
          WRITE(VTKIOUNIT,'(A8, A30, A10)')    'SCALARS ',   VARIABLE_NAME_X,     'double 1'
          CALL MID_DATA(VTKIOUNIT)
          DO I=1,MODEL_NUMP
            WRITE(VTKIOUNIT,'(I20)') VARIABLE((I-1)*3+1)
          END DO
          
          CALL START_DATA(VTKIOUNIT)
          WRITE(VTKIOUNIT,'(A8, A30, A10)')    'SCALARS ',   VARIABLE_NAME_Y,     'double 1'
          CALL MID_DATA(VTKIOUNIT)
          DO I=1,MODEL_NUMP
            WRITE(VTKIOUNIT,'(I20)') VARIABLE((I-1)*3+2)
          END DO
          
          CALL START_DATA(VTKIOUNIT)
          WRITE(VTKIOUNIT,'(A8, A30, A10)')    'SCALARS ',  VARIABLE_NAME_Z,     'double 1'
          CALL MID_DATA(50)
          DO I=1,MODEL_NUMP
            WRITE(VTKIOUNIT,'(I20)') VARIABLE((I-1)*3+3)
          END DO
      
      
      RETURN
      
      END SUBROUTINE
      
      
      
      SUBROUTINE VTK_WRITE_INT(VARIABLE_NAME, VARIABLE, MODEL_NUMP, VTKIOUNIT)
	  !
	  ! FUNCTION OF THIS LOGIC:
	  !
	  ! WRITE A DOUBLE TO THE VTK FILE
	  !
	  IMPLICIT NONE
      !
      !MODULAR INPUT:
      !
	  CHARACTER(*):: VARIABLE_NAME
      INTEGER, INTENT(IN):: MODEL_NUMP, VTKIOUNIT
      INTEGER, INTENT(IN):: VARIABLE(MODEL_NUMP)
      !
      !LOCAL VARIABLES
      !
      INTEGER:: I
      
      CALL START_DATA(VTKIOUNIT)
      
      WRITE(VTKIOUNIT,'(A8, A30, A10)')    'SCALARS ',   VARIABLE_NAME,     'double 1'

      CALL MID_DATA(VTKIOUNIT)
      !
      DO I=1,MODEL_NUMP
        WRITE(50,'(I20)') VARIABLE(I)
      END DO
      
      
      RETURN
      
      END SUBROUTINE
      
      
      
      SUBROUTINE START_DATA(VTKIOUNIT)
	  IMPLICIT NONE
      INTEGER, INTENT(IN):: VTKIOUNIT
      !
      WRITE(VTKIOUNIT,*) ' '
      !
      RETURN
      END SUBROUTINE
      
      
      
      SUBROUTINE MID_DATA(VTKIOUNIT)
	  IMPLICIT NONE
      INTEGER, INTENT(IN):: VTKIOUNIT
      !
      WRITE(VTKIOUNIT,'(A20)') 'LOOKUP_TABLE default'
      !
      RETURN
      END SUBROUTINE
      