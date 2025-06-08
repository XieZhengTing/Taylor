
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
    
    
       MODULE CONTROL
      !
    !$ACC DECLARE COPYIN(RK_DEGREE, RK_CONT, RK_IMPL, ITYPE_INT, RK_PSIZE, RK_DILA, TIME_OUTPUT, &
     !$ACC& IGRAVITY, QL_COEF, QL_LEN, QL, LFEM_OUTPUT, PERIDYNAMICS, KCONTACT, IKCONTACT, &
     !$ACC& HPC_SCHEME, NCORES, NCORES_INPUT, MAX_CORES, THIS_CORE_NO, LUSER_OMP_CORES, GHOST_BUFFER, &
     !$ACC& PDSEARCH, PDSTIME, AUTO_TS, DLT_FAC, TIME_SEARCH, LFINITE_STRAIN, LLAGRANGIAN, &
     !$ACC& STABILIZATION_CONTROL_COEF, USE_STAB_CONTROL, node_id_OUTPUT, displacement_OUTPUT, &
     !$ACC& velocity_OUTPUT, acceleration_OUTPUT, fint_OUTPUT, fixity_OUTPUT, bodyid_OUTPUT, &
     !$ACC& material_type_OUTPUT, initial_coordinate_OUTPUT, init_velocity_OUTPUT, &
     !$ACC& normalized_window_OUTPUT, nodal_volume_OUTPUT, nodal_mass_OUTPUT, Poissons_ratio_OUTPUT, &
     !$ACC& Youngs_modulus_OUTPUT, density_OUTPUT, physical_window_OUTPUT, BASIC_OUTPUT, ALL_OUTPUT, &
     !$ACC& UNF_OUTPUT, characteristic_distance_OUTPUT, max_wave_vel_OUTPUT, pre_disp_force_OUTPUT, &
     !$ACC& stress_OUTPUT, strain_OUTPUT, eps_OUTPUT, eps_shear_OUTPUT, eps_vol_OUTPUT, dam_OUTPUT, &
     !$ACC& IJKspace_OUTPUT, SHSUP)
      !
	  ! FUNCTION OF THIS LOGIC:
	  !
	  ! This module contains the "global" scalar-type parameters read into the control file. It is
      ! relatively harmless to include this module since all parameters are global and it should
      ! not interfere with the modular nature of the program
	  !
	  ! RK PARAMETERS
	  !
       INTEGER, SAVE::RK_DEGREE, RK_CONT, RK_IMPL, ITYPE_INT
       INTEGER, SAVE::RK_PSIZE
       DOUBLE PRECISION, SAVE:: RK_DILA, TIME_OUTPUT
       
	   !GRAVITY PARAMETERS
      ! LOGICAL::GRAVITY
       DOUBLE PRECISION, SAVE::IGRAVITY(3)
	   
       DOUBLE PRECISION, SAVE:: QL_COEF,QL_LEN
       LOGICAL, SAVE:: QL
	   
	   LOGICAL, SAVE::LFEM_OUTPUT ! MESH OUTPUT
	   
       !METHOD
       LOGICAL:: PERIDYNAMICS
	   LOGICAL:: KCONTACT
!       INTEGER, SAVE::IKCONTACT
	   !
	   !HPC VARIABLES
	   !
       INTEGER, SAVE::HPC_SCHEME
       INTEGER, SAVE::NCORES, NCORES_INPUT, MAX_CORES, THIS_CORE_NO
       LOGICAL:: LUSER_OMP_CORES
       INTEGER, SAVE:: GHOST_BUFFER
       INTEGER, SAVE:: PDSEARCH   ! IN SEMI LAGLANGIAN, WILL DO SOFT_SEARCH AFTER 'PDSEARCH' TIME STEPS, BY DEFAULT =1
!       DOUBLE PRECISION, SAVE::PDSTIME !PERIODIC SEARCH TIME
       !
	   ! SIMULATION PARAMETERS
	   !
       LOGICAL, SAVE:: AUTO_TS
       DOUBLE PRECISION, SAVE:: DLT_FAC,TIME_SEARCH
       LOGICAL, SAVE:: LFINITE_STRAIN,LLAGRANGIAN
       DOUBLE PRECISION, SAVE:: STABILIZATION_CONTROL_COEF
       LOGICAL, SAVE:: USE_STAB_CONTROL
	   !
	   ! OUTPUT PARAMETERS
	   !
       LOGICAL, SAVE:: node_id_OUTPUT
       LOGICAL, SAVE:: displacement_OUTPUT
       LOGICAL, SAVE:: velocity_OUTPUT
       LOGICAL, SAVE:: acceleration_OUTPUT
       LOGICAL, SAVE:: fint_OUTPUT
       LOGICAL, SAVE:: fixity_OUTPUT
       LOGICAL, SAVE:: bodyid_OUTPUT
       LOGICAL, SAVE:: material_type_OUTPUT
       LOGICAL, SAVE:: initial_coordinate_OUTPUT
       LOGICAL, SAVE:: init_velocity_OUTPUT
       LOGICAL, SAVE:: normalized_window_OUTPUT
       LOGICAL, SAVE:: nodal_volume_OUTPUT
       LOGICAL, SAVE:: nodal_mass_OUTPUT
       LOGICAL, SAVE:: Poissons_ratio_OUTPUT
       LOGICAL, SAVE:: Youngs_modulus_OUTPUT
       LOGICAL, SAVE:: density_OUTPUT
       LOGICAL, SAVE:: physical_window_OUTPUT
       LOGICAL, SAVE:: BASIC_OUTPUT
       LOGICAL, SAVE:: ALL_OUTPUT
       LOGICAL, SAVE:: UNF_OUTPUT
       LOGICAL, SAVE:: characteristic_distance_OUTPUT
       LOGICAL, SAVE:: max_wave_vel_OUTPUT
       LOGICAL, SAVE:: pre_disp_force_OUTPUT
       LOGICAL, SAVE:: stress_OUTPUT       
       LOGICAL, SAVE:: strain_OUTPUT 
	   
       LOGICAL, SAVE:: eps_OUTPUT   
       LOGICAL, SAVE:: eps_shear_OUTPUT       
       LOGICAL, SAVE:: eps_vol_OUTPUT   
       LOGICAL, SAVE:: dam_OUTPUT 

       LOGICAL:: SHSUP    !if ture, then using spherical support, =false by defult, then using box support
       INTEGER, SAVE::IKCONTACT
       DOUBLE PRECISION, SAVE::PDSTIME !PERIODIC SEARCH TIME
       LOGICAL, SAVE:: IJKspace_OUTPUT
       LOGICAL, SAVE, PARAMETER :: DETAILED_OUTPUT = .TRUE. ! For controlling detailed VTK output updates


       END MODULE
       