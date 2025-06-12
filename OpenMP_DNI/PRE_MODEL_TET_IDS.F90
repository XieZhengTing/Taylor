
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

    SUBROUTINE  GET_TET_IDS(T)
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! RETURN THE TET CONNECTIVITIES OF THE FACES
	  !
    
    IMPLICIT NONE
    
    INTEGER::T(12,4),J
    
    !
    ! TET DEFINATIONS
    !
    !(TET ID, NODE)
    J=1
    T(J,1) = 1
    T(J,2) = 2
    T(J,3) = 5
    
    J=2
    T(J,1) = 2
    T(J,2) = 5
    T(J,3) = 6
    
    J=3
    T(J,1) = 4
    T(J,2) = 8
    T(J,3) = 3
    
    J=4
    T(J,1) = 3
    T(J,2) = 7
    T(J,3) = 8
    
    J=5
    T(J,1) = 2
    T(J,2) = 3
    T(J,3) = 7
    
    J=6
    T(J,1) = 7
    T(J,2) = 2
    T(J,3) = 6
    
    J=7
    T(J,1) = 5
    T(J,2) = 8
    T(J,3) = 7
    
    J=8
    T(J,1) = 5
    T(J,2) = 7
    T(J,3) = 6
    
    J=9
    T(J,1) = 8
    T(J,2) = 4
    T(J,3) = 1
    
    J=10
    T(J,1) = 1
    T(J,2) = 8
    T(J,3) = 5 !3  !5?
    
    J=11
    T(J,1) = 1
    T(J,2) = 2
    T(J,3) = 3
    
    J=12
    T(J,1) = 1
    T(J,2) = 3
    T(J,3) = 4
    
    T(:,4) = 9
    
    RETURN
    
    END SUBROUTINE
       
       