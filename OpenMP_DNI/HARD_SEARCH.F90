
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
    
    SUBROUTINE HARD_SEARCH(GNUMP,GCOO,GN,GSTART,DIM_NN_LIST,GSTACK,GMAXN, GWIN, GXDIST_MAX, GYDIST_MAX, GZDIST_MAX)
    
    IMPLICIT NONE
      !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! PERFORM BASIC O(N^2) NODE SEARCH
	  !
      INTEGER, INTENT(INOUT)::GNUMP
      INTEGER, INTENT(INOUT)::GN(GNUMP)
      DOUBLE PRECISION, INTENT(IN)::GCOO(3,GNUMP)
      DOUBLE PRECISION, INTENT(IN)::GWIN(3,GNUMP)
	  INTEGER, INTENT(INOUT)::GSTART(GNUMP)              
	  INTEGER, INTENT(INOUT)::DIM_NN_LIST
	  INTEGER, INTENT(INOUT)::GSTACK(DIM_NN_LIST)
	  INTEGER, INTENT(INOUT):: GMAXN
      DOUBLE PRECISION, INTENT(IN)::GXDIST_MAX, GYDIST_MAX, GZDIST_MAX

    
    
    
    !LOCAL
    INTEGER:: I,J,K
    INTEGER:: DIM_NN_LIST_ACTUAL
    LOGICAL:: FIRST
    DOUBLE PRECISION:: XI, YI, ZI, XJ, YJ, ZJ, DIST_NORM
    
    DIM_NN_LIST_ACTUAL = 0
    
    GN = 0
    
    K = 0
    
    GMAXN = 0
    
    DO I=1,GNUMP
    !INTEGERATION POINT I
    
      XI=GCOO(1,I)
      YI=GCOO(2,I)
      ZI=GCOO(3,I)
    
      FIRST = .TRUE.
                  
      !NODES J
      DO J=1,GNUMP
        
        XJ=GCOO(1,J)
        YJ=GCOO(2,J)
        ZJ=GCOO(3,J)
        
        DIST_NORM = (ABS(XJ-XI)-GXDIST_MAX)/GWIN(1,J)
        !DIST_NORM = (ABS(XJ-XI)-GXDIST_MAX)/GWIN(1,I)        
        !DIST_NORM = ABS(XJ-XI)/GWIN(1,I)        
        
        IF (DIST_NORM.LE.(1.0d0)) THEN
        
          DIST_NORM = (ABS(YJ-YI)-GYDIST_MAX)/GWIN(2,J)
        
            IF (DIST_NORM.LE.(1.0d0)) THEN
        
              DIST_NORM = (ABS(ZJ-ZI)-GZDIST_MAX)/GWIN(3,J)
        
              IF (DIST_NORM.LE.(1.0d0)) THEN
        
                GN(I) = GN(I) + 1
                
                K = K + 1
                
                GSTACK(K) = J
                
                IF (FIRST) THEN
                  !
                  !SAVELOCATION IN STACK
                  !
                  GSTART(I) = K 
                  FIRST = .FALSE.
                END IF
              END IF
           END IF
         END IF
         
      END DO !J=1,GNUMP (NEIGHBOR NODES)
      
      IF (GN(I).GT.GMAXN) GMAXN = GN(I)
      
     END DO !I=1,GNUMP (INTEGRATION POINTS)
     
     DIM_NN_LIST_ACTUAL = K
                
     RETURN
     
     END SUBROUTINE
     
                
                
        