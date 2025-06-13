
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
    
      SUBROUTINE GET_SMOOTHING_LENGTHS(VOL,NUMP,COO,WIN, &   !IN
                                SM_LEN,SM_AREA,SM_VOL,XDIST_MAX, &
                                YDIST_MAX, ZDIST_MAX,&
                                X_MOM,Y_MOM,Z_MOM) !OUT
      
      IMPLICIT NONE
      
      !GLOBAL IN-OUT
      DOUBLE PRECISION, INTENT(IN):: VOL(NUMP)
      INTEGER, INTENT(IN):: NUMP
      DOUBLE PRECISION, INTENT(IN):: COO(3,NUMP)
      DOUBLE PRECISION, INTENT(IN):: WIN(3,NUMP)
      
      DOUBLE PRECISION, INTENT(OUT):: SM_LEN(6,NUMP)
      DOUBLE PRECISION, INTENT(OUT):: SM_AREA(3,NUMP)
      DOUBLE PRECISION, INTENT(OUT):: SM_VOL(NUMP)

      DOUBLE PRECISION, INTENT(OUT):: XDIST_MAX, YDIST_MAX,ZDIST_MAX   

      DOUBLE PRECISION, INTENT(OUT):: X_MOM(NUMP), Y_MOM(NUMP), Z_MOM(NUMP)       
      
      !LOCAL
      INTEGER::I
      DOUBLE PRECISION:: CRTV
      DOUBLE PRECISION:: THRD
      DATA THRD/0.333333333333333d0/
      DOUBLE PRECISION:: DIM_1, DIM_2, DIM_3, DIM_RATIO, VOL_WIN

      !GLOBAL MAX SMOOTHING LENGTH MEASURE FOR FINDING NEIGHBORS USE
      XDIST_MAX = 0.d0
      YDIST_MAX = 0.d0    
      ZDIST_MAX = 0.d0  
      
      DO I=1,NUMP
      
      
         DIM_1 = WIN(1,I)*2.0d0
         DIM_2 = WIN(2,I)*2.0d0
         DIM_3 = WIN(3,I)*2.0d0
         
         VOL_WIN = DIM_1*DIM_2*DIM_3
         
         DIM_RATIO = VOL(I)**THRD / VOL_WIN **THRD
         
         !
         ! COMPUTE THE MOMENT FOR NSNI
         !
         X_MOM(I) = (DIM_1*DIM_RATIO)**2.d0/12.d0
         Y_MOM(I) = (DIM_2*DIM_RATIO)**2.d0/12.d0
         Z_MOM(I) = (DIM_3*DIM_RATIO)**2.d0/12.d0
         
         DIM_RATIO = DIM_RATIO * 1.d0    !0.1d0
         !DIM_RATIO = DIM_RATIO * 0.1d0    !0.1d0         
         
         
         SM_LEN(1,I) = DIM_1 * DIM_RATIO /2.d0
         SM_LEN(2,I) = DIM_1 * DIM_RATIO /2.d0
         SM_LEN(3,I) = DIM_2 * DIM_RATIO /2.d0
         SM_LEN(4,I) = DIM_2 * DIM_RATIO /2.d0
         SM_LEN(5,I) = DIM_3 * DIM_RATIO /2.d0
         SM_LEN(6,I) = DIM_3 * DIM_RATIO /2.d0
         
         SM_VOL(I) = (SM_LEN(1,I) + SM_LEN(2,I)) * &
                     (SM_LEN(3,I) + SM_LEN(4,I)) * &
                     (SM_LEN(5,I) + SM_LEN(6,I)) 
         
         SM_AREA(1,I) = (SM_LEN(3,I) + SM_LEN(4,I)) * &
                        (SM_LEN(5,I) + SM_LEN(6,I)) 
                        
         SM_AREA(2,I) = (SM_LEN(1,I) + SM_LEN(2,I)) * &
                        (SM_LEN(5,I) + SM_LEN(6,I)) 
                        
         SM_AREA(3,I) = (SM_LEN(1,I) + SM_LEN(2,I)) * &
                        (SM_LEN(3,I) + SM_LEN(4,I)) 
                        
       !
       !
       !
       IF(SM_LEN(1,I) .GT. XDIST_MAX)  XDIST_MAX=SM_LEN(1,I)
       IF(SM_LEN(3,I) .GT. YDIST_MAX)  YDIST_MAX=SM_LEN(3,I)       
       IF(SM_LEN(6,I) .GT. ZDIST_MAX)  ZDIST_MAX=SM_LEN(6,I)                     
                        
      
      END DO
      
      RETURN
      END SUBROUTINE