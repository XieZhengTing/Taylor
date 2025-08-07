
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
    
      SUBROUTINE SECOND_MOMENT(VOL,NUMP,COO,WIN, &   !IN
                                NSNI_FAC) !OUT
      
      IMPLICIT NONE
      
      !GLOBAL IN-OUT
      DOUBLE PRECISION, INTENT(IN):: VOL(*)
      INTEGER, INTENT(IN):: NUMP
      DOUBLE PRECISION, INTENT(IN):: COO(3,*)
      DOUBLE PRECISION, INTENT(IN):: WIN(3,*)
      
      DOUBLE PRECISION, INTENT(OUT):: NSNI_FAC(3,*)
      
      !LOCAL
      INTEGER::I
      DOUBLE PRECISION:: THRD
      DATA THRD/0.333333333333333/
      DOUBLE PRECISION:: DIM_1, DIM_2, DIM_3, DIM_RATIO, VOL_WIN
      
      DO I=1,NUMP
      
      
         DIM_1 = WIN(1,I)*2.0d0
         DIM_2 = WIN(2,I)*2.0d0
         DIM_3 = WIN(3,I)*2.0d0
         
         VOL_WIN = DIM_1*DIM_2*DIM_3
         
         DIM_RATIO = VOL(I)**THRD / VOL_WIN **THRD
         
         !
         ! APPROXIMATE "H" IN EACH DIRECTION
         !
         DIM_1 = DIM_1*DIM_RATIO*0.5d0
         DIM_2 = DIM_2*DIM_RATIO*0.5d0
         DIM_3 = DIM_3*DIM_RATIO*0.5d0
         
         NSNI_FAC(1,I) = (DIM_1*2.0d0)**2/12.0d0
         NSNI_FAC(2,I) = (DIM_2*2.0d0)**2/12.0d0
         NSNI_FAC(3,I) = (DIM_3*2.0d0)**2/12.0d0
      
      END DO
      
      RETURN
      END SUBROUTINE