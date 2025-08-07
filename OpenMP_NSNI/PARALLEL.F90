
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
    
    
    
    
    SUBROUTINE ASSIGN_PARALLEL(HPC_SCHEME,MODEL_NUMP,LOCAL_NUMP)
                                
      IMPLICIT NONE
                                
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! OBTAIN THE SIZE OF LOCAL ARRAYS
	  !
    INTEGER, INTENT(IN):: HPC_SCHEME
    INTEGER, INTENT(IN):: MODEL_NUMP
    INTEGER, INTENT(OUT):: LOCAL_NUMP
    !
    ! LOCAL
    !
    INTEGER:: I
    
	  SELECT CASE (HPC_SCHEME)
	  
	  CASE(1)
	      !
	      ! SERIAL
	      !
      LOCAL_NUMP = MODEL_NUMP
      
	  CASE DEFAULT
	      !
	      !SOMETHING WENT WRONG
	      !
	      CALL EXIT_PROGRAM('INVALID HPC_SCHEME TYPE IN SUBROUTINE ASSEMBLER',0)
	      !
      END SELECT
		
      RETURN
		
      END SUBROUTINE
    
    
    
    
    SUBROUTINE ASSIGN_PARALLEL_MAP(HPC_SCHEME,MODEL_NUMP,LOCAL_NUMP,MODEL_MAP)
                                
      IMPLICIT NONE
                                
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! OBTAIN THE SIZE OF LOCAL ARRAYS, MAKE A MAP
	  !
    INTEGER, INTENT(IN):: HPC_SCHEME
    INTEGER, INTENT(IN):: LOCAL_NUMP
    INTEGER, INTENT(IN):: MODEL_NUMP
    INTEGER, INTENT(OUT):: MODEL_MAP(LOCAL_NUMP)
    !
    ! LOCAL
    !
    INTEGER:: I
    
	  SELECT CASE (HPC_SCHEME)
	  
	  CASE(1)
	      !
	      ! SERIAL: 1-to-1 mapping exactly
          !
      DO I=1,LOCAL_NUMP
      
         MODEL_MAP(I) = I
         
      END DO
      
      
	  CASE DEFAULT
	      !
	      !SOMETHING WENT WRONG
	      !
	      CALL EXIT_PROGRAM('INVALID HPC_SCHEME TYPE IN SUBROUTINE ASSEMBLER',0)
	      !
      END SELECT
		
      RETURN
		
      END SUBROUTINE
    
    
    
    SUBROUTINE PARALLEL_MODEL(MODEL_NUMP,   LOCAL_NUMP,   MODEL_MAP,                                     &
                              MODEL_SM_LEN, MODEL_SM_AREA, MODEL_SM_VOL,                                &
                              MODEL_WIN,    MODEL_VOL,     MODEL_NSNI_FAC, MODEL_VINIT, MODEL_MAT_TYPE, &
                              MODEL_PROP,   MODEL_BODY_ID,                                  &
                              LOCAL_SM_LEN, LOCAL_SM_AREA, LOCAL_SM_VOL,                                &
                              LOCAL_WIN,    LOCAL_VOL,     LOCAL_NSNI_FAC, LOCAL_VINIT, LOCAL_MAT_TYPE, &
                              LOCAL_PROP,   &
                              MODEL_XDIST_MAX, MODEL_YDIST_MAX, MODEL_ZDIST_MAX, &
                              LOCAL_XDIST_MAX, LOCAL_YDIST_MAX, LOCAL_ZDIST_MAX, &
                              MODEL_X_MOM, MODEL_Y_MOM, MODEL_Z_MOM , &
                              LOCAL_X_MOM, LOCAL_Y_MOM, LOCAL_Z_MOM , &   
                              LOCAL_BODY_ID)
                                
      IMPLICIT NONE
                                
	  !
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! ASSIGN THE GLOBAL VALUES TO THE LOCAL VALUES
	  !
                       
      INTEGER, INTENT(IN):: MODEL_NUMP
      DOUBLE PRECISION, INTENT(IN):: MODEL_SM_LEN(6,MODEL_NUMP)
      DOUBLE PRECISION, INTENT(IN):: MODEL_SM_AREA(3,MODEL_NUMP)
      DOUBLE PRECISION, INTENT(IN):: MODEL_SM_VOL(MODEL_NUMP)
      DOUBLE PRECISION, INTENT(IN):: MODEL_WIN(3,MODEL_NUMP)
      DOUBLE PRECISION, INTENT(IN):: MODEL_XDIST_MAX, MODEL_YDIST_MAX, MODEL_ZDIST_MAX 
      DOUBLE PRECISION, INTENT(IN):: MODEL_X_MOM(MODEL_NUMP)
      DOUBLE PRECISION, INTENT(IN):: MODEL_Y_MOM(MODEL_NUMP)
      DOUBLE PRECISION, INTENT(IN):: MODEL_Z_MOM(MODEL_NUMP)
      DOUBLE PRECISION, INTENT(IN):: MODEL_VOL(MODEL_NUMP)
      DOUBLE PRECISION, INTENT(IN):: MODEL_NSNI_FAC(3,MODEL_NUMP)
      DOUBLE PRECISION, INTENT(IN):: MODEL_VINIT(3,MODEL_NUMP)
      INTEGER, INTENT(IN):: MODEL_MAT_TYPE(MODEL_NUMP)
      DOUBLE PRECISION, INTENT(IN):: MODEL_PROP(30,MODEL_NUMP)
      INTEGER, INTENT(IN):: MODEL_BODY_ID(MODEL_NUMP)
      
      INTEGER, INTENT(IN):: LOCAL_NUMP
      INTEGER, INTENT(IN):: MODEL_MAP(LOCAL_NUMP)
      DOUBLE PRECISION, INTENT(OUT):: LOCAL_SM_LEN(6,LOCAL_NUMP)
      DOUBLE PRECISION, INTENT(OUT):: LOCAL_SM_AREA(3,LOCAL_NUMP)
      DOUBLE PRECISION, INTENT(OUT):: LOCAL_SM_VOL(LOCAL_NUMP)
      DOUBLE PRECISION, INTENT(OUT):: LOCAL_WIN(3,LOCAL_NUMP)
      DOUBLE PRECISION, INTENT(OUT):: LOCAL_XDIST_MAX, LOCAL_YDIST_MAX, LOCAL_ZDIST_MAX   
      DOUBLE PRECISION, INTENT(OUT):: LOCAL_X_MOM(LOCAL_NUMP)
      DOUBLE PRECISION, INTENT(OUT):: LOCAL_Y_MOM(LOCAL_NUMP)
      DOUBLE PRECISION, INTENT(OUT):: LOCAL_Z_MOM(LOCAL_NUMP)
      DOUBLE PRECISION, INTENT(OUT):: LOCAL_VOL(LOCAL_NUMP)
      DOUBLE PRECISION, INTENT(OUT):: LOCAL_NSNI_FAC(3,LOCAL_NUMP)
      DOUBLE PRECISION, INTENT(OUT):: LOCAL_VINIT(3,LOCAL_NUMP)
      INTEGER, INTENT(OUT):: LOCAL_MAT_TYPE(LOCAL_NUMP)
      DOUBLE PRECISION, INTENT(OUT):: LOCAL_PROP(30,LOCAL_NUMP)
      INTEGER, INTENT(OUT):: LOCAL_BODY_ID(LOCAL_NUMP)
      
      
      
      
      INTEGER:: I, II
      
      DO I = 1, LOCAL_NUMP
      
          II = MODEL_MAP(I)
      
          LOCAL_SM_LEN(1:6,I) = MODEL_SM_LEN(1:6,II)
          LOCAL_SM_AREA(1:3,I) = MODEL_SM_AREA(1:3,II)
          LOCAL_SM_VOL(I) = MODEL_SM_VOL(II)
          LOCAL_WIN(1:3,I) = MODEL_WIN(1:3,II)
          LOCAL_VOL(I) = MODEL_VOL(II)
          LOCAL_NSNI_FAC(1:3,I) = MODEL_NSNI_FAC(1:3,II)
          LOCAL_VINIT(1:3,I) = MODEL_VINIT(1:3,II)
          LOCAL_MAT_TYPE(I) = MODEL_MAT_TYPE(II)
          LOCAL_PROP(1:30,I) = MODEL_PROP(1:30,II)
          LOCAL_BODY_ID(I) = MODEL_BODY_ID(II)
          
          LOCAL_X_MOM(I) = MODEL_X_MOM(II) 
          LOCAL_Y_MOM(I) = MODEL_Y_MOM(II)           
          LOCAL_Z_MOM(I) = MODEL_Z_MOM(II)           
      
      END DO
      
      LOCAL_XDIST_MAX = MODEL_XDIST_MAX
      LOCAL_YDIST_MAX = MODEL_YDIST_MAX
      LOCAL_ZDIST_MAX = MODEL_ZDIST_MAX
      
      RETURN
      
      END SUBROUTINE
      