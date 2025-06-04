
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

	  SUBROUTINE HUGHES_WINGET(LMAT, & !IN
	                           ROT,STRAIN,D) !OUT
	  !$ACC ROUTINE SEQ
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! COMPUTE THE ROTATION AND STRAIN TENSORS USING
	  ! THE SO-CALLED HUGHES-WINGET ALGORITHM
	  !
      USE FINT_FUNCTIONS
      USE INVERSE_MOD
	  IMPLICIT NONE
	  !
	  !GLOBAL IN-OUT
	  DOUBLE PRECISION, INTENT(IN):: LMAT(3,3)
	  DOUBLE PRECISION, INTENT(OUT):: ROT(3,3),STRAIN(6)
      INTEGER :: ierr_inv
	  !
	  !LOCAL VARIABLES
	  DOUBLE PRECISION:: INV_LMAT(3,3)
	  DOUBLE PRECISION:: IDENT(3,3)
	  DOUBLE PRECISION:: A(3,3),IW(3,3),IW_INV(3,3)
	  DOUBLE PRECISION:: A_INV(3,3)
	  DOUBLE PRECISION:: G(3,3)
	  DOUBLE PRECISION:: G_T(3,3)
	  DOUBLE PRECISION:: E(3,3),W(3,3)
      DOUBLE PRECISION:: D(6)
	  DATA IDENT/ 1.0d0, 0.0d0, 0.0d0, &
	                  0.0d0, 1.0d0, 0.0d0, &
	                  0.0d0, 0.0d0, 1.0d0/
	  !
	  ! GET A=(I + 0.5*L)^-1
	  !
	  A = IDENT + 0.5d0*LMAT
	  !
      !CALL INVERSE(A, 3, A_INV)  
	  CALL INV3 (A, A_INV, ierr_inv)
      IF (ierr_inv /= 0) THEN 
          ROT = 0.0D0
          STRAIN = 0.0D0
          D = 0.0D0
          RETURN 
      END IF ! Basic error handling

	  !
	  ! GET G = L*A
	  !
	  G = MATMUL(LMAT,A_INV)
	  G_T=TRANSPOSE(G)
	  !
	  ! GET W = 1/2*(G-G^T)
	  ! GET E = 1/2*(G+G^T)
	  !
	  W =    0.5d0*(G - G_T)
      
      !ROT = I + (I-0.5D*W)^-1*W
      IW = IDENT-0.5D0*W
	  CALL INV3 (IW, IW_INV, ierr_inv)
      IF (ierr_inv /= 0) THEN 
          ROT = 0.0D0
          STRAIN = 0.0D0
          D = 0.0D0
          RETURN
      END IF ! Basic error handling

      !CALL INVERSE(IW, 3, IW_INV)       
      ROT = IDENT + MATMUL(IW_INV,W) 
      
      
	  E = 0.5d0*(G + G_T)
	  !
	  STRAIN = TENSOR_2_VTENSOR(E)
      D=STRAIN
	  !
      !TIMES FACT 2 FOR THE SHEAR COMPONENTS, TO BE CONSISTANT WITH THE STRAIN DEFINITION
      !
      STRAIN(4) = 2.D0*STRAIN(4) 
      STRAIN(5) = 2.D0*STRAIN(5)
      STRAIN(6) = 2.D0*STRAIN(6)      
      !WRITE(*,*)'FACT 2'
	  RETURN
	  END SUBROUTINE
	  

	  SUBROUTINE D_HUGHES_WINGET(LMAT,DLMAT, & !IN
	                           ROT,DSTRAIN) !OUT
	  !$ACC ROUTINE SEQ
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! COMPUTE THE ROTATION AND STRAIN TENSORS USING
	  ! THE SO-CALLED HUGHES-WINGET ALGORITHM
	  !
      USE FINT_FUNCTIONS
      USE INVERSE_MOD
      !
	  IMPLICIT NONE
	  !
	  !GLOBAL IN-OUT
	  DOUBLE PRECISION, INTENT(IN):: LMAT(3,3)
	  DOUBLE PRECISION, INTENT(IN):: DLMAT(3,3)
	  DOUBLE PRECISION, INTENT(OUT):: ROT(3,3),DSTRAIN(6)
      INTEGER :: ierr_inv
	  !
	  !LOCAL VARIABLES
	  DOUBLE PRECISION:: INV_LMAT(3,3)
	  DOUBLE PRECISION:: IDENT(3,3)
	  DOUBLE PRECISION:: A(3,3),IW(3,3),IW_INV(3,3),DA(3,3)
	  DOUBLE PRECISION:: A_INV(3,3), IDA(3,3), TEMP1(3,3), TEMP2(3,3)
	  DOUBLE PRECISION:: DG(3,3)
	  DOUBLE PRECISION:: DG_T(3,3)
	  DOUBLE PRECISION:: DE(3,3)
	  DATA IDENT/ 1.0d0, 0.0d0, 0.0d0, &
	                  0.0d0, 1.0d0, 0.0d0, &
	                  0.0d0, 0.0d0, 1.0d0/
	  !
	  ! GET A=(I + 0.5*L)^-1
	  !
	  A = IDENT + 0.5d0*LMAT
	  !
	  CALL INV3 (A, A_INV, ierr_inv)
      IF (ierr_inv /= 0) THEN 
          DSTRAIN = 0.0D0
          ROT = 0.0D0
          RETURN
      END IF ! Basic error handling

 
      
      !CALL INVERSE(A, 3, A_INV)
      
	  !
	  ! GET A,i
	  !
	  !!DA = 0.5d0*A
	  DA = 0.5d0*DLMAT      
      
	  !
	  ! GET INV(A,i)
	  !
	  IDA = MATMUL(-A_INV,DA)
	  IDA = MATMUL(IDA,A_INV)
	  !
	  ! GET TEMP MATS
	  !
	  !!TEMP1 = MATMUL(DLMAT,A)
	  !!TEMP2 = MATMUL(LMAT,DA)
      
	  TEMP1 = MATMUL(DLMAT,A_INV)
	  TEMP2 = MATMUL(LMAT,IDA)      
      
      
	  !
	  ! DG = DL*A + L*DA
	  !
	  DG = TEMP1 + TEMP2
	  DG_T=TRANSPOSE(DG)
	  !
	  ! GET DW = 1/2*(DG-DG^T)
	  ! GET DE = 1/2*(DG+DG^T)
	  !
	  DE = 0.5d0*(DG + DG_T)
	  !!DE = DG + DG_T      
      
	  !
	  DSTRAIN = TENSOR_2_VTENSOR(DE)
	  !
      !TIMES FACT 2 FOR THE SHEAR COMPONENTS, TO BE CONSISTANT WITH THE STRAIN DEFINITION
      !
      DSTRAIN(4) = 2.D0*DSTRAIN(4) 
      DSTRAIN(5) = 2.D0*DSTRAIN(5)
      DSTRAIN(6) = 2.D0*DSTRAIN(6)
      ROT = 0.0D0 ! Explicitly set ROT as it's an OUT parameter but not calculated
  
	  RETURN
	  END SUBROUTINE
	  