
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
    
    
	SUBROUTINE INVERSE(A,N, AINV)
      !$ACC ROUTINE SEQ	
        IMPLICIT NONE

        INTERFACE
            SUBROUTINE WARN(STRING_LINE)

                CHARACTER(*) :: STRING_LINE
            END SUBROUTINE WARN
        END INTERFACE

        INTEGER, INTENT(IN):: N
        DOUBLE PRECISION, INTENT(IN):: A(N,N)
        DOUBLE PRECISION, INTENT(OUT):: AINV(N,N)
	
	
	INTEGER:: IS(N),JS(N), L, K, I, J
	DOUBLE PRECISION T,D
	
		AINV = A
		
	L=1
	DO 100 K=1,N
		D=0.0D0
		DO 10 I=K,N
		DO 10 J=K,N
			IF(DABS(AINV(I,J)).GT.D) THEN
			D=DABS(AINV(I,J))
			IS(K)=I
			JS(K)=J
		END IF
10	CONTINUE
		IF(D+1.0D0 .EQ. 1.0D0) THEN
			L=0
			!RETURN
			GOTO 300
		END IF
20	FORMAT(1X,'ERR * * NOT INV')
		DO 30 J=1,N
			T=AINV(K,J)
			AINV(K,J)=AINV(IS(K),J)
			AINV(IS(K),J)=T
30	CONTINUE
		DO 40 I=1,N
			T=AINV(I,K)
			AINV(I,K)=AINV(I,JS(K))
			AINV(I,JS(K))=T
40	CONTINUE
		AINV(K,K)=1.0D0/AINV(K,K)
		DO 50 J=1,N
			IF(J.NE.K) THEN
				AINV(K,J)=AINV(K,J)*AINV(K,K)
			END IF
50	CONTINUE
		DO 70 I=1,N
			IF(I.NE.K) THEN
			DO 60 J=1,N
				IF(J.NE.K) THEN
				AINV(I,J)=AINV(I,J)-AINV(I,K)*AINV(K,J)
				END IF
60	CONTINUE
			END IF
70	CONTINUE
		DO 80 I=1,N
			IF(I.NE.K) THEN
				AINV(I,K)=-AINV(I,K)*AINV(K,K)
			END IF
80	CONTINUE
100	CONTINUE
		DO 130 K=N,1,-1
			DO 110 J=1,N
				T=AINV(K,J)
				AINV(K,J)=AINV(JS(K),J)
				AINV(JS(K),J)=T
110	CONTINUE
		DO 120 I=1,N
			T=AINV(I,K)
			AINV(I,K)=AINV(I,IS(K))
			AINV(I,IS(K))=T
120	CONTINUE
130	CONTINUE

300     CONTINUE

        IF (L.EQ.0) THEN
          CALL WARN('PROBLEM INVERTING MATRIX')
        END IF

		RETURN
		
		END SUBROUTINE
		

        
        

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
    
      SUBROUTINE INV3 (A, AINV)
	 !$ACC ROUTINE SEQ 
      ! THIS SUBROUTINE HAS BUG, COMMENTED, USE  INVERSE()
   
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! COMPUTE THE INVERSE OF A 3X3 MATRIX
	  !
      IMPLICIT NONE
   
      DOUBLE PRECISION, INTENT(IN):: A(3,3)
      DOUBLE PRECISION, INTENT(OUT):: AINV(3,3)
      DOUBLE PRECISION:: DET
   
      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)
   
      AINV(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2))*DET
      AINV(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))*DET
      AINV(3,1) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1))*DET
      AINV(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))*DET
      AINV(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1))*DET
      AINV(3,2) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))*DET
      AINV(1,3) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2))*DET
      AINV(2,3) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))*DET
      AINV(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1))*DET
   
   
      RETURN
   
      END SUBROUTINE
   
   !!!
   !!!
   !!!   SUBROUTINE DET3 (A, DET)
   !!!
	  !!!! FUNCTION OF THIS SUBROUTINE:
	  !!!!
	  !!!! COMPUTE THE INVERSE OF A 3X3 MATRIX
	  !!!!
   !!!   IMPLICIT NONE
   !!!
   !!!   DOUBLE PRECISION, INTENT(IN):: A(3,3)
   !!!   DOUBLE PRECISION, INTENT(OUT):: DET
   !!!
   !!!   DET =   A(1,1)*A(2,2)*A(3,3)  &
   !!!         - A(1,1)*A(2,3)*A(3,2)  &
   !!!         - A(1,2)*A(2,1)*A(3,3)  &
   !!!         + A(1,2)*A(2,3)*A(3,1)  &
   !!!         + A(1,3)*A(2,1)*A(3,2)  &
   !!!         - A(1,3)*A(2,2)*A(3,1)
   !!!
   !!!
   !!!   RETURN
   !!!
   !!!   END SUBROUTINE
	

	
	



      SUBROUTINE M44INV (A, AINV)
      !$ACC ROUTINE SEQ
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(4,4), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(4,4), INTENT(OUT) :: AINV

      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(4,4) :: COFACTOR


      DET =  A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)- &
             A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+ &
             A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)- &
             A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+ &
             A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))

      COFACTOR(1,1) = A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))
      COFACTOR(1,2) = A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))
      COFACTOR(1,3) = A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
      COFACTOR(1,4) = A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
      COFACTOR(2,1) = A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))
      COFACTOR(2,2) = A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))
      COFACTOR(2,3) = A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2))
      COFACTOR(2,4) = A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))
      COFACTOR(3,1) = A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2))
      COFACTOR(3,2) = A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3))
      COFACTOR(3,3) = A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1))
      COFACTOR(3,4) = A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2))
      COFACTOR(4,1) = A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3))
      COFACTOR(4,2) = A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(4,3) = A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2))
      COFACTOR(4,4) = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      RETURN

      END SUBROUTINE