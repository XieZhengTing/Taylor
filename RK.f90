
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

    SUBROUTINE RK1(X, DEG, MSIZE, CONT, IMPL, GCOO, GWIN, GNUMP, LSTACK, LN, LNMAX, EBCS,SELF_EBC, &
        QL, QL_COEF,QL_LEN,  &
        SHP, SHPD,SHSUP)

    IMPLICIT NONE 

    DOUBLE PRECISION, INTENT(IN):: X(3)
    INTEGER, INTENT(IN):: DEG, MSIZE, CONT, GNUMP, LN, LNMAX
    INTEGER, INTENT(IN):: LSTACK(LNMAX)
    LOGICAL, INTENT(IN):: EBCS(GNUMP)
    LOGICAL, INTENT(IN):: SELF_EBC

    DOUBLE PRECISION, INTENT(IN):: GCOO(3,GNUMP),GWIN(3,GNUMP)

    LOGICAL, INTENT(IN):: QL, SHSUP
    DOUBLE PRECISION, INTENT(IN):: QL_COEF,QL_LEN

    DOUBLE PRECISION, INTENT(OUT):: SHP(LNMAX)
    DOUBLE PRECISION, INTENT(OUT):: SHPD(3,LNMAX)

    !LOCAL
    DOUBLE PRECISION:: MINV(MSIZE,MSIZE)
    DOUBLE PRECISION:: H(MSIZE,1)
    DOUBLE PRECISION:: M_FULL(MSIZE,MSIZE)
    DOUBLE PRECISION:: M_FULL_STAR(MSIZE,MSIZE)
    DOUBLE PRECISION:: B(1,MSIZE),H0(1,MSIZE),C(1,1)
    DOUBLE PRECISION:: BX(1,MSIZE),H0X(1,MSIZE),CX(1,1)
    DOUBLE PRECISION:: BY(1,MSIZE),H0Y(1,MSIZE),CY(1,1)
    DOUBLE PRECISION:: BZ(1,MSIZE),H0Z(1,MSIZE),CZ(1,1)

    DOUBLE PRECISION:: H_FULL(MSIZE,1)
    DOUBLE PRECISION:: H_FULL_STAR(MSIZE,1)
    DOUBLE PRECISION:: H_FULL_TEMP(MSIZE,1)
    DOUBLE PRECISION:: C_STAR(1,1),WINDOW_MOD

    ! FOR DIRECT GRADIENT KC
    DOUBLE PRECISION:: H_X(MSIZE,1),H_Y(MSIZE,1),H_Z(MSIZE,1)
    DOUBLE PRECISION:: M_X(MSIZE,MSIZE),M_Y(MSIZE,MSIZE),M_Z(MSIZE,MSIZE)
    DOUBLE PRECISION:: MINV_X(MSIZE,MSIZE),MINV_Y(MSIZE,MSIZE),MINV_Z(MSIZE,MSIZE)
    DOUBLE PRECISION:: MINV_X1(MSIZE,MSIZE),MINV_Y1(MSIZE,MSIZE),MINV_Z1(MSIZE,MSIZE)
    DOUBLE PRECISION:: PHI_X(LNMAX), PHI_Y(LNMAX), PHI_Z(LNMAX), PHIX_X, PHIY_Y, PHIZ_Z
    DOUBLE PRECISION:: DRDX, DRDY, DRDZ
    INTEGER, INTENT(IN):: IMPL
    !
    INTEGER:: I,II, J, K, M

    DOUBLE PRECISION:: PHIX, PHIY, PHIZ
    DOUBLE PRECISION:: PHI(LNMAX)
    DOUBLE PRECISION:: XMXI_OA(LNMAX)
    DOUBLE PRECISION:: YMYI_OA(LNMAX)
    DOUBLE PRECISION:: ZMZI_OA(LNMAX)
    DOUBLE PRECISION:: DIA(LNMAX)

    DOUBLE PRECISION:: TEST

    DOUBLE PRECISION:: PHI_SUM,QL_PTS(3,6), QLX(3)
    DOUBLE PRECISION:: XM_QLX,YM_QLY,ZM_QLZ
    DOUBLE PRECISION:: DENOM
    LOGICAL:: ISZERO

    !
    ! WE NEED TO BE CAREFUL NOT TO EVALUATE THE SINGULAR KERNAL AT THE NODE
    ! ALSO, WE NEED TO OUTPUT PHYSICAL DISPLACEMENTS!
    ! #TODO
    !
    IF (SELF_EBC) THEN

        DO I=1,LN

            II = LSTACK(I)

            XMXI_OA(I) = (X(1) -  GCOO(1,II)) /GWIN(1,II)
            YMYI_OA(I) = (X(2) -  GCOO(2,II)) /GWIN(2,II)
            ZMZI_OA(I) = (X(3) -  GCOO(3,II)) /GWIN(3,II)

            TEST = DSQRT(XMXI_OA(I)**2 + YMYI_OA(I)**2 + ZMZI_OA(I)**2)
            IF (TEST.LT.(1.0d-13)) THEN
                SHP(I) = 1.0d0
            ELSE
                SHP(I) = 0.0d0
            END IF

        END DO
        RETURN
    END IF
    !
    ! GET THE MOMENT MATRIX
    !
    !  dMx,dMy,dMz KC

    PHI_SUM = 0.0d0
    M_FULL = 0.0d0
    M_FULL_STAR = 0.0d0

    M_X = 0.0d0
    M_Y = 0.0d0
    M_Z = 0.0d0

    DO I=1,LN

        II = LSTACK(I)
        
        
                IF (II.EQ.42) THEN
                CONTINUE
                END IF
        !TODO: MAKE NORMALIZED
        XMXI_OA(I) = (X(1) -  GCOO(1,II))
        YMYI_OA(I) = (X(2) -  GCOO(2,II))
        ZMZI_OA(I) = (X(3) -  GCOO(3,II))

        CALL FILL_H(XMXI_OA(I),YMYI_OA(I),ZMZI_OA(I),MSIZE,DEG,H_FULL)

        CALL DERIV_H(XMXI_OA(I),YMYI_OA(I),ZMZI_OA(I),MSIZE,DEG,H_X,H_Y,H_Z)


        XMXI_OA(I) = (X(1) -  GCOO(1,II)) /GWIN(1,II)
        YMYI_OA(I) = (X(2) -  GCOO(2,II)) /GWIN(2,II)
        ZMZI_OA(I) = (X(3) -  GCOO(3,II)) /GWIN(3,II)

        IF (SHSUP) THEN
 !           DENOM = (GWIN(1,II)*GWIN(2,II)*GWIN(3,II))**(1.D0/3.D0)
            DIA(I) = SQRT(XMXI_OA(I)**2+YMYI_OA(I)**2+ZMZI_OA(I)**2)

            CALL MLS_KERNEL0(DIA(I),GWIN(1,II),CONT,PHI(I),PHIX_X,ISZERO)
            
!            IF (IMPL.EQ.1) THEN
!
!            ELSE
                IF  (DIA(I).LE.(1.0d-13)) THEN
                    
                   ! IF  (YMYI_OA(I).LE.(1.0d-13)) THEN
            
                        !IF  (ZMZI_OA(I).LE.(1.0d-13)) THEN
                            DRDX = 0.0d0
                            DRDY = 0.0d0
                            DRDZ = 0.0d0
                       ! ENDIF
                    !ENDIF
                ELSE
                    DRDX = (X(1) -  GCOO(1,II))/GWIN(1,II)**2/DIA(I)
                    DRDY = (X(2) -  GCOO(2,II))/GWIN(2,II)**2/DIA(I)
                    DRDZ = (X(3) -  GCOO(3,II))/GWIN(3,II)**2/DIA(I)
                ENDIF
                
                PHI_X(I) = PHIX_X*DRDX
                PHI_Y(I) = PHIX_X*DRDY
                PHI_Z(I) = PHIX_X*DRDZ
            ENDIF
            
            PHI_SUM = PHI_SUM + PHI(I)

!        ELSE


            CALL MLS_KERNEL0(ABS(XMXI_OA(I)),GWIN(1,II),CONT,PHIX,PHIX_X,ISZERO)
            CALL MLS_KERNEL0(ABS(YMYI_OA(I)),GWIN(2,II),CONT,PHIY,PHIY_Y,ISZERO)
            CALL MLS_KERNEL0(ABS(ZMZI_OA(I)),GWIN(3,II),CONT,PHIZ,PHIZ_Z,ISZERO)

            !DENOM = GWIN(1,II)*GWIN(2,II)*GWIN(3,II)
            PHI(I) = PHIX*PHIY*PHIZ !/DENOM



!            IF (IMPL.EQ.1) THEN
!
!            ELSE
                IF (XMXI_OA(I).EQ.0) THEN
                    DRDX = 0.0d0
                ELSEIF (XMXI_OA(I).GE.0) THEN
                    DRDX = 1.0d0/GWIN(1,II)
                ELSEIF (XMXI_OA(I).LE.0) THEN
                    DRDX = -1.0d0/GWIN(1,II)
                ENDIF

                IF (YMYI_OA(I).EQ.0) THEN
                    DRDY = 0.0d0
                ELSEIF (YMYI_OA(I).GE.0) THEN
                    DRDY = 1.0d0/GWIN(1,II)
                ELSEIF (YMYI_OA(I).LE.0) THEN
                    DRDY = -1.0d0/GWIN(1,II)
                ENDIF

                IF (ZMZI_OA(I).EQ.0) THEN
                    DRDZ = 0.0d0
                ELSEIF (ZMZI_OA(I).GE.0) THEN
                    DRDZ = 1.0d0/GWIN(1,II)
                ELSEIF (ZMZI_OA(I).LE.0) THEN
                    DRDZ = -1.0d0/GWIN(1,II)
                ENDIF

                PHI_X(I) = PHIX_X*PHIY*PHIZ*DRDX
                PHI_Y(I) = PHIX*PHIY_Y*PHIZ*DRDY
                PHI_Z(I) = PHIX*PHIY*PHIZ_Z*DRDZ
!            ENDIF


            !
            ! SINGULAR KERNAL
            ! #TODO
            !
            IF (EBCS(II)) THEN

                WINDOW_MOD = 1.0d0/(DSQRT(XMXI_OA(I)**2 + YMYI_OA(I)**2 + ZMZI_OA(I)**2) + 1.0E-015)

                ! FOR DSHP
                PHI_X(I) = PHI_X(I)*WINDOW_MOD - PHI(I)*WINDOW_MOD*WINDOW_MOD*WINDOW_MOD*XMXI_OA(I)/GWIN(1,II)  ! x
                PHI_Y(I) = PHI_Y(I)*WINDOW_MOD - PHI(I)*WINDOW_MOD*WINDOW_MOD*WINDOW_MOD*YMYI_OA(I)/GWIN(2,II)  ! y
                PHI_Z(I) = PHI_Z(I)*WINDOW_MOD - PHI(I)*WINDOW_MOD*WINDOW_MOD*WINDOW_MOD*ZMZI_OA(I)/GWIN(3,II)  ! z
                ! FOR SHAP
                PHI(I) = PHI(I)*WINDOW_MOD

            END IF
            PHI_SUM = PHI_SUM + PHI(I)
!        ENDIF
        DO J = 1, MSIZE
            DO K = 1, MSIZE
                M_FULL(J,K) = M_FULL(J,K) + H_FULL(J,1)*H_FULL(K,1)*PHI(I)
!                IF (IMPL.EQ.1) THEN
!                ELSE
                    M_X(J,K) = M_X(J,K) + H_X(J,1)*H_FULL(K,1)*PHI(I)+H_FULL(J,1)*H_X(K,1)*PHI(I)+H_FULL(J,1)*H_FULL(K,1)*PHI_X(I)
                    M_Y(J,K) = M_Y(J,K) + H_Y(J,1)*H_FULL(K,1)*PHI(I)+H_FULL(J,1)*H_Y(K,1)*PHI(I)+H_FULL(J,1)*H_FULL(K,1)*PHI_Y(I)
                    M_Z(J,K) = M_Z(J,K) + H_Z(J,1)*H_FULL(K,1)*PHI(I)+H_FULL(J,1)*H_Z(K,1)*PHI(I)+H_FULL(J,1)*H_FULL(K,1)*PHI_Z(I)
!                END IF
!
            END DO
        END DO
        
        CONTINUE

    END DO

    !
    ! GET M* IF QUASI-LINEAR
    !
!    IF (QL) THEN
!
!        J=1
!        QL_PTS(1,J) = X(1) + QL_LEN
!        QL_PTS(2,J) = X(2)
!        QL_PTS(3,J) = X(3)
!
!        J=2
!        QL_PTS(1,J) = X(1) - QL_LEN
!        QL_PTS(2,J) = X(2)
!        QL_PTS(3,J) = X(3)
!
!        J=3
!        QL_PTS(1,J) = X(1)
!        QL_PTS(2,J) = X(2) + QL_LEN
!        QL_PTS(3,J) = X(3)
!
!        J=4
!        QL_PTS(1,J) = X(1)
!        QL_PTS(2,J) = X(2) - QL_LEN
!        QL_PTS(3,J) = X(3)
!
!        J=5
!        QL_PTS(1,J) = X(1)
!        QL_PTS(2,J) = X(2)
!        QL_PTS(3,J) = X(3) + QL_LEN
!
!        J=6
!        QL_PTS(1,J) = X(1)
!        QL_PTS(2,J) = X(2)
!        QL_PTS(3,J) = X(3) - QL_LEN
!
!        DO M = 1, 6
!
!            QLX(:) = QL_PTS(:,M)
!
!            XM_QLX = X(1)  - QLX(1)
!            YM_QLY = X(2)  - QLX(2)
!            ZM_QLZ = X(3)  - QLX(3)
!
!            CALL FILL_H(XM_QLX,YM_QLY,ZM_QLZ,MSIZE,DEG,H_FULL)
!
!            DO J = 1, MSIZE
!                DO K = 1, MSIZE
!                    M_FULL_STAR(J,K) = M_FULL_STAR(J,K) + H_FULL(J,1)*H_FULL(K,1)
!                END DO
!            END DO
!
!        END DO
!
!        M_FULL = M_FULL + M_FULL_STAR * PHI_SUM * QL_COEF
!
!    END IF

    IF (MSIZE.EQ.4) THEN
    CALL M44INV(M_FULL, MINV)
    ELSE
    CALL INVERSE(M_FULL, MSIZE, MINV)
    END IF
    
    H0 = 0.0d0
    H0(1,1) = 1.0d0

    B = MATMUL(H0,MINV)

    IF (IMPL.EQ.1) THEN

        H0X = 0.0d0
        H0X(1,2) = -1.0d0

        H0Y = 0.0d0
        H0Y(1,3) = -1.0d0

        H0Z = 0.0d0
        H0Z(1,4) = -1.0d0

        BX = MATMUL(H0X,MINV)
        BY = MATMUL(H0Y,MINV)
        BZ = MATMUL(H0Z,MINV)

    ELSE
        MINV_X1 = MATMUL(-MINV,M_X)
        MINV_Y1 = MATMUL(-MINV,M_Y)
        MINV_Z1 = MATMUL(-MINV,M_Z)

        MINV_X = MATMUL(MINV_X1,MINV)
        MINV_Y = MATMUL(MINV_Y1,MINV)
        MINV_Z = MATMUL(MINV_Z1,MINV)

        BX = MATMUL(H0,MINV_X)
        BY = MATMUL(H0,MINV_Y)
        BZ = MATMUL(H0,MINV_Z)

    ENDIF

    DO I=1, LN

        II = LSTACK(I)

        !TODO: MAKE NORMALIZED
        XMXI_OA(I) = (X(1) -  GCOO(1,II))
        YMYI_OA(I) = (X(2) -  GCOO(2,II))
        ZMZI_OA(I) = (X(3) -  GCOO(3,II))

        CALL FILL_H(XMXI_OA(I),YMYI_OA(I),ZMZI_OA(I),MSIZE,DEG,H_FULL)

!        IF (QL) THEN
!
!            H_FULL_STAR = 0.0d0
!
!            DO M = 1, 6
!
!                QLX(:) = QL_PTS(:,M)
!
!                XM_QLX = X(1)  - QLX(1)
!                YM_QLY = X(2)  - QLX(2)
!                ZM_QLZ = X(3)  - QLX(3)
!
!                CALL FILL_H(XM_QLX,YM_QLY,ZM_QLZ,MSIZE,DEG,H_FULL_TEMP)
!
!                H_FULL_STAR = H_FULL_STAR + H_FULL_TEMP
!
!            END DO
!
!            H_FULL = H_FULL + H_FULL_STAR * QL_COEF
!
!        END IF


        C = MATMUL(B,H_FULL)

        SHP(I) = C(1,1)*PHI(I)

        IF (IMPL.EQ.1) THEN
!
            CX = MATMUL(BX,H_FULL)
            CY = MATMUL(BY,H_FULL)
            CZ = MATMUL(BZ,H_FULL)
!
            SHPD(1,I) = CX(1,1)*PHI(I)
            SHPD(2,I) = CY(1,1)*PHI(I)
            SHPD(3,I) = CZ(1,1)*PHI(I)

        ELSE
            CX = MATMUL(BX,H_FULL)+MATMUL(B,H_X)
            CY = MATMUL(BY,H_FULL)+MATMUL(B,H_Y)
            CZ = MATMUL(BZ,H_FULL)+MATMUL(B,H_Z)

            SHPD(1,I) = CX(1,1)*PHI(I)+C(1,1)*PHI_X(I)
            SHPD(2,I) = CY(1,1)*PHI(I)+C(1,1)*PHI_Y(I)
            SHPD(3,I) = CZ(1,1)*PHI(I)+C(1,1)*PHI_Z(I)
        ENDIF


    END DO

    !IF(ABS(SUM(SHP(1:LN)-1.D0).GT. 1.0E-05)) THEN
    !    WRITE(*,*) 'SHP ERROR: SUM = ', SUM(SHP(1:LN))
    !    PAUSE
    !    !STOP
    !END IF
    !
    !
    CX = SUM(SHPD(1,1:LN))
    CY = 0.0D0
    DO I = 1, LN
        II = LSTACK(I)
        CY = CY +SHPD(1,I)*GCOO(1,II)
    ENDDO
    CY = CY-1.0D0
    !IF(ABS(SUM(SHPD(1,1:LN)-0.D0).GT. 1.0E-03)) THEN
    !WRITE(*,*) 'SHPD1 ERROR: SUM = ', SUM(SHPD(1,1:LN))
    !PAUSE
    !STOP
    !END IF
    !
    !IF(ABS(SUM(SHPD(2,1:LN)-0.D0).GT. 1.0E-03)) THEN
    !WRITE(*,*) 'SHPD2 ERROR: SUM = ', SUM(SHPD(2,1:LN))
    !PAUSE
    !!STOP
    !END IF
    !
    !IF(ABS(SUM(SHPD(3,1:LN)-0.D0).GT. 1.0E-03)) THEN
    !WRITE(*,*) 'SHPD3 ERROR: SUM = ', SUM(SHPD(3,1:LN))
    !PAUSE
    !!STOP
    !END IF


    !
    !DEBUG SUBROUTINE
    !
    !CALL TESTER(X,SHP,SHPD,LN,LSTACK,GCOO)
    !
    !CONTINUE
    !
    !CALL TESTERX(X,SHP,SHPD,LN,LSTACK,GCOO)
    !
    !CONTINUE
    !

    RETURN
    END SUBROUTINE








    SUBROUTINE FILL_H(XMXI_OA,YMYI_OA,ZMZI_OA,MSIZE,DEG,H_FULL)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)::XMXI_OA,YMYI_OA,ZMZI_OA
    INTEGER, INTENT(IN)::MSIZE,DEG
    DOUBLE PRECISION, INTENT(OUT)::H_FULL(MSIZE,1)

!    IF (DEG.EQ.0) THEN
!
!        H_FULL(1,1) = 1.0d0
!
    IF (DEG.EQ.1) THEN

        H_FULL(1,1) = 1.0d0
        H_FULL(2,1) = XMXI_OA
        H_FULL(3,1) = YMYI_OA
        H_FULL(4,1) = ZMZI_OA

!    ELSEIF (DEG.EQ.2) THEN
!
!        H_FULL(1,1) = 1.0d0
!        H_FULL(2,1) = XMXI_OA
!        H_FULL(3,1) = YMYI_OA
!        H_FULL(4,1) = ZMZI_OA
!
!        H_FULL(5,1) = XMXI_OA*XMXI_OA
!        H_FULL(6,1) = YMYI_OA*XMXI_OA
!        H_FULL(7,1) = ZMZI_OA*XMXI_OA
!
!        H_FULL(8,1) = YMYI_OA*YMYI_OA
!        H_FULL(9,1) = YMYI_OA*ZMZI_OA
!
!        H_FULL(10,1) = ZMZI_OA*ZMZI_OA
!
    END IF
!
!    RETURN
    END SUBROUTINE



    ! Derivative H function

    SUBROUTINE DERIV_H(XMXI_OA,YMYI_OA,ZMZI_OA,MSIZE,DEG,H_X,H_Y,H_Z)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)::XMXI_OA,YMYI_OA,ZMZI_OA
    INTEGER, INTENT(IN)::MSIZE,DEG
    DOUBLE PRECISION, INTENT(OUT)::H_X(MSIZE,1),H_Y(MSIZE,1),H_Z(MSIZE,1)

!    IF (DEG.EQ.0) THEN
!
!        H_X(1,1) = 0.0d0
!        H_Y(1,1) = 0.0d0
!        H_Z(1,1) = 0.0d0
!
    IF (DEG.EQ.1) THEN

        H_X(1,1) = 0.0d0
        H_X(2,1) = 1.0d0
        H_X(3,1) = 0.0d0
        H_X(4,1) = 0.0d0

        H_Y(1,1) = 0.0d0
        H_Y(2,1) = 0.0d0
        H_Y(3,1) = 1.0d0
        H_Y(4,1) = 0.0d0

        H_Z(1,1) = 0.0d0
        H_Z(2,1) = 0.0d0
        H_Z(3,1) = 0.0d0
        H_Z(4,1) = 1.0d0

!    ELSEIF (DEG.EQ.2) THEN
!
!        H_X(1,1) = 0.0d0
!        H_X(2,1) = 1.0d0
!        H_X(3,1) = 0.0d0
!        H_X(4,1) = 0.0d0
!        H_X(5,1) = 2*XMXI_OA
!        H_X(6,1) = YMYI_OA
!        H_X(7,1) = ZMZI_OA
!        H_X(8,1) = 0.0d0
!        H_X(9,1) = 0.0d0
!        H_X(10,1) = 0.0d0
!
!        H_Y(1,1) = 0.0d0
!        H_Y(2,1) = 0.0d0
!        H_Y(3,1) = 1.0d0
!        H_Y(4,1) = 0.0d0
!        H_Y(5,1) = 0.0d0
!        H_Y(6,1) = XMXI_OA
!        H_Y(7,1) = 0.0d0
!        H_Y(8,1) = 2*YMYI_OA
!        H_Y(9,1) = ZMZI_OA
!        H_Y(10,1) = 0.0d0
!
!        H_Z(1,1) = 0.0d0
!        H_Z(2,1) = 0.0d0
!        H_Z(3,1) = 0.0d0
!        H_Z(4,1) = 1.0d0
!        H_Z(5,1) = 0.0d0
!        H_Z(6,1) = 0.0d0
!        H_Z(7,1) = XMXI_OA
!        H_Z(8,1) = 0.0d0
!        H_Z(9,1) = YMYI_OA
!        H_Z(10,1) = 2*ZMZI_OA

    END IF

    RETURN
    END SUBROUTINE



!    SUBROUTINE UDFM_SHAPE_TENSOR(X, DEG, MSIZE, CONT, IMPL, GCOO, GVOL, GWIN, GNUMP, LSTACK, LN, LNMAX, EBCS,SELF_EBC, &
!        QL, QL_COEF,QL_LEN,  &
!        SHP, INVK_MATX)
!
!    !
!    ! THIS SUBROUTINE IS TO FORM THE UNDEFORMED SHAPE TENSOR FOR PERIDYNAMICS AT THE FIRST STEP
!    !
!
!    IMPLICIT NONE
!
!    DOUBLE PRECISION, INTENT(IN):: X(3)
!    INTEGER, INTENT(IN):: DEG, MSIZE, CONT, GNUMP, LN, LNMAX
!    INTEGER, INTENT(IN):: LSTACK(LNMAX)
!    LOGICAL, INTENT(IN):: EBCS(GNUMP)
!    LOGICAL, INTENT(IN):: SELF_EBC
!
!    DOUBLE PRECISION, INTENT(IN):: GCOO(3,GNUMP),GWIN(3,GNUMP),GVOL(GNUMP)
!
!    LOGICAL, INTENT(IN):: QL
!    DOUBLE PRECISION, INTENT(IN):: QL_COEF,QL_LEN
!
!    DOUBLE PRECISION, INTENT(OUT):: SHP(LNMAX)
!    DOUBLE PRECISION, INTENT(OUT):: INVK_MATX(MSIZE-1,MSIZE-1)
!
!    !LOCAL
!    DOUBLE PRECISION:: MINV(MSIZE,MSIZE)
!    DOUBLE PRECISION:: H(MSIZE,1)
!    DOUBLE PRECISION:: M_FULL(MSIZE,MSIZE)
!    DOUBLE PRECISION:: M_FULL_STAR(MSIZE,MSIZE)
!    DOUBLE PRECISION:: B(1,MSIZE),H0(1,MSIZE),C(1,1)
!    DOUBLE PRECISION:: BX(1,MSIZE),H0X(1,MSIZE),CX(1,1)
!    DOUBLE PRECISION:: BY(1,MSIZE),H0Y(1,MSIZE),CY(1,1)
!    DOUBLE PRECISION:: BZ(1,MSIZE),H0Z(1,MSIZE),CZ(1,1)
!
!    DOUBLE PRECISION:: H_FULL(MSIZE,1)
!    DOUBLE PRECISION:: H_FULL_STAR(MSIZE,1)
!    DOUBLE PRECISION:: H_FULL_TEMP(MSIZE,1)
!    DOUBLE PRECISION::C_STAR(1,1)
!
!    ! FOR DIRECT GRADIENT KC
!    DOUBLE PRECISION:: H_X(MSIZE,1),H_Y(MSIZE,1),H_Z(MSIZE,1)
!    DOUBLE PRECISION:: M_X(MSIZE,MSIZE),M_Y(MSIZE,MSIZE),M_Z(MSIZE,MSIZE)
!    DOUBLE PRECISION:: MINV_X(MSIZE,MSIZE),MINV_Y(MSIZE,MSIZE),MINV_Z(MSIZE,MSIZE)
!    DOUBLE PRECISION:: MINV_X1(MSIZE,MSIZE),MINV_Y1(MSIZE,MSIZE),MINV_Z1(MSIZE,MSIZE)
!    DOUBLE PRECISION:: PHI_X(LNMAX), PHI_Y(LNMAX), PHI_Z(LNMAX), PHIX_X, PHIY_Y, PHIZ_Z
!    DOUBLE PRECISION:: DRDX, DRDY, DRDZ
!    INTEGER, INTENT(IN):: IMPL
!    !
!    INTEGER:: I,II, J, K, M
!
!    DOUBLE PRECISION:: PHIX, PHIY, PHIZ
!    DOUBLE PRECISION:: PHI(LNMAX)
!    DOUBLE PRECISION:: XMXI_OA(LNMAX)
!    DOUBLE PRECISION:: YMYI_OA(LNMAX)
!    DOUBLE PRECISION:: ZMZI_OA(LNMAX)
!
!    DOUBLE PRECISION:: TEST
!
!    DOUBLE PRECISION:: PHI_SUM,QL_PTS(3,6), QLX(3)
!    DOUBLE PRECISION:: XM_QLX,YM_QLY,ZM_QLZ
!    DOUBLE PRECISION:: DENOM
!
!    DOUBLE PRECISION:: K_MATX(MSIZE-1,MSIZE-1)
!    !
!    LOGICAL:: ISZERO
!
!    !
!    ! KEEP SOME OF THE STATEMENT FOR LATER USE, SUCH AS SINGULAR KERNEL, QL
!    ! #TODO
!    !
!    !IF (SELF_EBC) THEN
!    IF (.FALSE.) THEN
!        DO I=1,LN
!
!            II = LSTACK(I)
!
!            XMXI_OA(I) = (X(1) -  GCOO(1,II)) /GWIN(1,II)
!            YMYI_OA(I) = (X(2) -  GCOO(2,II)) /GWIN(2,II)
!            ZMZI_OA(I) = (X(3) -  GCOO(3,II)) /GWIN(3,II)
!
!            TEST = DSQRT(XMXI_OA(I)**2 + YMYI_OA(I)**2 + ZMZI_OA(I)**2)
!            IF (TEST.LT.(1.0d-13)) THEN
!                SHP(I) = 1.0d0
!            ELSE
!                SHP(I) = 0.0d0
!            END IF
!
!        END DO
!        RETURN
!    END IF
!    !
!    ! GET THE MOMENT MATRIX
!    !
!    !  dMx,dMy,dMz KC
!
!    PHI_SUM = 0.0d0
!    M_FULL = 0.0d0
!    M_FULL_STAR = 0.0d0
!
!    M_X = 0.0d0
!    M_Y = 0.0d0
!    M_Z = 0.0d0
!
!    K_MATX = 0.d0
!
!    DO I=1,LN
!
!        II = LSTACK(I)
!
!        XMXI_OA(I) = -(X(1) -  GCOO(1,II))
!        YMYI_OA(I) = -(X(2) -  GCOO(2,II))
!        ZMZI_OA(I) = -(X(3) -  GCOO(3,II))
!
!        CALL FILL_H(XMXI_OA(I),YMYI_OA(I),ZMZI_OA(I),MSIZE,DEG,H_FULL)
!
!        !CALL DERIV_H(XMXI_OA(I),YMYI_OA(I),ZMZI_OA(I),MSIZE,DEG,H_X,H_Y,H_Z)
!
!        XMXI_OA(I) = -(X(1) -  GCOO(1,II)) /GWIN(1,II)
!        YMYI_OA(I) = -(X(2) -  GCOO(2,II)) /GWIN(2,II)
!        ZMZI_OA(I) = -(X(3) -  GCOO(3,II)) /GWIN(3,II)
!
!        CALL MLS_KERNEL0(ABS(XMXI_OA(I)),GWIN(1,II),CONT,PHIX,PHIX_X,ISZERO)
!        CALL MLS_KERNEL0(ABS(YMYI_OA(I)),GWIN(2,II),CONT,PHIY,PHIY_Y,ISZERO)
!        CALL MLS_KERNEL0(ABS(ZMZI_OA(I)),GWIN(3,II),CONT,PHIZ,PHIZ_Z,ISZERO)
!
!        !DENOM = GWIN(1,II)*GWIN(2,II)*GWIN(3,II)
!        PHI(I) = PHIX*PHIY*PHIZ !/DENOM
!
!        DO J = 1,MSIZE-1
!            DO K = 1,MSIZE-1
!                ! STANDART PERIDYNAMICS, LINEAR BASIS
!                K_MATX(J,K) = K_MATX(J,K) + H_FULL(1+J,1)*H_FULL(1+K,1)*PHI(I)* GVOL(II)
!            ENDDO
!        ENDDO
!
!        !
!        ! STORE INFLUENCE FUNCTION INTO SHP
!        !
!        SHP(I) = PHI(I)
!
!    END DO
!
!
!    CALL INVERSE(K_MATX, MSIZE-1, INVK_MATX)
!
!
!    RETURN
!    END SUBROUTINE



!    SUBROUTINE DFM_SHAPE_TENSOR(X_0,X_t, DEG, MSIZE, CONT, GCOO, GVOL, GWIN, GNUMP, LSTACK, LN, LNMAX,  &
!        LCOO_CUURENT, SHP, S_MATX)
!
!    !
!    ! THIS SUBROUTINE IS TO FORM THE DEFORMED SHAPE TENSOR FOR PERIDYNAMICS AT THE FIRST STEP
!    !
!
!    IMPLICIT NONE
!
!    DOUBLE PRECISION, INTENT(IN):: X_0(3),X_t(3)
!    INTEGER, INTENT(IN):: DEG, MSIZE, CONT, GNUMP, LN, LNMAX
!    INTEGER, INTENT(IN):: LSTACK(LNMAX)
!
!    DOUBLE PRECISION, INTENT(IN):: GCOO(3,GNUMP),GWIN(3,GNUMP),GVOL(GNUMP),LCOO_CUURENT(3,LNMAX)
!
!
!    DOUBLE PRECISION, INTENT(OUT):: SHP(LNMAX)
!
!    DOUBLE PRECISION:: H_FULL(MSIZE,1)
!    DOUBLE PRECISION:: H_FULL_STAR(MSIZE,1)
!
!
!    ! FOR DIRECT GRADIENT KC
!    DOUBLE PRECISION:: PHI_X(LNMAX), PHI_Y(LNMAX), PHI_Z(LNMAX), PHIX_X, PHIY_Y, PHIZ_Z
!
!    !
!    INTEGER:: I,II, J, K, M
!
!    DOUBLE PRECISION:: PHI(LNMAX)
!    DOUBLE PRECISION:: XMXI_OA(LNMAX)
!    DOUBLE PRECISION:: YMYI_OA(LNMAX)
!    DOUBLE PRECISION:: ZMZI_OA(LNMAX)
!
!
!    DOUBLE PRECISION:: S_MATX(MSIZE-1,MSIZE-1)
!
!    !
!    ! GET THE MOMENT MATRIX
!    !
!    !  dMx,dMy,dMz KC
!
!
!    H_FULL = 0.0d0
!    H_FULL_STAR = 0.0d0
!
!    S_MATX = 0.d0
!
!    DO I=1,LN
!
!        II = LSTACK(I)
!
!        !TODO: MAKE NORMALIZED
!        XMXI_OA(I) = GCOO(1,II) - X_0(1)
!        YMYI_OA(I) = GCOO(2,II) - X_0(2)
!        ZMZI_OA(I) = GCOO(3,II) - X_0(3)
!
!        CALL FILL_H(XMXI_OA(I),YMYI_OA(I),ZMZI_OA(I),MSIZE,DEG,H_FULL)
!
!
!        !TODO: MAKE NORMALIZED
!        XMXI_OA(I) = LCOO_CUURENT(1,I) - X_t(1)
!        YMYI_OA(I) = LCOO_CUURENT(2,I) - X_t(2)
!        ZMZI_OA(I) = LCOO_CUURENT(3,I) - X_t(3)
!        CALL FILL_H(XMXI_OA(I),YMYI_OA(I),ZMZI_OA(I),MSIZE,DEG,H_FULL_STAR)
!
!        PHI(I) = SHP(I)
!
!        DO J = 1,MSIZE-1
!            DO K = 1,MSIZE-1
!                ! STANDART PERIDYNAMICS, LINEAR BASIS
!                S_MATX(J,K) = S_MATX(J,K) + H_FULL_STAR(1+J,1)*H_FULL(1+K,1)*PHI(I) * GVOL(II)
!            ENDDO
!        ENDDO
!
!    END DO
!
!    RETURN
!    END SUBROUTINE
