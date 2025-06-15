
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
MODULE RK_PROCEDURES_MOD
    USE INVERSE_MOD
    USE FINT_FUNCTIONS
    USE DETERMINANT_MOD
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: RK1
    PUBLIC :: FILL_H
    PUBLIC :: DERIV_H
    PUBLIC :: UDFM_SHAPE_TENSOR
    PUBLIC :: DFM_SHAPE_TENSOR
    PUBLIC :: MLS_KERNEL0
    PUBLIC :: HUGHES_WINGET 
    PUBLIC :: D_HUGHES_WINGET 
    PUBLIC :: ROTATE_TENSOR 

CONTAINS
    SUBROUTINE RK1(X, DEG, MSIZE, CONT, IMPL, GCOO, GWIN, GNUMP, LSTACK, LN, LNMAX, EBCS,SELF_EBC, &
        QL, QL_COEF,QL_LEN,  &
        SHP, SHPD,SHSUP)
    !$ACC ROUTINE SEQ
!    USE INVERSE_MOD
!    IMPLICIT NONE 

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
    DOUBLE PRECISION:: DET
    INTEGER :: ierr_inv, ierr_mls
    INTEGER :: VALID_NEIGHBORS
    DOUBLE PRECISION :: MIN_DIST, MAX_DIST, AVG_DIST
    DOUBLE PRECISION :: DIST_TEMP
    DOUBLE PRECISION :: PHI_SUM_RAW
! 簡單的除錯輸出（避免複雜格式）
    LOGICAL, SAVE :: FIRST_RK1 = .TRUE.
IF (FIRST_RK1) THEN
    WRITE(*,*) 'DEBUG RK1: First call'
    WRITE(*,*) '  LN = ', LN
    WRITE(*,*) '  CONT = ', CONT
    WRITE(*,*) '  SHSUP = ', SHSUP
    WRITE(*,*) '  QL = ', QL
    FIRST_RK1 = .FALSE.
END IF

    SHP(:) = 0.0D0
    SHPD(:,:) = 0.0D0
    
    ! 初始化局部數組
    PHI(:) = 0.0D0
    PHI_X(:) = 0.0D0
    PHI_Y(:) = 0.0D0
    PHI_Z(:) = 0.0D0
    
    ! 初始化計數器
    VALID_NEIGHBORS = 0
    !
    ! WE NEED TO BE CAREFUL NOT TO EVALUATE THE SINGULAR KERNAL AT THE NODE
    ! ALSO, WE NEED TO OUTPUT PHYSICAL DISPLACEMENTS!
    ! #TODO
    !
    !IF (SELF_EBC) THEN
    IF (.FALSE.) THEN
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
        !RETURN
    END IF
    !
    ! GET THE MOMENT MATRIX
    !
    !  dMx,dMy,dMz KC
    WRITE(*,*) 'DEBUG RK1: Entry check'
    WRITE(*,*) '  LN (number of neighbors) = ', LN
    WRITE(*,*) '  CONT = ', CONT
    WRITE(*,*) '  QL = ', QL
    WRITE(*,*) '  SHSUP = ', SHSUP
    PHI_SUM = 0.0d0
    M_FULL = 0.0d0
    M_FULL_STAR = 0.0d0

    M_X = 0.0d0
    M_Y = 0.0d0
    M_Z = 0.0d0

    ! 診斷視窗大小與節點間距的關係
    IF (LN > 0) THEN
        II = LSTACK(1)
        WRITE(*,*) 'DEBUG: First neighbor - GWIN = ', GWIN(1,II), GWIN(2,II), GWIN(3,II)
        
        ! 計算節點間距統計
        MIN_DIST = 1.0D10
        MAX_DIST = 0.0D0
        AVG_DIST = 0.0D0
        DO I = 1, MIN(10, LN)  ! 檢查前10個鄰居
            II = LSTACK(I)
            DIST_TEMP = SQRT((X(1)-GCOO(1,II))**2 + (X(2)-GCOO(2,II))**2 + (X(3)-GCOO(3,II))**2)
            MIN_DIST = MIN(MIN_DIST, DIST_TEMP)
            MAX_DIST = MAX(MAX_DIST, DIST_TEMP)
            AVG_DIST = AVG_DIST + DIST_TEMP
        END DO
        AVG_DIST = AVG_DIST / MIN(10, LN)
        
        ! 報告診斷資訊
        WRITE(*,*) 'DEBUG: Node spacing statistics (first 10 neighbors):'
        WRITE(*,*) '  Min distance = ', MIN_DIST
        WRITE(*,*) '  Max distance = ', MAX_DIST
        WRITE(*,*) '  Avg distance = ', AVG_DIST
        WRITE(*,*) '  Window/Avg ratio = ', GWIN(1,II)/AVG_DIST, GWIN(2,II)/AVG_DIST, GWIN(3,II)/AVG_DIST
        
        ! 警告視窗大小可能不適當
        IF (GWIN(1,II)/AVG_DIST < 3.0D0) THEN
            WRITE(*,*) 'WARNING: Window size may be too small relative to node spacing!'
        END IF
    END IF
    
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


        ! 檢查視窗大小有效性並設定最小值
        ! 防止過小的視窗導致 PHI_SUM 異常
        IF (GWIN(1,II) <= 1.0D-3 .OR. GWIN(2,II) <= 1.0D-3 .OR. GWIN(3,II) <= 1.0D-3) THEN
            ! 警告但不跳過，使用最小視窗大小
            IF (I <= 5) THEN
                WRITE(*,*) 'WARNING: Small GWIN detected at node ', II
                WRITE(*,*) '  GWIN = ', GWIN(1,II), GWIN(2,II), GWIN(3,II)
            END IF
            ! 可以選擇設定最小值或跳過
            ! 選項1：跳過
            ! PHI(I) = 0.0D0; PHI_X(I) = 0.0D0; PHI_Y(I) = 0.0D0; PHI_Z(I) = 0.0D0
            ! CYCLE
            ! 選項2：繼續計算但發出警告
        END IF
        
        ! 歸一化座標用於核函數評估
        XMXI_OA(I) = (X(1) - GCOO(1,II)) / GWIN(1,II)
        YMYI_OA(I) = (X(2) - GCOO(2,II)) / GWIN(2,II)
        ZMZI_OA(I) = (X(3) - GCOO(3,II)) / GWIN(3,II)

IF (SHSUP) THEN
            ! 球形支撐域
            DIA(I) = SQRT(XMXI_OA(I)**2 + YMYI_OA(I)**2 + ZMZI_OA(I)**2)
            CALL MLS_KERNEL0(DIA(I), CONT, PHI(I), PHIX_X, ISZERO, ierr_mls)
            
            ! 球形支撐域（不進行體積歸一化）
            !PHI(I) = PHI(I) / (AVG_WIN**3)

            IF (ierr_mls /= 0) THEN
                SHP(I) = 0.0D0
                SHPD(:,I) = 0.0D0
                CYCLE
            END IF

            ! 簡單除錯輸出
            IF (I <= 3) THEN
                WRITE(*,*) 'DEBUG RK1 SHSUP: I=', I, ' DIA=', DIA(I)
                WRITE(*,*) '  GWIN=', GWIN(1,II), GWIN(2,II), GWIN(3,II)
                WRITE(*,*) '  PHI=', PHI(I)
            END IF
            
            IF (IMPL.EQ.1) THEN
                ! Implicit formulation (not implemented)
                ELSE
                    DRDX = (X(1) - GCOO(1,II))/GWIN(1,II)**2/DIA(I)
                    DRDY = (X(2) - GCOO(2,II))/GWIN(2,II)**2/DIA(I)
                    DRDZ = (X(3) - GCOO(3,II))/GWIN(3,II)**2/DIA(I)
                ENDIF
                
                PHI_X(I) = PHIX_X*DRDX
                PHI_Y(I) = PHIX_X*DRDY
                PHI_Z(I) = PHIX_X*DRDZ
     
ELSE
            ! 張量積支撐域（座標已在前面歸一化）
            CALL MLS_KERNEL0(ABS(XMXI_OA(I)), CONT, PHIX, PHIX_X, ISZERO, ierr_mls)
            IF (ierr_mls /= 0) THEN
                SHP(I) = 0.0D0; SHPD(:,I) = 0.0D0; CYCLE;
            END IF
            CALL MLS_KERNEL0(ABS(YMYI_OA(I)), CONT, PHIY, PHIY_Y, ISZERO, ierr_mls)
            IF (ierr_mls /= 0) THEN
                SHP(I) = 0.0D0; SHPD(:,I) = 0.0D0; CYCLE;
            END IF
            CALL MLS_KERNEL0(ABS(ZMZI_OA(I)), CONT, PHIZ, PHIZ_Z, ISZERO, ierr_mls)
            IF (ierr_mls /= 0) THEN
                SHP(I) = 0.0D0; SHPD(:,I) = 0.0D0; CYCLE;
            END IF

            ! 計算張量積（不除以體積，與 OpenMP 版本一致）
            !DENOM = GWIN(1,II)*GWIN(2,II)*GWIN(3,II)
            PHI(I) = PHIX*PHIY*PHIZ !/DENOM
            
            ! 檢查歸一化座標是否在合理範圍內
            IF (I <= 3 .OR. (PHI(I) > 0.1D0 .AND. I <= 10)) THEN
                WRITE(*,*) 'DEBUG RK1 TENSOR: I=', I, ' Node=', II
                WRITE(*,*) '  Raw distance: ', (X(1)-GCOO(1,II)), (X(2)-GCOO(2,II)), (X(3)-GCOO(3,II))
                WRITE(*,*) '  Window size: ', GWIN(1,II), GWIN(2,II), GWIN(3,II)
                WRITE(*,*) '  Normalized: ', ABS(XMXI_OA(I)), ABS(YMYI_OA(I)), ABS(ZMZI_OA(I))
                WRITE(*,*) '  PHIX=', PHIX, ' PHIY=', PHIY, ' PHIZ=', PHIZ
                WRITE(*,*) '  PHI=', PHI(I)
                ! 計算預期的貢獻
                IF (LN > 0) THEN
                    WRITE(*,*) '  Expected contribution to sum: PHI*LN ≈ ', PHI(I)*LN
                END IF
            END IF

            IF (XMXI_OA(I).EQ.0) THEN
                DRDX = 0.0d0
            ELSE
                DRDX = SIGN(1.0D0, XMXI_OA(I)) / GWIN(1,II)
            ENDIF

            IF (YMYI_OA(I).EQ.0) THEN
                DRDY = 0.0d0
            ELSE
                DRDY = SIGN(1.0D0, YMYI_OA(I)) / GWIN(2,II)
            ENDIF

            IF (ZMZI_OA(I).EQ.0) THEN
                DRDZ = 0.0d0
            ELSE
                DRDZ = SIGN(1.0D0, ZMZI_OA(I)) / GWIN(3,II)
            ENDIF

            ! 應用 chain rule
            PHI_X(I) = PHIX_X * PHIY * PHIZ * DRDX
            PHI_Y(I) = PHIX * PHIY_Y * PHIZ * DRDY
            PHI_Z(I) = PHIX * PHIY * PHIZ_Z * DRDZ

            !
            ! SINGULAR KERNAL
            !
            IF (EBCS(II)) THEN
                WINDOW_MOD = 1.0d0/(DSQRT(XMXI_OA(I)**2 + YMYI_OA(I)**2 + ZMZI_OA(I)**2) + 1.0E-015)

                ! FOR DSHP
                PHI_X(I) = PHI_X(I)*WINDOW_MOD - PHI(I)*WINDOW_MOD*WINDOW_MOD*WINDOW_MOD*XMXI_OA(I)/GWIN(1,II)
                PHI_Y(I) = PHI_Y(I)*WINDOW_MOD - PHI(I)*WINDOW_MOD*WINDOW_MOD*WINDOW_MOD*YMYI_OA(I)/GWIN(2,II)
                PHI_Z(I) = PHI_Z(I)*WINDOW_MOD - PHI(I)*WINDOW_MOD*WINDOW_MOD*WINDOW_MOD*ZMZI_OA(I)/GWIN(3,II)
                ! FOR SHAP
                PHI(I) = PHI(I)*WINDOW_MOD
            END IF
        ENDIF
        
        ! 檢查 PHI 值的合理性
        IF (PHI(I) < 0.0D0 .OR. PHI(I) > 1.1D0) THEN
            ! PHI 應該在 [0, 1] 範圍內
            PHI(I) = 0.0D0
        END IF
        
        ! 檢查這個鄰居是否在支撐域內（PHI > 0）
        IF (PHI(I) > 1.0D-12) THEN
            VALID_NEIGHBORS = VALID_NEIGHBORS + 1
        END IF
        
        PHI_SUM = PHI_SUM + PHI(I)
        DO J = 1, MSIZE
            DO K = 1, MSIZE
                M_FULL(J,K) = M_FULL(J,K) + H_FULL(J,1)*H_FULL(K,1)*PHI(I)
!                IF (IMPL.EQ.1) THEN
!                ELSE
                    M_X(J,K) = M_X(J,K) + H_X(J,1)*H_FULL(K,1)*PHI(I)+H_FULL(J,1)*H_X(K,1)*PHI(I)+H_FULL(J,1)*H_FULL(K,1)*PHI_X(I)
                    M_Y(J,K) = M_Y(J,K) + H_Y(J,1)*H_FULL(K,1)*PHI(I)+H_FULL(J,1)*H_Y(K,1)*PHI(I)+H_FULL(J,1)*H_FULL(K,1)*PHI_Y(I)
                    M_Z(J,K) = M_Z(J,K) + H_Z(J,1)*H_FULL(K,1)*PHI(I)+H_FULL(J,1)*H_Z(K,1)*PHI(I)+H_FULL(J,1)*H_FULL(K,1)*PHI_Z(I)
!                END IF

            END DO
        END DO
    !WRITE(*,*) 'DEBUG RK1: After neighbor loop - PHI_SUM = ', PHI_SUM, ' LN = ', LN
    ! 移除過多的除錯輸出，只保留關鍵診斷
    IF (ABS(PHI_SUM) < 1.0D-10 .OR. ABS(PHI_SUM - 1.0D0) > 0.1D0) THEN
        WRITE(*,*) 'WARNING RK1: Abnormal PHI_SUM = ', PHI_SUM, ' LN = ', LN
    END IF
        CONTINUE

    END DO
    
    ! 儲存原始PHI_SUM用於診斷
    DOUBLE PRECISION :: PHI_SUM_RAW
    PHI_SUM_RAW = PHI_SUM
    
    ! 檢查並修正異常的 PHI_SUM
    IF (PHI_SUM > 1000.0D0) THEN
        WRITE(*,*) 'ERROR: Raw PHI_SUM = ', PHI_SUM_RAW, ' at node with LN = ', LN
        WRITE(*,*) '  Forcing normalization to avoid numerical issues'
    END IF
    
    ! 歸一化PHI值
    IF (ABS(PHI_SUM-1.0D0) > 1.0D-8 .AND. PHI_SUM > 1.0D-12) THEN
        DO I = 1, LN
            PHI(I) = PHI(I) / PHI_SUM
        END DO
        PHI_SUM = 1.0D0
    END IF
    
    ! 報告原始和歸一化後的值
    IF (ABS(PHI_SUM_RAW - 1.0D0) > 0.1D0) THEN
        WRITE(*,*) 'WARNING RK1: Raw PHI_SUM = ', PHI_SUM_RAW, ' LN = ', LN
        WRITE(*,*) '  After normalization: PHI_SUM = ', PHI_SUM
    END IF
    WRITE(*,*) 'DEBUG: Valid neighbors = ', VALID_NEIGHBORS, ' out of ', LN

    IF (ABS(PHI_SUM - 1.0D0) > 0.1D0) THEN
        WRITE(*,*) 'WARNING: Poor partition of unity, PHI_SUM = ', PHI_SUM
    END IF
    !
    ! GET M* IF QUASI-LINEAR
    !
    IF (QL) THEN

        J=1
        QL_PTS(1,J) = X(1) + QL_LEN
        QL_PTS(2,J) = X(2)
        QL_PTS(3,J) = X(3)

        J=2
        QL_PTS(1,J) = X(1) - QL_LEN
        QL_PTS(2,J) = X(2)
        QL_PTS(3,J) = X(3)

        J=3
        QL_PTS(1,J) = X(1)
        QL_PTS(2,J) = X(2) + QL_LEN
        QL_PTS(3,J) = X(3)

        J=4
        QL_PTS(1,J) = X(1)
        QL_PTS(2,J) = X(2) - QL_LEN
        QL_PTS(3,J) = X(3)

        J=5
        QL_PTS(1,J) = X(1)
        QL_PTS(2,J) = X(2)
        QL_PTS(3,J) = X(3) + QL_LEN

        J=6
        QL_PTS(1,J) = X(1)
        QL_PTS(2,J) = X(2)
        QL_PTS(3,J) = X(3) - QL_LEN

        DO M = 1, 6

            QLX(:) = QL_PTS(:,M)

            XM_QLX = X(1)  - QLX(1)
            YM_QLY = X(2)  - QLX(2)
            ZM_QLZ = X(3)  - QLX(3)

            CALL FILL_H(XM_QLX,YM_QLY,ZM_QLZ,MSIZE,DEG,H_FULL)

            DO J = 1, MSIZE
                DO K = 1, MSIZE
                    M_FULL_STAR(J,K) = M_FULL_STAR(J,K) + H_FULL(J,1)*H_FULL(K,1)
                END DO
            END DO

        END DO

        M_FULL = M_FULL + M_FULL_STAR * PHI_SUM * QL_COEF

    END IF

! 新增除錯：檢查 M_FULL 矩陣
IF (LN > 0) THEN
    CALL DETERMINANT(M_FULL, DET)
    IF (ABS(DET) < 1.0D-12) THEN
        WRITE(*,*) 'WARNING: RK1 - M_FULL nearly singular, DET = ', DET
        WRITE(*,*) '  PHI_SUM = ', PHI_SUM
        WRITE(*,*) '  LN = ', LN
    END IF
END IF

IF (PHI_SUM < 1.0D-10) THEN
    ! PHI_SUM 太小，矩陣可能奇異
    WRITE(*,*) 'ERROR: PHI_SUM too small = ', PHI_SUM, ' LN = ', LN
    ierr_inv = 1
ELSE
    IF (MSIZE.EQ.4) THEN
        CALL M44INV(M_FULL, MINV, ierr_inv)
    ELSE
        CALL INVERSE(M_FULL, MSIZE, MINV, ierr_inv)
    END IF
END IF
    
    ! 錯誤處理 - 確保不會因為 GPU 限制而改變行為
    IF (ierr_inv /= 0) THEN 
        ! 設定錯誤值但不停止執行（與 OpenMP 版本一致）
        DO I=1,LNMAX
            SHP(I) = 0.0D0
            SHPD(1,I) = 0.0D0
            SHPD(2,I) = 0.0D0
            SHPD(3,I) = 0.0D0
        END DO
        RETURN
    END IF
    
    H0 = 0.0d0
    H0(1,1) = 1.0d0

    B = MATMUL(H0,MINV)

!    IF (IMPL.EQ.1) THEN

!        H0X = 0.0d0
!        H0X(1,2) = -1.0d0

!        H0Y = 0.0d0
!        H0Y(1,3) = -1.0d0

!        H0Z = 0.0d0
!        H0Z(1,4) = -1.0d0

!        BX = MATMUL(H0X,MINV)
!        BY = MATMUL(H0Y,MINV)
!        BZ = MATMUL(H0Z,MINV)

!    ELSE
        MINV_X1 = MATMUL(-MINV,M_X)
        MINV_Y1 = MATMUL(-MINV,M_Y)
        MINV_Z1 = MATMUL(-MINV,M_Z)

        MINV_X = MATMUL(MINV_X1,MINV)
        MINV_Y = MATMUL(MINV_Y1,MINV)
        MINV_Z = MATMUL(MINV_Z1,MINV)

        BX = MATMUL(H0,MINV_X)
        BY = MATMUL(H0,MINV_Y)
        BZ = MATMUL(H0,MINV_Z)

!    ENDIF

    DO I=1, LN

        II = LSTACK(I)

        !TODO: MAKE NORMALIZED
        XMXI_OA(I) = (X(1) -  GCOO(1,II))
        YMYI_OA(I) = (X(2) -  GCOO(2,II))
        ZMZI_OA(I) = (X(3) -  GCOO(3,II))

        CALL FILL_H(XMXI_OA(I),YMYI_OA(I),ZMZI_OA(I),MSIZE,DEG,H_FULL)

        IF (QL) THEN

            H_FULL_STAR = 0.0d0

            DO M = 1, 6

                QLX(:) = QL_PTS(:,M)

                XM_QLX = X(1)  - QLX(1)
                YM_QLY = X(2)  - QLX(2)
                ZM_QLZ = X(3)  - QLX(3)

                CALL FILL_H(XM_QLX,YM_QLY,ZM_QLZ,MSIZE,DEG,H_FULL_TEMP)

                H_FULL_STAR = H_FULL_STAR + H_FULL_TEMP

            END DO

            H_FULL = H_FULL + H_FULL_STAR * QL_COEF

        END IF


        C = MATMUL(B,H_FULL)

        SHP(I) = C(1,1)*PHI(I)


    IF (I <= 3) THEN
        WRITE(*,*) 'DEBUG RK1 SHP: I=', I
        WRITE(*,*) '  C11=', C(1,1), ' PHI=', PHI(I)
        WRITE(*,*) '  SHP=', SHP(I)
    END IF
!        IF (IMPL.EQ.1) THEN

!            CX = MATMUL(BX,H_FULL)
!            CY = MATMUL(BY,H_FULL)
!            CZ = MATMUL(BZ,H_FULL)

!            SHPD(1,I) = CX(1,1)*PHI(I)
!            SHPD(2,I) = CY(1,1)*PHI(I)
!            SHPD(3,I) = CZ(1,1)*PHI(I)

!        ELSE
            CX = MATMUL(BX,H_FULL)+MATMUL(B,H_X)
            CY = MATMUL(BY,H_FULL)+MATMUL(B,H_Y)
            CZ = MATMUL(BZ,H_FULL)+MATMUL(B,H_Z)

            SHPD(1,I) = CX(1,1)*PHI(I)+C(1,1)*PHI_X(I)
            SHPD(2,I) = CY(1,1)*PHI(I)+C(1,1)*PHI_Y(I)
            SHPD(3,I) = CZ(1,1)*PHI(I)+C(1,1)*PHI_Z(I)
!        ENDIF


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
    !$ACC ROUTINE SEQ
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)::XMXI_OA,YMYI_OA,ZMZI_OA
    INTEGER, INTENT(IN)::MSIZE,DEG
    DOUBLE PRECISION, INTENT(OUT)::H_FULL(MSIZE,1)

!    IF (DEG.EQ.0) THEN

!        H_FULL(1,1) = 1.0d0

    IF (DEG.EQ.1) THEN ! This is the active branch for DEG.EQ.1

        H_FULL(1,1) = 1.0d0
        H_FULL(2,1) = XMXI_OA
        H_FULL(3,1) = YMYI_OA
        H_FULL(4,1) = ZMZI_OA

!    ELSEIF (DEG.EQ.2) THEN

!        H_FULL(1,1) = 1.0d0
!        H_FULL(2,1) = XMXI_OA
!        H_FULL(3,1) = YMYI_OA
!        H_FULL(4,1) = ZMZI_OA

!        H_FULL(5,1) = XMXI_OA*XMXI_OA
!        H_FULL(6,1) = YMYI_OA*XMXI_OA
!        H_FULL(7,1) = ZMZI_OA*XMXI_OA

!        H_FULL(8,1) = YMYI_OA*YMYI_OA
!        H_FULL(9,1) = YMYI_OA*ZMZI_OA

!        H_FULL(10,1) = ZMZI_OA*ZMZI_OA

    END IF

    RETURN
    END SUBROUTINE



    ! Derivative H function

    SUBROUTINE DERIV_H(XMXI_OA,YMYI_OA,ZMZI_OA,MSIZE,DEG,H_X,H_Y,H_Z)
    !$ACC ROUTINE SEQ
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)::XMXI_OA,YMYI_OA,ZMZI_OA
    INTEGER, INTENT(IN)::MSIZE,DEG
    DOUBLE PRECISION, INTENT(OUT)::H_X(MSIZE,1),H_Y(MSIZE,1),H_Z(MSIZE,1)

!    IF (DEG.EQ.0) THEN

!        H_X(1,1) = 0.0d0
!        H_Y(1,1) = 0.0d0
!        H_Z(1,1) = 0.0d0

    IF (DEG.EQ.1) THEN ! This is the active branch for DEG.EQ.1

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



    SUBROUTINE UDFM_SHAPE_TENSOR(X, DEG, MSIZE, CONT, IMPL, GCOO, GVOL, GWIN, GNUMP, LSTACK, LN, LNMAX, EBCS,SELF_EBC, &
        QL, QL_COEF,QL_LEN,  &
        SHP, INVK_MATX)
    !$ACC ROUTINE SEQ
!    USE INVERSE_MOD
    !
    ! THIS SUBROUTINE IS TO FORM THE UNDEFORMED SHAPE TENSOR FOR PERIDYNAMICS AT THE FIRST STEP
    !

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: X(3)
    INTEGER, INTENT(IN):: DEG, MSIZE, CONT, GNUMP, LN, LNMAX
    INTEGER, INTENT(IN):: LSTACK(LNMAX)
    LOGICAL, INTENT(IN):: EBCS(GNUMP)
    LOGICAL, INTENT(IN):: SELF_EBC

    DOUBLE PRECISION, INTENT(IN):: GCOO(3,GNUMP),GWIN(3,GNUMP),GVOL(GNUMP)

    LOGICAL, INTENT(IN):: QL
    DOUBLE PRECISION, INTENT(IN):: QL_COEF,QL_LEN

    DOUBLE PRECISION, INTENT(OUT):: SHP(LNMAX)
    DOUBLE PRECISION, INTENT(OUT):: INVK_MATX(MSIZE-1,MSIZE-1)

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
    DOUBLE PRECISION::C_STAR(1,1)

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

    DOUBLE PRECISION:: K_MATX(MSIZE-1,MSIZE-1)
    !
    LOGICAL:: ISZERO
    INTEGER :: ierr_inv, ierr_mls
    INTEGER :: VALID_NEIGHBORS
    !
    ! KEEP SOME OF THE STATEMENT FOR LATER USE, SUCH AS SINGULAR KERNEL, QL
    ! #TODO
    ierr_inv = 0
    ierr_mls = 0
    VALID_NEIGHBORS = 0
    !IF (SELF_EBC) THEN
    IF (.FALSE.) THEN
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
    
    ! 初始化所有 PHI 相關陣列
    DO I = 1, LNMAX
        PHI(I) = 0.0D0
        PHI_X(I) = 0.0D0
        PHI_Y(I) = 0.0D0
        PHI_Z(I) = 0.0D0
        XMXI_OA(I) = 0.0D0
        YMYI_OA(I) = 0.0D0
        ZMZI_OA(I) = 0.0D0
        DIA(I) = 0.0D0
    END DO

    K_MATX = 0.d0

    DO I=1,LN

        II = LSTACK(I)

        XMXI_OA(I) = -(X(1) -  GCOO(1,II))
        YMYI_OA(I) = -(X(2) -  GCOO(2,II))
        ZMZI_OA(I) = -(X(3) -  GCOO(3,II))

        CALL FILL_H(XMXI_OA(I),YMYI_OA(I),ZMZI_OA(I),MSIZE,DEG,H_FULL)

        !CALL DERIV_H(XMXI_OA(I),YMYI_OA(I),ZMZI_OA(I),MSIZE,DEG,H_X,H_Y,H_Z)

        XMXI_OA(I) = -(X(1) -  GCOO(1,II)) /GWIN(1,II)
        YMYI_OA(I) = -(X(2) -  GCOO(2,II)) /GWIN(2,II)
        ZMZI_OA(I) = -(X(3) -  GCOO(3,II)) /GWIN(3,II)

        CALL MLS_KERNEL0(ABS(XMXI_OA(I)), CONT,PHIX,PHIX_X,ISZERO, ierr_mls)
        IF (ierr_mls /= 0) THEN
            SHP(I) = 0.0D0; CYCLE;
        END IF
        CALL MLS_KERNEL0(ABS(YMYI_OA(I)), CONT,PHIY,PHIY_Y,ISZERO, ierr_mls)
        IF (ierr_mls /= 0) THEN
            SHP(I) = 0.0D0; CYCLE;
        END IF
        CALL MLS_KERNEL0(ABS(ZMZI_OA(I)), CONT,PHIZ,PHIZ_Z,ISZERO, ierr_mls)
        IF (ierr_mls /= 0) THEN 
            ! Handle MLS_KERNEL0 error
            SHP(I) = HUGE(0.0D0); CYCLE;
        END IF

! 不除以體積（與 OpenMP 版本一致）
PHI(I) = PHIX*PHIY*PHIZ

            IF (I <= 5) THEN
                WRITE(*,*) 'DEBUG UDFM_SHAPE_TENSOR: I=', I, ' II=', II
                WRITE(*,*) '  XMXI_OA=', ABS(XMXI_OA(I)), ' PHIX=', PHIX
                WRITE(*,*) '  YMYI_OA=', ABS(YMYI_OA(I)), ' PHIY=', PHIY
                WRITE(*,*) '  ZMZI_OA=', ABS(ZMZI_OA(I)), ' PHIZ=', PHIZ
                WRITE(*,*) '  PHI=', PHI(I)
            END IF
        IF (PHI(I) > 1.0D-12) THEN
            VALID_NEIGHBORS = VALID_NEIGHBORS + 1
        END IF
            PHI_SUM = PHI_SUM + PHI(I)

        DO J = 1,MSIZE-1
            DO K = 1,MSIZE-1
                ! STANDART PERIDYNAMICS, LINEAR BASIS
                K_MATX(J,K) = K_MATX(J,K) + H_FULL(1+J,1)*H_FULL(1+K,1)*PHI(I)* GVOL(II)
            ENDDO
        ENDDO

        !
        ! STORE INFLUENCE FUNCTION INTO SHP
        !
        SHP(I) = PHI(I)

    END DO


    CALL INVERSE(K_MATX, MSIZE-1, INVK_MATX, ierr_inv)
    IF (ierr_inv /= 0) THEN 
        ! 處理 K_MATX 求逆失敗的情況 (e.g., INVK_MATX = HUGE(0.0D0); RETURN)
        INVK_MATX = HUGE(0.0D0); RETURN;
    END IF



    RETURN
    END SUBROUTINE



    SUBROUTINE DFM_SHAPE_TENSOR(X_0,X_t, DEG, MSIZE, CONT, GCOO, GVOL, GWIN, GNUMP, LSTACK, LN, LNMAX,  &
        LCOO_CUURENT, SHP, S_MATX)
    !$ACC ROUTINE SEQ
    !
    ! THIS SUBROUTINE IS TO FORM THE DEFORMED SHAPE TENSOR FOR PERIDYNAMICS AT THE FIRST STEP
    !

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: X_0(3),X_t(3)
    INTEGER, INTENT(IN):: DEG, MSIZE, CONT, GNUMP, LN, LNMAX
    INTEGER, INTENT(IN):: LSTACK(LNMAX)

    DOUBLE PRECISION, INTENT(IN):: GCOO(3,GNUMP),GWIN(3,GNUMP),GVOL(GNUMP),LCOO_CUURENT(3,LNMAX)


    DOUBLE PRECISION, INTENT(OUT):: SHP(LNMAX)

    DOUBLE PRECISION:: H_FULL(MSIZE,1)
    DOUBLE PRECISION:: H_FULL_STAR(MSIZE,1)


    ! FOR DIRECT GRADIENT KC
    DOUBLE PRECISION:: PHI_X(LNMAX), PHI_Y(LNMAX), PHI_Z(LNMAX), PHIX_X, PHIY_Y, PHIZ_Z

    !
    INTEGER:: I,II, J, K, M

    DOUBLE PRECISION:: PHI(LNMAX)
    DOUBLE PRECISION:: XMXI_OA(LNMAX)
    DOUBLE PRECISION:: YMYI_OA(LNMAX)
    DOUBLE PRECISION:: ZMZI_OA(LNMAX)


    DOUBLE PRECISION:: S_MATX(MSIZE-1,MSIZE-1)

    !
    ! GET THE MOMENT MATRIX
    !
    !  dMx,dMy,dMz KC


    H_FULL = 0.0d0
    H_FULL_STAR = 0.0d0

    S_MATX = 0.d0

    DO I=1,LN

        II = LSTACK(I)

        !TODO: MAKE NORMALIZED
        XMXI_OA(I) = GCOO(1,II) - X_0(1)
        YMYI_OA(I) = GCOO(2,II) - X_0(2)
        ZMZI_OA(I) = GCOO(3,II) - X_0(3)

        CALL FILL_H(XMXI_OA(I),YMYI_OA(I),ZMZI_OA(I),MSIZE,DEG,H_FULL)


        !TODO: MAKE NORMALIZED
        XMXI_OA(I) = LCOO_CUURENT(1,I) - X_t(1)
        YMYI_OA(I) = LCOO_CUURENT(2,I) - X_t(2)
        ZMZI_OA(I) = LCOO_CUURENT(3,I) - X_t(3)
        CALL FILL_H(XMXI_OA(I),YMYI_OA(I),ZMZI_OA(I),MSIZE,DEG,H_FULL_STAR)

        PHI(I) = SHP(I)

        DO J = 1,MSIZE-1
            DO K = 1,MSIZE-1
                ! STANDART PERIDYNAMICS, LINEAR BASIS
                S_MATX(J,K) = S_MATX(J,K) + H_FULL_STAR(1+J,1)*H_FULL(1+K,1)*PHI(I) * GVOL(II)
            ENDDO
        ENDDO

    END DO

    RETURN
    END SUBROUTINE
SUBROUTINE MLS_KERNEL0(XN,CONT,PHI,PHIX,ISZERO, ierr)
    !$ACC ROUTINE SEQ
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN):: XN  ! 已歸一化的距離
    INTEGER, INTENT(IN):: CONT
    DOUBLE PRECISION, INTENT(OUT):: PHI,PHIX
    LOGICAL, INTENT(OUT):: ISZERO
    INTEGER, INTENT(OUT) :: ierr
    
    DOUBLE PRECISION :: XSA
    
    ISZERO = .FALSE.
    ierr = 0
    
    XSA = XN  ! 輸入已經歸一化
    
    IF (CONT.EQ.3) THEN !CUBIC SPLINE
        IF (XSA.LE.0.5D0) THEN
            PHI = 2.0D0/3.0D0 - 4.0D0*XSA**2 + 4.0D0*XSA**3
            IF (XSA .GT. 1.0D-12) THEN
                PHIX = -8.0D0*XSA + 12.0D0*XSA**2
            ELSE
                PHIX = 0.0D0
            END IF
        ELSEIF (XSA.LE.1.0D0) THEN
            PHI = 4.0D0/3.0D0 - 4.0D0*XSA + 4.0D0*XSA**2 - 4.0D0/3.0D0*XSA**3
            PHIX = -4.0D0 + 8.0D0*XSA - 4.0D0*XSA**2
        ELSE
            PHI = 0.0D0
            PHIX = 0.0D0
            ISZERO = .TRUE.
        END IF
        
        ! 體積歸一化（如果需要）
        ! 對於3D張量積：體積 = (2*WIN)^3 = 8*WIN^3
        ! 對於球形支撐：體積 = 4/3*PI*WIN^3
        ! 目前暫時保留未歸一化，等待與OpenMP版本確認
        ! IF (NEED_VOLUME_NORMALIZATION) THEN
        !     DOUBLE PRECISION :: VOL_FACTOR
        !     VOL_FACTOR = 8.0D0  ! 張量積情況
        !     PHI = PHI / VOL_FACTOR
        !     PHIX = PHIX / VOL_FACTOR
        ! END IF
        
    ELSE
        ! 不支援的核函數類型
        PHI = 0.0D0
        PHIX = 0.0D0
        ierr = 1
    END IF
    
    RETURN
END SUBROUTINE MLS_KERNEL0


	  SUBROUTINE HUGHES_WINGET(LMAT, & !IN
                        ROT,STRAIN,D) !OUT
	  !$ACC ROUTINE SEQ
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! COMPUTE THE ROTATION AND STRAIN TENSORS USING
	  ! THE SO-CALLED HUGHES-WINGET ALGORITHM
	  !
      ! USE FINT_FUNCTIONS ! 已在模組層級 USE
      ! IMPLICIT NONE ! 已在模組層級定義
	  !
	  !GLOBAL IN-OUT
	  DOUBLE PRECISION, INTENT(IN):: LMAT(3,3)
	  DOUBLE PRECISION, INTENT(OUT):: ROT(3,3),STRAIN(6)
      INTEGER :: ierr_inv
	  !LOCAL VARIABLES
	  DOUBLE PRECISION:: INV_LMAT(3,3) ! 未被使用
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
	  CALL INV3 (A, A_INV, ierr_inv) ! INV3 來自 INVERSE_MOD
      IF (ierr_inv /= 0) THEN 
        ROT = 0.0D0; STRAIN = 0.0D0; D = 0.0D0; RETURN;
      END IF



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
        ROT = 0.0D0; STRAIN = 0.0D0; D = 0.0D0; RETURN;
      END IF


      !CALL INVERSE(IW, 3, IW_INV)
      ROT = IDENT + MATMUL(IW_INV,W)
      
      
	  E = 0.5d0*(G + G_T)
	  !
	  STRAIN = TENSOR_2_VTENSOR(E) ! TENSOR_2_VTENSOR 來自 FINT_FUNCTIONS
      D=STRAIN
	  !
      !TIMES FACT 2 FOR THE SHEAR COMPONENTS, TO BE CONSISTANT WITH THE STRAIN DEFINITION
      !
      STRAIN(4) = 2.D0*STRAIN(4)
      STRAIN(5) = 2.D0*STRAIN(5)
      STRAIN(6) = 2.D0*STRAIN(6)
      !WRITE(*,*)'FACT 2'
	  RETURN
	  END SUBROUTINE HUGHES_WINGET
	  
    
	  SUBROUTINE D_HUGHES_WINGET(LMAT,DLMAT, & !IN
	                           ROT,DSTRAIN) !OUT
	  !$ACC ROUTINE SEQ
	  ! FUNCTION OF THIS SUBROUTINE:
	  !
	  ! COMPUTE THE ROTATION AND STRAIN TENSORS USING
	  ! THE SO-CALLED HUGHES-WINGET ALGORITHM
	  !
      ! USE FINT_FUNCTIONS ! 已在模組層級 USE
      ! IMPLICIT NONE ! 已在模組層級定義
	  !
	  !GLOBAL IN-OUT
	  DOUBLE PRECISION, INTENT(IN):: LMAT(3,3)
	  DOUBLE PRECISION, INTENT(IN):: DLMAT(3,3)
	  DOUBLE PRECISION, INTENT(OUT):: ROT(3,3) ! ROT 未被賦值
	  DOUBLE PRECISION, INTENT(OUT):: DSTRAIN(6)
      INTEGER :: ierr_inv
	  !
	  !LOCAL VARIABLES
	  DOUBLE PRECISION:: INV_LMAT(3,3) ! 未被使用
	  DOUBLE PRECISION:: IDENT(3,3)
	  DOUBLE PRECISION:: A(3,3),DA(3,3) ! IW, IW_INV 在此未使用
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
	  CALL INV3 (A, A_INV, ierr_inv) ! INV3 來自 INVERSE_MOD
      IF (ierr_inv /= 0) THEN 
        DSTRAIN = 0.0D0; ROT = 0.0D0; RETURN;
      END IF


      
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
	  DSTRAIN = TENSOR_2_VTENSOR(DE) ! TENSOR_2_VTENSOR 來自 FINT_FUNCTIONS
	  !
      !TIMES FACT 2 FOR THE SHEAR COMPONENTS, TO BE CONSISTANT WITH THE STRAIN DEFINITION
      !
      DSTRAIN(4) = 2.D0*DSTRAIN(4)
      DSTRAIN(5) = 2.D0*DSTRAIN(5)
      DSTRAIN(6) = 2.D0*DSTRAIN(6)
      ! 注意：輸出參數 ROT 在此子常式中並未被賦值。
      ! 如果不需要輸出 ROT，應從參數列表中移除。
      ROT = 0.0D0
	  RETURN
	  END SUBROUTINE D_HUGHES_WINGET

      SUBROUTINE ROTATE_TENSOR(TRANSFORMATION_MATRIX_6X6, TENSOR_VOIGT_INOUT)
      !$ACC ROUTINE SEQ
      ! IMPLICIT NONE ! 已在模組層級定義
      ! FUNCTION OF THIS SUBROUTINE:
      ! ROTATE A (VOIGT NOTATION) TENSOR USING A GIVEN 6X6 TRANSFORMATION MATRIX
      !
      DOUBLE PRECISION, INTENT(IN)    :: TRANSFORMATION_MATRIX_6X6(6,6)
      DOUBLE PRECISION, INTENT(INOUT) :: TENSOR_VOIGT_INOUT(6)

      !LOCAL VARIABLES
      DOUBLE PRECISION :: TEMP_TENSOR_VOIGT(6)

      ! S_rotated_voigt = T_transformation * S_voigt_old
      TEMP_TENSOR_VOIGT = MATMUL(TRANSFORMATION_MATRIX_6X6, TENSOR_VOIGT_INOUT)
      TENSOR_VOIGT_INOUT = TEMP_TENSOR_VOIGT

      RETURN
      END SUBROUTINE ROTATE_TENSOR
END MODULE RK_PROCEDURES_MOD