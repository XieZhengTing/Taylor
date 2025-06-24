
    SUBROUTINE GET_BINS(NP,XYZ_POS,NODES_IN_BIN,MAX_NEIGH,NODELIST_IN_BIN, &
        NBINS,NBINSX,NBINSY,NBINSZ,ISPACE,JSPACE,KSPACE,XMIN,YMIN,ZMIN,XBIN,YBIN,ZBIN)

    IMPLICIT NONE
    INTEGER:: NP
    DOUBLE PRECISION:: XYZ_POS(3,NP)
    INTEGER:: ISPACE(NP),JSPACE(NP),KSPACE(NP)
    INTEGER:: NBINSX,NBINSY,NBINSZ,NBINS,MAX_NEIGH
    INTEGER::NODES_IN_BIN(NBINS)
    INTEGER::NODELIST_IN_BIN(NBINS,MAX_NEIGH)
    DOUBLE PRECISION:: XMIN,YMIN,ZMIN,XBIN,YBIN,ZBIN

    !LOCAL
    INTEGER:: I,ISX,ISY,ISZ,M
    DOUBLE PRECISION:: X,Y,Z

    !$ACC PARALLEL LOOP PRESENT_OR_CREATE(NODES_IN_BIN)
    DO I = 1, NBINS
        NODES_IN_BIN(I) = 0
    END DO
    !$ACC END PARALLEL LOOP

    DO I=1,NP

        X = XYZ_POS(1,I)
        Y = XYZ_POS(2,I)
        Z = XYZ_POS(3,I)

        !SAVE THE POSITION IN I-J SPACE

        ISX = FLOOR((X-XMIN)/XBIN)+1
        ISX=MAX(1,ISX)
        ISX=MIN(ISX,NBINSX)
        ISPACE(I) = ISX

        ISY = FLOOR((Y-YMIN)/YBIN)+1
        ISY=MAX(1,ISY)
        ISY=MIN(ISY,NBINSY)
        JSPACE(I) = ISY

        ISZ = FLOOR((Z-ZMIN)/ZBIN)+1
        ISZ=MAX(1,ISZ)
        ISZ=MIN(ISZ,NBINSZ)
        KSPACE(I) = ISZ


        !SAVE THE BIN LISTS OF NODES
        M = NBINSY*NBINSX*(KSPACE(I)-1) + NBINSX*(JSPACE(I)-1) + ISPACE(I)

        NODES_IN_BIN(M) = NODES_IN_BIN(M)+1
        IF (NODES_IN_BIN(M).GE.MAX_NEIGH) THEN
            WRITE(*,*) 'fatal error: (NODES_IN_BIN(M).GE.MAX_NEIGH), stopping'
            STOP

        END IF

        NODELIST_IN_BIN(M,NODES_IN_BIN(M)) = I

    END DO

    IF(SUM(NODES_IN_BIN).NE.NP)  THEN
        write(*,*)'MISS NODES',NP,SUM(NODES_IN_BIN)
        STOP
    ENDIF
    RETURN

    END SUBROUTINE




    SUBROUTINE BIN_INODE(IP,INODE,NODE_NO, &
        NP,NODES_IN_BIN,MAX_NEIGH,NODELIST_IN_BIN, &
        NBINS,NBINSX,NBINSY,NBINSZ,ISPACE,JSPACE,KSPACE)

    IMPLICIT NONE
    INTEGER:: NP
    INTEGER:: ISPACE(NP),JSPACE(NP), KSPACE(NP)
    INTEGER:: NBINSX,NBINSY,NBINSZ,NBINS,MAX_NEIGH
    INTEGER::NODES_IN_BIN(NBINS)
    INTEGER::NODELIST_IN_BIN(NBINS,MAX_NEIGH)
    INTEGER:: IP
    INTEGER:: INODE(10000)
    INTEGER:: NODE_NO

    !LOCAL
    INTEGER:: I,J,K,M
    INTEGER:: ISEARCH_DOWN,JSEARCH_DOWN,ISEARCH_UP,JSEARCH_UP,IX,IY,IZ
    INTEGER:: KSEARCH_DOWN,KSEARCH_UP
    INTEGER:: ISEARCH,JSEARCH,KSEARCH

    ISEARCH = 1 !CEILING(XGWINMAX/XBIN)
    JSEARCH = 1 !CEILING(YGWINMAX/YBIN)
    KSEARCH = 1 !CEILING(ZGWINMAX/ZBIN)


    ISEARCH_DOWN = MAX(ISPACE(NODE_NO)-ISEARCH,1)
    JSEARCH_DOWN = MAX(JSPACE(NODE_NO)-JSEARCH,1)
    KSEARCH_DOWN = MAX(KSPACE(NODE_NO)-KSEARCH,1)


    ISEARCH_UP = MIN(ISPACE(NODE_NO)+ISEARCH,NBINSX)
    JSEARCH_UP = MIN(JSPACE(NODE_NO)+JSEARCH,NBINSY)
    KSEARCH_UP = MIN(KSPACE(NODE_NO)+KSEARCH,NBINSZ)

    IP = 0

    DO IX = ISEARCH_DOWN,ISEARCH_UP
        DO IY = JSEARCH_DOWN,JSEARCH_UP
            DO IZ = KSEARCH_DOWN,KSEARCH_UP

                !THIS BIN NO
                !               M =(IZ-1)*NBINSX*NBINSY + (IY-1)*NBINSY + IX
                M =(IZ-1)*NBINSX*NBINSY + (IY-1)*NBINSX + IX

                DO J=1,NODES_IN_BIN(M)

                    IP = IP + 1

                    IF (IP.GT.(10000)) THEN
                        WRITE(*,*) 'fatal error, IP.GT.10000, stopping'
                        stop
                    end if

                    INODE(IP) = NODELIST_IN_BIN(M,J)

                END DO
            END DO
        END DO
    END DO



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

    SUBROUTINE SOFT_SEARCH(GNUMP,GCOO,GN,GSTART,DIM_NN_LIST,GSTACK,GMAXN, GWIN, &
        NODES_IN_BIN,MAX_NEIGH,NODELIST_IN_BIN, &
        NBINS,NBINSX,NBINSY,NBINSZ,ISPACE,JSPACE,KSPACE, &
        GXDIST_MAX, GYDIST_MAX, GZDIST_MAX,GSM_LEN)
        
    USE CONTROL
    IMPLICIT NONE
    !
    ! FUNCTION OF THIS SUBROUTINE:
    !
    ! PERFORM O(NlogN) NODE SEARCH
    !
    INTEGER, INTENT(IN)::GNUMP
    INTEGER, INTENT(INOUT)::GN(GNUMP)
    DOUBLE PRECISION, INTENT(IN)::GCOO(3,GNUMP)
    DOUBLE PRECISION, INTENT(IN)::GWIN(3,GNUMP)
    INTEGER, INTENT(INOUT)::GSTART(GNUMP)
    INTEGER, INTENT(INOUT)::DIM_NN_LIST
    INTEGER, INTENT(INOUT)::GSTACK(DIM_NN_LIST)
    INTEGER, INTENT(INOUT):: GMAXN
    DOUBLE PRECISION:: GSM_LEN(6,GNUMP)     !SMOOTHING LENGTHS FOR EACH NODE


    INTEGER:: ISPACE(GNUMP),JSPACE(GNUMP), KSPACE(GNUMP)
    INTEGER:: NBINSX,NBINSY,NBINSZ,NBINS,MAX_NEIGH
    INTEGER::NODES_IN_BIN(NBINS)
    INTEGER::NODELIST_IN_BIN(NBINS,MAX_NEIGH)
    INTEGER:: IP
    INTEGER:: INODE_BIN(10000)
    INTEGER:: IP_BIN

    DOUBLE PRECISION, INTENT(IN)::GXDIST_MAX, GYDIST_MAX, GZDIST_MAX

    !LOCAL
    INTEGER:: I,J,K,M,P, JJ
    INTEGER:: DIM_NN_LIST_ACTUAL
    LOGICAL:: FIRST
    DOUBLE PRECISION:: XI, YI, ZI, XJ, YJ, ZJ, DIST_NORM
    INTEGER:: NODE_NO
    LOGICAL:: LADD
    
    DOUBLE PRECISION:: LSM_PTS(3,7),LSM_LEN(6)
    INTEGER:: MMAX
                
    DIM_NN_LIST_ACTUAL = 0

    !$ACC PARALLEL LOOP PRESENT(GN)
    DO I = 1, GNUMP
        GN(I) = 0
    END DO
    !$ACC END PARALLEL LOOP

    K = 0

    GMAXN = 0
    IF (ITYPE_INT.EQ.0) THEN
      MMAX=1
    ELSE
      MMAX = 7
    END IF
            
        DO I=1,GNUMP
            !INTEGERATION POINT I

            XI=GCOO(1,I)
            YI=GCOO(2,I)
            ZI=GCOO(3,I)

            FIRST = .TRUE.

            NODE_NO = I

            CALL BIN_INODE(IP_BIN,INODE_BIN,NODE_NO, &
                GNUMP,NODES_IN_BIN,MAX_NEIGH,NODELIST_IN_BIN, &
                NBINS,NBINSX,NBINSY,NBINSZ,ISPACE,JSPACE,KSPACE)
                
            LSM_LEN(:) = GSM_LEN(:,I)
        
            DO JJ=1,IP_BIN

                J = INODE_BIN(JJ)
                
                LADD=.FALSE.
                ! GET THE SMOOTHING INFORMATION
                ! (1-6) = (+X, -X, +Y, -Y, +Z, -Z)
                !
                DO P = 1, 3
                    DO  M = 1, 7
                        LSM_PTS(P,M) = GCOO(P,J)
                    END DO
                END DO
                LSM_PTS(1,2) = GCOO(1,J) + LSM_LEN(1)
                LSM_PTS(1,3) = GCOO(1,J) - LSM_LEN(2)
                LSM_PTS(2,4) = GCOO(2,J) + LSM_LEN(3)
                LSM_PTS(2,5) = GCOO(2,J) - LSM_LEN(4)
                LSM_PTS(3,6) = GCOO(3,J) + LSM_LEN(5)
                LSM_PTS(3,7) = GCOO(3,J) - LSM_LEN(6)
                !
                !CHECK THE SMOOTHING POINTS
                !
                DO M = 1, MMAX

                    XJ = LSM_PTS(1,M)
                    YJ = LSM_PTS(2,M)
                    ZJ = LSM_PTS(3,M)
                    
                    IF (SHSUP) THEN !ELLIPSOIDAL SUPPORT
                
                        !DIST_NORM = SQRT(((ABS(XJ-XI)-GXDIST_MAX)/GWIN(1,J))**2 &
                        !                +((ABS(YJ-YI)-GYDIST_MAX)/GWIN(2,J))**2 &
                        !                +((ABS(ZJ-ZI)-GZDIST_MAX)/GWIN(3,J))**2)
                        
                        DIST_NORM = SQRT(((ABS(XJ-XI))/GWIN(1,J))**2 &
                                        +((ABS(YJ-YI))/GWIN(2,J))**2 &
                                        +((ABS(ZJ-ZI))/GWIN(3,J))**2)
                                        
                        IF (DIST_NORM.LT.(1.0d0)) LADD = .TRUE.


                    ELSE !CUBOID SUPPORT
        
                        DIST_NORM = ABS(XJ-XI)/GWIN(1,J)
                        !DIST_NORM = (ABS(XJ-XI)-GXDIST_MAX)/GWIN(1,J)

                        IF (DIST_NORM.LT.(1.0d0)) THEN

                            DIST_NORM = ABS(YJ-YI)/GWIN(2,J)
                            !DIST_NORM = (ABS(YJ-YI)-GYDIST_MAX)/GWIN(2,J)

                            IF (DIST_NORM.LT.(1.0d0)) THEN

                                DIST_NORM = ABS(ZJ-ZI)/GWIN(3,J)
                                !DIST_NORM = (ABS(ZJ-ZI)-GZDIST_MAX)/GWIN(3,J)

                                IF (DIST_NORM.LT.(1.0d0)) THEN
                            
                                LADD = .TRUE.
                    
                                END IF !DISTZ
                            END IF !DISTY
                        END IF !DISTZ
                    
                        END IF !TYPE OF SUPPORT
                    
                    IF (LADD) GOTO 10
                        
                    END DO !LOOP OF SMOOTHING POINTS
                    !
                    ! ADD THE NODE IF NEEDED
                    !
10                  IF (LADD) THEN
                        !
                        GN(I) = GN(I) + 1
                        !
                        K = K + 1
                        !
                        GSTACK(K) = J
                        !
                        IF (FIRST) THEN
                            !
                            !SAVE LOCATION IN STACK
                            !
                            GSTART(I) = K
                            FIRST = .FALSE.
                        END IF !FIRST
                        !
                    END IF !LADD
                    
                END DO !J=1,GNUMP (NEIGHBOR NODES)

          IF (GN(I).GT.GMAXN) GMAXN = GN(I)

        END DO !I=1,GNUMP (INTEGRATION POINTS)

        DIM_NN_LIST_ACTUAL = K


        RETURN

    END SUBROUTINE



