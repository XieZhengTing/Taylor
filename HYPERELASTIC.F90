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


    SUBROUTINE HYPERELASTIC(LPROP,LSTRESS,FTEN,LSTRAIN)

    !$ACC ROUTINE SEQ
    ! USE THIS SUBROUTINE TO GET THE 2ND PK STRESS BY DEFORMATION GRADIENT FOR COMPRESSIBLE NEO-HOOKEAN MATERAIL
    USE FINT_FUNCTIONS

    DOUBLE PRECISION, INTENT (IN):: LPROP(30)
    DOUBLE PRECISION:: LSTRESS(6),LSTRESSTEN(3,3),LSTRAIN(6)
    DOUBLE PRECISION:: FTEN(3,3),C(3,3),FT(3,3),CINV(3,3),DET
!F2PY INTENT(IN)::LPROP(30),FTEN(3,3)
!F2PY INTENT(OUT)::LSTRESS(6),LSTRAIN(6)
    DOUBLE PRECISION:: I1,I2,I3,I3T1,I3T2,JT
    DOUBLE PRECISION:: Q1,PTEMP
    INTEGER:: I,J
    DOUBLE PRECISION:: IDENT(3,3)
    DATA IDENT/ 1.0d0, 0.0d0, 0.0d0, &
        0.0d0, 1.0d0, 0.0d0, &
        0.0d0, 0.0d0, 1.0d0/


    FT=TRANSPOSE(FTEN)

    C=MATMUL(FT,FTEN)   !RIGHT CAUCHY GREEN TENSOR
    
    LSTRAIN=TENSOR_2_VTENSOR(C)

    CALL DETERMINANT(FTEN,DET)

    I1= C(1,1)+C(2,2)+C(3,3)

    CALL DETERMINANT(C,I3)
    
    JT=SQRT(I3)

    CALL INVERSE(C,3,CINV)

    I2=C(1,1)*C(2,2)+C(2,2)*C(3,3)+C(1,1)*C(3,3)-C(1,2)*C(2,1)-C(2,3)*C(3,2)-C(1,3)*C(3,1)

    SELECT CASE (INT(LPROP(1)))

!    CASE(1) !NEO_HOOKEAN(LPROP,LSTRESS,FMAT)
!
!        Q1=LPROP(4)
!        I3T1=Q1*I3**(-1.D0/3.D0)
!        PTEMP=LPROP(2)*(JT-1.0D0)*JT
!
!        !CASE(2)
!
!        !CALL MOONEY_RIVILING(LPROP,LSTRESS,DET,FMAT)

    CASE(3) !J.S. CHEN 96 SECOND EXAMPLE

        Q1=LPROP(4)+2*LPROP(5)*(I1-3)+3*LPROP(6)*(I1-3)**2
        I3T1=Q1*I3**(-1.D0/3.D0)
        PTEMP=LPROP(2)*(JT-1.0D0)*JT

    CASE DEFAULT

        CALL EXIT_PROGRAM('NOT A HYPERELASTIC MATERIAL',1)

    END SELECT

    DO I=1,3
        DO J=1,3
            LSTRESSTEN(I,J)=2.0D0*(I3T1*(IDENT(I,J)-(1.0D0/3.0D0)*I1*(CINV(I,J))))+PTEMP*CINV(I,J)
        END DO
    END DO

    CALL SPK2CAUCHY(LSTRESSTEN,DET,FTEN,FT)

    LSTRESS=TENSOR_2_VTENSOR(LSTRESSTEN)



    END SUBROUTINE

    !*********************************************************

    SUBROUTINE SPK2CAUCHY(LSTRESSTEN,DET,FTEN,FT)
    !$ACC ROUTINE SEQ

    DOUBLE PRECISION:: LSTRESSTEN(3,3),LSTEMP(3,3)
    DOUBLE PRECISION:: FTEN(3,3),FT(3,3),DET

    LSTEMP=MATMUL(FTEN,LSTRESSTEN)

    LSTRESSTEN=(1.0D0/DET)*MATMUL(LSTEMP,FT)

    END SUBROUTINE