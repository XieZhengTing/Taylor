    SUBROUTINE TESTER(X,SHP,SHPD,LN,LSTACK,GCOO)
    
    DOUBLE PRECISION, INTENT(IN):: X(3)
    DOUBLE PRECISION, INTENT(IN):: SHP(*)
    DOUBLE PRECISION, INTENT(IN):: SHPD(3,*)
    DOUBLE PRECISION, INTENT(IN):: GCOO(3,*)
    INTEGER, INTENT(IN):: LN,LSTACK(*)
      
    WRITE(*,*) '-------------------------------'
    WRITE(*,*) '-------------------------------'
    WRITE(*,*) '-------errors in SHP:----------'
    WRITE(*,*) '-------------------------------'
    WRITE(*,*) '-------------------------------'
    
    
    WRITE(*,*) '-----------------------'
    WRITE(*,*) '-----0ST ORDER---------'
    WRITE(*,*) '-----------------------'
    
      TEST = SUM(SHP(1:LN))
      TEST = TEST - 1.0d0
    WRITE(*,'(E25.15)') TEST
    
    WRITE(*,*) '-----------------'
!    CONTINUE
!      TEST = SUM(SHPD(1,1:LN))
!    WRITE(*,'(E25.15)') TEST
!    
!    CONTINUE
!      TEST = SUM(SHPD(2,1:LN))
!    WRITE(*,'(E25.15)') TEST
!    
!    CONTINUE
!      TEST = SUM(SHPD(3,1:LN))
!    WRITE(*,'(E25.15)') TEST
!    WRITE(*,*) '-----------------'
!    
!    CONTINUE
      
    WRITE(*,*) '-----------------------'
    WRITE(*,*) '-----1ST ORDER---------'
    WRITE(*,*) '-----------------------'
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHP(I)*GCOO(1,II)
    END DO
    TEST = TEST - X(1)
    WRITE(*,'(E25.15)') TEST
    CONTINUE
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHP(I)*GCOO(2,II)
    END DO
    TEST = TEST - X(2)
    WRITE(*,'(E25.15)') TEST
    
    CONTINUE
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHP(I)*GCOO(3,II)
    END DO
    TEST = TEST - X(3)
    WRITE(*,'(E25.15)') TEST
    
    CONTINUE
    
    WRITE(*,*) '-----------------------'
    WRITE(*,*) '-----2ND ORDER---------'
    WRITE(*,*) '-----------------------'
    
    
    
    
    CONTINUE
      
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHP(I)*GCOO(1,II)**2
    END DO
    TEST = TEST - X(1)**2
    WRITE(*,'(E25.15)') TEST
    CONTINUE
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHP(I)*GCOO(2,II)**2
    END DO
    TEST = TEST - X(2)**2
    WRITE(*,'(E25.15)') TEST
    
    CONTINUE
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHP(I)*GCOO(3,II)**2
    END DO
    TEST = TEST - X(3)**2
    WRITE(*,'(E25.15)') TEST
    
    CONTINUE
    
    
      
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHP(I)*GCOO(1,II)*GCOO(2,II)
    END DO
    TEST = TEST - X(1)*X(2)
    WRITE(*,'(E25.15)') TEST
    CONTINUE
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHP(I)*GCOO(1,II)*GCOO(3,II)
    END DO
    TEST = TEST - X(1)*X(3)
    WRITE(*,'(E25.15)') TEST
    
    CONTINUE
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHP(I)*GCOO(2,II)*GCOO(2,II)
    END DO
    TEST = TEST - X(2)*X(2)
    WRITE(*,'(E25.15)') TEST
    
    CONTINUE
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHP(I)*GCOO(2,II)*GCOO(3,II)
    END DO
    TEST = TEST - X(2)*X(3)
    WRITE(*,'(E25.15)') TEST
    
    CONTINUE
    
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHP(I)*GCOO(3,II)*GCOO(3,II)
    END DO
    TEST = TEST - X(3)*X(3)
    WRITE(*,'(E25.15)') TEST
    
    CONTINUE
    
    RETURN
    
    END SUBROUTINE
    
    
    
      

    SUBROUTINE TESTERX(X,SHP,SHPD,LN,LSTACK,GCOO)
    
    DOUBLE PRECISION, INTENT(IN):: X(3)
    DOUBLE PRECISION, INTENT(IN):: SHP(*)
    DOUBLE PRECISION, INTENT(IN):: SHPD(3,*)
    DOUBLE PRECISION, INTENT(IN):: GCOO(3,*)
    INTEGER, INTENT(IN):: LN,LSTACK(*)
      
    WRITE(*,*) '-------------------------------'
    WRITE(*,*) '-------------------------------'
    WRITE(*,*) '-------errors in DSHP:----------'
    WRITE(*,*) '-------------------------------'
    WRITE(*,*) '-------------------------------'
    
    
    WRITE(*,*) '-----------------------'
    WRITE(*,*) '-----0ST ORDER---------'
    WRITE(*,*) '-----------------------'
    
    
    CONTINUE
      TEST = SUM(SHPD(1,1:LN))
    WRITE(*,'(E25.15)') TEST
    
    CONTINUE
      TEST = SUM(SHPD(2,1:LN))
    WRITE(*,'(E25.15)') TEST
    
    CONTINUE
      TEST = SUM(SHPD(3,1:LN))
    WRITE(*,'(E25.15)') TEST
    WRITE(*,*) '-----------------'
    
    CONTINUE
      
    WRITE(*,*) '-----------------------'
    WRITE(*,*) '-----1ST ORDER---------'
    WRITE(*,*) '-----------------------'
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHPD(1,I)*GCOO(1,II)
    END DO
    TEST = TEST - 1.0d0
    WRITE(*,'(E25.15)') TEST
    CONTINUE
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHPD(2,I)*GCOO(1,II)
    END DO
    TEST = TEST 
    WRITE(*,'(E25.15)') TEST
    CONTINUE
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHPD(3,I)*GCOO(1,II)
    END DO
    TEST = TEST 
    WRITE(*,'(E25.15)') TEST
    CONTINUE
    
    !----
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHPD(2,I)*GCOO(2,II)
    END DO
    TEST = TEST - 1.0d0
    WRITE(*,'(E25.15)') TEST
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHPD(2,I)*GCOO(1,II)
    END DO
    TEST = TEST 
    WRITE(*,'(E25.15)') TEST
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHPD(2,I)*GCOO(3,II)
    END DO
    TEST = TEST 
    WRITE(*,'(E25.15)') TEST
    
    CONTINUE
    
    !----
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHPD(3,I)*GCOO(1,II)
    END DO
    TEST = TEST 
    WRITE(*,'(E25.15)') TEST
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHPD(3,I)*GCOO(2,II)
    END DO
    TEST = TEST 
    WRITE(*,'(E25.15)') TEST
    
      TEST = 0.0d0
    DO I=1, LN
      II = LSTACK(I)
      TEST = TEST + SHPD(3,I)*GCOO(3,II)
    END DO
    TEST = TEST - 1.0d0
    WRITE(*,'(E25.15)') TEST
    
    
    !PAUSE
    
    
    
    CONTINUE
    
    RETURN
    
    END SUBROUTINE
    

    
     !TEST HUGHS-WINDET ROTATION ALGORITHM, LATER SHOULD BE REMOVED
    SUBROUTINE ROTATION_TEST(LOCAL_DSP,LOCAL_COO,LOCAL_NUMP,TIME_COUNTER,DLT)
    !   
    INTEGER, INTENT(IN):: LOCAL_NUMP    
    DOUBLE PRECISION, INTENT(OUT)::LOCAL_DSP(3*LOCAL_NUMP)
    DOUBLE PRECISION, INTENT(IN)::LOCAL_COO(3,LOCAL_NUMP) 
    DOUBLE PRECISION ::TIME_COUNTER,DLT  
    DOUBLE PRECISION :: W,X,Y,T
    
    W = 3.1415926*100.d0 ! 3.1415926*0.1
    T = TIME_COUNTER
    
    DO I=1, LOCAL_NUMP
        X = LOCAL_COO(1,I)
        Y = LOCAL_COO(2,I)
        LOCAL_DSP((I-1)*3+1) = DLT*W*(-SIN(W*T)*X - COS(W*T)*Y) 
        LOCAL_DSP((I-1)*3+2) = DLT*W*( COS(W*T)*X - SIN(W*T)*Y)   
        LOCAL_DSP((I-1)*3+3) = 0.D0
    
    END DO
    
    RETURN
    
    END SUBROUTINE