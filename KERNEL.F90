
         
       SUBROUTINE MLS_KERNEL0(XSA,AJ,ISPLINE,      & !INPUT
                             PHI,PHI_X,ISZERO)   !OUTPUT
                             
      IMPLICIT NONE
       
      !GLOBAL VARIABLES
      DOUBLE PRECISION:: XS
      INTEGER:: ISPLINE
      DOUBLE PRECISION:: AJ,PHI,PHI_X
      
      !LOCAL VARIABLES
      DOUBLE PRECISION::XSA
      INTEGER::I,J,K,II,JJ,KK
      
      INTEGER::CL,L,FL,F2LPO
      LOGICAL:: ISZERO
      
      !XSA = XSA/AJ
      ISZERO=.FALSE.
      
        IF (ISPLINE.GE.5) THEN
        
        IF (XSA.LE.(1.0D0)) THEN
            L=ISPLINE+1
            CALL FACTI(L,FL)
            CALL FACTI(2*L+1,F2LPO)
            !CL=FACTORIAL(2*L+1)/(2**(2*L+1)*FACTORIAL(L)**2)
            CL = F2LPO / (2**(2*L+1)*FL**2)

             PHI = CL*(1-XSA**2)**L
          
        ELSE
            PHI = 0.0D0
            ISZERO=.TRUE.
        END IF


        ELSEIF (ISPLINE.EQ.4) THEN

        IF (XSA.LE.(1.0D0/3.0D0)) THEN
            PHI = 11.0D0/20.0D0 - 9.0D0/2.0D0*XSA**2  +      81.0D0/4.0D0*XSA**4                -81.0D0/4.0D0*XSA**5
            PHI_X = - 9.0D0*XSA  +      81.0D0*XSA**3                -81.0D0/4.0D0*5.0D0*XSA**4
        ELSEIF (XSA.LE.(2.0D0/3.0D0)) THEN
            PHI = 17.0D0/40.0D0 + 15.0D0/8.0D0*XSA - 63.0D0/4.0D0*XSA**2       +135.0D0/4.0D0*XSA**3 -     243.0D0/8.0D0*XSA**4 +   81.0D0/8.0D0*XSA**5
            PHI_X = 15.0D0/8.0D0 - 63.0D0/2.0D0*XSA       +135.0D0/4.0D0*3.0D0*XSA**2 -     243.0D0/2.0D0*XSA**3 +   81.0D0/8.0D0*5.0D0*XSA**4
        ELSEIF (XSA.LE.(1.0D0)) THEN
            PHI = 81.0D0/40.0D0 - 81.0D0/8.0D0*XSA + 81.0D0/4.0D0*XSA**2  -     81.0D0/4.0D0*XSA**3    +     81.0D0/8.0D0*XSA**4 -   81.0D0/40.0D0*XSA**5
            PHI_X = - 81.0D0/8.0D0 + 81.0D0/2.0D0*XSA  -     81.0D0/4.0D0*3.0D0*XSA**2    +     81.0D0/2.0D0*XSA**3 -   81.0D0/8.0D0*XSA**4
        ELSE
            PHI = 0.0D0
            PHI_X = 0.0D0
            ISZERO=.TRUE.
        END IF

        ELSEIF (ISPLINE.EQ.3) THEN
            
            ! C3 CONTINUETY
            
            IF (XSA.LE.(0.5D0)) THEN
                PHI=  2.0D0/3.0D0    - 4.0D0*XSA**2  +  4.0D0*XSA**3
                PHI_X=  - 8.0D0*XSA  +  12.0D0*XSA**2
            ELSEIF (XSA.LE.(1.0D0)) THEN
                PHI=  4.0D0/3.0D0 -4.0D0*XSA +  4.0D0*XSA**2  - 4.0D0/3.0D0*XSA**3
                PHI_X= -4.0D0 +  8.0D0*XSA  - 4.0D0*XSA**2
            ELSE   
                PHI= 0.0D0
                PHI_X= 0.0D0
            ISZERO=.TRUE.
            END IF

        ELSEIF (ISPLINE.EQ.2) THEN
            
            ! C2 CONTINUETY
            
            IF (XSA.LE.(1.0D0)) THEN
                
                PHI =   1.0D0 -  6.0D0*XSA**2      +  8.0D0*XSA**3      - 3.0D0*XSA**4
                PHI_X =   -  12.0D0*XSA      +  8.0D0*3.0D0*XSA**2      - 3.0D0*4.0D0*XSA**3
            ELSE
                PHI= 0.0D0
                PHI_X= 0.0D0
                ISZERO=.TRUE.
        END IF
            
        ELSEIF (ISPLINE.EQ.1) THEN
            
            XSA=XSA*1.5
            IF (XSA.LE.(0.5D0)) THEN
                PHI =  3.0D0/4.0D0-XSA**2
                PHI_X =  3.0D0/2.0D0-XSA
            ELSEIF (XSA.LE.(1.5D0)) THEN
                PHI =   9.0D0/8.0D0-3.0D0/2.0D0*XSA+0.5*XSA**2
                PHI_X =   -3.0D0/2.0D0+0.5*2.0D0*XSA
            ELSE
                PHI= 0.0D0
                PHI_X= 0.0D0
                ISZERO=.TRUE.
        END IF
            
        ELSEIF (ISPLINE.EQ.0) THEN
            
            ! C0 CONTINUETY
            
            IF (XSA.LE.(1.0D0)) THEN
                
                  PHI=  1.0d0  -XSA
                  PHI_X=  -1.0d0
            ELSE
                  PHI=  0.0D0
                  PHI_X=  0.0D0
                  ISZERO=.TRUE.
        END IF

        ELSEIF (ISPLINE.EQ.-1) THEN
            
            ! C-1 CONTINUETY
            
            IF (XSA.LE.(1.0D0)) THEN
                
                PHI=1.0D0
                PHI_X=0.0D0
            ELSE
                PHI=0.0D0
                PHI_X=0.0D0
                ISZERO=.TRUE.
        END IF

        ELSE
            CALL EXIT_PROGRAM('KERNEL TYPE NOT SUPPORTED',1)
            PAUSE
        END IF
        
        !PHI = PHI / AJ
      !  PHI_X = PHI_X / AJ
      
      RETURN
      
    END SUBROUTINE
    
 
      

      
      
      SUBROUTINE FACTI(INTEG,FACT_INTEG)
      IMPLICIT NONE
      
      INTEGER::INTEG,FACT_INTEG
      INTEGER:: I
      
      FACT_INTEG = 1
      DO I=1, INTEG
        FACT_INTEG = FACT_INTEG * I
      END DO
      
      RETURN
      
      END SUBROUTINE
      
    