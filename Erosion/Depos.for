      SUBROUTINE DEPOS(XB,XE,CDEP,A,B,C,PHI,THETA,ILAST,DL,LDLAST,
     1    XINPUT,KTRATO,DETACH,LOAD,TC,QOSTAR)
C
C     CALCULATES DEPOSITION IN EACH SEGMENT OF THE HILLSLOPE
C
C     MODULE ADAPTED FROM WEPP VERSION 2004.7 AND CALLED FROM 
C     SUBROUTINE ROUTE
C
C     AUTHOR(S): D.C. FLANAGAN AND J.C. ASCOUGH II
C     DATE LAST MODIFIED: 4-1-2005
C
C     + + + ARGUMENT DECLARATIONS + + +
C
      REAL XB, XE, CDEP, A, B, C, PHI, THETA, DL, LDLAST, XTERM
      REAL XINPUT(101), KTRATO, DETACH(101), LOAD(101), TC(101)
      REAL QOSTAR
      INTEGER ILAST
C
C     + + + ARGUMENT DEFINITIONS + + +
C
C     XB     - N.D. DISTANCE WHERE DEPOSITION BEGINS
C     XE     - N.D. DISTANCE WHERE DEPOSITION ENDS
C     CDEP   - PORTION OF SOLUTION TO DEPOSITION EQUATION
C     A      - SHEAR STRESS EQUATION COEFFICIENT
C     B      - SHEAR STRESS EQUATION COEFFICIENT
C     C      - SHEAR STRESS EQUATION COEFFICIENT
C     PHI    - N.D. DEPOSITION PARAMETER
C     THETA  - N.D. INTERILL DETACHMENT PARAMETER
C     ILAST  - INDEX COUNTER AT LAST POINT WHERE LOAD COMPUTED
C     DL     - N.D. DEPOSITION RATE AT DISTANCE X=XE
C     LDLAST - N.D. SEDIMENT LOAD AT DISTANCE X=XE
C     XTERM  -
C     XINPUT -
C     KTRATO - N.D. SEDIMENT TRANSPORT EQUATION COEFFICIENT
C     DETACH - N.D. DETACHMENT FOR EACH POINT DOWN OFE
C     LOAD   - N.D. SEDIMENT LOAD AT A POINT
C     TC     - SEDIMENT TRANSPORT CAPACITY AT EACH POINT (KG/S/M)
C     QOSTAR - NONDIMENSIONAL FLOW DISCHARGE ONTO OFE
C
C     + + + LOCAL VARIABLES + + +
C
      INTEGER IBEG, I, LOOPFG
      REAL TCLAST
C
C     + + + LOCAL VARIABLE DEFINITIONS + + +
C
C     IBEG   - COUNTER VARIABLE VALUE AT FIRST DEPOSITION POINT
C     TCLAST - N.D. TRANSPORT CAPACITY AT DISTANCE X=XE
C     LOOPFG - FLAG.  1 = EXIT L3 LOOP.
C
C     + + + SUBROUTINES CALLED + + +
C
C     DEPEQS
C
C     BEGIN SUBROUTINE DEPOS
C
      IBEG = ILAST + 1
C     
      IF (IBEG.LT.102) THEN
C        
         IF (XINPUT(IBEG).GT.XE) THEN
C           
            IF (QOSTAR.LE.-1.0.OR.QOSTAR.GE.0.0.OR.XE.LE.-QOSTAR) THEN
               CALL DEPEQS(XB,CDEP,A,B,PHI,THETA,XE,DL,KTRATO,QOSTAR)
               XTERM = A * XE ** 2 + B * XE + C
               TCLAST = XTERM * KTRATO
               IF (TCLAST.LE.0.0) TCLAST = 0.0
               LDLAST = TCLAST - DL * (XE+QOSTAR) / PHI
            ELSE
               TCLAST = 0.0
               LDLAST = 0.0
            END IF
         ELSE
C           
            I = ILAST
            LOOPFG = 0
C           
   10       CONTINUE
C           
            I = I + 1
C           
            IF (XINPUT(I).LE.XE) THEN
C              
C              CHECK IF POINT IS PAST END OF RUNOFF ON A CASE 4 PLANE
C              
               IF (QOSTAR.LE.-1.0.OR.QOSTAR.GE.0.0.OR.XINPUT(I)
     1             .LE.-QOSTAR) THEN
C                 
                  CALL DEPEQS(XB,CDEP,A,B,PHI,THETA,XINPUT(I),
     1                DETACH(I),KTRATO,QOSTAR)
                  XTERM = A * XINPUT(I) ** 2 + B * 
     1                XINPUT(I) + C
                  TC(I) = XTERM * KTRATO
                  IF (TC(I).LT.0.0) TC(I) = 0.0
C                 
                  LOAD(I) = TC(I) - DETACH(I)*(XINPUT(I)+QOSTAR) / PHI
C                 
C                 ADDED TO PREVENT ERRONEOUS CALCULATION OF DETACHMENT 
C                 BY DEPOSITION EQUATION (FOR CASE 4 PLANE)
C                 
                  IF (THETA.LE.0.0.AND.I.GT.1.AND.LOAD(I).GT.LOAD(I-1))
     1                LOAD(I) = LOAD(I-1)
               ELSE
                  LOAD(I) = 0.0
                  TC(I) = 0.0
               END IF
C              
               IF (LOAD(I).LT.0.0) LOAD(I) = 0.0
               ILAST = I
               IF (XINPUT(I).GE.1.0) LOOPFG = 1
            ELSE
               LOOPFG = 1
            END IF
C           
            IF (I.LT.101.AND.LOOPFG.EQ.0) GO TO 10
C           
C           CORRECTIONS MADE TO PREVENT BOMBING FOR A CASE 4 
C           PLANE WHERE XE IS GREATER THAN -QOSTAR
C           
C           CASE 1 - PLANE WHERE FLOW DOES NOT END
C           
            IF (QOSTAR.GE.0.0.OR.QOSTAR.LE.-1.0) THEN
               CALL DEPEQS(XB,CDEP,A,B,PHI,THETA,XE,DL,KTRATO,QOSTAR)
               XTERM = A * XE ** 2 + B * XE + C
               TCLAST = XTERM * KTRATO
               IF (TCLAST.LT.0.0) TCLAST = 0.0
               LDLAST = TCLAST - DL * (XE+QOSTAR) / PHI
            ELSE
C              
C              CASE 2 - PLANE WHERE FLOW ENDS BUT ON THE CURRENT 
C              SLOPE SEGMENT
C              
               IF (XE.LT.-QOSTAR) THEN
                  CALL DEPEQS(XB,CDEP,A,B,PHI,THETA,XE,DL,KTRATO,QOSTAR)
                  XTERM = A * XE ** 2 + B * XE + C
                  TCLAST = XTERM * KTRATO
                  IF (TCLAST.LT.0.0) TCLAST = 0.0
                  LDLAST = TCLAST - DL * (XE+QOSTAR) / PHI
               ELSE
C                 
C                 CASE 3 - PLANE WHERE FLOW ENDS ON THE CURRENT SLOPE 
C                 SEGMENT
C                 
                  TCLAST = 0.0
                  LDLAST = 0.0
                  DL = 0.0
C              
               END IF
            END IF
C           
            IF (LDLAST.LT.0.0) LDLAST = 0.0
            IF (TCLAST.LT.0.0) TCLAST = 0.0
C        
         END IF ! IF (XINPUT(IBEG).GT.XE) THEN
      END IF    ! IF (IBEG.LT.102) THEN
C     
      RETURN
      END
