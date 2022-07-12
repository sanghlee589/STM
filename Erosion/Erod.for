      SUBROUTINE EROD(XB,XE,A,B,C,ATC,BTC,CTC,EATA,TAUC,THETA,PHI,ILAST,
     1    DL,LDLAST,XDBEG,NDEP,XINPUT,KTRATO,LOAD,TC,DETACH,QOSTAR)
C
C     SUBROUTINE EROD CALCULATES SOIL PARTICLE DETACHMENT IN EACH 
C     SLOPE SEGMENT ON FLOW PLANE 
C
C     MODULE ADAPTED FROM WEPP VERSION 2004.7 AND CALLED FROM 
C     SUBROUTINE ROUTE
C
C     AUTHOR(S): D.C. FLANAGAN AND J.C. ASCOUGH II
C     DATE LAST MODIFIED: 4-1-2005
C
C     + + + ARGUMENT DECLARATIONS + + +
C
      REAL XB, XE, A, B, C, EATA, TAUC, THETA, DL, XDBEG, LDLAST
      REAL ATC, BTC, CTC, PHI, QOSTAR
      REAL XINPUT(101), KTRATO, LOAD(101), TC(101), DETACH(101)
      INTEGER ILAST, NDEP
C
C     + + + ARGUMENT DEFINITIONS + + +
C
C     XB     - DISTANCE AT BEGINNING OF DETACHMENT REGION
C     XE     - DISTANCE AT END OF DETATCHMENT REGION
C     A      - SHEAR STRESS EQUATION COEFFICIENT
C     B      - SHEAR STRESS EQUATION COEFFICIENT
C     C      - SHEAR STRESS EQUATION COEFFICIENT
C     ATC    - TRANSPORT EQ. COEF
C     BTC    - TRANSPORT EQ. COEF
C     CTC    - TRANSPORT EQ. COEF
C     EATA   - N.D. RILL ERODIBILITY PARAMETER
C     TAUC   - N.D. CRITICAL SHEAR STRESS PARAMETER
C     THETA  - N.D. INTERRILL ERODIBILITY PARAMETER
C     PHI    - N.D. DEPOSITION PARAMETER
C     ILAST  - COUNTER VALUE FOR LAST POINT WHERE LOAD COMPUTED
C     DL     - N.D. DEPOSITION RATE
C     LDLAST - N.D. SEDIMENT LOAD CALCULATED AT I=ILAST
C     XDBEG  - N.D. DISTANCE WHERE DEPOSITION BEGINS
C     NDEP   - FLAG INDICATING IF DEPOSITION BEGINS IN
C              DETACHMENT SEGMENT
C     XINPUT - DISTANCE ON OFE WHERE DETACHMENT IS TO BE CALCULATED
C     KTRATO - N.D. SEDIMENT TRANSPORT EQUTION COEFFICIENT
C     LOAD   - N.D. SEDIMENT LOAD AT A POINT
C     TC     - SEDIMENT TRANSPORT CAPACITY AT A POINT (KG/S/M)
C     DETACH - DETACHMENT AT A POINT ON OFE
C     QOSTAR - N.D. FLOW DISCHARGE ONTO TOP OF AN OFE.
C
C     + + + LOCAL VARIABLES + + +
C
      REAL LDTRY, LDRAT, DX, XLAST, TCAP, XTRY, XFRT, DETFRT, DETTRY, 
     1    DETLST, LDRAT2, TCLAST, XTERM, XTRMTC, XX
      INTEGER IBEG, KFLAG, I, LOOPFG, CURRPT
      DATA XX /0.0/
C
C     + + + LOCAL VARIABLE DEFINITIONS + + +
C
C     LDTRY  - N.D. LOAD RETURNED FROM RUNGE-KUTTA PROCEDURE - IT
C              IS USED AS LONG AS LOAD IS LESS THAN TCAP
C     IBEG   - BEGINNING POINT COUNTER
C     DX     - DELTA X VALUE SENT TO RUNGE-KUTTA PROCEDURE
C     XLAST  - N.D. DISTANCE AT LAST POINT WHERE LOAD WAS COMPUTED
C     DCAP   - DETACHMENT CAPACITY AT A POINT
C     TCAP   - TRANSPORT CAPACITY AT A POINT
C     XTRY   - ITERATIVE SOLUTION FOR DISTANCE TO FIND POINT
C              WHERE DEPOSITION BEGINS IN A DETACHMENT SEGMENT
C     XFRT   - N.D. DISTANCE AT FRONT POINT WHERE DEPOSITION IS
C              PREDICTED TO OCCUR
C     DETFRT - RATIO AT FRONT POINT USED TO DETERMINE IF DETACH.
C              OR DEPOSITION CONDITIONS OCCUR AT X=XFRT
C     DETTRY - RATIO RELATING SEDIMENT LOAD AND TRANSPORT CAP.
C              AT X=XTRY
C     DETLST - RATIO RELATING SEDIMENT LOAD AND TRANSPORT CAP.
C              AT LAST POINT (X=XLAST)
C     LDRAT  - RATIO RELATING SED. LOAD AND TRANS. CAP. AT POINT
C              USED WHEN TRAN. CAP. > 0, AND SED. LOAD = 0
C     LDRAT2 - RATIO RELATING SED. LOAD AND TRANS. CAP. AT POINT
C              USED WHEN SEDIMENT LOAD IS GREATER THAN 0
C     TCLAST - TRANSPORT CAPACITY AT X=XLAST
C     KFLAG  - FLAG INDICATING TYPE OF HYDROLOGIC PLANE AND
C              WHICH RATIO TYPE TO USE
C     LOOPFG - FLAG. 1 = EXIT THE L2 & N2 LOOPS.
C     CURRPT - NUMBER OF THE CURRENT POINT WHEN L2 LOOP EXITED.
C
C     + + + SAVES + + +
C
      SAVE XFRT, DETLST, DETFRT
      SAVE XX, EATAX, TAUCX, SHR, DCAP
C
C     + + + SUBROUTINES CALLED + + +
C
C     RUNGE
C
C     + + + FUNCTION DECLARATIONS + + +
C
      REAL CROSS
C
C     BEGIN SUBROUTINE EROD
C
      LDRAT = 0.0
      LDRAT2 = 0.0
      NDEP = 0
      IBEG = ILAST + 1
C     
C     VERIFY THAT THE BEGINNING POINT HAS NOT EXCEEDED THE END OF
C     THE SLOPE SEGMENT
C     
      IF (IBEG.LT.102) THEN
         IF (XINPUT(IBEG).LE.XE) THEN
C           
C           LOOPFG IS SET TO 1 TO FORCE AN EXIT FROM THE LOOP
C           
            LOOPFG = 0
            I = ILAST
C           
   10       CONTINUE
C           
            I = I + 1
C           
            IF (XINPUT(I).LE.XE) THEN
C              
C              DETERMINE THE NONDIMENSIONAL HORIZONTAL DISTANCE BETWEEN
C              THE BEGINNING OF THE CURRENT SLOPE SEGMENT AND NEXT OF 
C              THE 101 POINTS OF THE PLANE
C              
C              ALSO DETERMINE TRANSPORT CAPACITY AT THE BEGINNING OF 
C              THE SLOPE SEGMENT [EQUATION 11.4.6]
C              
               IF (I.LE.IBEG) THEN
                  DX = XINPUT(I) - XB
                  XLAST = XB
                  XTERM = A * XB ** 2 + B * XB + C
                  XTRMTC = ATC * XB ** 2 + BTC * XB + CTC
                  TCLAST = XTRMTC * KTRATO
                  IF (TCLAST.LT.0.0) TCLAST = 0.0
C              
C              FOR ALL OTHER POINTS IN THE DETACHMENT SEGMENT, SET 
C              THE NONDIMENSIONAL HORIZONTAL DISTANCE INCREMENT TO 
C              THE DEFAULT OF 0.01 
C              
C              SET DISTANCE AND OBTAIN LOAD AND TRANSPORT CAPACITY 
C              FOR THE PREVIOUS POINT
               ELSE
                  DX = 0.01
                  XLAST = XINPUT(I-1)
                  LDLAST = LOAD(I-1)
                  TCLAST = TC(I-1)
               END IF
C              
C              CALCULATE DIMENSIONLESS SHEAR STRESS AT THE CURRENT
C              POINT  [EQUATION 11.4.1]
C              
               XTERM = A * XINPUT(I) ** 2 + B * XINPUT(I) + C
               XTRMTC = ATC * XINPUT(I) ** 2 + BTC * XINPUT(I) + CTC
C              
               IF (XTERM.NE.XX) THEN
                  IF (XTERM.GT.0.0) THEN
                     SHR = EXP(0.666667*LOG(XTERM))
                  ELSE
                     SHR = 0.0
                  END IF
C                 
                  XX = XTERM
C                 
C                 CALCULATE DETACHMENT CAPACITY AT THE CURRENT POINT
C                 [EQUATIONS 11.2.3, 11.3.7, 11.3.8, AND 11.4.1]
C                 
                  DCAP = EATA * (SHR-TAUC)
                  IF (DCAP.LT.0.0) DCAP = 0.0
                  EATAX = EATA
                  TAUCX = TAUC
               ELSE IF (EATAX.NE.EATA.OR.TAUCX.NE.TAUC) THEN
                  DCAP = EATA * (SHR-TAUC)
                  IF (DCAP.LT.0.0) DCAP = 0.0
                  EATAX = EATA
                  TAUCX = TAUC
               END IF
C              
C              CALCULATE DIMENSIONLESS TRANSPORT CAPACITY AT THE CURRENT
C              POINT [EQUATION 11.4.6]
C              
               TCAP = XTRMTC * KTRATO
               IF (TCAP.LT.0.0) TCAP = 0.0
               TC(I) = TCAP
C              
C              CHECK WHETHER ON A CASE 4 PLANE PAST WHERE RUNOFF ENDS
C              IF PAST POINT, SET LOAD EQUAL TO ZERO, AND FLAG THE 
C              PLANE (KFLAG = 4)
C              
               IF (QOSTAR.GT.-1.0.AND.QOSTAR.LT.0.0.AND.
     1             XINPUT(I).GT.-QOSTAR) THEN
                  LOAD(I) = 0.0
                  KFLAG = 4
                  NDEP = 0
C              
C              USE RUNGE-KUTTA NUMERICAL PROCEDURE TO SOLVE FOR 
C              SEDIMENT LOAD AT THE CURRENT POINT 
C              [NSERL REPORT 10; EQUATION 11.3.13]
C              
               ELSE
C                 
                  CALL RUNGE(A,B,C,ATC,BTC,CTC,EATA,TAUC,THETA,DX,XLAST,
     1                LDLAST,LOAD(I),XX,EATAX,TAUCX,SHR,DCAP,KTRATO)
C                 
C                 IF TRANSPORT CAPACITY AT THE CURRENT POINT IS GREATER
C                 THAN ZERO, CALCULATE RATIO VALUES USED TO TEST IF
C                 CURRENT POINT IS IN DEPOSITION
C                 
C                 KFLG = 1  INDICATES  TC > ZERO
C                 KFLG = 2  INDICATES  LOAD > ZERO
C                 
                  IF (TCAP.GT.0.0) THEN
                     LDRAT = 1.0 - (LOAD(I)/TCAP)
                     KFLAG = 1
                     DETACH(I) = DCAP * LDRAT
C                    
                     IF (LOAD(I).GT.0.0) THEN
                        LDRAT2 = (TCAP/LOAD(I)) - 1.0
                        KFLAG = 2
                     END IF
C                 
C                 WHEN TRANSPORT CAPACITY AT THE CURRENT POINT <= 
C                 ZERO - IF LOAD EXCEEDS ZERO MUST USE LDRAT2 RATIO
C                 
                  ELSE
                     IF (LOAD(I).GT.0.0) THEN
                        LDRAT2 = (TCAP/LOAD(I)) - 1.0
                        KFLAG = 2
C                    
C                    WHEN BOTH LOAD AND TRANSPORT CAPACITY AT POINT 
C                    ARE ZERO (KFLG = 3)
                     ELSE
                        LOAD(I) = 0.0
                        KFLAG = 3
                     END IF
                  END IF        ! IF (TCAP.GT.0.0) THEN
               END IF   ! IF (QOSTAR.GT.-1.0.AND.QOSTAR.LT.0.0.AND...)
C              
C              IF DEPOSITION IS PREDICTED AT CURRENT POINT, SET FLAG
C              NDEP = 1 AND EXIT DETACHMENT CALCULATIONS
C              
               IF ((KFLAG.EQ.2.AND.LDRAT2.LT.0.0).OR.(KFLAG.EQ.1.AND.
     1             LDRAT.LT.0.0)) THEN
                  NDEP = 1
                  LOOPFG = 1
               ELSE
                  ILAST = I
               END IF
C              
               IF (XINPUT(I).GE.1.0) LOOPFG = 1
C           
            ELSE
               LOOPFG = 1
            END IF     ! IF (XINPUT(I).LE.XE) THEN
C           
C           IF END OF SEGMENT HAS NOT BEEN REACHED AND HAVE NOT
C           ENCOUNTERED DEPOSITION THEN GO BACK THROUGH THE 
C           "IF (XINPUT(IBEG.LE.XE) THEN" LOOP
C           
            IF (LOOPFG.EQ.0.AND.I.LT.102) GO TO 10
C           
C           REMEMBER NUMBER OF CURRENT POINT
C           
            CURRPT = I
C           
C           ON THE LAST OF THE 101 OFE POINTS IN THIS SEGMENT, IF
C           DEPOSITION IS NOT OCCURRING, COMPUTE LOAD AT THE END OF
C           THE SEGMENT
C           
            IF (NDEP.EQ.0) THEN
C              
C              ON A SEGMENT WHERE FLOW IS PRESENT (NOT A CASE 4),
C              USE RUNGE-KUTTA SOLUTION [PAGE 11.6, SECTION 11.3.5]
C              
               IF (KFLAG.NE.4) THEN
                  IF (XE.NE.XINPUT(ILAST)) THEN
                     DX = XE - XINPUT(ILAST)
                     CALL RUNGE(A,B,C,ATC,BTC,CTC,EATA,TAUC,THETA,DX,
     1                   XINPUT(ILAST),LOAD(ILAST),LDLAST,XX,
     1                   EATAX,TAUCX,SHR,DCAP,KTRATO)
                     XLAST = XE
                  ELSE
                     LDLAST = LOAD(ILAST)
                     XLAST = XINPUT(ILAST)
                  END IF
C              
C              ON A SEGMENT WHERE FLOW IS NOT PRESENT (CASE 4),
C              PAST WHERE RUNOFF ENDS, SET LOAD TO ZERO AND RETURN 
C              TO ROUTE
               ELSE
                  LDLAST = 0.0
                  XLAST = XE
                  DL = 0.0
                  RETURN
               END IF
C              
               XTERM = A * XLAST ** 2 + B * XLAST + C
               XTRMTC = ATC * XLAST ** 2 + BTC * XLAST + CTC
C              
               IF (XTERM.NE.XX) THEN
C                 
C                 CALCULATE SHEAR STRESS AT END OF SEGMENT (X=XE)
C                 [EQUATION 11.4.1]
C                 
                  IF (XTERM.GT.0.0) THEN
                     SHR = EXP(0.666667*LOG(XTERM))
                  ELSE
                     SHR = 0.0
                  END IF
C                 
                  XX = XTERM
C                 
C                 CALCULATE DETACHMENT CAPACITY AT END OF SEGMENT 
C                 (X=XE) [EQUATION 11.2.3, AND OTHERS]
C                 
                  DCAP = EATA * (SHR-TAUC)
                  IF (DCAP.LT.0.0) DCAP = 0.0
                  EATAX = EATA
                  TAUCX = TAUC
               ELSE IF (EATAX.NE.EATA.OR.TAUCX.NE.TAUC) THEN
                  DCAP = EATA * (SHR-TAUC)
                  IF (DCAP.LT.0.0) DCAP = 0.0
                  EATAX = EATA
                  TAUCX = TAUC
               END IF
C              
C              CALCULATE TRANSPORT CAPACITY AT END OF SEGMENT (X=XE)
C              [EQUATION 11.4.6]
C              
               TCAP = XTRMTC * KTRATO
               IF (TCAP.LT.0.0) TCAP = 0.0
C              
C              TRANSPORT CAPACITY GREATER THAN ZERO
C              
               IF (TCAP.GT.0.0) THEN
                  LDRAT = 1.0 - (LDLAST/TCAP)
                  DL = DCAP * LDRAT
                  KFLAG = 1
C                 
C                 IF LOAD IS LESS THAN TRANSPORT CAPACITY AT END OF 
C                 SEGMENT (STILL IN DETACHMENT CONDITION) THEN RETURN 
C                 TO ROUTE
C                 
                  IF (LDRAT.GE.0.0) RETURN
C              
C              TRANSPORT CAPACITY IS ZERO
C              
               ELSE
C                 
C                 IF LOAD AT END OF SEGMENT IS ALSO ZERO THEN RETURN 
C                 TO ROUTE
C                 
                  IF (LDLAST.LE.0.0) THEN
                     LDLAST = 0.0
                     DL = 0.0
                     RETURN
                  END IF
C              
               END IF
C              
C              SET UP THE LAST POINT (XLAST) AND FRONT POINT (XFRT) 
C              RATIOS USED TO DETERMINE WHERE DEPOSITION BEGINS 
C              (X=XDBEG)
C              
               LDRAT2 = (TCAP/LDLAST) - 1.0
               KFLAG = 2
               DETFRT = LDRAT2
               IF (LOAD(ILAST).GT.0.0) DETLST = (TC(ILAST)/LOAD(ILAST))
     1             - 1.0
               NDEP = 1
               XFRT = XLAST
C              
               IF (XINPUT(ILAST).EQ.XFRT) THEN
                  XLAST = XINPUT(ILAST-1)
                  IF (DETFRT.EQ.LDRAT2) THEN
                     IF (LOAD(ILAST-1).GT.0.0) DETLST = (TC(ILAST-1)/
     1                   LOAD(ILAST-1)) - 1.0
                  ELSE
                     IF (TC(ILAST-1).GT.0.0) DETLST = 1.0 - (
     1                   LOAD(ILAST-1)/TC(ILAST-1))
                  END IF
               ELSE
                  XLAST = XINPUT(ILAST)
               END IF
C           
            ELSE
C              
C              ON THE LAST OF THE 101 POINTS ON THIS SEGMENT, 
C              DEPOSITION IS OCCURRING - COMPUTE WHERE DEPOSITION 
C              BEGINS (X=XDBEG)
C              
               XFRT = XINPUT(CURRPT)
C              
               IF (XLAST.LE.0.0.AND.TCLAST.LE.0.0.AND.LDLAST.LE.0.0) 
     1             THEN
                  KFLAG = 5
                  DETLST = DL
                  DETFRT = (PHI/(PHI+1.0)) * (KTRATO*(ATC*XFRT*XFRT+BTC*
     1                XFRT+CTC)-THETA)
               END IF
C              
               IF (KFLAG.EQ.1) THEN
                  DETFRT = LDRAT
                  IF (TCLAST.GT.0.0) THEN
                     DETLST = 1.0 - (LDLAST/TCLAST)
                  ELSE
                     DETLST = 0.0
                  END IF
               ELSE IF (KFLAG.EQ.2) THEN
                  DETFRT = LDRAT2
                  IF (LDLAST.GT.0.0) THEN
                     DETLST = (TCLAST/LDLAST) - 1.0
                  ELSE
                     DETLST = 0.0
                  END IF
               END IF
C              
C              IF AT TOP OF FLOW PLANE AND RATIOS ARE NOT POSITIVE AT 
C              BOTH THE BEGINNING AND END OF THE SEGMENT, ASSUME
C              DEPOSITION BEGINS AT TOP OF THE PLANE
C              
               IF (DETFRT.LE.0.0.AND.DETLST.LE.0.0.AND.XLAST.LE.0.0) 
     1             THEN
                  XDBEG = 0.0
                  RETURN
               END IF
C              
C              PREVENT A NEGATIVE VALUE FOR DETACHMENT AT THE LAST POINT
C              
               IF (DETLST.LT.0.0) DETLST = 0.0
C           
            END IF      ! IF (NDEP.EQ.0) THEN
C           
C           ITERATIVE PROCEEDURE TO FIND POINT WHERE DEPOSITION BEGINS
C           (X=XDBEG)
C           
            I = 0
C           
   20       I = I + 1
C           
C           USE CROSS FUNCTION TO SOLVE FOR POINT (X=XTRY) WHERE
C           TCAP = LDTRY, OR WHERE TEST RATIOS EQUAL ZERO
C           
            XTRY = CROSS(XLAST,DETLST,XFRT,DETFRT)
            DX = XTRY - XLAST
C           
C           USE RUNGE-KUTTA PROCEDURE TO ESTIMATE LOAD AT X=XTRY
C           
            CALL RUNGE(A,B,C,ATC,BTC,CTC,EATA,TAUC,THETA,DX,XLAST,
     1          LDLAST,LDTRY,XX,EATAX,TAUCX,SHR,DCAP,KTRATO)
C           
            TCAP = (ATC*XTRY**2+BTC*XTRY+CTC) * KTRATO
            IF (TCAP.LT.0.0) TCAP = 0.0
C           
            LOOPFG = 0
C           
            IF (KFLAG.EQ.2) THEN
               IF (LDTRY.LE.0.0) LDTRY = 0.00001
               IF (ABS((TCAP-LDTRY)/LDTRY).LT.0.001) THEN
                  LOOPFG = 1
               ELSE
                  DETTRY = (TCAP/LDTRY) - 1.0
               END IF
            ELSE IF (KFLAG.EQ.1) THEN
               IF (TCAP.LE.0.0) TCAP = 0.00001
               IF (ABS((LDTRY-TCAP)/TCAP).LT.0.001) THEN
                  LOOPFG = 1
               ELSE
                  DETTRY = 1.0 - (LDTRY/TCAP)
               END IF
            ELSE IF (KFLAG.EQ.5) THEN
               DETTRY = (PHI/(PHI+1.0)) * (TCAP-THETA)
            END IF
C           
            IF (LOOPFG.EQ.0) THEN
               IF (DETTRY.LE.0.0) THEN
                  DETFRT = DETTRY
                  XFRT = XTRY
               ELSE
                  XLAST = XTRY
                  DETLST = DETTRY
                  LDLAST = LDTRY
               END IF
            END IF
C           
            IF (I.LT.10.AND.LOOPFG.EQ.0) GO TO 20
C           
C           
            XDBEG = XTRY
            DL = 0.0
            LDLAST = LDTRY
C        
         END IF ! IF (XINPUT(IBEG).LE.XE) THEN
      END IF    ! IF (IBEG.LT.102) THEN
C     
      RETURN
      END
