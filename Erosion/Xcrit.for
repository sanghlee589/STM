      SUBROUTINE XCRIT(A,B,C,TAUC,XB,XE,XC1,XC2,MSHEAR)
C
C     SUBROUTINE XCRIT DETERMINES WHETHER SHEAR STRESS EXCEEDS 
C     CRITICAL SHEAR STRESS FOR A CERTAIN SEGMENT AND RETURNS A FLAG
C
C     MODULE ADAPTED FROM WEPP VERSION 2004.7 AND CALLED FROM 
C     SUBROUTINE ROUTE
C
C     AUTHOR(S): D.C. FLANAGAN AND J.C. ASCOUGH II
C     DATE LAST MODIFIED: 9-30-2004
C
C     + + + ARGUMENT DECLARATIONS + + +
C
      REAL A, B, C, TAUC, XB, XE, XC1, XC2
      INTEGER MSHEAR
C
C     + + + ARGUMENT DEFINITIONS + + +
C
C     A     - SHEAR STRESS EQUATION COEFFICIENT
C     B     - SHEAR STRESS EQUATION COEFFICIENT
C     C     - SHEAR STRESS EQUATION COEFFICIENT
C     TAUC  - N.D. CRITICAL SHEAR FOR THE OVERLAND FLOW ELEMENT
C     XB    - BEGINNING DISTANCE FOR SEGMENT OF OFE
C     XE    - ENDING DISTANCE FOR SEGMENT OF OFE
C     XC1   - 1ST POINT WHERE SHEAR = CRITICAL SHEAR (IF ONE EXISTS)
C     XC2   - 2ND POINT WHERE SHEAR = CRITICAL SHEAR (IF ONE EXISTS)
C     MSHEAR- FLAG INDICATING WHAT SHEAR CONDITIONS EXIST ON SEGMENT
C
C     MSHEAR |      MEANING
C     ----------------------------------------------------------------
C       1    |  SHEAR LESS THAN CRITICAL THROUGHOUT SEGMENT
C       2    |  SHEAR GREATER THAN CRITICAL THROUGHOUT SEGMENT
C       3    |  SHEAR EQUAL CRITICAL SHEAR WITHIN A SEGMENT
C            |      UPSLOPE  SHEAR < CRITICAL SHEAR
C            |      DOWNSLOPE SHEAR > CRITICAL SHEAR
C       4    |  SHEAR EQUAL CRITICAL SHEAR WITHIN A SEGMENT
C            |      SHEAR DECREASING WITHIN SEGMENT
C            |      UPSLOPE  SHEAR > CRITICAL SHEAR
C            |      DOWNSLOPE SHEAR < CRITICAL SHEAR
C       5    |  SHEAR EQUAL CRITICAL SHEAR TWO PLACES IN A SEGMENT
C            |      SHEAR INCREASES THEN DECREASES
C     ----------------------------------------------------------------
C
C     + + + LOCAL VARIABLES + + +
C
      DOUBLE PRECISION TAUCHK, X1, X2
C
C     + + + LOCAL VARIABLE DEFINITIONS + + +
C
C     TAUCHK - PARTIAL SOLUTION TO QUADRATIC EQUATION USED HERE
C     TAUB   - N.D. SHEAR STRESS CALCULATED AT THE BEG. OF SEGMENT
C     TAUE   - N.D. SHEAR STRESS CALCULATED AT THE END OF THE SEGMENT
C     PART   - PARTIAL SOLUTION TO QUADRATIC EQUATION USED HERE
C
C     + + + SUBROUTINES CALLED + + +
C
C     ROOT
C
C     + + + FUNCTION DECLARATIONS + + +
C
      REAL SHEAR
C
C     BEGIN SUBROUTINE XCRIT
C
      TAUCHK = TAUC ** 1.5 - C
      TAUB = SHEAR(A,B,C,XB)
      TAUE = SHEAR(A,B,C,XE)
C     
      IF (A.EQ.0.0) THEN
C        
C        UNIFORM SLOPE SEGMENT - DETERMINE LOCATION WHERE CRITICAL 
C        SHEAR STRESS IS EXCEEDED
C        
         IF (B.NE.0.0) THEN
            XC1 = TAUCHK / B
         ELSE
            XC1 = 1000.
         END IF
C        
         IF (TAUE.GT.TAUB) THEN
            MSHEAR = 3
            IF (XC1.LE.XB) MSHEAR = 2
            IF (XC1.GE.XE) MSHEAR = 1
         ELSE
            MSHEAR = 4
            IF (XC1.GE.XE) MSHEAR = 2
            IF (XC1.LE.XB) MSHEAR = 1
         END IF
C     
      ELSE IF (A.GT.0.0.AND.TAUE.GT.TAUB) THEN
C        
C        CONVEX SEGMENT ON AN PLANE ON WHICH SHEAR INCREASES DOWNSLOPE
C        
C        IF SHEAR STRESS AT THE BEGINNING OF THE CONVEX SEGMENT 
C        EXCEEDS CRITICAL, THEN THE ENTIRE SEGMENT EXCEEDS CRITICAL 
C        SHEAR STRESS
C        
         IF (TAUB.GE.TAUC) THEN
            MSHEAR = 2
         ELSE
C           
C           IF SHEAR STRESS AT THE END OF THE CONVEX SEGMENT IS LESS
C           THAN CRITICAL, THEN THE ENTIRE SEGMENT IS BELOW CRITICAL 
C           SHEAR STRESS
C           
            IF (TAUE.LE.TAUC) THEN
               MSHEAR = 1
            ELSE
               CALL ROOT(A,B,TAUCHK,X1,X2)
               MSHEAR = 3
C              
C              DETERMINE THE POINT WHERE SHEAR STRESS EXCEEDS
C              CRITICAL SHEAR STRESS
C              
               IF (X1.GE.XB.AND.X1.LE.XE) THEN
                  XC1 = X1
               ELSE
                  IF (X2.GE.XB.AND.X2.LE.XE) XC1 = X2
               END IF
            END IF
         END IF
C     
      ELSE
C        
C        ANY OTHER TYPE OF SEGMENT:
C        1) CONVEX WITH SHEAR DECREASING DOWN SEGMENT
C        2) CONCAVE WITH SHEAR INCREASING OR DECREASING DOWN SEGMENT
C        3) UNIFORM WITH SHEAR INCREASING OR DECREASING DOWN SEGMENT
C        
         IF (TAUE.GE.TAUC.AND.TAUB.GE.TAUC) THEN
C           
C           IF SHEAR STRESS EXCEEDS CRITICAL AT BOTH ENDS OF THE 
C           SEGMENT, IT EXCEEDS CRITICAL ALL ALONG THE SEGMENT
C           
            MSHEAR = 2
C        
         ELSE
C           
C           DETERMINE IF SHEAR EXCEEDS CRITICAL SHEAR SOMEWHERE
C           ALONG THE SEGMENT
C           
            PART = B ** 2 + 4.0 * A * TAUCHK
C           
C           IF SOLUTION OF QUADRATIC EQUATION HAS NO REAL ROOTS, THEN
C           CRITICAL SHEAR STRESS IS NOT EXCEEDED ALONG ENTIRE SEGMENT
C           
            IF (PART.LE.0.0) THEN
               MSHEAR = 1
            ELSE
C              
C              DETERMINE WHERE SHEAR STRESS EQUALS CRITICAL SHEAR
C              STRESS, BY SOLVING FOR ROOTS OF QUADRATIC EQUATION
C              
               CALL ROOT(A,B,TAUCHK,X1,X2)
C              
C              IF SHEAR STRESS INCREASES ON SEGMENT, SET MSHEAR 
C              FLAG = 3 AND RETURN LOCATION WHERE SHEAR = CRITICAL
C              
               IF (TAUB.LE.TAUC.AND.TAUE.GE.TAUC) THEN
                  MSHEAR = 3
                  IF (X1.LE.XB.OR.X1.GE.XE) THEN
                     XC1 = X2
                  ELSE
                     XC1 = X1
                  END IF
               ELSE
C                 
C                 IF SHEAR STRESS DECREASES ON SEGMENT, SET MSHEAR 
C                 FLAG = 4 AND RETURN LOCATION WHERE SHEAR = CRITICAL
C                 
                  IF (TAUB.GE.TAUC.AND.TAUE.LE.TAUC) THEN
                     MSHEAR = 4
                     IF (X1.LE.XB.OR.X1.GE.XE) THEN
                        XC1 = X2
                     ELSE
                        XC1 = X1
                     END IF
                  ELSE
C                    
C                    IF SHEAR AT BOTH TOP AND BOTTOM OF SEGMENT IS 
C                    BELOW CRITICAL, RETURN TWO LOCATIONS WHERE IT 
C                    EQUALS CRITICAL
C                    
                     IF (TAUB.LE.TAUC.AND.TAUE.LE.TAUC) THEN
                        MSHEAR = 5
                        XC1 = X1
                        XC2 = X2
C                       
C                       CHECK TO MAKE SURE THAT FOR A CASE 4, BOTH 
C                       POINTS FALL BETWEEN XB AND XE, AND THAT THEY 
C                       ARE NOT THE SAME POINT
C                       
                        IF (X1.LT.XB.OR.X1.GT.XE.OR.X2.LT.XB.OR.X2.GT.XE
     1                      .OR.X1.EQ.X2) MSHEAR = 1
                     END IF     ! IF (TAUB.LE.TAUC.AND.TAUE.LE.TAUC) THE
C                 
                  END IF        ! IF (TAUB.GE.TAUC.AND.TAUE.LE.TAUC) THE
C              
               END IF   ! IF (TAUB.LE.TAUC.AND.TAUE.GE.TAUC) THEN
C           
            END IF      ! IF (PART.LE.0.0) THEN
C        
         END IF ! IF (TAUE.GE.TAUC.AND.TAUB.GE.TAUC) THEN
C     
      END IF    ! IF (A.EQ.0.0) THEN 
C     
      RETURN
      END
