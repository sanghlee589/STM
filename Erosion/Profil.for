      SUBROUTINE PROFIL(A,B,AVGSLP,NSLPTS,SLPLEN,XINPUT,XU,XL,Y,TOTLEN)
C
C     SUBROUTINE PROFIL CALCULATES SLOPE INPUT COEFFICIENTS
C
C     MODULE ADAPTED FROM WEPP VERSION 2004.7 AND CALLED FROM THE
C     MAIN PROGRAM
C
C     AUTHOR(S): D.C FLANAGAN AND J.C. ASCOUGH II
C     DATE LAST MODIFIED: 4-1-2005
C
C     + + + PARAMETER DECLARATIONS + + + 
C
      INTEGER MXSLP
      PARAMETER (MXSLP = 40)
C
C     + + + ARGUMENT DECLARATIONS + + +
C
      REAL A(MXSLP), B(MXSLP), AVGSLP, SLPLEN, XINPUT(101), XU(MXSLP), 
     1    XL(MXSLP), Y(101), TOTLEN
      INTEGER NSLPTS
C      
C     + + + ARGUMENT DEFINITIONS + + +
C
C
C     + + + LOCAL DECLARATIONS + + + 
C                              
      REAL SLPINP(MXSLP), SLEN, SSTAR(MXSLP), XSTAR(MXSLP),
     1     YL(MXSLP), YU(MXSLP), C(MXSLP) 
C
C     + + + LOCAL VARIABLE DEFINITIONS + + + 
C
C     BEGIN SUBROUTINE PROFIL
C       
      READ (7,*) NSLPTS, SLPLEN
      READ (7,*) (XINPUT(J),SLPINP(J),J = 1,NSLPTS)
C     
      SLEN = XINPUT(NSLPTS)
      Y(NSLPTS) = 0.0
C     
      DO K = 1, NSLPTS - 1
         KM = NSLPTS - K
         Y(KM) = Y(KM+1) + (XINPUT(KM+1) - XINPUT(KM)) * 
     1           (SLPINP(KM) + SLPINP(KM+1)) / 2.0
      END DO
C     
      AVGSLP = Y(1) / SLEN
C     
      IF (AVGSLP.LE.0.0) AVGSLP = 0.000001
C     
      DO K = 1, NSLPTS
         SSTAR(K) = SLPINP(K) / AVGSLP
         XSTAR(K) = XINPUT(K) / SLEN
      END DO
C     
      DO K = 2, NSLPTS
         A(K) = (SSTAR(K)-SSTAR(K-1)) / (XSTAR(K)-XSTAR(K-1))
         B(K) = SSTAR(K-1) - A(K) * XSTAR(K-1)
      END DO
C     
      YL(1) = 1.0
      XL(1) = 0.0
C     
      DO K = 2, NSLPTS
         YU(K) = YL(K-1)
         XU(K) = XL(K-1)
         C(K) = YU(K) + A(K) * XSTAR(K-1) ** 2 / 2.0 + B(K) * XSTAR(K-1)
         YL(K) = -A(K) * XSTAR(K) ** 2 / 2.0 - B(K) * XSTAR(K) + C(K)
         XL(K) = XSTAR(K)
      END DO
C     
      K = 2
      Y(1) = 1.0
C     
      DO L = 2, 101
         XINPUT(L) = FLOAT(L-1) * 0.01
   10    IF (XINPUT(L).GT.XSTAR(K)) THEN
            K = K + 1
            GO TO 10
         END IF
         Y(L) = -A(K) * XINPUT(L) ** 2 / 2.0 - B(K) * XINPUT(L) + C(K)
      END DO
C     
      TOTLEN = SLPLEN
C
      RETURN
      END
