      FUNCTION SHIELD(REYN)
C
C     + + + PURPOSE + + +
C
C     FUNCTION SHIELD GENERATES PARAMETERS (SHIELD PARAMETERS)
C     BY INTERPOLATING VALUES FROM A TABLE (SHIELD DIAGRAM)
C     FOR GIVEN REYNOLDS NUMBERS.
C
C     CALLED FROM: SRS TRNCAP, YALIN
C     AUTHOR(S): ASCOUGH II, R. VAN DER ZWEEP, V. LOPES
C     REFERENCE IN USER GUIDE:
C
C     VERSION:
C     DATE RECODED:
C     RECODED BY: JIM ASCOUGH II
C
C     + + + KEYWORDS + + +
C
C     + + + PARAMETERS + + +
C
C     + + + ARGUMENT DECLARATIONS + + +
C
      REAL REYN
C
C     + + + ARGUMENT DEFINITIONS + + +
C
C     REYN -
C
C     + + + COMMON BLOCKS + + +
C
C     + + + LOCAL VARIABLES + + +
C
      REAL Y(8), R(8), SHIELD, SLOPE, YCR
      INTEGER I
C
C     + + + LOCAL DEFINITIONS + + +
C
C     REAL VARIABLES
C
C     Y(8)   -
C     R(8)   -
C     SHIELD -
C     SLOPE  -
C     YCR    -
C
C     INTEGER VARIABLES
C
C     I -
C
C     + + + SAVES + + +
C
      SAVE
C
C     + + + SUBROUTINES CALLED + + +
C
C     + + + DATA INITIALIZATIONS + + +
C
      DATA Y /0.0772, 0.0579, 0.04, 0.035, 0.034, 0.045, 0.055, 0.057/
      DATA R /1.0, 2.0, 4.0, 8.0, 12.0, 100.0, 400.0, 1000.0/
C
C     + + + END SPECIFICATIONS + + +
C
C
      IF (REYN.LT.R(1)) THEN
         I = 2
         SLOPE = (ALOG(Y(I))-ALOG(Y(I-1))) / (ALOG(R(I))-ALOG(R(I-1)))
         YCR = ALOG(Y(1)) - SLOPE * (ALOG(R(1))-ALOG(REYN))
      ELSE IF (REYN.GT.R(8)) THEN
         I = 8
         SLOPE = (ALOG(Y(I))-ALOG(Y(I-1))) / (ALOG(R(I))-ALOG(R(I-1)))
         YCR = Y(8) + SLOPE * (ALOG(REYN)-ALOG(R(8)))
C     
      ELSE
C        
         DO 10 I = 2, 8
C           
            IF (REYN.GE.R(I-1).AND.REYN.LE.R(I)) THEN
               SLOPE = (ALOG(Y(I))-ALOG(Y(I-1))) / (ALOG(R(I))-
     1             ALOG(R(I-1)))
               YCR = ALOG(Y(I-1)) + SLOPE * (ALOG(REYN)-ALOG(R(I-1)))
               GO TO 20
            END IF
C        
   10    CONTINUE
C     
      END IF
C     
   20 SHIELD = EXP(YCR)
C     
      RETURN
      END
