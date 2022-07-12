      SUBROUTINE YALIN(EFFSH,TOTTC,SAND,DIA,SPG,TCF1,NPART,FRAC)
C                                                                   
C      THIS ROUTINE IS CALLED FROM FN TRCOEF TO COMPUTE SEDIMENT    
C      TRANSPORT CAPACITY USING THE YALIN EQUATION. IT IS CALLED    
C      FN SHIELD.                                                   
C                                                                   
C
C     + + + PARAMETER DECLARATIONS + + +
C
      INTEGER MXPART, MXNSL
      PARAMETER (MXPART = 10,MXNSL = 10)
C
      REAL EFFSH, TOTTC
      REAL DIA(MXPART), SPG(MXPART), SAND(MXNSL), TCF1(MXPART),
     1     FRAC(MXPART)
      INTEGER NPART
C                                                                   
C      ARGUMENTS                                                    
C         EFFSH - EFFECTIVE SHEER STRESS                            
C         TOTTC - TOTAL SEDIMENT TRANSPORT CAPACITY                 
C                                                                   
C    LOCAL VARIABLES:                                               *
C       WS     : SEDIMENT TRANSPORT CAPACITY FOR PARTICLE CLASS K   *
C                (KG/M*S)                                           *
C       COEF   : PORTION OF YALIN EQUATION SOLUTION                 *
C       YCRIT  : ORDINATE FROM SHIELDS DIAGRAM FOR DIMENSIONLESS    *
C                CRITICAL SHEAR FOR TRANSPORT                       *
C       DELTA  : PORTION OF YALIN EQUATION SOLUTION                 *
C       SIGMA  : PORTION OF YALIN EQUATION SOLUTION                 *
C       P      : SEDIMENT TRANSPORT CAPACITY FOR PARTICLE CLASS     *
C                K (NONDIMENSIONAL)                                 *
C       DLTRAT : PORTION OF YALIN EQUATION SOLUTION                 *
C       TOTTC  : TOTAL SEDIMENT TRANSPORT CAPACITY (KG/M*S)         *
C                                                                   *
C********************************************************************
C
      SAVE
      REAL WS(MXPART), COEF(MXPART), YCRIT(MXPART), DELTA(MXPART), 
     1    SIGMA(MXPART), P(MXPART), DLTRAT(MXPART), REYN, SHIELD, T, 
     1    VSTAR, YALCON, OLDTOT, ADJTC, MSDENS, KINVIS
C
C INITIALIZATION:
C
      YALCON = 0.635
      T = 0.0
      TOTTC = 0.0
C     
      MSDENS = 1000.0
      KINVIS = 1.0E-06
      ACCGAV = 9.807
C     
C     THE CONSTANT 0.635 WAS DERIVED EMPIRICALLY BY YALIN
C     
C     COMPUTE SHEAR VELOCITY (VSTAR):
C     
      VSTAR = SQRT(EFFSH/MSDENS)
C     
C     COMPUTE COEFFICIENT COEF=VSTAR*MSDENS*DIA*SPG FOR EACH
C     PARTICLE CLASSES:
C     
C     THIS IS GOOFY - CHECK OUT
C     
      COEF(NPART) = VSTAR * MSDENS
      DO 10 K = 1, NPART
         COEF(K) = COEF(NPART) * DIA(K) * SPG(K)
   10 CONTINUE
C     
C     COMPUTE REYNOLD'S NUMBER (REYN), DIMENSIONLESS CRITICAL SHEAR
C     PARAMETER FROM THE SHIELDS DIAGRAM (YCRIT), PARAMETERS DELTA AND
C     SIGMA, AND THE DIMENSIONLESS SEDIMENT TRANSPORT CAPACITY (P) FOR
C     EACH PARTICLE CLASS:
C     
      DO 20 K = 1, NPART
         REYN = VSTAR * DIA(K) / KINVIS
         YCRIT(K) = SHIELD(REYN)
         DELTA(K) = (VSTAR**2/(SPG(K)-1.0)/ACCGAV/DIA(K)/
     1               YCRIT(K)) - 1.0
C        
         IF (DELTA(K).GT.0.0) THEN
            SIGMA(K) = DELTA(K) * 2.45 * SPG(K) ** (-0.4) * 
     1          SQRT(YCRIT(K))
            P(K) = YALCON * DELTA(K) * (1.0-1.0/SIGMA(K)*
     1          ALOG(1.0+SIGMA(K)))
            T = T + DELTA(K)
         ELSE
            DELTA(K) = 0.0
            P(K) = 0.0
         END IF
   20 CONTINUE  
C     
C     COMPUTE THE TRANSPORT CAPACITY (MASS PER UNIT WIDTH PER UNIT TIME)
C     WS FOR EACH PARTICLE CLASS:
C     
      IF (T.EQ.0.0) T = 1000.0
      DO 30 K = 1, NPART
         DLTRAT(K) = DELTA(K) / T
         WS(K) = P(K) * DLTRAT(K) * COEF(K)
C        
C        RISSE 9/20/94 ADD WEIGHTING SCHEME TO TRANSPORT CAPACITY TO ACCOUNT
C        FOR THE AMOUNT OF EACH SEDIMENT CLASS BEING TRANSPORTED
C        NOTE: THIS WILL INCREASE TC FOR CLAYS AND SILTS AND DECREASE FOR SANDS
C        
C        WS(K)=WS(K)*(FRAC(K,IPLANE)/0.2)
C        REPLACED ABOVE EQUATION WITH FOLLOWING TO TAKE CARE OF SITUATION
C        WHERE OTHER THAN 5 PARTICLE SIZE CLASSES HAVE BEEN USED.
C        DCF 11/18/94
C        XXX ADDITIONAL QUESTION ON THIS IS WHETHER WE WANT TO USE
C        XXX "FRAC(K,IPLANE)"  OR  "FRCFLW(K,IPLANE)" FOR THE WEIGHTING?????
C        
         WS(K) = WS(K) * (FRAC(K)*FLOAT(NPART))
C        
         TOTTC = TOTTC + WS(K)
   30 CONTINUE
C     
C     
C     XXX ADD CHANGES TO INCLUDE NEARING ALTERATION TO TC THAT WAS
C     XXX PREVIOUSLY INCLUDED IN TCEND CALCULATION IN PARAM.FOR.  DCF
      OLDTOT = TOTTC
      IF (SAND(1).GT.0.5) THEN
         ADJTC = 0.3 + 0.7 * EXP(-12.52*(SAND(1)-0.5))
         IF (ADJTC.LT.0.30) ADJTC = 0.30
         TOTTC = TOTTC * ADJTC
      END IF
C     
      DO 40 K = 1, NPART
         IF (OLDTOT.GT.0.0) WS(K) = (WS(K)/OLDTOT) * TOTTC
   40 CONTINUE
C     
      DO 50 K = 1, NPART
C        
C        LINES ADDED 10/16/89 TO PREVENT DIVIDE BY ZERO
C        IF TOTTC IS ZERO.
C        
         IF (TOTTC.GT.0.0) THEN
            TCF1(K) = WS(K) / TOTTC
         ELSE
            TCF1(K) = 0.0
         END IF
   50 CONTINUE
      RETURN
      END
