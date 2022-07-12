      SUBROUTINE TRCOEFF(TRCOEF,SHRSOL,SAND,DIA,SPG,TCF1,NPART,FRAC)
C
C************************************************************
C                                                           *
C   THIS FUNCTION IS CALLED FROM SR PARAM TO COMPUTE THE    *
C   SEDIMENT TRANSPORT COEFFICIENT FOR RILL FLOW. IT CALLS  *
C   SR YALIN.                                               *
C                                                           *
C************************************************************
C                                                           *
C   ARGUMENT                                                *
C      SHRSOL  : FLOW SHEAR STRESS                          *
C                                                           *
C************************************************************
C
C     + + + PARAMETER DECLARATIONS + + +
C
      INTEGER MXPART, MXNSL
      PARAMETER (MXPART = 10,MXNSL = 10)
C
      REAL SAND(MXNSL), DIA(MXPART), SPG(MXPART), TCF1(MXPART),
     1     FRAC(MXPART)
      INTEGER NPART
C      
      CALL YALIN(SHRSOL,TOTTC,SAND,DIA,SPG,TCF1,NPART,FRAC)
      TRCOEF = TOTTC / SHRSOL ** 1.5
      IF (TRCOEF.EQ.0.0) TRCOEF = 0.000000001
C     
      RETURN
      END
