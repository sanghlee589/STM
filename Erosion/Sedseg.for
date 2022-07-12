      SUBROUTINE SEDSEG(DSLOST,JUN,IYEAR,NOOUT,DSTOT,STDIST,IRDGDX,
     1    YSDIST,AVGSLP,SLPLEN,Y)
C
C     + + + PURPOSE + + +
C     THIS SUBROUTINE BREAKS THE HILLSLOPE PROFILE INTO DETACHMENT
C     AND DEPOSITION SEGMENTS.  IT CALLS SR SEDIST AND SR SEDSTA
C
C     CALLED FROM SUBROUTINE SEDOUT
C     AUTHOR(S): M. NEARING, D. FLANAGAN
C
C     VERSION: THIS MODULE MODIFIED FROM WEPP V2004.7 CODE
C     DATE MODIFIED: 10-2004
C     MODIFIED BY: D. FLANAGAN
C
C     + + + ARGUMENT DECLARATIONS + + +
C
      INTEGER JUN, IYEAR, NOOUT
      REAL DSLOST(100), DSTOT(1000), STDIST(1000), IRDGDX, 
     1    YSDIST(1000), AVGSLP, SLPLEN, Y(101)
C
C     + + + ARGUMENT DEFINITIONS + + +
C
C     DSLOST  - ARRAY CONTAINING THE NET SOIL LOSS/GAIN AT EACH
C               OF THE 100 POINTS ON EACH OFE
C     JUN     - UNIT NUMBER OF FILE TO WRITE OUTPUT TO
C     IYEAR   - FLAG FOR PRINTING OUT ANNUAL SOIL LOSS OUTPUT
C     NOOUT   - FLAG INDICATING PRINTING OF EVENT BY EVENT
C               SUMMARY FILES IS DESIRED (SEE SEDOUT)
C     DSTOT   -
C     STDIST  -
C     IRDGDX  - AVERAGE INTERRILL DETACHMENT RATE ON AN OFE
C
C
C     + + + LOCAL VARIABLES + + +
C
      REAL DETDIS(100), DEPDIS(100), DETPT1(100), DETPT2(100), 
     1    DTAVLS(100)
      REAL DEPPT1(100), DEPPT2(100), DPAVLS(100)
      REAL SUM1, SUM2, FILOSS, FIDEP, TOTMAX, PDTMX, PDPMX
      REAL DEPMAX(100), DEPMIN(100), DEPSTD(100)
      REAL DETMAX(100), DETMIN(100), DETSTD(100)
      REAL PDPMAX(100), PDPMIN(100), PDTMAX(100), PDTMIN(100)
      REAL DELXX
C     REAL MTF,KGMTPA
      INTEGER IDTSIN(100), IDPSIN(100), JADET, JADEP, JDET, ITSIN, IPSIN
      INTEGER LEND, I, ICNT, J, KBEG, KK, LBEG
      INTEGER JFLAG(100), IBEGIN, LSEG
      INTEGER IMODEL, IOUTPT, IOUTSS
C     CHARACTER*7 UNIT(9)
C
C     + + + LOCAL DEFINITIONS + + +
C
C     DEPPT1 : DISTANCE WHERE DEPOSITION SECTION BEGINS (M)
C     DEPPT2 : DISTANCE WHERE DEPOSITION SECTION ENDS (M)
C     DPAVLS : AVERAGE DEPOSITION IN SECTION (KG/M**2)
C     DETPT1 : DISTANCE WHERE DETACHMENT SECTION BEGINS (M)
C     DETPT2 : DISTANCE WHERE DETACHMENT SECTION ENDS (M)
C     DTAVLS : AVERAGE DETACHMENT IN SECTION (KG/M**2)
C     DETSTD : STANDARD DEVIATION OF DETACHMENT IN SECT.(KG/M**2)
C     DETMAX : MAXIMUM DETACHMENT IN SECTION (KG/M**2)
C     PDTMAX : POINT OF MAXIMUM DETACHMENT (M)
C     DETMIN : MINIMUM DETACHMENT IN SECTION (KG/M**2)
C     PDTMIN : POINT OF MINIMUM DETACHMENT (M)
C     DEPSTD : STANDARD DEVIATION OF DEPOSITION IN SECT.(KG/M**2)
C     DEPMAX : MAXIMUM DEPOSITION IN SECTION (KG/M**2)
C     PDPMAX : POINT OF MAXIMUM DEPOSITION IN SECTION (M)
C     DEPMIN : MINIMUM DEPOSITION IN SECTION (KG/M**2)
C     PDPMIN : POINT OF MINIMUM DEPOSITION IN SECTION (M)
C     DELXX  - DELTA X INCREMENTS BETWEEN EACH POINT DOWN OFE (M)
C
C     DATA UNIT /'     MM', '  KG/M2','    IN.','    T/A','    FT.'
C    1,' LBS/FT','    LBS','   KG/M','      M'/
C
C
C     + + + END SPECIFICATIONS + + +
C
C     AVEDEP = 0.0
C     AVEDET = 0.0
C     MAXDEP = 0.0
C     MAXDET = 0.0
C     PTDEP = 0.0
C     PTDET = 0.0
C     TDEP(IHILL) = 0.0
C     TDET(IHILL) = 0.0
      IMODEL = 2
      IOUTPT = 1
      IOUTSS = 0
C     
      WRITE (JUN,1000)
C     
      CALL SEDIST(DSLOST,DSTOT,STDIST,DELXX,SLPLEN,AVGSLP,Y,
     1    YSDIST)
C     
      LSEG = 1
      JADEP = 0
      JADET = 0
      SUM1 = 0.0
      SUM2 = 0.0
      JDET = 0
      JDEP = 0
      ITSIN = 0
      IPSIN = 0
      LEND = 100
      KBEG = 0
      LBEG = 1
      ICNT = 0
C     
      DO 10 J = 1, LEND
         IF (DSTOT(J).NE.0.0) THEN
            IF (KBEG.EQ.0) THEN
               LBEG = J
               KBEG = 1
            END IF
            ICNT = J
         END IF
   10 CONTINUE
C     
      IBEGIN = LBEG
      IF (ICNT.LT.LEND) THEN
         LEND = ICNT + 1
      END IF
C     
      IF (DSTOT(LBEG).GT.0.0) JFLAG(LSEG) = 1
      IF (DSTOT(LBEG).LT.0.0) JFLAG(LSEG) = 0
C     
C     
C     ADDED BY DCF 4/16/90 TO COVER POSSIBILITY IF DSTOT = 0
C     
      IF (DSTOT(LBEG).EQ.0.0) JFLAG(LSEG) = 2
C     
      DO 20 I = LBEG + 1, LEND
         IF ((JFLAG(LSEG).EQ.1.AND.DSTOT(I).LE.0.0).OR.(I.EQ.LEND.AND.
     1       DSTOT(I).GT.0.0)) THEN
C           
            JADET = JADET + 1
            IEND = I - 1
            IF (I.EQ.LEND.AND.DSTOT(I).GT.0.0) IEND = I
C           
C           IF THE BEGINNING OF THE DETACHMENT IS THE FIRST POINT ON
C           THE SLOPE SET THE POINT TO ZERO OTHERWISE AVERAGE I 
C           WITH THE POINT BEFORE IT
C           
            IF (IBEGIN.EQ.1) THEN
               DETPT1(JADET) = 0.0
            ELSE
               DETPT1(JADET) = STDIST(IBEGIN-1)
            END IF
C           
            DETPT2(JADET) = STDIST(IEND)
C           
            DETDIS(JADET) = DETPT2(JADET) - DETPT1(JADET)
C           
            IF (I.EQ.LEND.AND.DSTOT(I).LT.0.0) THEN
               JADEP = JADEP + 1
               IDPSIN(JADEP) = 1
               DPAVLS(JADEP) = DSTOT(LEND)
               DEPSTD(JADEP) = 0.0
               DEPPT1(JADEP) = STDIST(LEND-1)
               DEPPT2(JADEP) = STDIST(LEND)
               DEPDIS(JADEP) = DEPPT2(JADEP) - DEPPT1(JADEP)
               DEPMAX(JADEP) = DSTOT(LEND)
               PDPMAX(JADEP) = STDIST(LEND)
               DEPMIN(JADEP) = DSTOT(LEND)
               PDPMIN(JADEP) = STDIST(LEND)
               IPSIN = 1
               JDEP = 1
            END IF
            IF (IBEGIN.EQ.IEND) THEN
               IDTSIN(JADET) = 1
               DETPT1(JADET) = STDIST(IBEGIN-1)
               DETPT2(JADET) = STDIST(IBEGIN)
               DETDIS(JADET) = DETPT2(JADET) - DETPT1(JADET)
               DTAVLS(JADET) = DSTOT(IBEGIN)
               DETSTD(JADET) = 0.0
               DETMAX(JADET) = DSTOT(IBEGIN)
               PDTMAX(JADET) = STDIST(IBEGIN)
               DETMIN(JADET) = DSTOT(IBEGIN)
               PDTMIN(JADET) = STDIST(IBEGIN)
               ITSIN = 1
               JDET = 1
            ELSE
               IDTSIN(JADET) = 0
               CALL SEDSTA(JADET,DTAVLS,DETSTD,DETMAX,PDTMAX,DETMIN,
     1             PDTMIN,IBEGIN,IEND,JFLAG,LSEG,DSTOT,STDIST,DELXX)
               JDET = 1
            END IF
C           
            IBEGIN = IEND + 1
            LSEG = LSEG + 1
            IF (DSTOT(I).EQ.0.0) JFLAG(LSEG) = 2
            IF (DSTOT(I).LT.0.0) JFLAG(LSEG) = 0
C        
         ELSE IF ((JFLAG(LSEG).EQ.0.AND.DSTOT(I).GE.0.0).OR.(I.EQ.LEND
     1       .AND.DSTOT(I).LT.0.0)) THEN
C           
            JADEP = JADEP + 1
            IEND = I - 1
            IF (I.EQ.LEND.AND.DSTOT(I).LT.0.0) IEND = I
            IF (IBEGIN.EQ.1) THEN
               DEPPT1(JADEP) = 0.0
            ELSE
               DEPPT1(JADEP) = STDIST(IBEGIN-1)
            END IF
C           
            DEPPT2(JADEP) = STDIST(IEND)
            DEPDIS(JADEP) = DEPPT2(JADEP) - DEPPT1(JADEP)
C           
            IF (I.EQ.LEND.AND.DSTOT(I).GT.0.0) THEN
               JADET = JADET + 1
               IDTSIN(JADET) = 1
               DTAVLS(JADET) = DSTOT(LEND)
               DETSTD(JADET) = 0.0
               DETPT1(JADET) = STDIST(LEND-1)
               DETPT2(JADET) = STDIST(LEND)
               DETDIS(JADET) = DETPT2(JADET) - DETPT1(JADET)
               DETMAX(JADET) = DSTOT(LEND)
               DETMIN(JADET) = DSTOT(LEND)
               PDTMAX(JADET) = STDIST(LEND)
               PDTMIN(JADET) = STDIST(LEND)
               ITSIN = 1
               JDET = 1
            END IF
            IF (IBEGIN.EQ.IEND) THEN
               IDPSIN(JADEP) = 1
               DEPPT1(JADEP) = STDIST(IBEGIN-1)
               DEPPT2(JADEP) = STDIST(IBEGIN)
               DEPDIS(JADEP) = DEPPT2(JADEP) - DEPPT1(JADEP)
               DPAVLS(JADEP) = DSTOT(IBEGIN)
               DEPMIN(JADEP) = DSTOT(IBEGIN)
               DEPMAX(JADEP) = DSTOT(IBEGIN)
               DEPSTD(JADEP) = 0.0
               PDPMAX(JADEP) = STDIST(IBEGIN)
               PDPMIN(JADEP) = STDIST(IBEGIN)
               IPSIN = 1
               JDEP = 1
            ELSE
               IDPSIN(JADEP) = 0
               CALL SEDSTA(JADEP,DPAVLS,DEPSTD,DEPMAX,PDPMAX,DEPMIN,
     1             PDPMIN,IBEGIN,IEND,JFLAG,LSEG,DSTOT,STDIST,DELXX)
               JDEP = 1
            END IF
C           
            IBEGIN = IEND + 1
            LSEG = LSEG + 1
            IF (DSTOT(I).EQ.0.0) JFLAG(LSEG) = 2
            IF (DSTOT(I).GT.0.0) JFLAG(LSEG) = 1
C        
         ELSE IF (JFLAG(LSEG).EQ.2.AND.DSTOT(I).NE.0.0) THEN
            IEND = I - 1
            IF (I.EQ.LEND) THEN
               IF (JFLAG(LSEG).EQ.1) THEN
                  JADET = JADET + 1
                  IDTSIN(JADET) = 1
                  DETPT1(JADET) = STDIST(LEND-1)
                  DETPT2(JADET) = STDIST(LEND)
                  DETDIS(JADET) = DETPT2(JADET) - DETPT1(JADET)
                  DTAVLS(JADET) = DSTOT(LEND)
                  DETSTD(JADET) = 0.0
                  DETMAX(JADET) = DSTOT(LEND)
                  DETMIN(JADET) = DSTOT(LEND)
                  PDTMIN(JADET) = STDIST(LEND)
                  PDTMAX(JADET) = STDIST(LEND)
                  ITSIN = 1
                  JDET = 1
               ELSE IF (JFLAG(LSEG).EQ.0) THEN
                  JADEP = JADEP + 1
                  IDPSIN(JADEP) = 1
                  DEPPT1(JADEP) = LEND
                  DEPPT2(JADEP) = LEND
                  DEPDIS(JADEP) = DEPPT2(JADEP) - DEPPT1(JADEP)
                  DPAVLS(JADEP) = DSTOT(LEND)
                  DEPSTD(JADEP) = 0.0
                  DEPMAX(JADEP) = DSTOT(LEND)
                  DEPMIN(JADEP) = DSTOT(LEND)
                  PDPMIN(JADEP) = STDIST(LEND)
                  PDPMAX(JADEP) = STDIST(LEND)
                  IPSIN = 1
                  JDEP = 1
               END IF
            END IF
            LSEG = LSEG + 1
            IBEGIN = IEND + 1
            IF (DSTOT(I).GT.0.0) JFLAG(LSEG) = 1
            IF (DSTOT(I).LT.0.0) JFLAG(LSEG) = 0
         END IF
C     
   20 CONTINUE
C     
C     
      IF (JADET.GT.0) THEN
         IF (NOOUT.LE.1) WRITE (JUN,1100)
         IF (JDET.GT.0) THEN
            TOTMAX = DETMAX(1)
            PDTMX = PDTMAX(1)
            DO 30 KK = 1, JADET
               SUM1 = SUM1 + (DTAVLS(KK)*DETDIS(KK))
               SUM2 = SUM2 + DETDIS(KK)
C              
               IF (DETMAX(KK).GT.TOTMAX) THEN
                  TOTMAX = DETMAX(KK)
                  PDTMX = PDTMAX(KK)
               END IF
C           
   30       CONTINUE
C           
            IF (SUM2.NE.0.0) THEN
               FILOSS = SUM1 / SUM2
C              TDET(IHILL) = SUM2*FWIDTH*FILOSS
C              KG/M2>T/A
C              KGMTPA=4.4605
C              M>FT
C              MTF=3.2808
               IF (NOOUT.LE.1) THEN
C                 IF(OUTOPT.EQ.1.AND.UNITS.EQ.1)THEN
C                 ABBREVIATED ENGLISH UNITS
C                 WRITE (JUN,2700) FILOSS*KGMTPA,UNIT(4)
C                 WRITE (JUN,2800) TOTMAX*KGMTPA,UNIT(4),
C                 1               PDTMX*MTF,UNIT(5)
C                 ELSE
C                 METRIC UNITS
                  WRITE (JUN,1200) FILOSS
                  WRITE (JUN,1300) TOTMAX, PDTMX
               END IF
               IF (IMODEL.EQ.2.OR.(IMODEL.EQ.1.AND.IYEAR.NE.1.AND.IOUTPT
     1             .EQ.1)) THEN
C                 DO 40 III = 1, NPLANE
C                    IF(OUTOPT.EQ.1.AND.UNITS.EQ.1)THEN
C                    ABBREVIATED ENGLISH UNITS
C                    WRITE (JUN,2900) IRDGDX(III)*KGMTPA,UNIT(4), III
C                    ELSE
                     WRITE (JUN,1400) IRDGDX
C                 END IF
C  40             CONTINUE
               END IF
            END IF
C        AVEDET = FILOSS
C        MAXDET = TOTMAX
C        PTDET = PDTMX
         END IF
         IF (NOOUT.LE.1) THEN
C           IF(OUTOPT.EQ.1.AND.UNITS.EQ.1)THEN
C           ABBREVIATED ENGLISH UNITS
C           IF (IOUTSS.NE.2) WRITE (JUN,3500)UNIT(5),UNIT(4),
C           1              UNIT(4),UNIT(4),UNIT(5),UNIT(4),UNIT(5)
C           ELSE
            IF (IOUTSS.NE.2) WRITE (JUN,2000)
C        END IF
         END IF
C        
         DO 50 KK = 1, JADET
            IF (IOUTSS.NE.2) THEN
               IF (NOOUT.LE.1) THEN
C                 IF(OUTOPT.EQ.1.AND.UNITS.EQ.1)THEN
C                 ABBREVIATED ENGLISH UNITS
C                 WRITE (JUN,3300) DETPT1(KK)*MTF,DETPT2(KK)*MTF,
C                 1               DTAVLS(KK)*KGMTPA,DETSTD(KK)*KGMTPA,
C                 1               DETMAX(KK)*KGMTPA, PDTMAX(KK)*MTF,
C                 1               DETMIN(KK)*KGMTPA, PDTMIN(KK)*MTF
C                 ELSE
C                 METRIC UNITS
                  WRITE (JUN,1800) DETPT1(KK), DETPT2(KK), DTAVLS(KK), 
     1                DETSTD(KK), DETMAX(KK), PDTMAX(KK), DETMIN(KK), 
     1                PDTMIN(KK)
C              END IF
               END IF
            END IF
C        
C        ENDIF
C        
   50    CONTINUE
      END IF
C     
C     OUTPUT SUMMARY OF EROSION AND DETATCHMENT TO ABBREV.RAW
C     
C     
      IF (ITSIN.GT.0) THEN
         DO 60 KK = 1, JADET
            IF (IDTSIN(KK).EQ.1) THEN
               IF (IOUTSS.NE.2) THEN
C                 IF(OUTOPT.EQ.1.AND.UNITS.EQ.1)THEN
C                 ABBREVIATED ENGLISH UNITS
C                 IF (NOOUT.LE.1) WRITE (JUN,3400)UNIT(5),UNIT(4)
C                 ELSE
                  IF (NOOUT.LE.1) WRITE (JUN,1900)
C                 END IF
C                 
                  GO TO 70
               END IF
            END IF
C        
   60    CONTINUE
C        
   70    DO 80 KK = 1, JADET
            IF (IDTSIN(KK).EQ.1) THEN
               IF (IOUTSS.NE.2) THEN
C                 IF(OUTOPT.EQ.1.AND.UNITS.EQ.1)THEN
C                 ABBREVIATED ENGLISH UNITS
C                 IF (NOOUT.LE.1) WRITE (JUN,3600) DETPT2(KK)*MTF,
C                 1              DTAVLS(KK)*KGMTPA
C                 ELSE
                  IF (NOOUT.LE.1) WRITE (JUN,2100) DETPT2(KK), 
     1                DTAVLS(KK)
C              END IF
               END IF
            END IF
   80    CONTINUE
      END IF
C     END IF  !  NEED TO CHECK THE NESTING HERE
C     
C     
      IF (JADEP.GT.0) THEN
         IF (NOOUT.LE.1) WRITE (JUN,1500)
         IF (JDEP.GT.0) THEN
            SUM1 = 0.0
            SUM2 = 0.0
            TOTMAX = DEPMAX(1)
            PDPMX = PDPMAX(1)
C           
            DO 90 KK = 1, JADEP
               SUM1 = SUM1 + (DPAVLS(KK)*DEPDIS(KK))
               SUM2 = SUM2 + DEPDIS(KK)
               IF (DEPMAX(KK).LT.TOTMAX) THEN
                  TOTMAX = DEPMAX(KK)
                  PDPMX = PDPMAX(KK)
               END IF
C           
   90       CONTINUE
C           
            IF (SUM2.NE.0.0) THEN
               FIDEP = SUM1 / SUM2
C              TDEP(IHILL)=SUM2*FWIDTH*FIDEP
               IF (NOOUT.LE.1) THEN
C                 IF(OUTOPT.EQ.1.AND.UNITS.EQ.1)THEN
C                 ABBREVIATED ENGLISH UNITS
C                 WRITE (JUN,3100) FIDEP*KGMTPA,UNIT(4)
C                 WRITE (JUN,3200) TOTMAX*KGMTPA,UNIT(4),PDPMX*MTF,UNIT(5)
C                 ELSE
                  WRITE (JUN,1600) FIDEP
                  WRITE (JUN,1700) TOTMAX, PDPMX
C              END IF
               END IF
C           AVEDEP = FIDEP
C           MAXDEP = TOTMAX
C           PTDEP = PDPMX
            END IF
            IF (NOOUT.LE.1) THEN
C              IF(OUTOPT.EQ.1.AND.UNITS.EQ.1)THEN
C              ABBREVIATED ENGLISH UNITS
C              IF (IOUTSS.NE.2) WRITE (JUN,3700)UNIT(5),UNIT(4),UNIT(4),
C              1           UNIT(4),UNIT(5),UNIT(4),UNIT(5)
C              ELSE
               IF (IOUTSS.NE.2) WRITE (JUN,2200)
C           END IF
            END IF
C           
            DO 100 KK = 1, JADEP
C              IF(IDPSIN(KK).NE.1) THEN
               IF (IOUTSS.NE.2) THEN
                  IF (NOOUT.LE.1) THEN
C                    IF(OUTOPT.EQ.1.AND.UNITS.EQ.1)THEN
C                    ABBREVIATED ENGLISH UNITS
C                    
C                    WRITE (JUN,3300) DEPPT1(KK)*MTF, DEPPT2(KK)*MTF,
C                    1              DPAVLS(KK)*KGMTPA, DEPSTD(KK)*KGMTPA,
C                    1              DEPMAX(KK)*KGMTPA, PDPMAX(KK)*MTF,
C                    1              DEPMIN(KK)*KGMTPA, PDPMIN(KK)*MTF
C                    ELSE
                     WRITE (JUN,1800) DEPPT1(KK), DEPPT2(KK), 
     1                   DPAVLS(KK), DEPSTD(KK), DEPMAX(KK), 
     1                   PDPMAX(KK), DEPMIN(KK), PDPMIN(KK)
C                 END IF
                  END IF
               END IF
C           ENDIF
  100       CONTINUE
         END IF
         IF (IPSIN.GT.0) THEN
            DO 110 KK = 1, JADEP
               IF (IDPSIN(KK).EQ.1) THEN
                  IF (IOUTSS.NE.2) THEN
C                    IF(OUTOPT.EQ.1.AND.UNITS.EQ.1)THEN
C                    ABBREVIATED ENGLISH UNITS
C                    
C                    IF (NOOUT.LE.1) WRITE (JUN,3800)UNIT(5),UNIT(4)
C                    ELSE
                     IF (NOOUT.LE.1) WRITE (JUN,2300)
C                    END IF
                     GO TO 120
                  END IF
               END IF
  110       CONTINUE
C           
  120       DO 130 KK = 1, JADEP
               IF (IDPSIN(KK).EQ.1) THEN
                  IF (IOUTSS.NE.2) THEN
C                    IF(OUTOPT.EQ.1.AND.UNITS.EQ.1)THEN
C                    ABBREVIATED ENGLISH UNITS
C                    
C                    IF (NOOUT.LE.1) WRITE (JUN,3600) DEPPT2(KK)*MTF,
C                    1                DPAVLS(KK)*KGMTPA
C                    ELSE
                     IF (NOOUT.LE.1) WRITE (JUN,2100) DEPPT2(KK), 
     1                   DPAVLS(KK)
C                 END IF
C                 
                  END IF
               END IF
  130       CONTINUE
         END IF
      END IF
      RETURN
C     
 1000 FORMAT (//'II.  ON SITE EFFECTS  ON SITE EFFECTS',
     1    '  ON SITE EFFECTS',/,5X,(3(15('-'),2X))/)
 1100 FORMAT (/2X,'A.  AREA OF NET SOIL LOSS')
C2700 FORMAT (/6X,'** SOIL LOSS (AVG. OF NET DETACHMENT',' AREAS) = ',F8
C     1    .1,A,' **')
 1200 FORMAT (/6X,'** SOIL LOSS (AVG. OF NET DETACHMENT',' AREAS) = ',F8
     1    .3,' KG/M2 **')
C2800 FORMAT (6X,'** MAXIMUM SOIL LOSS  = ',F8.1,A,' AT ',F7.1,
C     1    A,' **'/)
 1300 FORMAT (6X,'** MAXIMUM SOIL LOSS  = ',F8.3,' KG/M2 AT ',F7.2,
     1    ' METERS **'/)
C2900 FORMAT (6X,'** INTERRILL CONTRIBUTION = ',F8.1,A,
C     1    ' FOR OFE #',I2)
 1400 FORMAT (6X,'** INTERRILL CONTRIBUTION = ',F8.3,' KG/M2 ')
 1500 FORMAT (/2X,'B.  AREA OF SOIL DEPOSITION')
C3100 FORMAT (/6X,'** SOIL DEPOSITION (AVG. OF NET DEPOSITION',
C     1    ' AREAS) = ',F9.1,A,' **')
 1600 FORMAT (/6X,'** SOIL DEPOSITION (AVG. OF NET DEPOSITION',
     1    ' AREAS) = ',F9.3,' KG/M2 **')
C3200 FORMAT (6X,'** MAXIMUM SOIL DEPOSITION  = ',F9.1,A,' AT ',F7
C     1    .1,A,' **'/)
 1700 FORMAT (6X,'** MAXIMUM SOIL DEPOSITION  = ',F9.3,' KG/M2 AT ',F7
     1    .2,' METERS **'/)
C3300 FORMAT (F7.1,'-',F7.1,1X,F8.1,2X,F8.1,2X,F9.1,1X,F7.1,2X,F8.1,2X,
C     1    F7.1)
 1800 FORMAT (F7.2,'-',F7.2,1X,F8.3,2X,F8.3,2X,F9.3,1X,F7.2,2X,F8.3,2X,
     1    F7.2)
C3400 FORMAT (/6X,'SINGLE POINT',5X,'SINGLE POINT',/,7X,'SOIL AREA',10X,
C     1    'LOSS',/,9X,A,13X,A,/,6X,30('-'))
 1900 FORMAT (/6X,'SINGLE POINT',5X,'SINGLE POINT',/,7X,'SOIL AREA',10X,
     1    'LOSS',/,9X,'(M)',13X,'(KG/M2/)',/,6X,30('-'))
C3500 FORMAT (/,6X,'AREA OF',4X,'SOIL LOSS',3X,'SOIL LOSS',3X,'MAX',3X,
C     1    'MAX LOSS',3X,'MIN',3X,'MIN LOSS',/,6X,'NET LOSS',6X,'MEAN',6
C     1    X,'STDEV',6X,'LOSS',4X,'POINT',4X,'LOSS',3X,'POINT',/,4X,
C     1    A,5X,A,4X,A,3X,A,2X,A,1X,A,1X,A,/,72('-'))
 2000 FORMAT (/,6X,'AREA OF',4X,'SOIL LOSS',3X,'SOIL LOSS',3X,'MAX',3X,
     1    'MAX LOSS',3X,'MIN',3X,'MIN LOSS',/,6X,'NET LOSS',6X,'MEAN',6
     1    X,'STDEV',6X,'LOSS',4X,'POINT',4X,'LOSS',3X,'POINT',/,8X,
     1    '(M)',7X,'(KG/M2)',5X,'(KG/M2)',2X,'(KG/M2)',4X,'(M)',4X,
     1    '(KG/M2)',2X,'(M)',/,72('-'))
C3600 FORMAT (6X,F7.1,9X,F9.1)
 2100 FORMAT (6X,F7.2,9X,F9.3)
C3700 FORMAT (6X,'AREA OF',4X,'SOIL DEP',3X,'SOIL DEP',5X,'MAX',4X,
C     1    'MAX DEP',3X,'MIN',3X,'MIN DEP',/,6X,'NET DEP',7X,'MEAN',5X,
C     1    'STDEV',7X,'DEP',5X,'POINT',4X,'DEP',4X,'POINT',/,4X,
C     1    A,5X,A,4X,A,3X,A,2X,A,1X,A,1X,A,/,72('-'))
 2200 FORMAT (6X,'AREA OF',4X,'SOIL DEP',3X,'SOIL DEP',5X,'MAX',4X,
     1    'MAX DEP',3X,'MIN',3X,'MIN DEP',/,6X,'NET DEP',7X,'MEAN',5X,
     1    'STDEV',7X,'DEP',5X,'POINT',4X,'DEP',4X,'POINT',/,8X,'(M)',7X,
     1    '(KG/M2)',4X,'(KG/M2)',3X,'(KG/M2)',4X,'(M)',4X,'(KG/M2)',2X,
     1    '(M)',/,72('-'))
C3800 FORMAT (/6X,'SINGLE POINT',5X,'SINGLE POINT',/,7X,'SOIL AREA',10X,
C     1    'DEP',/,9X,A,13X,A,/,6X,30('-'))
 2300 FORMAT (/6X,'SINGLE POINT',5X,'SINGLE POINT',/,7X,'SOIL AREA',10X,
     1    'DEP',/,9X,'(M)',13X,'(KG/M2/)',/,6X,30('-'))
C     
      END
