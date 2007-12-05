      SUBROUTINE METRIC (IBLOCK, ILVL, IDIM, JDIM, KDIM, R, S, VOL, 
     1                   IERRCD)
C
C Routine to calculate metrics of the grid 
C
C This routine is organized as follows:
C     1. Calculate and normalize I metrics
C     2. Calculate and normalize J metrics
C     3. Calculate and normalize K metrics
C     4. Calculate cell volumes
C     5. Check for zero or negative metrics
C     6. Check for zero or negative volumes
C
C IBLOCK         : Current block
C ILVL           : Current grid level
C IDIM,JDIM,KDIM : Dimensions of current block
C R              : Position vector for the grid
C S              : Metrics
C VOL            : Cell volume
C IERRCD         : Error code
C
      include '../header/common.h'
C
      DIMENSION R     (0:JDIM+1,0:KDIM+1,0:IDIM+1,3),
     1          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     2          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
C 0.  Control printout based on grid level
C
      IF (ILVL .EQ. 1) THEN
         MAXPRT = 50
      ELSE
         MAXPRT =  0
      ENDIF
C
C 1.  Calculate xsi(x), xsi(y), xsi(z). 
C     Calculate S(J,K,I,4,1) = |del(xsi)|.
C
      DO 200 J = 1, JDIM+1
         DO 180 K = 1, KDIM+1
            DO 160 I = 0, IDIM+1
               S(J,K,I,1,1) = 0.5E0 * ( (R(J,K,I,2)-R(J-1,K-1,I,2))*
     1                                  (R(J-1,K,I,3)-R(J,K-1,I,3)) -
     2                                  (R(J,K,I,3)-R(J-1,K-1,I,3))*
     3                                  (R(J-1,K,I,2)-R(J,K-1,I,2)) )
               S(J,K,I,2,1) = 0.5E0 * ( (R(J,K,I,3)-R(J-1,K-1,I,3))*
     1                                  (R(J-1,K,I,1)-R(J,K-1,I,1)) -
     2                                  (R(J,K,I,1)-R(J-1,K-1,I,1))*
     3                                  (R(J-1,K,I,3)-R(J,K-1,I,3)) )
               S(J,K,I,3,1) = 0.5E0 * ( (R(J,K,I,1)-R(J-1,K-1,I,1))*
     1                                  (R(J-1,K,I,2)-R(J,K-1,I,2)) -
     2                                  (R(J,K,I,2)-R(J-1,K-1,I,2))*
     3                                  (R(J-1,K,I,1)-R(J,K-1,I,1)) )
               S(J,K,I,4,1) = SQRT ( S(J,K,I,1,1)*S(J,K,I,1,1) + 
     1                               S(J,K,I,2,1)*S(J,K,I,2,1) + 
     2                               S(J,K,I,3,1)*S(J,K,I,3,1) )
 160        CONTINUE
 180     CONTINUE
 200  CONTINUE
C
C     Normalize and store S(J,K,I,1,1) = xsi(x) / |del(xsi)|.
C     Normalize and store S(J,K,I,2,1) = xsi(y) / |del(xsi)|.
C     Normalize and store S(J,K,I,3,1) = xsi(z) / |del(xsi)|.
C
      DO 230 J = 1, JDIM+1
         DO 220 K = 1, KDIM+1
            DO 210 I = 0, IDIM+1
               IF (S(J,K,I,4,1) .GT. 0.E0) THEN
                  S(J,K,I,1,1) = S(J,K,I,1,1) / S(J,K,I,4,1)
                  S(J,K,I,2,1) = S(J,K,I,2,1) / S(J,K,I,4,1)
                  S(J,K,I,3,1) = S(J,K,I,3,1) / S(J,K,I,4,1)
               ELSE
                  S(J,K,I,1,1) = SHATX
                  S(J,K,I,2,1) = SHATY
                  S(J,K,I,3,1) = SHATZ
               ENDIF
  210       CONTINUE
  220    CONTINUE
  230 CONTINUE
C
C 2.  Calculate eta(x), eta(y), eta(z). 
C     Calculate S(J,K,I,4,2) = |del(eta)|.
C
      DO 400 I = 1, IDIM+1
         DO 380 K = 1, KDIM+1
            DO 360 J = 0, JDIM+1
               S(J,K,I,1,2) = 0.5E0 * ( (R(J,K,I,3)-R(J,K-1,I-1,3))* 
     1                                  (R(J,K,I-1,2)-R(J,K-1,I,2)) -
     2                                  (R(J,K,I,2)-R(J,K-1,I-1,2))*
     3                                  (R(J,K,I-1,3)-R(J,K-1,I,3)) )
               S(J,K,I,2,2) = 0.5E0 * ( (R(J,K,I,1)-R(J,K-1,I-1,1))* 
     1                                  (R(J,K,I-1,3)-R(J,K-1,I,3)) -
     2                                  (R(J,K,I,3)-R(J,K-1,I-1,3))*
     3                                  (R(J,K,I-1,1)-R(J,K-1,I,1)) )
               S(J,K,I,3,2) = 0.5E0 * ( (R(J,K,I,2)-R(J,K-1,I-1,2))* 
     1                                  (R(J,K,I-1,1)-R(J,K-1,I,1)) -
     2                                  (R(J,K,I,1)-R(J,K-1,I-1,1))*
     3                                  (R(J,K,I-1,2)-R(J,K-1,I,2)) )
               S(J,K,I,4,2) = SQRT ( S(J,K,I,1,2)*S(J,K,I,1,2) + 
     1                               S(J,K,I,2,2)*S(J,K,I,2,2) + 
     2                               S(J,K,I,3,2)*S(J,K,I,3,2) )
 360        CONTINUE
 380     CONTINUE
 400  CONTINUE
C
C     Normalize and store S(J,K,I,1,2) = eta(x) / |del(eta)|.
C     Normalize and store S(J,K,I,2,2) = eta(y) / |del(eta)|.
C     Normalize and store S(J,K,I,3,2) = eta(z) / |del(eta)|.
C
      DO 430 I = 1, IDIM+1
         DO 420 K = 1, KDIM+1
            DO 410 J = 0, JDIM+1
               IF (S(J,K,I,4,2) .GT. 0.E0) THEN
                  S(J,K,I,1,2) = S(J,K,I,1,2) / S(J,K,I,4,2)
                  S(J,K,I,2,2) = S(J,K,I,2,2) / S(J,K,I,4,2)
                  S(J,K,I,3,2) = S(J,K,I,3,2) / S(J,K,I,4,2)
               ELSE
                  S(J,K,I,1,2) = SHATX
                  S(J,K,I,2,2) = SHATY
                  S(J,K,I,3,2) = SHATZ
               ENDIF
  410       CONTINUE
  420    CONTINUE
  430 CONTINUE
C
C 3.  Calculate zeta(x), zeta(y), zeta(z). 
C     Calculate S(J,K,I,4,3) = |del(zeta)|.
C
      DO 600 I = 1, IDIM+1
         DO 580 J = 1, JDIM+1
            DO 560 K = 0, KDIM+1
               S(J,K,I,1,3) = 0.5E0 * ( (R(J,K,I,2)-R(J-1,K,I-1,2))*
     1                                  (R(J,K,I-1,3)-R(J-1,K,I,3)) -
     2                                  (R(J,K,I,3)-R(J-1,K,I-1,3))*
     3                                  (R(J,K,I-1,2)-R(J-1,K,I,2)) )
               S(J,K,I,2,3) = 0.5E0 * ( (R(J,K,I,3)-R(J-1,K,I-1,3))*
     1                                  (R(J,K,I-1,1)-R(J-1,K,I,1)) -
     2                                  (R(J,K,I,1)-R(J-1,K,I-1,1))*
     3                                  (R(J,K,I-1,3)-R(J-1,K,I,3)) )
               S(J,K,I,3,3) = 0.5E0 * ( (R(J,K,I,1)-R(J-1,K,I-1,1))*
     1                                  (R(J,K,I-1,2)-R(J-1,K,I,2)) -
     2                                  (R(J,K,I,2)-R(J-1,K,I-1,2))*
     3                                  (R(J,K,I-1,1)-R(J-1,K,I,1)) )
               S(J,K,I,4,3) = SQRT ( S(J,K,I,1,3)*S(J,K,I,1,3) + 
     1                               S(J,K,I,2,3)*S(J,K,I,2,3) + 
     2                               S(J,K,I,3,3)*S(J,K,I,3,3) )
 560        CONTINUE
 580     CONTINUE
 600  CONTINUE
C
C     Normalize and store S(J,K,I,1,3) = zeta(x) / |del(zeta)|.
C     Normalize and store S(J,K,I,2,3) = zeta(y) / |del(zeta)|.
C     Normalize and store S(J,K,I,3,3) = zeta(z) / |del(zeta)|.
C
      DO 630 I = 1, IDIM+1
         DO 620 J = 1, JDIM+1
            DO 610 K = 0, KDIM+1
               IF (S(J,K,I,4,3) .GT. 0.E0) THEN
                  S(J,K,I,1,3) = S(J,K,I,1,3) / S(J,K,I,4,3)
                  S(J,K,I,2,3) = S(J,K,I,2,3) / S(J,K,I,4,3)
                  S(J,K,I,3,3) = S(J,K,I,3,3) / S(J,K,I,4,3)
               ELSE
                  S(J,K,I,1,3) = SHATX
                  S(J,K,I,2,3) = SHATY
                  S(J,K,I,3,3) = SHATZ
               ENDIF
  610       CONTINUE
  620    CONTINUE
  630 CONTINUE
C
C 4.  Calculate the cell volume from eqn. 2c of AIAA J. Vol. 21
C     No. 6 pp. 917-918 by Kordulla and Vinokur
C
      DO 800 K = 1, KDIM+1
         DO 780 J = 1, JDIM+1
            DO 760 I = 1, IDIM+1
               DX = ( S(J,K,I-1,1,1)*S(J,K,I-1,4,1) +
     1                S(J-1,K,I,1,2)*S(J-1,K,I,4,2) +
     2                S(J,K-1,I,1,3)*S(J,K-1,I,4,3) ) *
     3              ( R(J,K,I,1) - R(J-1,K-1,I-1,1) )
               DY = ( S(J,K,I-1,2,1)*S(J,K,I-1,4,1) +
     1                S(J-1,K,I,2,2)*S(J-1,K,I,4,2) +
     2                S(J,K-1,I,2,3)*S(J,K-1,I,4,3) ) *
     3              ( R(J,K,I,2) - R(J-1,K-1,I-1,2) )
               DZ = ( S(J,K,I-1,3,1)*S(J,K,I-1,4,1) +
     1                S(J-1,K,I,3,2)*S(J-1,K,I,4,2) +
     2                S(J,K-1,I,3,3)*S(J,K-1,I,4,3) ) *
     3              ( R(J,K,I,3) - R(J-1,K-1,I-1,3) ) 
               VOL(J,K,I) = (DX + DY + DZ) / 3.E0
 760        CONTINUE
 780     CONTINUE
 800  CONTINUE
C
C 5.  Check for negative or zero areas
C
C I faces
C
      ISNEG  = 0
      DO 830 J = 2, JDIM
         DO 820 K = 2, KDIM
            DO 810 I = 0, IDIM+1
               IF (S(J,K,I,4,1) .LE. 0.E0) THEN
                  ISNEG = ISNEG + 1
                  IF (ISNEG .LE. MAXPRT) THEN
                     WRITE (IOUT,1100) 'I', S(J,K,I,4,1), IBLOCK, ILVL,
     1                                  I, J, K
                  ENDIF
               ENDIF
  810       CONTINUE
  820    CONTINUE
  830 CONTINUE
      IF (ISNEG .GT. 0) THEN
         WRITE (IOUT,1110) ISNEG, 'I', IBLOCK, ILVL
      ENDIF
C
C J faces
C
      ISNEG  = 0
      DO 860 I = 2, IDIM
         DO 850 K = 2, KDIM
            DO 840 J = 0, JDIM+1
               IF (S(J,K,I,4,2) .LE. 0.E0) THEN
                  ISNEG = ISNEG + 1
                  IF (ISNEG .LE. MAXPRT) THEN
                     WRITE (IOUT,1100) 'J', S(J,K,I,4,2), IBLOCK, ILVL,
     1                                  I, J, K
                  ENDIF
               ENDIF
  840       CONTINUE
  850    CONTINUE
  860 CONTINUE
      IF (ISNEG .GT. 0) THEN
         WRITE (IOUT,1110) ISNEG, 'J', IBLOCK, ILVL
      ENDIF
C
C K faces
C
      ISNEG  = 0
      DO 890 I = 2, IDIM
         DO 880 J = 2, JDIM
            DO 870 K = 0, KDIM+1
               IF (S(J,K,I,4,3) .LE. 0.E0) THEN
                  ISNEG = ISNEG + 1
                  IF (ISNEG .LE. MAXPRT) THEN
                     WRITE (IOUT,1100) 'K', S(J,K,I,4,3), IBLOCK, ILVL,
     1                                  I, J, K
                  ENDIF
               ENDIF
  870       CONTINUE
  880    CONTINUE
  890 CONTINUE
      IF (ISNEG .GT. 0) THEN
         WRITE (IOUT,1110) ISNEG, 'K', IBLOCK, ILVL
      ENDIF
C
C 6. Check for negative or zero volumes - left-handed coordinate
C    systems or collapsed cells
C
      IERROR = 0
      DO 930 I = 2, IDIM
         DO 920 K = 2, KDIM
            DO 910 J = 2, JDIM
               IF (VOL(J,K,I) .LE. 0.E0) THEN
                  IERROR = IERROR + 1
                  IF (IERROR .LE. MAXPRT) THEN
                     WRITE (IOUT,1010) IBLOCK, ILVL, I, J, K, VOL(J,K,I)
                  ENDIF
               ENDIF
  910       CONTINUE
  920    CONTINUE
  930 CONTINUE
      IF (IERROR .GT. 0) THEN
         WRITE (IOUT,1015) IERROR, IBLOCK, ILVL, IDIM, JDIM, KDIM
         IERRCD = IYES
      ENDIF
C
C Finished with the metric calculations.
C
      RETURN
C
 1010 FORMAT (' ','METRIC: ERROR-> Volume <= 0 at BLOCK,LEVEL,I,J,K = ',
     1            I3,I3,3I4,1PE10.3)
 1015 FORMAT (' ','METRIC: ERROR-> ',I9,' negative or zero volumes ',
     1            'in grid ',
     2       /' ',16X,' BLOCK: LEVEL: IDIM, JDIM, KDIM = ',
     3            I3,':',I2,':',3I5/)
 1100 FORMAT (' ','METRIC: ',A,' face non-positive ',1PE9.2,
     1            ' at BLOCK:LEVEL:I,J,K ',I3,':',I2,':',3I4)
 1110 FORMAT (' ','METRIC: ',I9,' non-positive ',A,' face areas ',
     1            'for BLOCK ',I3,' LEVEL ',I2/)
C
      END
