C
C     This file is part of NuWTun, see <http://nuwtun.berlios.de>, and was
C     originally taken from ISAAC Version 4.2, release date October 2001. 
C     This file may have been modified; for a list of changes, see the 
C     changes.txt file in the docs directory and the subversion log.
C
C     Portions Copyright (C) 2001 Joseph H. Morrison
C
C     This code is part of ISAAC.
C
C     This program is distributed under the terms of the ISAAC Public Source
C     License. This program is distributed WITHOUT ANY WARRANTY; without
C     even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
C     PURPOSE. 
C
C     You should have received a copy of the ISAAC Public Source License
C     with this program. If you did not, you may get a copy of the license
C     at <http://isaac-cfd.sourceforge.net>
C
      SUBROUTINE OUTPUT (IDIM, JDIM, KDIM, Q, RMU, RMUT, YPLUS, ITURB, 
     1                   RC, IBLOCK, NPRINT, IPRINT)
C
C     Routine to write output
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     Q              : Primitive variables at cell centers
C     RMU            : Molecular viscosity at cell centers
C     RMUT           : Turbulent eddy viscosity at cell centers
C     YPLUS          : Y+ at cell centers
C     ITURB          : Turbulence model
C     RC             : Calculated position of grid at cell centers
C     IBLOCK         : Block number to identify output
C     NPRINT         : Number of print statements for current block
C     IPRINT         : Print control vector for current block
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          RMU   (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     2          RMUT  (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     3          YPLUS (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     4          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      DIMENSION IPRINT(10,NPRINT)
C
      DIMENSION ISTRT(3), IEND(3), INC(3)
C
      CHARACTER*10 TITLE(3)
      DATA TITLE / 'J, K =    ',
     1             'K, I =    ',
     2             'I, J =    ' /
C
C     Loop over all print statements
C
      DO 800 IPR = 1, NPRINT
         IDIR     = IPRINT(1,IPR)
         ISTRT(1) = IPRINT(2,IPR)
         ISTRT(2) = IPRINT(3,IPR)
         ISTRT(3) = IPRINT(4,IPR)
         IEND(1)  = IPRINT(5,IPR)
         IEND(2)  = IPRINT(6,IPR)
         IEND(3)  = IPRINT(7,IPR)
         INC(1)   = IPRINT(8,IPR)
         INC(2)   = IPRINT(9,IPR)
         INC(3)   = IPRINT(10,IPR)
         DO 530 I3 = ISTRT(ICY(IDIR,3)), IEND(ICY(IDIR,3)), 
     1               INC(ICY(IDIR,3))
            DO 520 I2 = ISTRT(ICY(IDIR,2)), IEND(ICY(IDIR,2)), 
     1                  INC(ICY(IDIR,2))
C
C     Output Header
C
               IF (ITURB .EQ. ITLMNR) THEN
                  WRITE (IOUT,1003) IBLOCK, TITLE(IDIR), I2, I3
               ELSE IF (ITURB .EQ. ITBLMX) THEN
                  WRITE (IOUT,1005) IBLOCK, TITLE(IDIR), I2, I3
               ELSE IF (ITURB .EQ. ITKE) THEN
                  WRITE (IOUT,1010) IBLOCK, TITLE(IDIR), I2, I3
               ELSE IF (ITURB .EQ. ITKW) THEN
                  WRITE (IOUT,1015) IBLOCK, TITLE(IDIR), I2, I3
               ELSE IF (ITURB .EQ. ITRS) THEN
                  WRITE (IOUT,1025) IBLOCK, TITLE(IDIR), I2, I3
               ELSE IF (ITURB .EQ. ITLES) THEN
                  WRITE (IOUT,1030) IBLOCK, TITLE(IDIR), I2, I3
               ELSE
                  WRITE (IOUT,1000) IBLOCK, TITLE(IDIR), I2, I3
               ENDIF
C
C     Inner Loop to Print out
C
               DO 510 I1 = ISTRT(ICY(IDIR,1)), IEND(ICY(IDIR,1)), 
     1                     INC(ICY(IDIR,1))
                  II = I1*IKD(ICY(IDIR,1),1) + I2*IKD(ICY(IDIR,2),1) +
     1                 I3*IKD(ICY(IDIR,3),1)
                  JJ = I1*IKD(ICY(IDIR,1),2) + I2*IKD(ICY(IDIR,2),2) +
     1                 I3*IKD(ICY(IDIR,3),2)
                  KK = I1*IKD(ICY(IDIR,1),3) + I2*IKD(ICY(IDIR,2),3) +
     1                 I3*IKD(ICY(IDIR,3),3)
                  PPINF = GAMMA * Q(JJ,KK,II,5)
                  TTINF = GAMMA * Q(JJ,KK,II,5) / Q(JJ,KK,II,1)
                  IF (ITURB .EQ. ITLMNR) THEN
                     WRITE (IOUT,1103) II, JJ, KK,
     1                  RC(JJ,KK,II,1), RC(JJ,KK,II,2), RC(JJ,KK,II,3),
     2                  Q(JJ,KK,II,1),  Q(JJ,KK,II,2),  Q(JJ,KK,II,3), 
     3                  Q(JJ,KK,II,4),  PPINF, TTINF, RMU(JJ,KK,II)
                  ELSE IF (ITURB .EQ. ITBLMX) THEN
                     WRITE (IOUT,1105) II, JJ, KK,
     1                  RC(JJ,KK,II,1), RC(JJ,KK,II,2), RC(JJ,KK,II,3),
     2                  Q(JJ,KK,II,1),  Q(JJ,KK,II,2),  Q(JJ,KK,II,3), 
     3                  Q(JJ,KK,II,4),  PPINF, TTINF, RMU(JJ,KK,II),
     4                  RMUT(JJ,KK,II)
                  ELSE IF (ITURB. EQ. ITKE) THEN
                     RHO   = Q(JJ,KK,II,1)
                     P     = Q(JJ,KK,II,5)
                     TKE   = Q(JJ,KK,II,6)
                     EPSLN = Q(JJ,KK,II,7)
                     ASQ   = GAMMA * P / RHO
                     RMT   = SQRT (2.E0 * TKE / ASQ)
                     RET   = RHO*TKE*TKE/(EPSLN*RMU(JJ,KK,II)+RSMALL) *
     1                       RE / FSMACH
                     WRITE (IOUT,1110) II, JJ, KK,
     1                  RC(JJ,KK,II,1), RC(JJ,KK,II,2), RC(JJ,KK,II,3),
     2                  Q(JJ,KK,II,1),  Q(JJ,KK,II,2),  Q(JJ,KK,II,3), 
     3                  Q(JJ,KK,II,4),  PPINF, TTINF, RMU(JJ,KK,II),
     4                  RMUT(JJ,KK,II), Q(JJ,KK,II,6),  Q(JJ,KK,II,7), 
     5                  YPLUS(JJ,KK,II), RMT, RET
                  ELSE IF (ITURB .EQ. ITKW) THEN
                     RHO   = Q(JJ,KK,II,1)
                     P     = Q(JJ,KK,II,5)
                     TKE   = Q(JJ,KK,II,6)
                     EPSLN = TKE * Q(JJ,KK,II,7) 
                     ASQ   = GAMMA * P / RHO
                     RMT   = SQRT (2.E0 * TKE / ASQ)
                     RET   = RHO*TKE*TKE/(EPSLN*RMU(JJ,KK,II)+RSMALL) *
     1                       RE / FSMACH
                     WRITE (IOUT,1110) II, JJ, KK,
     1                  RC(JJ,KK,II,1), RC(JJ,KK,II,2), RC(JJ,KK,II,3),
     2                  Q(JJ,KK,II,1),  Q(JJ,KK,II,2),  Q(JJ,KK,II,3), 
     3                  Q(JJ,KK,II,4),  PPINF, TTINF, RMU(JJ,KK,II),
     4                  RMUT(JJ,KK,II), Q(JJ,KK,II,6),  Q(JJ,KK,II,7), 
     5                  YPLUS(JJ,KK,II), RMT, RET
                  ELSE IF (ITURB. EQ. ITRS) THEN
                     WRITE (IOUT,1125) II, JJ, KK,
     1                  RC(JJ,KK,II,1), RC(JJ,KK,II,2), RC(JJ,KK,II,3),
     2                  Q(JJ,KK,II,1),  Q(JJ,KK,II,2),  Q(JJ,KK,II,3), 
     3                  Q(JJ,KK,II,4),  PPINF, TTINF, RMU(JJ,KK,II),
     4                  RMUT(JJ,KK,II), Q(JJ,KK,II,6),  
     5                  Q(JJ,KK,II,7),  Q(JJ,KK,II,8),  Q(JJ,KK,II,9),  
     6                  Q(JJ,KK,II,10), Q(JJ,KK,II,11), Q(JJ,KK,II,12), 
     7                  YPLUS(JJ,KK,II)
                  ELSE IF (ITURB .EQ. ITLES) THEN
                     WRITE (IOUT,1130) II, JJ, KK,
     1                  RC(JJ,KK,II,1), RC(JJ,KK,II,2), RC(JJ,KK,II,3),
     2                  Q(JJ,KK,II,1),  Q(JJ,KK,II,2),  Q(JJ,KK,II,3), 
     3                  Q(JJ,KK,II,4),  PPINF, TTINF, RMU(JJ,KK,II),
     4                  RMUT(JJ,KK,II)
                  ELSE
                     WRITE (IOUT,1100) II, JJ, KK,
     1                  RC(JJ,KK,II,1), RC(JJ,KK,II,2), RC(JJ,KK,II,3),
     2                  Q(JJ,KK,II,1),  Q(JJ,KK,II,2),  Q(JJ,KK,II,3), 
     3                  Q(JJ,KK,II,4),  PPINF, TTINF
                  ENDIF
  510          CONTINUE
  520       CONTINUE
  530    CONTINUE
  800 CONTINUE
C
      RETURN
C
 1000 FORMAT (//' ','BLOCK = ',I5,10X,A10,2I5/
     1        ' ','   I   J   K',9X,'X',16X,'Y',16X,'Z',
     2            15X,'RHO',15X,'U',16X,'V',
     3            16X,'W',13X,'P/Pinf',12X,'T/Tinf')
 1003 FORMAT (//' ','BLOCK = ',I5,10X,A10,2I5/
     1        ' ','   I   J   K',7X,'X',12X,'Y',12X,'Z',
     2            11X,'RHO',11X,'U',12X,'V',
     3            12X,'W',9X,'P/Pinf',8X,'T/Tinf',9X,'RMU')
 1005 FORMAT (//' ','BLOCK = ',I5,10X,A10,2I5/
     1        ' ','   I   J   K',7X,'X',12X,'Y',12X,'Z',
     2            11X,'RHO',11X,'U',12X,'V',
     3            12X,'W',9X,'P/Pinf',8X,'T/Tinf',9X,'RMU',
     4            8X,'RMUT')
 1010 FORMAT (//' ','BLOCK = ',I5,10X,A10,2I5/
     1        ' ','   I   J   K',7X,'X',12X,'Y',12X,'Z',
     2            11X,'RHO',11X,'U',12X,'V',
     3            12X,'W',9X,'P/Pinf',8X,'T/Tinf',9X,'RMU',
     4            8X,'RMUT',9X,'k',10X,'EPSILON',8X,'Yplus',
     5            7X,'Mach_t',7X,'Re_t')
 1015 FORMAT (//' ','BLOCK = ',I5,10X,A10,2I5/
     1        ' ','   I   J   K',7X,'X',12X,'Y',12X,'Z',
     2            11X,'RHO',11X,'U',12X,'V',
     3            12X,'W',9X,'P/Pinf',8X,'T/Tinf',9X,'RMU',
     4            8X,'RMUT',9X,'k',11X,'OMEGA',8X,'Yplus',
     5            7X,'Mach_t',7X,'Re_t')
 1025 FORMAT (//' ','BLOCK = ',I5,10X,A10,2I5/
     1        ' ','   I   J   K',7X,'X',12X,'Y',12X,'Z',
     2            11X,'RHO',11X,'U',12X,'V',
     3            12X,'W',9X,'P/Pinf',8X,'T/Tinf',9X,'RMU',
     4            8X,'RMUT',9X,'TAUXX',8X,'TAUYY',8X,'TAUZZ',8X,
     5            'TAUXY',8X,'TAUXZ',8X,'TAUYZ',7X,'EPSILON',8X,'Yplus')
 1030 FORMAT (//' ','BLOCK = ',I5,10X,A10,2I5/
     1        ' ','   I   J   K',7X,'X',12X,'Y',12X,'Z',
     2            11X,'RHO',11X,'U',12X,'V',
     3            12X,'W',9X,'P/Pinf',8X,'T/Tinf',9X,'RMU',
     4            8X,'RMUT')
 1100 FORMAT (' ',3I4,10(1PE17.9))
 1103 FORMAT (' ',3I4, 3(1PE13.5), 8(1PE13.5))
 1105 FORMAT (' ',3I4, 3(1PE13.5), 9(1PE13.5))
 1110 FORMAT (' ',3I4, 3(1PE13.5),14(1PE13.5))
 1125 FORMAT (' ',3I4, 3(1PE13.5),17(1PE13.5))
 1130 FORMAT (' ',3I4, 3(1PE13.5), 9(1PE13.5))
C
      END
