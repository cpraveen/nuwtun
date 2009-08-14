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
      SUBROUTINE LIMITR (NPTS, DQP, DQM, LIMIT, RKAP, LBEG, LEND, 
     1                   IERROR)
C
C The second order limited correction to the first order interpolation
C is calculated by this routine and stored in DQP for the q(+) term and 
C DQM for the q(-) term.  Delta(q(i)) = q(i+1) - q(i) is stored in DQP(i) 
C on entry.
C
C N.B. The values of DQP(NPTS) and DQM(1) are not defined properly on exit 
C as this leads to better vectorization and they are unnecessary.
C
C This routine is organized as follows:
C     0. Set constants for the interpolation
C     1. Venkatakrishnan's limiter
C     2. MinMod limiter
C     3. Smooth limiter
C     4. Unlimited scheme
C     5. Error
C
C NPTS           : Length of array to limit successive values of delta(q)
C DQP,DQM        : Storage for delta q = q(i) - q(i-1) on entry.
C                  Storage of the higher order correction to the 
C                  first order scheme on exit
C LIMIT          : Limiter to use for MUSCL interpolation
C                  ILNONE [0] => Unlimited Scheme
C                  ILSMTH [1] => Smooth Limiter
C                  ILMNMD [2] => Minmod Limiter
C                  ILVNKT [3] => Smooth limiter (fixed at kappa = 1/3) 
C                           from Venkatakrishnan, AIAA J. Vol. 29 No. 7
C                           also Thomas, J. L., 5th Copper Mtn. Multigrid
C RKAP           : Kappa to use in MUSCL interpolation
C LBEG,LEND      : Starting,stopping equation index for limiting
C IERROR         : Error code
C
      include '../header/common.h'
C
      DIMENSION DQP   (NPTS,NQ),
     1          DQM   (NPTS,NQ)
C
C     Initialize error
C
      IERROR = INO
C
C 0.  Set constants for the interpolation.
C
      OMKAP = 1.E0 - RKAP
      OPKAP = 1.E0 + RKAP
C
C 1.  Hemker-Koren Limiter
C     This is wrongly called Venkat's Smooth Limiter in ISAAC
C
      IF (LIMIT .EQ. ILVNKT) THEN
         DO L = LBEG, LEND
            DO I = 2, NPTS
               CALL SMILIM(DQP(I-1,L), DQP(I,L), DQM(I,L))
            END DO
         END DO
C
C 2.  MinMod Limiter
C
      ELSE IF (LIMIT .EQ. ILMNMD) THEN
         B = (3.E0 - RKAP) / (1.E0 - RKAP)
         DO L = LBEG, LEND
            DO I = 2, NPTS
               CALL RMMLIM(DQP(I-1,L), DQP(I,L), DQM(I,L), 
     1                     OMKAP, OPKAP, B)
            END DO
         END DO
C
C 3.  Smooth Limiter
C
      ELSE IF (LIMIT .EQ. ILSMTH) THEN
         DO L = LBEG, LEND
            DO I = 2, NPTS
               CALL SMTLIM(DQP(I-1,L), DQP(I,L), DQM(I,L), RKAP)
            END DO
         END DO
C
C 4.  Unlimited Scheme
C
      ELSE IF (LIMIT .EQ. ILNONE) THEN
         DO L = LBEG, LEND
            DO I = 2, NPTS
               CALL RNOLIM(DQP(I-1,L), DQP(I,L), DQM(I,L), OMKAP, OPKAP)
            END DO
         END DO
C
C 5.  Error
C
      ELSE
         WRITE (IOUT,1000) LIMIT
         IERROR = IYES
         STOP
      ENDIF
C
C Finished with limiting routine.
C
      RETURN
 1000 FORMAT(' ','LIMITR  : Invalid value of LIMIT = ',I10,
     1      /' ','Run Aborting.')
      END
