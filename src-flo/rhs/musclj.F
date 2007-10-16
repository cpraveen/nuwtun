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
      SUBROUTINE MUSCLJ (IDIM, JDIM, KDIM, I,
     1                   Q, DQP, DQM, QP, QM,
     2                   IBCDAT, NBCS, IORDER, LIMIT, RKAP)
C
C Routine to do MUSCL interpolation in the J direction. This routine 
C interpolates for Q(+,-) on the J faces for j = 1, jdim by k = 2, kdim 
C on the current I plane. This allows for greater vectorization and also 
C minimizes intermediate storage to a cross section.  
C
C This routine is organized as follows:
C     0. Set up beginning/ending indices for first and second order sections
C     1. Calculate Q(+,-)
C        a. First order
C        b. Second order
C           i.   Calculate delta(Q)
C           ii.  Limit delta(Q)
C           iii. Construct Q(+,-)
C     2. Correct wall values of Q(+,-)
C
C IDIM,JDIM,KDIM : Dimensions of current block
C I              : I Plane to do interpolation on
C Q              : Primitive variables at cell centers
C DQP,DQM        : Storage for delta q = q(j) - q(j-1) and
C                  for the higher order interpolant to be added
C                  to the first order scheme to create QP,QM
C QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C IBCDAT         : Data controlling boundary conditions:
C                  IBCDAT(1,ibc) = Boundary Condition Type
C                         2      = Direction (1=i, 2=j, 3=k)
C                         3      = j(1 or jdim)
C                         4      = istrt
C                         5      = iend
C                         6      = jstrt
C                         7      = jend
C                         8      = kstrt
C                         9      = kend
C                         10     = profile identifier
C NBCS           : Number of boundary conditions for current block
C IORDER         : Order of MUSCL interpolation for I direction 
C                  IORDER = 1 first order
C                         = 2 second order (or 'third' order)
C LIMIT          : Limiter to use for MUSCL interpolation
C RKAP           : Kappa to use in MUSCL interpolation
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION DQP   (0:JDIM+1,2:KDIM,NQ),
     1          DQM   (0:JDIM+1,2:KDIM,NQ),
     2          QP    (JDIM,2:KDIM,NQ), 
     3          QM    (JDIM,2:KDIM,NQ) 
C
      DIMENSION IBCDAT(10,NBCS)
C
C Local Storage
C
      LOGICAL   DOFRST,
     1          DOSCND
C
C 0.  Set up beginning/ending indices for first and second order sections
C
      LFOBEG = 0
      LFOEND = -1
      LSOBEG = 1
      LSOEND = NQ
      DOFRST = .FALSE.
      DOSCND = .TRUE.
C
      IF (IORDER .EQ. 1) THEN
         LFOBEG = 1
         LFOEND = NQ
         LSOBEG = 0
         LSOEND = -1
         DOFRST = .TRUE.
         DOSCND = .FALSE.
      ELSE IF (FOTURB) THEN
         LFOBEG = 6
         LFOEND = NQ
         LSOBEG = 1
         LSOEND = 5
         DOFRST = .TRUE.
         DOSCND = .TRUE.
      END IF
C
C 1.  Calculate Q(+,-)
C     Loop over j,k formulating q(extrapolated). J is stored first so that 
C     successive values of j (j=1,2,3...) are stored together for limiting.
C
C 1-a. First Order
C
      IF (DOFRST) THEN
         DO 130 L = LFOBEG, LFOEND
            DO 120 K = 2, KDIM
               DO 110 J = 1, JDIM
                  QP(J,K,L) = Q(J+1,K,I,L) 
                  QM(J,K,L) = Q(J,K,I,L)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
      ENDIF
C
C 1-b. Second Order or Third Order
C 1-b-i. Calculate delta(Q)
C        Calculate delta(q(j)) = q(j+1) - q(j) to formulate q(extrapolated) 
C        using limiting.
C
      IF (DOSCND) THEN
         DO 230 L = LSOBEG, LSOEND
            DO 220 K = 2, KDIM
               DO 200 J = 0, JDIM + 1
                  DQP(J,K,L) = Q(J+1,K,I,L) - Q(J,K,I,L)
  200          CONTINUE
  220       CONTINUE
  230    CONTINUE
C
C     Set delta(q) to 2*delta(q) to adjust for the wall boundary condition 
C     being on the wall rather than at a cell center.
C     Set delta(q) = 0 on wall for wall function to adjust for grid spacing.
C
         DO 300 IBC = 1, NBCS
            IBCTYP = IBCDAT(1,IBC)
            IDIR   = IBCDAT(2,IBC)
            IF ( IDIR .EQ. 2 .AND.
     1           ( IBCTYP .EQ. ITAN  .OR. IBCTYP .EQ. IWALL .OR.
     2             IBCTYP .EQ. IWALFN ) ) THEN
               INROUT = IBCDAT(3,IBC)
               ISTRT  = IBCDAT(4,IBC)
               IEND   = IBCDAT(5,IBC)
               KSTRT  = IBCDAT(8,IBC)
               KEND   = IBCDAT(9,IBC)
               IF (INROUT .EQ. 1) THEN
                  JBCPT = 1
               ELSE
                  JBCPT = JDIM
               ENDIF
               IF (IBCTYP .EQ. IWALFN) THEN
                  FACTR = 0.E0
               ELSE
                  FACTR  = 2.E0
               ENDIF
               IF ( I .GE. ISTRT .AND. I .LE. IEND ) THEN
                  DO 290 L = LSOBEG, LSOEND
                     DO 280 K = KSTRT, KEND
                        DQP(JBCPT,K,L) = FACTR * DQP(JBCPT,K,L)
  280                CONTINUE
  290             CONTINUE
               ENDIF
            ENDIF
  300    CONTINUE
C
C 1-b-ii.  Limit delta(Q)
C     Limit from j = 0 to j = jdim+1 by k = 2 to k = kdim times.
C     DQP is used to pass in the values of delta q.
C     The values of the limited deltas are stored in DQP and DQM.
C     These values not valid at the endpoints of each j line, but are
C     calculated here to maximize vectorization on a jxk plane.
C
         NPTS = ( JDIM + 2 ) * ( KDIM - 1 )
         CALL LIMITR (NPTS, DQP, DQM, LIMIT, RKAP, LSOBEG, LSOEND, 
     1                IERROR)
C
C     Calculate q(extrapolated,j+1/2) from q and the 
C     limited values of delta(q).
C
         DO 380 L = LSOBEG, LSOEND
            DO 370 K = 2, KDIM
               DO 360 J = 1, JDIM
                  QP(J,K,L) = Q(J+1,K,I,L) + DQP(J,K,L)
                  QM(J,K,L) = Q(J,K,I,L)   + DQM(J,K,L)
  360          CONTINUE
  370       CONTINUE
  380    CONTINUE
      ENDIF
C
C 2.  Correct wall values of Q(+,-)
C     Fix boundary values of q(extrapolated) due to boundary conditions
C     for either the first order or higher order scheme.
C     WALL boundary condition.  Reset QP(wall)=QM(wall)=Q(wall)
C
      DO 800 IBC = 1, NBCS
         IBCTYP = IBCDAT(1,IBC)
         IDIR   = IBCDAT(2,IBC)
         IF ( IDIR .EQ. 2 .AND.
     1        ( IBCTYP .EQ. ITAN  .OR. IBCTYP .EQ. IWALL .OR.
     2          IBCTYP .EQ. IWALFN ) ) THEN
            INROUT = IBCDAT(3,IBC)
            ISTRT  = IBCDAT(4,IBC)
            IEND   = IBCDAT(5,IBC)
            KSTRT  = IBCDAT(8,IBC)
            KEND   = IBCDAT(9,IBC)
            IF (INROUT .EQ. 1) THEN
               JQPM  = 1
               JBCPT = 1
            ELSE
               JQPM  = JDIM
               JBCPT = JDIM + 1
            ENDIF
            IF ( I .GE. ISTRT .AND. I .LE. IEND ) THEN
               DO 790 L = 1, NQ
                  DO 780 K = KSTRT, KEND
                     QP(JQPM,K,L) = Q(JBCPT,K,I,L)
                     QM(JQPM,K,L) = Q(JBCPT,K,I,L)
  780             CONTINUE
  790          CONTINUE
            ENDIF
         ENDIF
  800 CONTINUE
C
C Finished with Q(+/-) calculation for MUSCL scheme in j.
C
      RETURN
      END

