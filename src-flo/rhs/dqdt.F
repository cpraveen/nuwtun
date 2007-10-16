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
      SUBROUTINE DQDT   (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, QN, VOL, DTJ, RES, QCP, TKE, IMODEL)
C
C Routine to calculate the dQ/dt term to add to residual.
C     t-ts: Physical Time Subiteration (Rumsey et al. Computers&Fluids, 1996)
C           dQ/dt = [(1 + phi) Q^(n+1) - (1 + 2 phi) Q^(n) + phi Q^(n-1)] / dt
C
C     tau-ts: Pseudo Time Subiteration (Rumsey et al. Computers&Fluids, 1996)
C           dQ/dt = [(1 + phi) Q^(n+1) - (1 + 2 phi) Q^(n) + phi Q^(n-1)] / dt
C                 + [phi' (Q^(m) - Q^(m-1))] / dtau
C
C     TDTHET : Theta for Time Dependent calculation (physical time term)
C     TDPHI  : Phi   for Time Dependent calculation (physical time term)
C     TDPHIP : Phi'  for Time Dependent calculation (pseudo   time term)
C
C +-------+-----------------+-------------+-------------------+---------------+
C |Scheme | Order           | Theta  Phi  |  Phi' store dQ^(m)| TDPDTJ TDPDEL |
C +-------+-----------------+-------------+-------------------+---------------+
C |t-ts   | 1st             |   1     0   |   -         No    | 1+Phi     0   |
C |       | 2nd             |   1    1/2  |   -         No    | 1+Phi     0   |
C |tau-ts | 2nd, 1st pseudo |   1    1/2  |   0         No    | 1+Phi'  1+Phi |
C |       | 2nd, 2nd pseudo |   1    1/2  |  1/2       Yes    | 1+Phi'  1+Phi |
C +-------+-----------------+-------------+-------------------+---------------+
C
C N.B. This routine assumes that the previous time level solutions are
C      passed as conserved variables rather than primitive variables.
C
C This routine is organized as follows:
C     1. Copy Q to work array
C     2. Convert Q from primitive to conserved variables
C     3. Calculate dQ/dt and add to residual array
C     4. Q^(m) - Q^(m-1) term for tau-ts 
C        a. Add Q^(m) - Q^(m-1) term for tau-ts with 2nd order psuedo 
C           time stepping
C        b. Shift Q^(m) to Q^(m-1) for next subiteration
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                  (used in dimensioning temporary space)
C Q              : Primitive variables at iterative time level p
C QN             : Conserved variables at time levels n, n-1, etc.
C VOL            : Cell volumes
C DTJ            : Time step divided by cell volume
C RES            : Residual
C QCP            : Work array for conserved Q at level p
C TKE            : Work array for turbulent kinetic energy
C IMODEL         : Turbulence model
C
      include '../header/common.h'
C
      INTEGER   IDIM, JDIM, KDIM
C
      DIMENSION Q   (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QN  (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ,NTMLVL),
     2          VOL (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     3          DTJ (JDIM,KDIM,IDIM),
     4          RES (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF)
C
      DIMENSION QCP (2:JDIM,2:KDIM,NQ),
     1          TKE (2:JDIM,2:KDIM)
C
      DIMENSION IMODEL(NMDL)
C
C Local variables
C
      INTEGER   I, J, K, L, NPTS
C
C     Turbulence model
C
      ITURB = IMODEL(1)
C
C     Constants for temporal accuracy 
C     (see e.g. Pulliam AIAA 93-3360, Rai AIAA 95-0583, 
C               Rumsey et al. Computers & Fluids 1996)
C          CONP constant for level p, CONN constant for time n, 
C          CONNM1 constant for time n-1
C
      CONP   = 1.E0 + TDPHI
      CONN   = 1.E0 + 2.E0 * TDPHI
      CONNM1 = TDPHI
C
C      Loop over each cross section
C
      DO 900 I = IPSTRT, IPEND
C
C 1.   Copy Q to work array
C
         DO 130 L = 1, NQ
            DO 120 K = 2, KDIM
               DO 110 J = 2, JDIM
                  QCP(J,K,L) = Q(J,K,I,L)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
C
C 2.   Convert Q from primitive to conserved
C
         NPTS = (JDIM-1) * (KDIM-1)
         CALL CONP2C (ITURB, NPTS, QCP, TKE)
C
C 3.   Calculate dQ/dt
C
         DO 330 L = 1, NQ
            DO 320 K = 2, KDIM
               DO 310 J = 2, JDIM
                  DQDTP = ( CONP   * QCP(J,K,L) - CONN * QN(J,K,I,L,1) 
     1                    + CONNM1 * QN(J,K,I,L,2) ) * VOL(J,K,I)/TIMEDT
                  RES(J,K,I,L) = RES(J,K,I,L) + DQDTP
  310          CONTINUE
  320       CONTINUE
  330    CONTINUE
C
C 4.   Q^(m) - Q^(m-1) term for tau-ts 
C
         IF (ITMTYP .EQ. ITAUTS .AND. TDPHIP .NE. 0.E0) THEN
            IT = NTIME + 1
C
C 4-a. Add Q^(m) - Q^(m-1) term for tau-ts with 2nd order psuedo time stepping
C
            DO 430 L = 1, NQ
               DO 420 K = 2, KDIM
                  DO 410 J = 2, JDIM
                     DQMDT = TDPHIP * (QCP(J,K,L) - QN(J,K,I,L,IT))
     1                       / DTJ(J,K,I)
                     RES(J,K,I,L) = RES(J,K,I,L) + DQMDT
  410             CONTINUE
  420          CONTINUE
  430       CONTINUE
C
C 4-b. Shift Q^(m) to Q^(m-1) for next subiteration
C
            DO 460 L = 1, NQ
               DO 450 K = 2, KDIM
                  DO 440 J = 2, JDIM
                     QN(J,K,I,L,IT) = QCP(J,K,L)
  440             CONTINUE
  450          CONTINUE
  460       CONTINUE
         ENDIF
C
C      Next cross section
C
  900 CONTINUE
C
C      Finished calculating dQ/dt
C
      RETURN
      END
