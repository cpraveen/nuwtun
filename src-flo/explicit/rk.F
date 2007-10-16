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

      SUBROUTINE RK (IBLK, IDIM, JDIM, KDIM, MXSECT, 
     1               NSTAGE, ISTAGE, RKALPH,
     2               Q, QN, NQNPTS, PROPS, NPPTS, S, VOL, RC, DTJ, RES,
     3               FHAT, DQP, DQM, QP, QM, TAU, F, FM, STEMP, WORKX,
     4               IBCDAT, NBCS, IORDER, LIMIT, RKAP, IFLUX, VISCOS,
     5               IMODEL, IFDELQ, W,
     6               QMIN, NFIXQ, NRELIZ, NPRLIM, NPRNEG, IERRCD)
C
C Routine for a single stage of multistage Runge-Kutta time integration
C
C IBLK           : Current block
C IDIM,JDIM,KDIM : Dimensions of current block
C MXSECT         : Maximum size of any (ixj, jxk or ixk) 
C                  (used in dimensioning of temporary space used 
C                   for calculating q(+,-) and f.
C NSTAGE         : Number of stages of multistage scheme
C ISTAGE         : Current stage of multistage scheme
C RKALPH         : Coefficient for current stage
C Q              : Primitive variables at cell centers
C QN             : Conserved variables at time levels n, n-1, etc.
C NQNPTS         : Number of QN points. This allows avoiding storing QN when
C                  running time accurate.
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity 
C                  PROPS(3) = Y+
C NPPTS          : Number of points at which PROPS are stored.  This allows
C                  the storage of PROPS for viscous flows and the avoidance
C                  of the additional storage for inviscid flows.
C S              : Metrics
C VOL            : Cell volumes
C RC             : Position vector at cell centers
C DTJ            : Time step divided by cell volume 
C RES            : Residual
C FHAT           : Storage for fluxes for 4th order
C DQP,DQM        : Storage for delta q = q(i) - q(i-1)
C QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C TAU            : Storage for Reynolds stresses
C F,FM           : Fluxes at cell interfaces - 
C                  if Flux Difference Splitting then F = Flux
C                  if Flux Vector Splitting then F = Fplus, FM = Fminus
C STEMP          : Storage for metrics passed to flux evaluation routine.
C WORKX          : Work array for cross-section (used for damping, etc.)
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C IORDER         : Order of MUSCL interpolation 
C                  IORDER = 1 first order
C                         = 2 second order (or 'third' order)
C LIMIT          : Limiter to use for MUSCL interpolation
C RKAP           : Kappa to use in MUSCL interpolation
C IFLUX          : Flux evaluation scheme
C VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C IMODEL         : Turbulence modeling scheme array
C IFDELQ         : Array controlling calculation of delq
C W              : Temporary storage for Q in the Runge-Kutta scheme
C QMIN           : Control array for minimum value of Q
C NFIXQ          : Number of times reset Q
C NRELIZ         : Number of times violate realizability
C NPRLIM         : Number of times reset production
C NPRNEG         : Number of times invoke positivity preservation for production
C IERRCD         : Error code
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QN    (NQNPTS,NQ,NTMLVL),
     2          PROPS (NPPTS,NP),
     3          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     4          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1), 
     5          RC    (JDIM+1,KDIM+1,IDIM+1,3), 
     6          DTJ   (JDIM,KDIM,IDIM)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IDIM+1,NF),
     1          FHAT  (JDIM+1,KDIM+1,IDIM+1,NF),
     2          DQP   (MXSECT,NQ),
     3          DQM   (MXSECT,NQ),
     4          QP    (MXSECT,NQ),
     5          QM    (MXSECT,NQ),
     6          TAU   (MXSECT,6),
     7          F     (MXSECT,NF),
     8          FM    (MXSECT,NF),
     9          STEMP (MXSECT,3),
     A          WORKX (MXSECT,NW)
C
      DIMENSION IMODEL(NMDL),
     1          IFDELQ(NQ)
C
      DIMENSION W     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IORDER(3), RKAP(3), LIMIT(3)
C
      LOGICAL   VISCOS(3)
C
      DIMENSION QMIN  (NQ)
      DIMENSION NFIXQ (NQ),
     1          NRELIZ(NQ)
C
      ITURB  = IMODEL(1)
C
C     IPSTRT         : I Plane (cell center) to start calculations at 
C                      (Minimum Value: 2)
C     IPEND          : I Plane (cell center) to end caclulations at
C                      (Maximum Value: IDIM)
C
      IPSTRT = 2
      IPEND  = IDIM 
C
C     For first stage, Set W = Q(t=n)
C
      IF (ISTAGE .EQ. 1) THEN
         DO 140 L = 1, NQ
            DO 130 I = 0, IDIM + 2
               DO 120 K = 0, KDIM + 2
                  DO 110 J = 0, JDIM + 2
                     W(J,K,I,L) = Q(J,K,I,L) 
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
      ENDIF
C
C     Perform one stage of Runge-Kutta integration
C
C     Set Q = Q(t=n) 
C
      DO 540 L = 1, NQ
         DO 530 I = 2, IDIM
            DO 520 K = 2, KDIM
               DO 510 J = 2, JDIM
                  Q(J,K,I,L)   = W(J,K,I,L) 
  510          CONTINUE
  520       CONTINUE
  530    CONTINUE
  540 CONTINUE
C
C     Multiply residual by Runge-Kutta coefficient and DTJ to get Delta(Q)
C
      DO 640 L = 1, NF
         DO 630 I = 2, IDIM
            DO 620 K = 2, KDIM
               DO 610 J = 2, JDIM
                  RES(J,K,I,L) = - RKALPH * DTJ(J,K,I) * RES(J,K,I,L)
  610          CONTINUE
  620       CONTINUE
  630    CONTINUE
  640 CONTINUE
C
C     Update Conserved Variables
C
      CALL UPDATE (IDIM, JDIM, KDIM, IPSTRT, IPEND, Q, RES, ITURB,
     1             QP, QM, QMIN, NFIXQ, NRELIZ, FHAT, IERRCD)
C
C     Finished with one stage of Runge-Kutta time step 
C
      RETURN
      END
