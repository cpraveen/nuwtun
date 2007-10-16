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
      SUBROUTINE FLUXKA (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   RES, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                   IBCDAT, NBCS, 
     4                   IORDER, LIMIT, RKAP, IFLUX, VISCOS, IMODEL)
C
C Routine to calculate the axisymmetric flux contriubtion to the residual 
C from the k direction.  This routine is used for either the time dependent 
C calculation or the marching calculation.  
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at 
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                  (used in dimensioning temporary space used for 
C                   calculating q(+,-) and f.
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity 
C                  PROPS(3) = Y+
C NPPTS          : Number of points at which PROPS are stored.  This allows
C                  the storage of PROPS for viscous flows and the avoidance
C                  of the additional storage for inviscid flows.
C S              : Metrics
C VOL            : Cell volumes
C RES            : Residual
C DQP,DQM        : Storage for delta q = q(k) - q(k-1)
C QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C TAU            : Storage for Reynolds stresses at interfaces
C F,FM           : Fluxes at cell interfaces -
C                  if Flux Difference Splitting then F = Flux
C                  if Flux Vector Splitting then F = Fplus, FM = Fminus
C STEMP          : Storage for metrics passed to flux evaluation routine.
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C IORDER         : Order of MUSCL interpolation for K direction 
C                  IORDER = 1 first order
C                         = 2 second order (or 'third' order)
C LIMIT          : Limiter to use for MUSCL interpolation
C RKAP           : Kappa to use in MUSCL interpolation
C IFLUX          : Flux evaluation scheme
C VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C IMODEL         : Turbulence model array
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          DQP   (MXSECT,NQ), 
     2          DQM   (MXSECT,NQ),
     3          QP    (KDIM,2:JDIM,NQ), 
     4          QM    (MXSECT,NQ),
     5          TAU   (MXSECT,6),
     6          F     (KDIM,2:JDIM,NF), 
     7          FM    (KDIM,2:JDIM,NF), 
     8          STEMP (KDIM,2:JDIM,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IMODEL(NMDL)
C
      LOGICAL   VISCOS
C
      ITURB = IMODEL(1)
C
C     Loop over I calculating axisymmetric fluxes and updating the residual.
C
      DO 500 I = IPSTRT, IPEND
C
C 1.  Inviscid fluxes: Pressure only on the v- and w-momentum equations
C
         DO 110 J = 2, JDIM
            P            = Q(J,2,I,5)
C
            F(1,J,3)     = P * S(J,1,I,2,3)
            F(2,J,3)     = P * S(J,2,I,2,3)
            RES(J,2,I,3) = RES(J,2,I,3) + 
     1                     ( F(2,J,3) * S(J,2,I,4,3) -
     2                       F(1,J,3) * S(J,1,I,4,3) )
C
            F(1,J,4)     = P * S(J,1,I,3,3)
            F(2,J,4)     = P * S(J,2,I,3,3)
            RES(J,2,I,4) = RES(J,2,I,4) + 
     1                     ( F(2,J,4) * S(J,2,I,4,3) -
     2                       F(1,J,4) * S(J,1,I,4,3) )
  110    CONTINUE
C
C 2.  Viscous Residuals
C
         IF ( VISCOS ) THEN
C
C     Calculate the derivatives and averages of stored variables at the
C     interfaces.
C     Temporary variables passed to the viscous flux evaluation routine:
C         QP    for Q      (variables at interface)
C         FM    for PROPSI (properties at interface)
C               NOTE: This works properly only as long as NP <= NQ
C         QM    for DQDX   (derivatives of Q in x direction at interface)
C         DQP   for DQDY   (derivatives of Q in y direction at interface)
C         DQM   for DQDZ   (derivatives of Q in z direction at interface)
C     Temporary variables used internally in DELQ:
C         F(2-NQ) for SK (STEMP is dimensioned incorrectly to be used for
C                     the averaged metrics in the J direction.  F is large
C                     enough to provide the necessary space.  Do not use
C                     F for anything else in this call).
C         F(1)    for VOLAUX 
C
            CALL DELQK (IDIM, JDIM, KDIM, I,
     1                  Q, PROPS, S, VOL, IBCDAT, NBCS,
     2                  F(1,2,2), F(1,2,1), QP, FM, QM, DQP, DQM)
C
C     Metrics from the current direction.
C
            DO 320 K = 1, KDIM
               DO 310 J = 2, JDIM
                  STEMP(K,J,1) = S(J,K,I,1,3)
                  STEMP(K,J,2) = S(J,K,I,2,3)
                  STEMP(K,J,3) = S(J,K,I,3,3)
  310          CONTINUE
  320       CONTINUE
            NPTS = KDIM * (JDIM - 1)
C
C     Viscous Terms 
C
            CALL FVIS(IMODEL, NPTS, QP, FM, QM, DQP, DQM, STEMP, TAU, F)
C
C     Add fluxes to residual
C
            DO 380 L = 1, NF
               DO 360 K = 2, KDIM
                  DO 340 J = 2, JDIM
                     RES(J,K,I,L) = RES(J,K,I,L) - 
     1                              ( F(K,J,L)   * S(J,K,I,4,3) -
     2                                F(K-1,J,L) * S(J,K-1,I,4,3) )
  340             CONTINUE
  360          CONTINUE
  380       CONTINUE
C
         ENDIF
C
  500 CONTINUE
C
C     Finished with k direction axisymmetric residual computation.
C
      RETURN
      END
