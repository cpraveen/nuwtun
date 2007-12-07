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
      SUBROUTINE FLUXJ_X (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                  Q, PROPS, NPPTS, S, SB, VOL, RES, RESB, 
     2                  DQP, DQM, QP, QM, TAU, F, FB, FM, STEMP, STEMPB,
     3                  IBCDAT, NBCS, IORDER, LIMIT, RKAP, IFLUX,
     4                  VISCOS, IMODEL)
C
C     Routine to calculate the contribution to the residual from the
C     j direction.  This routine is used for either the 
C     time dependent calculation or the marching calculation.  
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     IPSTRT         : I Plane (cell center) to start calculations at 
C                      (Minimum Value: 2)
C     IPEND          : I Plane (cell center) to end caclulations at
C                      (Maximum Value: IDIM)
C     MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                      (used in dimensioning temporary space used for 
C                       calculating q(+,-) and f.
C     Q              : Primitive variables at cell centers
C     PROPS          : Properties at cell centers
C                      PROPS(1) = molecular viscosity
C                      PROPS(2) = turbulent eddy viscosity 
C                      PROPS(3) = Y+
C     NPPTS          : Number of points at which PROPS are stored.  This allows
C                      the storage of PROPS for viscous flows and the avoidance
C                      of the additional storage for inviscid flows.
C     S              : Metrics
C     VOL            : Cell volumes
C     RES            : Residual
C     DQP,DQM        : Storage for delta q = q(j) - q(j-1)
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     TAU            : Storage for Reynolds stresses at interfaces
C     F,FM           : Fluxes at cell interfaces -
C                      if Flux Difference Splitting then F = Flux
C                      if Flux Vector Splitting then F = Fplus, FM = Fminus
C     STEMP          : Storage for metrics passed to flux evaluation routine.
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IORDER         : Order of MUSCL interpolation for J direction 
C                      IORDER = 1 first order
C                             = 2 second order (or 'third' order)
C     LIMIT          : Limiter to use for MUSCL interpolation
C     RKAP           : Kappa to use in MUSCL interpolation
C     IFLUX          : Flux evaluation scheme
C     VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C     IMODEL         : Turbulence model array
C
      include '../../src-flo/header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          DQP   (MXSECT,NQ), 
     2          DQM   (MXSECT,NQ),
     3          QP    (MXSECT,NQ), 
     4          QM    (MXSECT,NQ),
     5          TAU   (MXSECT,6),
     6          F     (JDIM,2:KDIM,NF), 
     7          FM    (JDIM,2:KDIM,NF), 
     8          STEMP (JDIM,2:KDIM,3)
C
      DIMENSION RESB  (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     2          SB    (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     6          FB    (JDIM,2:KDIM,NF), 
     8          STEMPB(JDIM,2:KDIM,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IMODEL(NMDL)
C
      LOGICAL   VISCOS
C
      ITURB = IMODEL(1)
C
C     Loop over I calculating Q(+,-), F(Q(+,-)) and updating the Residual.
C
      DO 500 I = IPSTRT, IPEND
C
C     Inviscid Residuals
C
         CALL MUSCLJ (IDIM, JDIM, KDIM, I,
     1                Q, DQP, DQM, QP, QM,
     2                IBCDAT, NBCS, IORDER, LIMIT, RKAP)
C
C     Calculate the fluxes from q(plus) and q(minus).  Pass only
C     the appropriate metrics from the current direction.
C
         DO 120 K = 2, KDIM
            DO 110 J = 1, JDIM
               DO 100 L = 1, 3
                  STEMP(J,K,L) = S(J,K,I,L,2)
  100          CONTINUE
  110       CONTINUE
  120    CONTINUE
         NPTS = JDIM * (KDIM - 1)
C
C     Flux Difference Splitting
C
         IF (IFCHAR(IFLUX,1) .EQ. IFDS) THEN
C
C     Add the flux contributions to the residual.
C
            CALL FDS (IFLUX, ITURB, NPTS, QP, QM, STEMP, F)
         ENDIF
C
C     Backward sweep
C
         FB     (:,:,:) = 0.0
         STEMPB (:,:,:) = 0.0
C
C     Flux Difference Splitting
C
         IF (IFCHAR(IFLUX,1) .EQ. IFDS) THEN
C
C     Add the flux contributions to the residual.
C
            DO 180 L = 1, NF
               DO 160 K = 2, KDIM
                  DO 140 J = 2, JDIM
                     FB(J,K,L)       = FB(J,K,L)     + 
     1                                 RESB(J,K,I,L) * S(J,K,I,4,2)
                     FB(J-1,K,L)     = FB(J-1,K,L)   - 
     1                                 RESB(J,K,I,L) * S(J-1,K,I,4,2)
                     SB(J,K,I,4,2)   = SB(J,K,I,4,2)   +
     1                                 F(J,K,L) * RESB(J,K,I,L)
                     SB(J-1,K,I,4,2) = SB(J-1,K,I,4,2) -
     1                                 F(J-1,K,L) * RESB(J,K,I,L)
  140             CONTINUE
  160          CONTINUE
  180       CONTINUE
            CALL FDS_X(IFLUX, ITURB, NPTS, QP, QM, STEMP, STEMPB, F, FB)
         ENDIF
C
C     Add contributions from temporary array STEMPB to SB
C
         DO K = 2, KDIM
            DO J = 1, JDIM
               DO L = 1, 3
                  SB(J,K,I,L,2) = SB(J,K,I,L,2) + STEMPB(J,K,L)
               ENDDO
            ENDDO
         ENDDO
C
C     Viscous Residuals - To be continued
C
  500 CONTINUE
C
C     Finished with j direction residual computation.
C
      RETURN
      END
