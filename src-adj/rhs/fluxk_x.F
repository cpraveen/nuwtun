      SUBROUTINE FLUXK_X (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                  Q, PROPS, NPPTS, S, SB, VOL, RES, RESB,
     2                  DQP, DQM, QP, QM, TAU, F, FB, FM, STEMP, STEMPB,
     3                  IBCDAT, NBCS, IORDER, LIMIT, RKAP, IFLUX,
     4                  VISCOS, IMODEL)
C
C     Routine to calculate the contribution to the residual from the
C     k direction.  This routine is used for either the 
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
C     DQP,DQM        : Storage for delta q = q(k) - q(k-1)
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     TAU            : Storage for Reynolds stresses at interfaces
C     F,FM           : Fluxes at cell interfaces -
C                      if Flux Difference Splitting then F = Flux
C                      if Flux Vector Splitting then F = Fplus, FM = Fminus
C     STEMP          : Storage for metrics passed to flux evaluation routine.
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IORDER         : Order of MUSCL interpolation for K direction 
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
     6          F     (KDIM,2:JDIM,NF), 
     7          FM    (KDIM,2:JDIM,NF), 
     8          STEMP (KDIM,2:JDIM,3)
C
      DIMENSION RESB  (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          SB    (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     2          FB    (KDIM,2:JDIM,NF), 
     3          STEMPB(KDIM,2:JDIM,3)
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
         CALL MUSCLK (IDIM, JDIM, KDIM, I,
     1                Q, DQP, DQM, QP, QM,
     2                IBCDAT, NBCS, IORDER, LIMIT, RKAP)
C
C     Calculate the fluxes from q(plus) and q(minus).  Pass only
C     the appropriate metrics from the current direction.
C
         DO 120 J = 2, JDIM
            DO 110 K = 1, KDIM
               DO 100 L = 1, 3
                  STEMP(K,J,L) = S(J,K,I,L,3)
  100          CONTINUE
  110       CONTINUE
  120    CONTINUE
         NPTS = KDIM * (JDIM - 1)
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
         FB    (:,:,:) = 0.0E0
         STEMPB(:,:,:) = 0.0E0
C
C     Flux Difference Splitting
C
         IF (IFCHAR(IFLUX,1) .EQ. IFDS) THEN
C
C     Add the flux contributions to the residual.
C
            DO 180 L = 1, NF
               DO 160 J = 2, JDIM
                  DO 140 K = 2, KDIM
                     FB(K,  J,L)     = FB(K,J,L) +
     1                                 RESB(J,K,I,L) * S(J,K,I,4,3)
                     FB(K-1,J,L)     = FB(K-1,J,L) -
     1                                 RESB(J,K,I,L) * S(J,K-1,I,4,3)
                     SB(J,K,I,4,3)   = SB(J,K,I,4,3) +
     1                                 F(K,J,L) * RESB(J,K,I,L)
                     SB(J,K-1,I,4,3) = SB(J,K-1,I,4,3) -
     1                                 F(K-1,J,L) * RESB(J,K,I,L)
  140             CONTINUE
  160          CONTINUE
  180       CONTINUE
            CALL FDS_X(IFLUX, ITURB, NPTS, QP, QM, STEMP, STEMPB, F, FB)
         ENDIF
C
C     Add contributions from temporary array STEMPB to SB
C
         DO J = 2, JDIM
            DO K = 1, KDIM
               DO L = 1, 3
                  SB(J,K,I,L,3) = SB(J,K,I,L,3) + STEMPB(K,J,L)
               ENDDO
            ENDDO
         ENDDO
C
C     Viscous Residuals - TO BE CONTINUED
C
  500 CONTINUE
C
C     Finished with k direction residual computation.
C
      RETURN
      END
