      SUBROUTINE RESID_X (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT, Q,
     1                   QN, NQNPTS, PROPS, NPPTS, S, SB, VOL, RC, DTJ,
     2                   RES, RESB, FHAT, DQP, DQM, QP, QM, TAU, F, FB, 
     3                   FM, STEMP, STEMPB, WORKX, IBCDAT, 
     4                   NBCS, IORDER, LIMIT, RKAP, IFLUX,
     5                   VISCOS, IMODEL, INITRS, IFDELQ, MRCHNG, 
     6                   VGNRN, SIGMA, NPRLIM, NPRNEG, IFL2DQ, R2DQDT)
C
C Routine to calculate the residual throughout the flowfield for either 
C the time dependent calculation or the marching calculation.  
C The input logical INITRS determines whether to initialize the residual
C to zero or to not initialize the residual. INITRS = false is used when
C the residual is initialized in the multigrid cycle to the coarse grid
C forcing function.
C
C Time dependent calculation: IPSTRT = 2, IPEND = IDIM
C Marching calculation      : IPSTRT = IPEND = I plane presently at
C
C This routine is organized as follows:
C     1. If INITRS then zero residuals
C     2. If iterative time dependent then add dQ/dt term
C     3. Residual in I direction - Add dF/dx term
C     4. Residual in J direction - Add dG/dy term
C     5. If 3D or axisymmetric then residual in K direction - add dH/dz term
C     6. Calculate Source Terms
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at 
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end calculations at
C                  (Maximum Value: IDIM)
C MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                  (used in dimensioning temporary space used for 
C                   calculating q(+,-) and f.
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
C TAU            : Storage for Reynolds stresses to improve modularity
C F,FM           : Fluxes at cell interfaces -
C                  if Flux Difference Splitting then F = Flux
C                  if Flux Vector Splitting then F = Fplus, FM = Fminus
C STEMP          : Temporary storage for metrics
C WORKX          : Work array for cross-section (used for damping, etc.)
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C IORDER         : Order of MUSCL interpolation for current direction 
C                  IORDER = 1 first order
C                         = 2 second order (or 'third' order)
C LIMIT          : Limiter to use for MUSCL interpolation
C RKAP           : Kappa to use in MUSCL interpolation
C IFLUX          : Flux evaluation scheme
C VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C IMODEL         : Turbulence modeling scheme array
C INITRS         : Logical controlling initializing RES to zero for MG
C IFDELQ         : Array controlling calculation of delq
C MRCHNG         : Logical for Marching
C VGNRN          : Logical for Vigneron approach in marching flux evaluation
C SIGMA          : Safety factor for Vigneron technique when marching
C NPRLIM         : Number of times reset production
C NPRNEG         : Number of times invoke positivity preservation for production
C IFL2DQ         : Logical to calculate L2 norm of dQ/dt term for unsteady
C R2DQDT         : L2 norm of dQ/dt term for unsteady
C
      include '../../src-flo/header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QN    (NQNPTS,NQ,NTMLVL),
     2          PROPS (NPPTS,NP),
     3          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     4          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     5          RC    (JDIM+1,KDIM+1,IDIM+1,3),
     6          DTJ   (JDIM,KDIM,IDIM)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          FHAT  (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
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
C     Adjoint variables
C
      DIMENSION RESB  (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          FB    (MXSECT,NF),
     2          SB    (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          STEMPB(MXSECT,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IORDER(3), LIMIT(3), RKAP(3)
C
      DIMENSION IMODEL(NMDL),
     1          IFDELQ(NQ)
C
      DIMENSION R2DQDT(NF+1)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),
     2          IEND  (3)
C
      LOGICAL   VISCOS(3), INITRS, MRCHNG, VGNRN, IFL2DQ
C
      ITURB = IMODEL(1)
C
C 1.  If INITRS then zero residuals
C
C       INITRS = .true.
C       IF (INITRS) THEN
C          SB  (:,:,:,:,:) = 0.0E0
C       ENDIF
C
C 3.  Residual in I direction - Add dF/dx term
C
      CALL FLUXI_X (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1              Q, PROPS, NPPTS, S, SB, VOL, 
     2              RES, RESB, DQP, DQM, QP,
     3              QM, TAU, F, FB, FM, STEMP, STEMPB,
     4              IBCDAT, NBCS, IORDER(1), LIMIT(1), RKAP(1),
     5              IFLUX, VISCOS(1), IMODEL)
C
C 4.  Residual in J direction - Add dG/dy term
C
      CALL FLUXJ_X  (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1              Q, PROPS, NPPTS, S, SB, VOL, 
     2              RES, RESB, DQP, DQM, QP,
     3              QM, TAU, F, FB, FM, STEMP, STEMPB,
     3              IBCDAT, NBCS, IORDER(2), LIMIT(2), RKAP(2), 
     4              IFLUX, VISCOS(2), IMODEL)
C
C 5.  If 3D or axisymmetric then Residual in K direction - add dH/dz term
C
      IF ( THREED ) THEN
         CALL FLUXK_X(IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                Q, PROPS, NPPTS, S, SB, VOL, 
     2                RES, RESB, DQP, DQM, QP, 
     3                QM, TAU, F, FB, FM, STEMP, STEMPB,
     3                IBCDAT, NBCS, IORDER(3), LIMIT(3), RKAP(3), 
     4                IFLUX, VISCOS(3), IMODEL)
      ELSE IF ( AXISYM ) THEN
         STOP "RESID_Q: FLUXKA NOT DONE"
         CALL FLUXKA (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                Q, PROPS, NPPTS, S, VOL, 
     2                RES, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                IBCDAT, NBCS, IORDER(3), LIMIT(3), RKAP(3), 
     4                IFLUX, VISCOS(3), IMODEL)
      ENDIF
C
C     Finished with residual computation.
C
      RETURN
      END
