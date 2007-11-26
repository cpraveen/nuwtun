      SUBROUTINE RK (IDIM, JDIM, KDIM, ISTAGE, RKALPH,
     1               QB, DTJ, QA, W)
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
      include '../../src-flo/header/common.h'
C
      DIMENSION QB    (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          DTJ   (JDIM,KDIM,IDIM),
     2          QA    (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     3          W     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
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
                     W(J,K,I,L) = QA(J,K,I,L) 
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
c-    DO 540 L = 1, NQ
c-       DO 530 I = 2, IDIM
c-          DO 520 K = 2, KDIM
c-             DO 510 J = 2, JDIM
c-                QA(J,K,I,L)   = W(J,K,I,L) 
c-510          CONTINUE
c-520       CONTINUE
c-530    CONTINUE
c-540 CONTINUE
C
C     Update adjoint Variables
C
      DO 640 L = 1, NQ
         DO 630 I = 2, IDIM
            DO 620 K = 2, KDIM
               DO 610 J = 2, JDIM
c-                DRESB       = - RKALPH * DTJ(J,K,I) * QB(J,K,I,L)
c-                QA(J,K,I,L) = QA(J,K,I,L) + DRESB
                  DRESB       = - DTJ(J,K,I) * QB(J,K,I,L)
                  QA(J,K,I,L) = RKALPH*W(J,K,I,L) + 
     1                          (1.0-RKALPH)*(QA(J,K,I,L) + DRESB)
  610          CONTINUE
  620       CONTINUE
  630    CONTINUE
  640 CONTINUE
C
C     Finished with one stage of Runge-Kutta time step 
C
      RETURN
      END
