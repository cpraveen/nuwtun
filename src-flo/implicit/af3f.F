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
      SUBROUTINE AF3F (IBLK, IDIM, JDIM, KDIM, MXSECT, MXABC,
     1                 Q, QN, NQNPTS, PROPS, NPPTS, S, VOL, RC, DTJ, 
     2                 RES, FHAT, DQP, DQM, QP, QM,TAU,F,FM,STEMP,WORKX,
     3                 IBCDAT, NBCS, 
     4                 IORDER, LIMIT, RKAP, IFLUX, VISCOS, DIAG, 
     5                 IMODEL, IFDELQ, ABC, RTEMP, AT, AI,
     6                 QMIN, NFIXQ, NRELIZ, NPRLIM, NPRNEG, IERRCD)
C
C Routine to do spatially split Three Factor Approximate Factorisation
C using block inversions with implicit terms of the form 
C     df/dq = T |Lambda| T(-1).
C For inviscid, laminar, algebraic turbulence model:
C     (I + DTJ A) (I + DTJ B) (I + DTJ C) = - DTJ RES
C for implicit source terms (including 2 eqn turbulence models)
C     (I + DTJ D) (I + DTJ A) (I + DTJ B) (I + DTJ C) = - DTJ RES
C where
C     A = dF/dQ, B = dG/dQ, C = dH/dQ and D = dS/dQ (S = source terms).
C IJKMAX is a parameter passed from the main storage holding the
C value of the maximum dimension anywhere in the grid which is
C used to set up the temporary storage.
C
C IBLK           : Current block
C IDIM,JDIM,KDIM : Dimensions of current block
C MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                  (used in dimensioning of temporary space used
C                   for calculating q(+,-) and f.
C MXABC          : Amount of memory to use to store the tridiagonal
C                  systems to invert.  Must be at least
C                  max (IDIM*JDIM, JDIM*KDIM) * number of planes to
C                  vectorize inversion on
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
C DIAG           : Control for diagonalized scheme (TRUE) or block (FALSE)
C IMODEL         : Turbulence modeling scheme array
C IFDELQ         : Array controlling calculation of delq
C ABC            : Storage for the tridiagonal system to invert
C RTEMP          : Storage for the residual to be put in the proper
C                  form to use with the inversion routine
C AT             : Storage for the Jacobian calculations
C AI             : Identity Matrix
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
     9          STEMP (MXSECT,4),
     A          WORKX (MXSECT,NW)
C
      DIMENSION IMODEL(NMDL),
     1          IFDELQ(NQ)
C
      DIMENSION ABC   (MXABC*NRANK*NRANK,3),
     1          RTEMP (MXABC*NF),
     2          AT    (MXSECT*NRANK*NRANK,3)
C
      DIMENSION AI    (NRANK,NRANK)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IORDER(3), RKAP(3), LIMIT(3)
C
      LOGICAL   VISCOS(3), DIAG(3)
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
C     Set RES = - RES / ( (1+phi')/(dtau*J) + (1+phi)/(dt*J) ) for inversions
C
      DO 140 L = 1, NF
         DO 130 I = 2, IDIM 
            DO 120 K = 2, KDIM 
               DO 110 J = 2, JDIM 
                  DTTERM       = 1.0E0 / 
     1                           ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                           + TDPDEL * VOL(J,K,I) / TIMEDT )
                  RES(J,K,I,L) = - DTTERM * RES(J,K,I,L)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
C
C     Perform the (I + DTJ D) inversion for the case with implicit source 
C     terms -> IMPSRC not equal INO and 2 equation or higher turbulence model
C
      IF (IMPSRC .NE. INO) THEN
         IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW .OR.
     1       ITURB .EQ. ITRS) THEN
C
C     Calculate INPL to be the number of planes that will fit
C     in the memory allocated.  Actually vectorize over NPL 
C     which takes into account if the last loop would be longer
C     then IDIM.  K is the primary vectorization direction and
C     I is the secondary vectorization direction.
C     Make sure that INPL is no larger than the maximum size of the grid.
C     Temporary storage: Use WORKX(1)    for TKE ,     WORKX(2)     for ETILDE,
C                            WORKX(3)    for XSIW,     WORKX(4-5)   for FDAMP,
C                            WORKX(6)    for EPSC,     WORKX(7)     for PDIL,
C                            WORKX(8-13) for SIJ ,     WORKX(14-19) for WIJ,
C                            F           for SRC
C
            INPL = MXABC / ( (KDIM-1) * (JDIM-1) )
            INPL = MIN (INPL, (IDIM-1))
            NPL  = INPL
            DO 190 I = 2, IDIM, INPL
               IF ( (I+INPL) .GT. IDIM ) NPL = IDIM - I + 1
               CALL PISRC (IDIM, JDIM, KDIM, NPL, I, IMODEL, IFDELQ,
     1                     Q, PROPS, NPPTS, S, VOL, RC, DTJ, RES,
     2                     MXSECT, DQP, DQM, QP, QM, TAU, FM, STEMP,
     3                     WORKX(1,1),WORKX(1,2),WORKX(1,3),WORKX(1,4), 
     4                     WORKX(1,6),WORKX(1,7),WORKX(1,8),WORKX(1,14),
     5                     F, IBCDAT, NBCS, ABC(1,2), RTEMP, AI,AT(1,1))
  190       CONTINUE
         ENDIF
      ENDIF
C
C     I implicit sweep
C
C     I Diagonalized Inversion:
C     Use F(1-3), FM(1-3) for lhat, mhat.  Use F(5) for TKE.  FM(5) for VOLJ.
C     Use QM for PROPSJ.  AT(3) for EIGENP, DQP for EIGENM, DQM for EIGEN.
C
      IF ( DIAG(1) ) THEN
         DO 200 K = 2, KDIM
            CALL AF3DI  (IDIM, JDIM, KDIM, K,
     1                   Q, PROPS, NPPTS, S, VOL, DTJ, RES,
     2                   QP, QM, STEMP, FM(1,5), F(1,1), FM(1,1),F(1,5),
     3                   IBCDAT, NBCS, IFLUX, VISCOS(1), ITURB,
     4                   ABC(1,1), ABC(1,2), ABC(1,3), RTEMP,
     5                   AT(1,1), AT(1,2), AT(1,3), DQP, DQM)
  200    CONTINUE
C
C     I Block Inversion:
C     Calculate INPL to be the number of planes that will fit
C     in the memory allocated.  Actually vectorize over NPL 
C     which takes into account if the last loop would be longer
C     then KDIM.  J is the primary vectorization direction and
C     K is the secondary vectorization direction.
C     Make sure that INPL is no larger than the maximum size of the grid.
C
      ELSE
         INPL = MXABC / ( (IDIM-1) * (JDIM-1) )
         INPL = MIN (INPL, (KDIM-1))
         NPL  = INPL
         DO 210 K = 2, KDIM, INPL
            IF ( (K+INPL) .GT. KDIM ) NPL = KDIM - K + 1
            CALL AF3SWI (IDIM, JDIM, KDIM, NPL, K,
     1                   Q, PROPS, NPPTS, S, VOL, DTJ, RES, QP, QM, 
     2                   STEMP, IBCDAT, NBCS, IFLUX, VISCOS(1), ITURB,
     3                   ABC(1,1), ABC(1,2), ABC(1,3), RTEMP,
     4                   AT(1,1), AT(1,2), AT(1,3), AI, DQP)
  210    CONTINUE
      ENDIF
C
C     J implicit sweep
C
C     J Diagonal Inversion:
C     Use F(1-3), FM(1-3) for lhat, mhat.  Use F(5) for TKE.  FM(5) for VOLJ.
C     Use QM for PROPSJ.  AT(3) for EIGENP, DQP for EIGENM, DQM for EIGEN.
C
      IF ( DIAG(2) ) THEN
         DO 300 K = 2, KDIM
            CALL AF3DJ  (IDIM, JDIM, KDIM, K,
     1                   Q, PROPS, NPPTS, S, VOL, DTJ, RES, 
     2                   QP, QM, STEMP, FM(1,5), F(1,1), FM(1,1),F(1,5), 
     3                   IBCDAT, NBCS, IFLUX, VISCOS(2), ITURB,
     4                   ABC(1,1), ABC(1,2), ABC(1,3), RTEMP,
     5                   AT(1,1), AT(1,2), AT(1,3), DQP, DQM) 
  300    CONTINUE
C
C     J Block Inversion:
C     Calculate INPL to be the number of planes that will fit
C     in the memory allocated.  Actually vectorize over NPL 
C     which takes into account if the last loop would be longer
C     then IDIM.  K is the primary vectorization direction and
C     I is the secondary vectorization direction.
C     Make sure that INPL is no larger than the maximum size of the grid.
C
      ELSE
         INPL = MXABC / ( (KDIM-1) * (JDIM-1) )
         INPL = MIN (INPL, (IDIM-1))
         NPL  = INPL
         DO 310 I = 2, IDIM, INPL
            IF ( (I+INPL) .GT. IDIM ) NPL = IDIM - I + 1
            CALL AF3SWJ (IDIM, JDIM, KDIM, NPL, I,
     1                   Q, PROPS, NPPTS, S, VOL, DTJ, RES, QP, QM, 
     2                   STEMP, IBCDAT, NBCS, 
     3                   IFLUX, VISCOS(2), ITURB,
     4                   ABC(1,1), ABC(1,2), ABC(1,3), RTEMP,
     5                   AT(1,1), AT(1,2), AT(1,3), AI, DQP)
  310    CONTINUE
      ENDIF
C
C     K implicit sweep
C
      IF ( THREED ) THEN
C
C     K Diagonal Inversion:
C     Use F(1-3), FM(1-3) for lhat, mhat.  Use F(5) for TKE.  FM(5) for VOLJ.
C     Use QM for PROPSJ.  AT(3) for EIGENP, DQP for EIGENM, DQM for EIGEN.
C
         IF ( DIAG(3) ) THEN
            DO 400 I = 2, IDIM
               CALL AF3DK  (IDIM, JDIM, KDIM, I,
     1                   Q, PROPS, NPPTS, S, VOL, DTJ, RES, 
     2                   QP, QM, STEMP, FM(1,5), F(1,1), FM(1,1),F(1,5), 
     3                   IBCDAT, NBCS, IFLUX, VISCOS(3), ITURB,
     4                   ABC(1,1), ABC(1,2), ABC(1,3), RTEMP,
     5                   AT(1,1), AT(1,2), AT(1,3), DQP, DQM) 
  400       CONTINUE
C
C     K Block Inversion:
C     Calculate INPL to be the number of planes that will fit
C     in the memory allocated.  Actually vectorize over NPL 
C     which takes into account if the last loop would be longer
C     then IDIM.  J is the primary vectorization direction and
C     I is the secondary vectorization direction.
C     Make sure that INPL is no larger than the maximum size of the grid.
C
         ELSE
            INPL = MXABC / ( (JDIM-1) * (KDIM-1) )
            INPL = MIN (INPL, (IDIM-1))
            NPL  = INPL
            DO 410 I = 2, IDIM, INPL
               IF ( (I+INPL) .GT. IDIM ) NPL = IDIM - I + 1
               CALL AF3SWK (IDIM, JDIM, KDIM, NPL, I,
     1                      Q, PROPS, NPPTS, S, VOL, DTJ, RES, QP, QM, 
     2                      STEMP, IBCDAT, NBCS, IFLUX, VISCOS(3),ITURB,
     3                      ABC(1,1), ABC(1,2), ABC(1,3), RTEMP,
     5                      AT(1,1), AT(1,2), AT(1,3), AI, DQP)
  410       CONTINUE
         ENDIF
      ENDIF
C
C     Conservative Variable Update
C
      CALL UPDATE (IDIM, JDIM, KDIM, IPSTRT, IPEND, Q, RES, ITURB,
     1             QP, QM, QMIN, NFIXQ, NRELIZ, FHAT, IERRCD)
C
C     Finished with 3 Factor AF
C
      RETURN
      END






      SUBROUTINE AF3SWI (IDIM, JDIM, KDIM, NPL, KSTART, 
     1                   Q, PROPS, NPPTS, S, VOL, DTJ, RES, QP, QM, 
     2                   STEMP, IBCDAT, NBCS, IFLUX, VISCOS, ITURB,
     3                   A, B, C, RTEMP, APT, AMT, ART, AI, WORK)
C
C     Routine to do the I implicit sweep of three factor
C     Approximate Factorisation
C
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     NPL            : Number of IxJ planes to invert simultaneously
C     KSTART         : K to begin the multiple plane work on
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
C     DTJ            : Time step divided by cell volume
C     RES            : Residual - must be over the entire field for
C                      this time integration scheme
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     STEMP          : Storage for metrics passed to Jacobian routine
C                      including the area term
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IFLUX          : Flux evaluation scheme
C     VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C     ITURB          : Turbulence modeling scheme
C     A,B,C          : Tridiagonal systems to invert
C                      A*DQ(I-1) + B*DQ(I) + C*DQ(I+1) = RTEMP
C                      Implicit in the I direction.  Vectorization
C                      by inverting multiple systems at one time to
C                      avoid the recursion of the inversion.  The primary
C                      vectorization direction is J and the secondary 
C                      (for multiple plane vectorization) is K.  These
C                      are stored in the first index for a unit stride
C                      in the inversion.
C     RTEMP          : Storage for the residual array in the form 
C                      appropriate for the inversion routine
C     APT,AMT,ART    : Storage on a plane (IxJ) to calculate the
C                      Jacobians
C     AI             : Identity matrix
C     WORK           : Temporary array of size NPTS*4 for Jacobian calculation
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          DTJ   (JDIM,KDIM,IDIM)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IDIM+1,NF),
     1          QP    (IDIM,2:JDIM,NQ), 
     2          QM    (IDIM,2:JDIM,NQ),
     3          STEMP (IDIM,2:JDIM,4),
     4          WORK  (IDIM,2:JDIM,4)
C
C     Storage for the implicit system to invert
C     Store the system only in the internal portion of
C     the field and use explicit boundary conditions
C     of the form Delta(Q) boundaries = 0.
C
      DIMENSION A     (2:JDIM,NPL,2:IDIM,NRANK,NRANK),
     1          B     (2:JDIM,NPL,2:IDIM,NRANK,NRANK),
     2          C     (2:JDIM,NPL,2:IDIM,NRANK,NRANK),
     3          RTEMP (2:JDIM,NPL,2:IDIM,NF)
C
C     Storage for the Jacobians used to calculate A,B and C
C
      DIMENSION APT   (IDIM,2:JDIM,NRANK,NRANK),
     1          AMT   (IDIM,2:JDIM,NRANK,NRANK),
     2          ART   (IDIM,2:JDIM,NRANK,NRANK)
C
      DIMENSION AI    (NRANK,NRANK)
C
      DIMENSION IBCDAT(10,NBCS)
C
      LOGICAL   VISCOS
C
C     Loop over NPL planes of IxJ, forming the tridiagonal system
C     on on IxJ plane at a time.  Invert the system after all NPL
C     planes have been created.  This reduces the amount of temporary
C     storage required to store QP,QM and the Jacobian matrices
C     but also allows the vectorization of the inversion over 
C     NPL*(JDIM-1) as memory allows.
C
      DO 500 KPL = 1, NPL
         K = KSTART + KPL - 1
C
C     Calculate QP,QM using first order
C
         DO 130 L = 1, NQ
            DO 120 J = 2, JDIM
               DO 110 I = 1, IDIM
                  QP(I,J,L) = Q(J,K,I+1,L)
                  QM(I,J,L) = Q(J,K,I,L)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
C
C     Store the metrics for the Jacobian calculation 
C     including the area term
C
         DO 190 L = 1, 4
            DO 180 J = 2, JDIM
               DO 170 I = 1, IDIM
                  STEMP(I,J,L) = S(J,K,I,L,1)
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
C
C     Calculate Jacobians and set up matrices to invert
C
         NPTS = IDIM * (JDIM - 1)
C
C     Roe's scheme
C          APT = dF/dQ(QP) 
C          AMT = dF/dQ(QM) 
C
         IF (IFLUX .EQ. IFROE) THEN
            CALL DFROE (ITURB, NPTS, QP, QM, STEMP, APT, AMT, ART, WORK)
C
            DO 214 L = 1, NRANK
               DO 213 M = 1, NRANK
                  DO 212 I = 2, IDIM
                     DO 211 J = 2, JDIM
                        DTTERM         = 0.5E0 / 
     1                                 ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                                 + TDPDEL * VOL(J,K,I) / TIMEDT )
                        A(J,KPL,I,L,M) = - DTTERM *
     1                                   (AMT(I-1,J,L,M)+ART(I-1,J,L,M))
                        B(J,KPL,I,L,M) =   AI(L,M) + DTTERM *
     1                                   ( AMT(I,J,L,M)-APT(I-1,J,L,M)
     2                                   + ART(I,J,L,M)+ART(I-1,J,L,M) )
                        C(J,KPL,I,L,M) =   DTTERM *
     1                                   ( APT(I,J,L,M)-ART(I,J,L,M) )
  211                CONTINUE
  212             CONTINUE
  213          CONTINUE
  214       CONTINUE
C
         ELSE
            write (iout,'(a)') ' *** Invalid Flux Function ***'
            stop
         ENDIF
C
C     Implicit Viscous Terms
C          PROPS(2) is RMUT
C          AP = DF(i+1/2)/DQ(i+1)
C          AM = DF(i+1/2)/DQ(i)
C          PROPS(1) is RMU and PROPS(2) is RMUT
C          Use QM as temporary space to store U,V,W,RMUI,RMUTI
C          Use QP as temporary space to store T - being
C          careful as T is dimensioned as IDIM+1, therefore
C          store T in QP(1) and leave QP(2) unused in this call.
C
c        IF ( VISCOS ) THEN
c           CALL DFVTLI (IDIM, JDIM, KDIM, K, 2, IDIM,
c    1                   Q, PROPS(1,1), PROPS(1,2), S, VOL, 
c    2                   IBCDAT, NBCS, 
c    3                   QP(1,2,1), QM(1,2,1), QM(1,2,2), QM(1,2,3), 
c    4                   QM(1,2,4), QM(1,2,5), APT, AMT)
C
c           DO 324 L = 1, NRANK
c              DO 323 M = 1, NRANK
c                 DO 322 I = 2, IDIM
c                    DO 321 J = 2, JDIM
c                       DTTERM       = 1.0E0 / 
c    1                                 ( TDPDTJ / (TDTHET*DTJ(J,K,I))
c    2                                 + TDPDEL * VOL(J,K,I) / TIMEDT )
c                       A(J,KPL,I,L,M) = A(J,KPL,I,L,M)
c    1                    + DTTERM * AMT(I-1,J,L,M)
c                       B(J,KPL,I,L,M) = B(J,KPL,I,L,M)
c    1                    + DTTERM * 
c    2                            ( APT(I,J,L,M)-AMT(I-1,J,L,M) )
c                       C(J,KPL,I,L,M) = C(J,KPL,I,L,M) 
c    1                    - DTTERM * APT(I,J,L,M)
c 321                CONTINUE
c 322             CONTINUE
c 323          CONTINUE
c 324       CONTINUE
c        ENDIF
C
C     Store Residual array for inversion
C
         DO 490 L = 1, NF
            DO 480 I = 2, IDIM
               DO 470 J = 2, JDIM
                  RTEMP(J,KPL,I,L) = RES(J,K,I,L)
  470          CONTINUE
  480       CONTINUE
  490    CONTINUE
C
  500 CONTINUE
C
C     Invert the block tridiagonal system:
C     A*DQ(I-1) + B*DQ(I) + C*DQ(I+1) = RTEMP
C
      NPTS = IDIM - 1
      NSIM = NPL * (JDIM - 1)
      CALL VBTRI (NPTS, NSIM, NRANK, A, B, C, RTEMP)
C
C     Update the residual array
C
      DO 890 KPL = 1, NPL
         K = KSTART + KPL - 1
         DO 880 L = 1, NF
            DO 870 I = 2, IDIM
               DO 860 J = 2, JDIM
                  RES(J,K,I,L) = RTEMP(J,KPL,I,L)
  860          CONTINUE
  870       CONTINUE
  880    CONTINUE
  890 CONTINUE
C
C     Finished with the I implicit sweep of three factor AF
C
      RETURN
      END






      SUBROUTINE AF3SWJ (IDIM, JDIM, KDIM, NPL, ISTART, 
     1                   Q, PROPS, NPPTS, S, VOL, DTJ, RES, QP, QM, 
     2                   STEMP, IBCDAT, NBCS, IFLUX, VISCOS, ITURB,
     3                   A, B, C, RTEMP, APT, AMT, ART, AI, WORK)
C
C     Routine to do the J implicit sweep of three factor
C     Approximate Factorisation
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     NPL            : Number of JxK planes to invert simultaneously
C     ISTART         : I to begin the multiple plane work on
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
C     DTJ            : Time step divided by cell volume
C     RES            : Residual - must be over the entire field for
C                      this time integration scheme
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     STEMP          : Storage for metrics passed to Jacobian routine
C                      including the area term
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IFLUX          : Flux evaluation scheme
C     VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C     ITURB          : Turbulence modeling scheme
C     A,B,C          : Tridiagonal systems to invert
C                      A*DQ(J-1) + B*DQ(J) + C*DQ(J+1) = RTEMP
C                      Implicit in the J direction.  Vectorization
C                      by inverting multiple systems at one time to
C                      avoid the recursion of the inversion.  The primary
C                      vectorization direction is I and the secondary 
C                      (for multiple plane vectorization) is K.  These
C                      are stored in the first index for a unit stride
C                      in the inversion.
C     RTEMP          : Storage for the residual array in the form 
C                      appropriate for the inversion routine
C     APT,AMT,ART    : Storage on a plane (IxJ) to calculate the
C                      Jacobians
C     AI             : Identity matrix 
C     WORK           : Temporary array of size NPTS*4 for Jacobian calculation
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          DTJ   (JDIM,KDIM,IDIM)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IDIM+1,NF),
     1          QP    (JDIM,2:KDIM,NQ), 
     2          QM    (JDIM,2:KDIM,NQ),
     3          STEMP (JDIM,2:KDIM,4),
     4          WORK  (JDIM,2:KDIM,4)
C
C     Storage for the implicit system to invert
C     Store the system only in the internal portion of
C     the field and use explicit boundary conditions
C     of the form Delta(Q) boundaries = 0.
C
      DIMENSION A     (2:KDIM,NPL,2:JDIM,NRANK,NRANK),
     1          B     (2:KDIM,NPL,2:JDIM,NRANK,NRANK),
     2          C     (2:KDIM,NPL,2:JDIM,NRANK,NRANK),
     3          RTEMP (2:KDIM,NPL,2:JDIM,NF)
C
C     Storage for the Jacobians used to calculate A,B and C
C
      DIMENSION APT   (JDIM,2:KDIM,NRANK,NRANK),
     1          AMT   (JDIM,2:KDIM,NRANK,NRANK),
     2          ART   (JDIM,2:KDIM,NRANK,NRANK)
C
      DIMENSION AI    (NRANK,NRANK)
C
      DIMENSION IBCDAT(10,NBCS)
C
      LOGICAL   VISCOS
C
C     Loop over NPL planes of JxK, forming the tridiagonal system
C     on on JxK plane at a time.  Invert the system after all NPL
C     planes have been created.  This reduces the amount of temporary
C     storage required to store QP,QM and the Jacobian matrices
C     but also allows the vectorization of the inversion over 
C     NPL*(KDIM-1) as memory allows.
C
      DO 500 IPL = 1, NPL
         I = ISTART + IPL - 1
C
C     Calculate QP,QM using first order
C
         DO 130 L = 1, NQ
            DO 120 K = 2, KDIM
               DO 110 J = 1, JDIM
                  QP(J,K,L) = Q(J+1,K,I,L)
                  QM(J,K,L) = Q(J,K,I,L)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
C
C     Store the metrics for the Jacobian calculation
C     including the area term
C
         DO 190 L = 1, 4
            DO 180 K = 2, KDIM
               DO 170 J = 1, JDIM
                  STEMP(J,K,L) = S(J,K,I,L,2)
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
C
C     Calculate Jacobians and set up matrices to invert
C
         NPTS = JDIM * (KDIM - 1)
C
C     Roe's scheme
C          APT = dF/dQ(QP)
C          AMT = dF/dQ(QM)
C
         IF (IFLUX .EQ. IFROE) THEN
            CALL DFROE (ITURB, NPTS, QP, QM, STEMP, APT, AMT, ART, WORK)
C
            DO 214 L = 1, NRANK
               DO 213 M = 1, NRANK
                  DO 212 K = 2, KDIM
                     DO 211 J = 2, JDIM
                        DTTERM         = 0.5E0 / 
     1                                 ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                                 + TDPDEL * VOL(J,K,I) / TIMEDT )
                        A(K,IPL,J,L,M) = - DTTERM *
     1                                   (AMT(J-1,K,L,M)+ART(J-1,K,L,M))
                        B(K,IPL,J,L,M) =   AI(L,M) + DTTERM *
     1                                   ( AMT(J,K,L,M)-APT(J-1,K,L,M)
     2                                   + ART(J,K,L,M)+ART(J-1,K,L,M) )
                        C(K,IPL,J,L,M) =   DTTERM *
     1                                   ( APT(J,K,L,M)-ART(J,K,L,M) )
  211                CONTINUE
  212             CONTINUE
  213          CONTINUE
  214       CONTINUE
C
         ELSE
            write (iout,'(a)') ' *** Invalid Flux Function ***'
            stop
         ENDIF
C
C     Implicit Viscous Terms
C
C          PROPS(2) is RMUT
C          AP = DF(j+1/2)/DQ(j+1)
C          AM = DF(j+1/2)/DQ(j)
C          PROPS(1) is RMU and PROPS(2) is RMUT
C          Use QM as temporary space to store U,V,W,RMUJ,RMUTJ
C          Use QP as temporary space to store T - being
C          careful as T is dimensioned as JDIM+1, therefore
C          store T in QP(1) and leave QP(2) unused in this call.
C
         IF ( VISCOS ) THEN
            CALL DFVISJ (IDIM, JDIM, KDIM, I, ITURB,
     1                   Q, PROPS(1,1), PROPS(1,2), S, VOL, 
     2                   IBCDAT, NBCS,
     3                   QP(1,2,1), QM(1,2,1), QM(1,2,2), QM(1,2,3),
     4                   QM(1,2,4), QM(1,2,5), APT, AMT)
C
            DO 324 L = 1, NRANK
               DO 323 M = 1, NRANK
                  DO 322 K = 2, KDIM
                     DO 321 J = 2, JDIM
                        DTTERM         = 1.0E0 / 
     1                                 ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                                 + TDPDEL * VOL(J,K,I) / TIMEDT )
                        A(K,IPL,J,L,M) = A(K,IPL,J,L,M)
     1                    + DTTERM * AMT(J-1,K,L,M)
                        B(K,IPL,J,L,M) = B(K,IPL,J,L,M)
     1                    + DTTERM * ( APT(J,K,L,M)-AMT(J-1,K,L,M) )
                        C(K,IPL,J,L,M) = C(K,IPL,J,L,M) 
     1                    - DTTERM * APT(J,K,L,M)
  321                CONTINUE
  322             CONTINUE
  323          CONTINUE
  324       CONTINUE
         ENDIF
C
C     Store Residual array for inversion
C
         DO 490 L = 1, NF
            DO 480 K = 2, KDIM
               DO 470 J = 2, JDIM
                  RTEMP(K,IPL,J,L) = RES(J,K,I,L)
  470          CONTINUE
  480       CONTINUE
  490    CONTINUE
C
  500 CONTINUE
C
C     Invert the block tridiagonal system:
C     A*DQ(J-1) + B*DQ(J) + C*DQ(J+1) = RTEMP
C
      NPTS = JDIM - 1
      NSIM = NPL * (KDIM - 1)
      CALL VBTRI (NPTS, NSIM, NRANK, A, B, C, RTEMP)
C
C     Update the residual array
C
      DO 890 IPL = 1, NPL
         I = ISTART + IPL - 1
         DO 880 L = 1, NF
            DO 870 K = 2, KDIM
               DO 860 J = 2, JDIM
                  RES(J,K,I,L) = RTEMP(K,IPL,J,L)
  860          CONTINUE
  870       CONTINUE
  880    CONTINUE
  890 CONTINUE
C
C     Finished with the J implicit sweep of three factor AF
C
      RETURN
      END






      SUBROUTINE AF3SWK (IDIM, JDIM, KDIM, NPL, ISTART, 
     1                   Q, PROPS, NPPTS, S, VOL, DTJ, RES, QP, QM, 
     2                   STEMP, IBCDAT, NBCS, IFLUX, VISCOS, ITURB,
     3                   A, B, C, RTEMP, APT, AMT, ART, AI, WORK)
C
C     Routine to do the K implicit sweep of three factor
C     Approximate Factorisation
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     NPL            : Number of KxJ planes to invert simultaneously
C     ISTART         : I to begin the multiple plane work on
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
C     DTJ            : Time step divided by cell volume
C     RES            : Residual - must be over the entire field for
C                      this time integration scheme
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     STEMP          : Storage for metrics passed to Jacobian routine
C                      including the area term
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IFLUX          : Flux evaluation scheme
C     VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C     ITURB          : Turbulence modeling scheme
C     A,B,C          : Tridiagonal systems to invert
C                      A*DQ(K-1) + B*DQ(K) + C*DQ(K+1) = RTEMP
C                      Implicit in the I direction.  Vectorization
C                      by inverting multiple systems at one time to
C                      avoid the recursion of the inversion.  The primary
C                      vectorization direction is J and the secondary 
C                      (for multiple plane vectorization) is I.  These
C                      are stored in the first index for a unit stride
C                      in the inversion.
C     RTEMP          : Storage for the residual array in the form 
C                      appropriate for the inversion routine
C     APT,AMT,ART    : Storage on a plane (KxJ) to calculate the
C                      Jacobians
C     AI             : Identity matrix 
C     WORK           : Temporary array of size NPTS*4 for Jacobian calculation
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          DTJ   (JDIM,KDIM,IDIM)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IDIM+1,NF),
     1          QP    (KDIM,2:JDIM,NQ), 
     2          QM    (KDIM,2:JDIM,NQ),
     3          STEMP (KDIM,2:JDIM,4),
     4          WORK  (KDIM,2:JDIM,4)
C
C     Storage for the implicit system to invert
C     Store the system only in the internal portion of
C     the field and use explicit boundary conditions
C     of the form Delta(Q) boundaries = 0.
C
      DIMENSION A     (2:JDIM,NPL,2:KDIM,NRANK,NRANK),
     1          B     (2:JDIM,NPL,2:KDIM,NRANK,NRANK),
     2          C     (2:JDIM,NPL,2:KDIM,NRANK,NRANK),
     3          RTEMP (2:JDIM,NPL,2:KDIM,NF)
C
C     Storage for the Jacobians used to calculate A,B and C
C
      DIMENSION APT   (KDIM,2:JDIM,NRANK,NRANK),
     1          AMT   (KDIM,2:JDIM,NRANK,NRANK),
     2          ART   (KDIM,2:JDIM,NRANK,NRANK)
C
      DIMENSION AI    (NRANK,NRANK)
C
      DIMENSION IBCDAT(10,NBCS)
C
      LOGICAL   VISCOS
C
C     Loop over NPL planes of KxJ, forming the tridiagonal system
C     on on KxJ plane at a time.  Invert the system after all NPL
C     planes have been created.  This reduces the amount of temporary
C     storage required to store QP,QM and the Jacobian matrices
C     but also allows the vectorization of the inversion over 
C     NPL*(JDIM-1) as memory allows.
C
      DO 500 IPL = 1, NPL
         I = ISTART + IPL - 1
C
C     Calculate QP,QM using first order
C
         DO 130 L = 1, NQ
            DO 120 J = 2, JDIM
               DO 110 K = 1, KDIM
                  QP(K,J,L) = Q(J,K+1,I,L)
                  QM(K,J,L) = Q(J,K,I,L)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
C
C     Store the metrics for the Jacobian calculation
C     including the area term
C
         DO 190 L = 1, 4
            DO 180 J = 2, JDIM
               DO 170 K = 1, KDIM
                  STEMP(K,J,L) = S(J,K,I,L,3)
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
C
C     Calculate Jacobians and set up matrices to invert
C
         NPTS = KDIM * (JDIM - 1)
C
C     Roe's scheme
C          APT = dF/dQ(QP)
C          AMT = dF/dQ(QM)
C
         IF (IFLUX .EQ. IFROE) THEN
            CALL DFROE (ITURB, NPTS, QP, QM, STEMP, APT, AMT, ART, WORK)
C
            DO 214 L = 1, NRANK
               DO 213 M = 1, NRANK
                  DO 212 K = 2, KDIM
                     DO 211 J = 2, JDIM
                        DTTERM         = 0.5E0 / 
     1                                 ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                                 + TDPDEL * VOL(J,K,I) / TIMEDT )
                        A(J,IPL,K,L,M) = - DTTERM *
     1                                   (AMT(K-1,J,L,M)+ART(K-1,J,L,M))
                        B(J,IPL,K,L,M) =   AI(L,M) + DTTERM *
     1                                   ( AMT(K,J,L,M)-APT(K-1,J,L,M)
     2                                   + ART(K,J,L,M)+ART(K-1,J,L,M) )
                        C(J,IPL,K,L,M) =   DTTERM *
     1                                   ( APT(K,J,L,M)-ART(K,J,L,M) )
  211                CONTINUE
  212             CONTINUE
  213          CONTINUE
  214       CONTINUE
C
         ELSE
            write (iout,'(a)') ' *** Invalid Flux Function ***'
            stop
         ENDIF
C
C     Implicit Viscous Terms
C
C          PROPS(2) is RMUT
C          AP = DF(i+1/2)/DQ(i+1)
C          AM = DF(i+1/2)/DQ(i)
C          PROPS(1) is RMU and PROPS(2) is RMUT
C          Use QM as temporary space to store U,V,W,RMUK,RMUTK
C          Use QP as temporary space to store T - being
C          careful as T is dimensioned as KDIM+1, therefore
C          store T in QP(1) and leave QP(2) unused in this call.
C
         IF ( VISCOS ) THEN
            CALL DFVISK (IDIM, JDIM, KDIM, I, ITURB,
     1                   Q, PROPS(1,1), PROPS(1,2), S, VOL, 
     2                   IBCDAT, NBCS,
     3                   QP(1,2,1), QM(1,2,1), QM(1,2,2), QM(1,2,3),
     4                   QM(1,2,4), QM(1,2,5), APT, AMT)
C
            DO 324 L = 1, NRANK
               DO 323 M = 1, NRANK
                  DO 322 K = 2, KDIM
                     DO 321 J = 2, JDIM
                        DTTERM         = 1.0E0 / 
     1                                 ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                                 + TDPDEL * VOL(J,K,I) / TIMEDT )
                        A(J,IPL,K,L,M) = A(J,IPL,K,L,M)
     1                    + DTTERM * AMT(K-1,J,L,M)
                        B(J,IPL,K,L,M) = B(J,IPL,K,L,M)
     1                    + DTTERM * ( APT(K,J,L,M)-AMT(K-1,J,L,M) )
                        C(J,IPL,K,L,M) = C(J,IPL,K,L,M) 
     1                    - DTTERM * APT(K,J,L,M)
  321                CONTINUE
  322             CONTINUE
  323          CONTINUE
  324       CONTINUE
         ENDIF
C
C     Store Residual array for inversion
C
         DO 490 L = 1, NF
            DO 480 K = 2, KDIM
               DO 470 J = 2, JDIM
                  RTEMP(J,IPL,K,L) = RES(J,K,I,L)
  470          CONTINUE
  480       CONTINUE
  490    CONTINUE
C
  500 CONTINUE
C
C     Invert the block tridiagonal system:
C     A*DQ(K-1) + B*DQ(K) + C*DQ(K+1) = RTEMP
C
      NPTS = KDIM - 1
      NSIM = NPL * (JDIM - 1)
      CALL VBTRI (NPTS, NSIM, NRANK, A, B, C, RTEMP)
C
C     Update the residual array
C
      DO 890 IPL = 1, NPL
         I = ISTART + IPL - 1
         DO 880 L = 1, NF
            DO 870 K = 2, KDIM
               DO 860 J = 2, JDIM
                  RES(J,K,I,L) = RTEMP(J,IPL,K,L)
  860          CONTINUE
  870       CONTINUE
  880    CONTINUE
  890 CONTINUE
C
C     Finished with the K implicit sweep of three factor AF
C
      RETURN
      END





      SUBROUTINE PISRC (IDIM, JDIM, KDIM, NPL, ISTART, IMODEL, IFDELQ,
     1                  Q, PROPS, NPPTS, S, VOL, RC, DTJ, RES,
     2                  MXSECT, DQP, DQM, QP, QM, TAU, FM, STEMP,
     3                  TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ,
     4                  SRC, IBCDAT, NBCS, B, RTEMP, AI, AT)
C
C Routine to perform the Point Implicit Source term treatment
C This routine is organized as follows:
C     1. Diagonal implicit source term treatment
C     2. Block    implicit source term treatment
C     3. Error
C
C IDIM,JDIM,KDIM : Dimensions of current block
C NPL            : Number of JxK planes to invert simultaneously
C ISTART         : I to begin the multiple plane work on
C IMODEL         : Turbulence model array
C IFDELQ         : Array controlling calculation of delq
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
C RC             : Position vector at cell centers
C DTJ            : Time step divided by cell volume
C RES            : Residual - must be over the entire field for
C                  this time integration scheme
C MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                  (used in dimensioning of temporary space used
C                   for calculating q(+,-) and f.
C DQP,DQM        : Storage for delta q = q(j) - q(j-1)
C QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C TAU            : Storage for Reynolds stresses
C FM             : Fluxes at cell interfaces - currently temp storage
C STEMP          : Temporary storage for metrics - currently temp stor.
C TKE            : Turbulent kinetic energy
C ETILDE         : Epsilon^tilde
C XSIW           :
C FDAMP          : Damping functions F1, F2
C EPSC           :
C PDIL           :
C SIJ            :
C WIJ            :
C SRC            : Source terms at cell (J,K)
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C B              : Diagonal term = I + DTJ D (D = dS/dQ)
C AI             : Identity matrix
C AT             : Temporary storage of implicit source terms
C
      include '../header/common.h'
C
      DIMENSION IMODEL(NMDL),
     1          IFDELQ(NQ)
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          RC    (JDIM+1,KDIM+1,IDIM+1,3),
     5          DTJ   (JDIM,KDIM,IDIM)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IDIM+1,NF)
C
      DIMENSION DQP   (MXSECT,NQ),
     1          DQM   (MXSECT,NQ),
     2          QP    (MXSECT,NQ),
     3          QM    (MXSECT,NQ),
     4          TAU   (MXSECT,6),
     5          FM    (MXSECT,NQ),
     6          STEMP (MXSECT,4),
     7          TKE   (MXSECT),
     8          ETILDE(MXSECT),
     9          XSIW  (MXSECT),
     A          FDAMP (MXSECT,2),
     B          EPSC  (MXSECT),
     C          PDIL  (MXSECT),
     D          SIJ   (MXSECT,6),
     E          WIJ   (MXSECT,6),
     F          SRC   (2:JDIM,2:KDIM,NF)
C
      DIMENSION B     (2:KDIM,NPL,2:JDIM,NRANK,NRANK),
     1          RTEMP (2:KDIM,NPL,2:JDIM,NF)
C
      DIMENSION AI    (NRANK,NRANK)
C
      DIMENSION AT    (2:JDIM,2:KDIM,NRANK,NRANK)
C
      DIMENSION IBCDAT(10,NBCS)
C
C Local variable:
C
      LOGICAL   ISAVPR
C
      ITURB  = IMODEL(1)
C
C Do not store CMUSTAR from local array back into global array when doing
C implicit source terms.  This has been done in RESID.
C
      ISAVPR = .FALSE.
C
C 1.  Diagonal implicit source term treatment
C
c     impsrc = 0
      IF (IMPSRC .EQ. ISDIAG) THEN
         NPRLIM = 0
         NPRNEG = 0
         DO 200 IPL = 1, NPL
            I = ISTART + IPL - 1
C
C     Temporary variables passed to the source term evaluation routine:
C         FM(2-) for PROPSC (properties in one-dimensional array - this is
C                            passed in using FM.  Must always be sure that
C                            NP+1 <= NQ or else will need to fix this)
C         STEMP      for SN   (Wall normals - needed in Reynolds stress model)
C         QM         for DQDX (derivatives of Q in x direction at cell centers)
C         DQP        for DQDY (derivatives of Q in y direction at cell centers)
C         DQM        for DQDZ (derivatives of Q in z direction at cell centers)
C         FM(1)      for DEL2K (Laplacian of turbulent kinetic energy)
C
            CALL SOURCE (IDIM, JDIM, KDIM, I, NPRLIM, NPRNEG,
     1                   Q, PROPS, NPPTS, S, VOL, RC,
     2                   QP, FM(1,2), STEMP, QM, DQP, DQM, FM(1,1),
     3                   TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ,
     4                   TAU, SRC, IBCDAT, NBCS, 
     5                   IFLUX, IMODEL, IFDELQ, ISAVPR)
C
C     Divide residual by 1 - delta(t)*min(source/(rho phi), 0) as the
C     diagonal inversion of the point implicit treatment (after work
C     of Huang and Spallart,Allmaras)
C
            LSTRT = 6
            DO 190 L = LSTRT, NQ
               DO 180 K = 2, KDIM
                  DO 170 J = 2, JDIM
                     DTTERM = 1.0E0 / 
     1                        ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                        + TDPDEL * VOL(J,K,I) / TIMEDT )
                     PHI    = Q(J,K,I,L)
                     SGNPHI = SIGN (1.E0, PHI)
                     RHOPHI = Q(J,K,I,1) * PHI + SGNPHI * RSMALL
c--                  BDIAG  = MIN (SRC(J,K,L)*VOL(J,K,I)/RHOPHI, 0.E0)
                     BDIAG  = - ABS ( SRC(J,K,L)*VOL(J,K,I)/RHOPHI )
                     DENOM  = 1.E0 - DTTERM * BDIAG
                     RES(J,K,I,L) = RES(J,K,I,L) / DENOM
  170             CONTINUE
  180          CONTINUE
  190       CONTINUE
  200    CONTINUE
C
C 2.  Block implicit source terms
C     Loop over NPL planes of JxK
C
      ELSE IF (IMPSRC .EQ. ISBLOC) THEN
         DO 500 IPL = 1, NPL
            I = ISTART + IPL - 1
C
C     Temporary variables passed to the source term evaluation routine:
C         FM(2-) for PROPSC (properties in one-dimensional array - this is
C                            passed in using FM.  Must always be sure that
C                            NP+1 <= NQ or else will need to fix this)
C         STEMP      for SN   (Wall normals - needed in Reynolds stress model)
C         QM         for DQDX (derivatives of Q in x direction at cell centers)
C         DQP        for DQDY (derivatives of Q in y direction at cell centers)
C         DQM        for DQDZ (derivatives of Q in z direction at cell centers)
C         FM(1)      for DEL2K (Laplacian of turbulent kinetic energy)
C
            CALL SOURCE (IDIM, JDIM, KDIM, I, NPRLIM, NPRNEG,
     1                   Q, PROPS, NPPTS, S, VOL, RC,
     2                   QP, FM(1,2), STEMP, QM, DQP, DQM, FM(1,1),
     3                   TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ,
     4                   TAU, SRC, IBCDAT, NBCS, 
     5                   IFLUX, IMODEL, IFDELQ, ISAVPR)
C
C     Add the implicit source terms to the identity matrix 
C     Temporary variables passed to the source term evaluation routine:
C         QP     for QC (variables in one-dimensional array)
C         FM(2-) for PROPSC (properties in one-dimensional array - this is
C                            passed in using FM.  Must always be sure that
C                            NP+1 <= NQ or else will need to fix this)
C         QM         for DQDX (derivatives of Q in x direction at cell centers)
C         DQP        for DQDY (derivatives of Q in y direction at cell centers)
C         DQM        for DQDZ (derivatives of Q in z direction at cell centers)
C         FM(1)      for DEL2K (Laplacian of turbulent kinetic energy)
C
            CALL DSRCDQ (IDIM, JDIM, KDIM, I, IMODEL, 
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   QP, FM(1,2), QM, DQP, DQM, FM(1,1), 
     3                   TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ,
     4                   TAU, SRC, IBCDAT, NBCS, AT)
C
            DO 340 L = 1, NRANK
               DO 330 M = 1, NRANK
                  DO 320 K = 2, KDIM 
                     DO 310 J = 2, JDIM 
                        DTTERM = 1.0E0 / 
     1                           ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                           + TDPDEL * VOL(J,K,I) / TIMEDT )
                        B(K,IPL,J,L,M) = AI(L,M) - DTTERM * AT(J,K,L,M) 
  310                CONTINUE
  320             CONTINUE
  330          CONTINUE
  340       CONTINUE
C
C     Store residual array for inversion
C
            DO 490 L = 1, NF
               DO 480 K = 2, KDIM
                  DO 470 J = 2, JDIM
                     RTEMP(K,IPL,J,L) = RES(J,K,I,L)
  470             CONTINUE
  480          CONTINUE
  490       CONTINUE
C
  500    CONTINUE
C
C     Finished forming B = I - DTJ D(S)
C
C     Invert the diagonal system:
C          B * DQ = RTEMP
C
         NPTS = JDIM - 1
         NSIM = NPL * (KDIM - 1)
         CALL VBDIA (NPTS, NSIM, NRANK, B, RTEMP)
C
C     Update the residual array
C
         DO 890 IPL = 1, NPL
            I = ISTART + IPL - 1 
            DO 880 L = 1, NF
               DO 870 K = 2, KDIM
                  DO 860 J = 2, JDIM
                     RES(J,K,I,L) = RTEMP(K,IPL,J,L)
  860             CONTINUE
  870          CONTINUE
  880       CONTINUE
  890    CONTINUE
C
C 3.  Error
C
      ELSE
         WRITE (IOUT,1000) IMPSRC
         STOP
      ENDIF
C
C     Finished with the implicit source term inversion
C
      RETURN
 1000 FORMAT (' ','PISRC : ERROR-> Invalid implicit source terms ',
     1            'specified = ',I10,
     2       /' ', 8X,'RUN ABORTING!')
      END
