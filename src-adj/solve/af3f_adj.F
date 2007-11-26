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
      SUBROUTINE AF3F_ADJ (IBLK, IDIM, JDIM, KDIM, MXSECT, MXABC,
     1                 Q, QN, QA, NQNPTS, PROPS, NPPTS, S, VOL, RC, DTJ,
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
      include '../../src-flo/header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QN    (NQNPTS,NQ,NTMLVL),
     2          PROPS (NPPTS,NP),
     3          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     4          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     5          RC    (JDIM+1,KDIM+1,IDIM+1,3),
     6          DTJ   (JDIM,KDIM,IDIM),
     7          QA    (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION RES   (0:JDIM+2,0:KDIM+2,0:IDIM+2,NF),
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
            STOP "PISRC not finished for adjoint"
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
            CALL AF3DI_ADJ  (IDIM, JDIM, KDIM, K,
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
         STOP "AF3SWI not finished for adjoint"
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
            CALL AF3DJ_ADJ  (IDIM, JDIM, KDIM, K,
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
         STOP "AF3SWJ not finished for adjoint"
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
               CALL AF3DK_ADJ  (IDIM, JDIM, KDIM, I,
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
            STOP "AF3SWK not finished for adjoint"
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
C     Adjoint Variable Update
C
      DO J=2,JDIM
         DO K=2,KDIM
            DO I=2,IDIM
               QA(J,K,I,:) = QA(J,K,I,:) + RES(J,K,I,:)
            ENDDO
         ENDDO
      ENDDO
C
C     Finished with 3 Factor AF
C
      RETURN
      END
