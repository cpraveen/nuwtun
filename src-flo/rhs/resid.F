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
      SUBROUTINE RESID (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                  Q, QN, NQNPTS, PROPS, NPPTS, S, VOL, RC, DTJ,
     2                  RES, FHAT, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                  WORKX, IBCDAT, NBCS, IORDER, LIMIT, RKAP, IFLUX,
     4                  VISCOS, IMODEL, INITRS, IFDELQ, MRCHNG, 
     5                  VGNRN, SIGMA, NPRLIM, NPRNEG, IFL2DQ, R2DQDT)
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
      IF (INITRS) THEN
         DO 140 L = 1, NF
            DO 130 I = IPSTRT - 1, IPEND + 1
               DO 120 K = 1, KDIM + 1
                  DO 110 J = 1, JDIM + 1
                     RES(J,K,I,L) = 0.E0
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
      ENDIF
C
C 2.  If iterative time dependent then add dQ/dt term and
C     calculate L2 Norm of dQ/dt residual
C
      IF (ITIMED .EQ. IUNSTD) THEN
         CALL DQDT   (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                Q, QN, VOL, DTJ, RES, QP, QM, IMODEL)
         IF (IFL2DQ) THEN
            CALL L2NORM (IDIM, JDIM, KDIM, IPSTRT, IPEND,
     1                   RES, VOL, R2DQDT)

         ENDIF
      ENDIF
C
C 3.  Residual in I direction - Add dF/dx term
C         FLUXM if marching calculation
C         FLUXI if time-dependent calculation
C
      IF ( MRCHNG ) THEN
         CALL FLUXM (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1               Q, PROPS, NPPTS, S, VOL, 
     2               RES, DQP, DQM, QP, QM, F, FM, STEMP,
     3               IBCDAT, NBCS, IORDER(1), LIMIT(1), 
     4               VGNRN, SIGMA, IMODEL)
      ELSE
         IF (FOURTH) THEN
            CALL FLUXI4 (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   RES, FHAT, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                   IBCDAT, NBCS, IORDER(1), LIMIT(1), RKAP(1), 
     4                   IFLUX, VISCOS(1), IMODEL)
         ELSE
            CALL FLUXI  (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   RES, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                   IBCDAT, NBCS, IORDER(1), LIMIT(1), RKAP(1), 
     4                   IFLUX, VISCOS(1), IMODEL)
         ENDIF
      ENDIF
C
C 4.  Residual in J direction - Add dG/dy term
C
      IF (FOURTH) THEN
         CALL FLUXJ4 (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                Q, PROPS, NPPTS, S, VOL, 
     2                RES, FHAT, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                IBCDAT, NBCS, IORDER(2), LIMIT(2), RKAP(2), 
     4                IFLUX, VISCOS(2), IMODEL)
      ELSE
         CALL FLUXJ  (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                Q, PROPS, NPPTS, S, VOL, 
     2                RES, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                IBCDAT, NBCS, IORDER(2), LIMIT(2), RKAP(2), 
     4                IFLUX, VISCOS(2), IMODEL)
      ENDIF
C
C 5.  If 3D or axisymmetric then Residual in K direction - add dH/dz term
C
      IF ( THREED ) THEN
         IF (FOURTH) THEN
            CALL FLUXK4 (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   RES, FHAT, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                   IBCDAT, NBCS, IORDER(3), LIMIT(3), RKAP(3), 
     4                   IFLUX, VISCOS(3), IMODEL)
         ELSE
            CALL FLUXK  (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   RES, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                   IBCDAT, NBCS, IORDER(3), LIMIT(3), RKAP(3), 
     4                   IFLUX, VISCOS(3), IMODEL)
         ENDIF
      ELSE IF ( AXISYM ) THEN
         CALL FLUXKA (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                Q, PROPS, NPPTS, S, VOL, 
     2                RES, DQP, DQM, QP, QM, TAU, F, FM, STEMP,
     3                IBCDAT, NBCS, IORDER(3), LIMIT(3), RKAP(3), 
     4                IFLUX, VISCOS(3), IMODEL)
      ENDIF
C
C 6.  Calculate Source Terms
C     Temporary variables passed to the source term evaluation routine:
C         QP     for QC (variables in one-dimensional array)
C         FM(2-) for PROPSC (properties in one-dimensional array - this is 
C                            passed in using FM.  Must always be sure that 
C                            NP+1 <= NF or else will need to fix this)
C         STEMP        for SN   (Wall normals - needed in Reynolds stress model)
C         QM           for DQDX (derivatives of Q wrt x at cell centers)
C         DQP          for DQDY (derivatives of Q wrt y at cell centers)
C         DQM          for DQDZ (derivatives of Q wrt z at cell centers)
C         FM(1)        for DEL2K (Laplacian of turbulent kinetic energy)
C         WORKX(1)     for TKE
C         WORKX(2)     for ETILDE
C         WORKX(3)     for XSIW
C         WORKX(4-5)   for FDAMP
C         WORKX(6)     for EPSC
C         WORKX(7)     for PDIL
C         WORKX( 8-13) for SIJ
C         WORKX(14-19) for WIJ
C
      IF ( ISRCE ) THEN
         IF (ITURB.EQ.ITKE .OR. ITURB.EQ.ITKW .OR. ITURB.EQ.ITRS) THEN
            IF (NF .LT. (NP + 1)) THEN
               WRITE (IOUT,1000) NP+1, NF
               STOP
            ENDIF
C
            CALL SRCRES (IDIM, JDIM, KDIM, IPSTRT, IPEND,
     1                   Q, PROPS, NPPTS, S, VOL, RC, RES,
     2                   QP, FM(1,2), STEMP, QM, DQP, DQM, FM(1,1), 
     3                   WORKX(1,1), WORKX(1,2), WORKX(1,3),WORKX(1,4),
     4                   WORKX(1,6), WORKX(1,7), WORKX(1,8),WORKX(1,14),
     5                   TAU, F, IBCDAT, NBCS, IFLUX, IMODEL, IFDELQ,
     6                   NPRLIM, NPRNEG)
         ENDIF
      ENDIF
C
C Wall Functions: reset residuals to zero on the turbulence quantities that
C are set from the wall functions (epsilon is set from wall functions, all
C other turbulence quantities, tau_ij or k, are calculated at cell center).
C The assumption in the wall functions is that we replace the calculated
C value of the turbulence quantities with an 'exact' value.  This
C then implies that the residual at these points is zero.  Necessary as
C the implicit time integration couples the mean and turbulent variables.
C
      IF ( IFWALF ) THEN
C
C Set up equation indices to reset residual to zero
C
         IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
            LSTRT = 6
            LEND  = IEPSLN
         ELSE IF (ITURB .EQ. ITRS) THEN
            LSTRT = IEPSLN
            LEND  = IEPSLN
         ELSE
            LSTRT = 1
            LEND  = 0
         ENDIF
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
         IJKDIM(1) = IDIM
         IJKDIM(2) = JDIM
         IJKDIM(3) = KDIM
C
C Loop over all boundary conditions
C
         DO 490 IBC = 1, NBCS
            IBCTYP = IBCDAT(1,IBC)
            IF (IBCTYP .EQ. IWALFN) THEN
               IDIR      = IBCDAT(2,IBC)
               INROUT    = IBCDAT(3,IBC)
               ISTRT(1)  = MAX (IBCDAT(4,IBC), IPSTRT)
               IEND (1)  = MIN (IBCDAT(5,IBC), IPEND)
               ISTRT(2)  = IBCDAT(6,IBC)
               IEND (2)  = IBCDAT(7,IBC)
               ISTRT(3)  = IBCDAT(8,IBC)
               IEND (3)  = IBCDAT(9,IBC)
               IF (INROUT .EQ. 1) THEN
                  IBCN = 2
               ELSE
                  IBCN = IJKDIM(IDIR)
               ENDIF
C
               IC1 = ICY (IDIR, 1)
               IC2 = ICY (IDIR, 2)
               IC3 = ICY (IDIR, 3)
C
               DO 430 L = LSTRT, LEND
                  DO 420 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
                     DO 410 I2 = ISTRT(IC2), IEND(IC2)
C
                        IN = IBCN*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                        JN = IBCN*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                        KN = IBCN*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                        RES(JN,KN,IN,L) = 0.E0
C
  410                CONTINUE
  420             CONTINUE
  430          CONTINUE
            ENDIF
  490    CONTINUE
      END IF
C
C     Finished with residual computation.
C
      RETURN
 1000 FORMAT (' ','***** FATAL ERROR - STOP in routine RESID ',
     1       /' ','Temporary storage required to evaluate source terms',
     2       /' ','is insufficient.  Contact ISAAC developer to fix.',
     3       /' ','NP+1 = ', I5, ' NF = ', I5)
      END



      SUBROUTINE SAVEPR (IDIM, JDIM, KDIM, I, PROPSC, PROPS)
C
C Routine to store CMUSTAR calculated in TAUIJ from local array in TAUIJ
C back into global array for the rest of the program.
C
C IDIM,JDIM,KDIM : Dimensions of block
C I              : Current I = const. plane
C PROPSC         : Local array for Properties (stored on an I = const. plane)
C PROPS          : Global array for Properties
C
      include '../header/common.h'
C
      DIMENSION PROPSC(2:JDIM,2:KDIM,NP),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP)
C
C Copy data from temporary storage into permanent PROPS array
C
      DO 120 K = 2, KDIM
         DO 110 J = 2, JDIM
            PROPS(J,K,I,ICMUST) = PROPSC(J,K,ICMUST)
  110    CONTINUE
  120 CONTINUE
C
C Do boundary conditions so that Cmu* is consistently defined on boundaries.
C
      DO 210 K = 2, KDIM
         PROPS(1     ,K,I,ICMUST) = PROPS(2   ,K,I,ICMUST)
         PROPS(JDIM+1,K,I,ICMUST) = PROPS(JDIM,K,I,ICMUST)
  210 CONTINUE
C
      DO 310 J = 1, JDIM+1
         PROPS(J,1     ,I,ICMUST) = PROPS(J,2   ,I,ICMUST)
         PROPS(J,KDIM+1,I,ICMUST) = PROPS(J,KDIM,I,ICMUST)
  310 CONTINUE
C
C If I = 2   , then do I = 1      boundary
C If I = IDIM, then do I = IDIM+1 boundary
C
      IF (I .EQ. 2) THEN
         DO 420 K = 1, KDIM+1
            DO 410 J = 1, JDIM+1
               PROPS(J,K,1     ,ICMUST) = PROPS(J,K,2   ,ICMUST)
  410       CONTINUE
  420    CONTINUE
C
      ELSE IF (I .EQ. IDIM) THEN
         DO 440 K = 1, KDIM+1
            DO 430 J = 1, JDIM+1
               PROPS(J,K,IDIM+1,ICMUST) = PROPS(J,K,IDIM,ICMUST)
  430       CONTINUE
  440    CONTINUE
      ENDIF
C
      RETURN
      END





      SUBROUTINE VSCSTY (NPTS, T, RMU)
C
C     Routine to calculate the molecular viscosity.
C
C
C     NPTS           : Number of points at which to calculate the molecular
C                      viscosity
C     T              : Temperature (T/Tinf)
C     RMU            : Molecular Viscosity
C
      include '../header/common.h'
C
      DIMENSION T     (NPTS),
     1          RMU   (NPTS)
C
C     Sutherland's Law
C
C     Constants for Sutherland's Law
C     Sref = 198.6 degrees Rankine.  If using SI units, correct to Kelvins.
C
      SREF  = 198.6E0
      IF ( SIUNIT ) SREF = SREF / 1.8E0
      STINF = SREF / TINF
      SNUM  = 1.E0 + STINF
C
C     Calculate molecular viscosity
C
      DO 100 I = 1, NPTS
         RMU(I) = T(I) * SQRT( T(I) ) * SNUM / ( T(I) + STINF )
  100 CONTINUE
C
C     Finished calculating molecular viscosity
C
      RETURN
      END





      SUBROUTINE YPLUSS (IDIM, JDIM, KDIM, IPSTRT, IPEND, ITURB,
     1                   Q, PROPS, S, VOL, IBCDAT, NBCS)
C
C Routine to calculate Y+ in the flowfield.
C This routine is used for either the time dependent calculation 
C or the marching calculation.  
C
C This routine is organized as follows:
C     1. Initialize Y+ = infinity 
C     2. Loop over all boundary conditions - for wall/wall functions do
C        a. Calculate Y+ = min (Y+, Y+ calculated from current wall)
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT,IPEND   : I plane to begin,end calculations at
C ITURB          : Turbulence model
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
C S              : Metrics
C VOL            : Cell volumes
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION IBCDAT(10,NBCS)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),
     2          IEND  (3)
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
      IJKDIM(1) = IDIM
      IJKDIM(2) = JDIM
      IJKDIM(3) = KDIM
C
C 1. Initialize Y+ to <infinity> in case there are no walls
C
      DO 130 I = IPSTRT, IPEND
         DO 120 K = 1, KDIM + 1
            DO 110 J = 1, JDIM + 1
               PROPS(J,K,I,3) = RMAX
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C 2. Loop over WALL boundary conditions
C    Calculate Y+ = MIN (y+, Y+) to allow for multiple walls in
C    the same coordinate direction.
C
      DO 500 LBC = 1, NBCS
         IBCTYP = IBCDAT(1,LBC)
         IF (IBCTYP .EQ. IWALL .OR. IBCTYP .EQ. IWALFN) THEN
            IDIR      = IBCDAT(2,LBC)
            INROUT    = IBCDAT(3,LBC)
            ISTRT(1)  = MAX (IBCDAT(4,LBC), IPSTRT)
            IEND (1)  = MIN (IBCDAT(5,LBC), IPEND)
            ISTRT(2)  = IBCDAT(6,LBC)
            IEND (2)  = IBCDAT(7,LBC)
            ISTRT(3)  = IBCDAT(8,LBC)
            IEND (3)  = IBCDAT(9,LBC)
C
C Modify beginning, ending indices to extend boundary condition to edge/corner
C
            DO 210 L = 1, 3
               IF (L .NE. IDIR) THEN
                  IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
                  IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
               END IF
  210       CONTINUE
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCPT2 is the boundary condition location outside the block from IBCPT
C      IBCN   is the point to the inside of the block from IBCPT
C      IBCM   is the location of the metrics
C
C
            IF (INROUT .EQ. 1) THEN
               IBCPT  = 1
               IBCPT2 = 0
               IBCN   = 2
               IBCM   = 1
            ELSE
               IBCPT  = IJKDIM(IDIR) + 1
               IBCPT2 = IJKDIM(IDIR) + 2
               IBCN   = IJKDIM(IDIR)
               IBCM   = IJKDIM(IDIR)
            ENDIF
C
            IC1 = ICY (IDIR, 1)
            IC2 = ICY (IDIR, 2)
            IC3 = ICY (IDIR, 3)
C
            DO 320 I3 = ISTRT(IC3), IEND(IC3)
               DO 310 I2 = ISTRT(IC2), IEND(IC2)
C
                  IBC  = IBCPT *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JBC  = IBCPT *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KBC  = IBCPT *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  IBC2 = IBCPT2*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JBC2 = IBCPT2*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KBC2 = IBCPT2*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  IN   = IBCN  *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JN   = IBCN  *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KN   = IBCN  *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  IM   = IBCM  *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JM   = IBCM  *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KM   = IBCM  *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  RHOW   = Q    (JBC,KBC,IBC,1)
                  RMUW   = PROPS(JBC,KBC,IBC,1)
                  YNORMW = ABS (PROPS(JBC,KBC,IBC,4))
C
                  YNORMN = ABS (PROPS(JN ,KN ,IN ,4))
C
C Calculate wall shear or use value from wall functions
C
                  IF ( IFWALF ) THEN
                     TAUW   = ABS (PROPS(JBC2,KBC2,IBC2,ITQTAU))
                  ELSE
C
C Calculate the wall value of du/dn
C
                     UBARN = S(JM,KM,IM,1,IDIR) * Q(JN,KN,IN,2)
     1                     + S(JM,KM,IM,2,IDIR) * Q(JN,KN,IN,3)
     2                     + S(JM,KM,IM,3,IDIR) * Q(JN,KN,IN,4)
                     VELN = SQRT ( 
     1                     (Q(JN,KN,IN,2)-UBARN*S(JM,KM,IM,1,IDIR))**2
     2                   + (Q(JN,KN,IN,3)-UBARN*S(JM,KM,IM,2,IDIR))**2
     3                   + (Q(JN,KN,IN,4)-UBARN*S(JM,KM,IM,3,IDIR))**2 )
C
C Calculate first derivative term as du/dn = delta(u_parallel) / delta(normal)
C
                     DVELDN = VELN / (ABS (YNORMN - YNORMW) + RSMALL)
C
C TAUW = mu * du/dy |(wall)
C
                     TAUW   = RMUW * DVELDN
                  ENDIF
C
C CTAUW = constants for non-dimensionalization * ustar
C         ustar = sqrt (Tau_wall / rho_wall)
C
                  CTAUW  = SQRT ( RE / FSMACH * TAUW / RHOW )
C
C Loop normal to the wall setting Y+
C
CDIR$ IVDEP
                  DO 300 I1 = 1, IJKDIM(IDIR) + 1
C
                     II = I1*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
                     JJ = I1*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
                     KK = I1*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
C Yplus = y ustar / nu 
C
                     RHO               = Q    (JJ,KK,II,1)
                     RMU               = PROPS(JJ,KK,II,1)
                     YNORML            = ABS (PROPS(JJ,KK,II,4))
                     YPL               = CTAUW * YNORML * RHO / RMU
                     PROPS(JJ,KK,II,3) = MIN (YPL, PROPS(JJ,KK,II,3))
  300             CONTINUE
  310          CONTINUE
  320       CONTINUE
C
         ENDIF
C
C Next boundary condition
C
  500 CONTINUE
C
C Finished calculating Y+
C
      RETURN
      END
