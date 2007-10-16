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
      SUBROUTINE PRPRTY (IDIM, JDIM, KDIM, IPSTRT, IPEND, NXSECT,
     1                   Q, PROPS, S, VOL, IBCDAT, NBCS, IMODEL, IFDELQ,
     2                   T, QC, PROPSC, VOLC, DQDX, DQDY, DQDZ, WORKX)
C
C Routine to calculate the PROPERTIES of the flowfield.
C This routine is used for either the time dependent calculation 
C or the marching calculation.  
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM+1)
C NXSECT         : Maximum size of any (ixj, jxk or ixk)
C                  (used in dimensioning temporary space used for 
C                   calculating q(+,-) and f.
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity 
C                  PROPS(3) = Y+
C S              : Metrics
C VOL            : Cell volumes
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C IMODEL         : Turbulence model array
C IFDELQ         : Array to control calculation of delq
C T              : Temporary storage for the Temperature
C QC,PROPSC,VOLC : Q, PROPS and VOL in one-dimensional array
C DQDX,DQDY,DQDZ : Derivatives of Q at cell centers
C WORKX          : Work array
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION T     (0:JDIM+2,0:KDIM+2),
     1          QC    (1:JDIM+1,1:KDIM+1,NQ),
     2          PROPSC(1:JDIM+1,1:KDIM+1,NP),
     3          VOLC  (1:JDIM+1,1:KDIM+1),
     4          DQDX  (1:JDIM+1,1:KDIM+1,NQ),
     5          DQDY  (1:JDIM+1,1:KDIM+1,NQ),
     6          DQDZ  (1:JDIM+1,1:KDIM+1,NQ)
C
      DIMENSION WORKX (NXSECT,NW)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IMODEL(NMDL)
C
      ITURB = IMODEL(1)
      ISGS  = IMODEL(6)
C
C     Loop over I 
C
      DO 300 I = IPSTRT, IPEND
C
C     Calculating the molecular viscosity
C
         DO 220 K = 0, KDIM + 2
            DO 210 J = 0, JDIM + 2
               T(J,K) = GAMMA * Q(J,K,I,5) / Q(J,K,I,1)
  210       CONTINUE
  220    CONTINUE
C
C     Calculate the molecular viscosity and store it into the PROPS array
C
         NPTS = (JDIM + 3) * (KDIM + 3)
         CALL VSCSTY (NPTS, T, PROPS(0,0,I,1))
  300 CONTINUE
C
C     Calculate the values of Y+
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW .OR. 
     1    ITURB .EQ. ITRS .OR. ITURB .EQ. ITLES) THEN
         CALL YPLUSS (IDIM, JDIM, KDIM, IPSTRT, IPEND, ITURB,
     1                Q, PROPS, S, VOL, IBCDAT, NBCS)
      ENDIF
C
C     Calculate the Turbulent Eddy Viscosity
C
C     Loop over I
C
      DO 600 I = IPSTRT, IPEND
C
C     Baldwin-Lomax Turbulence Model
C     Update turbulent eddy viscosity with algebraic model if algebraic model 
C     is specified or if NITALG < ITER such that the algebraic model is used
C     in the initial transient to provide a better convergence for the
C     higher order models.
C
C     Only call Baldwin-Lomax for 2 < i < idim-1
C     PROPS(2) is RMUT
C     Temporary storage: Use WORKX(1)     for VORT
C                            WORKX(2-5)   for TVIS
C                            WORKX(6-9)   for DIST
C                            WORKX(10)    for RHO
C                            WORKX(11)    for VORTJ
C                            WORKX(12)    for RMUTJ
C                            WORKX(13)    for Y
C                            WORKX(14)    for FUN
C
         IF (MUTALG .OR. ITURB .EQ. ITBLMX) THEN
            IF (I .GE. 2 .AND. I .LE. IDIM) THEN
               MAXDIM = MAX (JDIM, KDIM)
               CALL TURB (IDIM, JDIM, KDIM, MAXDIM, I, ITURB,
     1                    Q, PROPS(0,0,0,2), S, VOL, IBCDAT, NBCS,
     2                    WORKX(1, 1), WORKX(1, 2), WORKX(1, 6),
     3                    WORKX(1,10), WORKX(1,11), WORKX(1,12), 
     4                    WORKX(1,13), WORKX(1,14))
            ENDIF
C
C     k-epsilon, Reynolds stress turbulence model
C     Temporary storage: Use WORKX(1)  for TKE
C                            WORKX(2)  for FDAMP
C
         ELSE IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITRS) THEN
C
            DO 330 L = 1, NQ
               DO 320 K = 1, KDIM+1
                  DO 310 J = 1, JDIM+1
                     QC(J,K,L) = Q(J,K,I,L)
  310             CONTINUE
  320          CONTINUE
  330       CONTINUE
C
            DO 360 L = 1, NP
               DO 350 K = 1, KDIM+1
                  DO 340 J = 1, JDIM+1
                     PROPSC(J,K,L) = PROPS(J,K,I,L)
  340             CONTINUE
  350          CONTINUE
  360       CONTINUE
C
            NPTS  = (JDIM + 1) * (KDIM + 1)
            CALL GETTKE (ITURB,  NPTS, QC, WORKX(1,1))
            CALL RMUTKE (IMODEL, NPTS, QC, PROPSC,WORKX(1,1),WORKX(1,2))
C
            DO 380 K = 1, KDIM + 1
               DO 370 J = 1, JDIM + 1
                  PROPS(J,K,I,2) = PROPSC(J,K,2)
  370          CONTINUE
  380       CONTINUE
C
C     k-omega turbulence model
C     PROPS(1) is RMU and PROPS(2) is RMUT
C
         ELSE IF (ITURB .EQ. ITKW) THEN
            CALL RMUTKW (IDIM, JDIM, KDIM, I, Q, PROPS(0,0,0,2)) 
C
C     Large-Eddy Simulation
C     Temporary storage: WORKX    for FDAMP
C
         ELSE IF (ITURB .EQ. ITLES) THEN
C
            DO 470 K = 1, KDIM+1
               DO 465 J = 1, JDIM+1
                  VOLC(J,K) = VOL(J,K,I)
  465          CONTINUE
  470       CONTINUE
C
            NPTS  = (JDIM + 1) * (KDIM + 1)
C
C     Temporary variables used internally in DELQC:
C         WORKX(1-2) for QJ
C         WORKX(3-4) for QK
C         WORKX(5-6) for QPI
C         WORKX(7-8) for QMI
C     QJ is dimensioned one larger than normal routines therefore assign it
C     2 planes of WORKX. Similarly for QK. QPI and QMI require only on plane
C     as they are dimensioned normally - but they are allocated one additional
C     plane as a precaution.
C
C     Calculate derivatives at interior cells and boundaries.
C
            NBOUND = 1
            CALL DELQC (IDIM, JDIM, KDIM, I, IFDELQ, NBOUND,
     1                  Q, PROPS, S, VOL, IBCDAT, NBCS,
     2                  WORKX(1,1), WORKX(1,3), WORKX(1,5), WORKX(1,7),
     3                  QC, PROPSC, DQDX, DQDY, DQDZ)
C
C     Calculate Structure function for SF model
C     Use WORKX(1) to store SF
C
            IF (ISGS .EQ. ISGSFM) THEN
               CALL GETSF2 (IMODEL, IDIM, JDIM, KDIM, I, Q, PROPS, S, 
     1                      VOL, IBCDAT, NBCS, WORKX(1,1), IERROR)
            ENDIF
C
C     Calculate SGS coefficient
C
            CALL RMUSGS (IMODEL, NPTS, QC, PROPSC, VOLC, DQDX, DQDY, 
     1                   DQDZ, WORKX(1,1), WORKX(1,2), IERROR)
            IF (IERROR .NE. INO) THEN
               STOP
            ENDIF
C
            DO 490 K = 1, KDIM + 1
               DO 480 J = 1, JDIM + 1
                  PROPS(J,K,I,2) = PROPSC(J,K,2)
  480          CONTINUE
  490       CONTINUE
C
C     No turbulence - set RMUT to zero
C
         ELSE
            DO 520 K = 0, KDIM + 2
               DO 510 J = 0, JDIM + 2
                  PROPS(J,K,I,2) = 0.E0
  510          CONTINUE
  520       CONTINUE
         ENDIF
C
  600 CONTINUE
C
C     Finished with the PROPERTIES
C
      RETURN
      END
