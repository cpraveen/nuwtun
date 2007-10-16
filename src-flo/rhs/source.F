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
      SUBROUTINE SOURCE (IDIM, JDIM, KDIM, I, NPRLIM, NPRNEG,
     1                   Q, PROPS, NPPTS, S, VOL, RC,
     2                   QC, PROPSC, SN, DQDX, DQDY, DQDZ, DEL2K, 
     3                   TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ, 
     4                   TAU, SRC, IBCDAT, NBCS, 
     5                   IFLUX, IMODEL, IFDELQ, ISAVPR)
C
C Routine to calculate the contriubtion the source terms.  This routine is 
C used for either the time dependent calculation or the marching calculation.  
C This routine is used for explicit and implicit source term treatment.
C
C This routine is organized as follows:
C      1. Calculate the Laplacian of the turbulent kinetic energy
C      2. Calculate derivatives of Q at cell centers
C      3. Calculate turbulent kinetic energy in 1D array
C      4. Calculate the Reynold's stresses at the cell centers.  
C      5. Store CMU* from the local array back into the global storage
C      6. Calculate S_ij, W_ij
C      7. Call routines to calculate turbulence model source terms
C         a. K-epsilon
C         b. K-omega
C         c. Reynolds stress
C
C IDIM,JDIM,KDIM : Dimensions of current block
C I              : I Plane (cell center) to calculate source terms at
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
C RES            : Residual
C QC             : Storage for Q in one-dimensional array
C PROPSC         : Storage for PROPS array in one-dimensional array
C SN             : Storage for wall normals in one-dimensional array
C DQDX,DQDY,DQDZ : Storage for derivatives of Q at cell centers
C DEL2K          : Laplacian of turbulent kinetic energy
C TKE            :
C ETILDE         :
C XSIW           :
C FDAMP          :
C EPSC           :
C PDIL           :
C SIJ            :
C WIJ            :
C TAU            : Storage for Reynolds stress terms
C SRC            : Source terms for cell (J,K)
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C IFLUX          : Flux evaluation scheme
C IMODEL         : Turbulence model array
C IFDELQ         : Array controlling calculation of delq
C ISAVPR         : Logical controlling storage of CMUSTAR from local
C                  array where it is calculated in TAUIJ back into
C                  global array.  This is done only when called from RESID.
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      DIMENSION QC    (2:JDIM,2:KDIM,NQ), 
     1          PROPSC(2:JDIM,2:KDIM,NP),
     2          SN    (2:JDIM,2:KDIM,3),
     3          DQDX  (2:JDIM,2:KDIM,NQ),
     4          DQDY  (2:JDIM,2:KDIM,NQ), 
     5          DQDZ  (2:JDIM,2:KDIM,NQ),
     6          DEL2K (2:JDIM,2:KDIM),
     7          TKE   (2:JDIM,2:KDIM),
     8          ETILDE(2:JDIM,2:KDIM),
     9          XSIW  (2:JDIM,2:KDIM),
     A          FDAMP (2:JDIM,2:KDIM,2),
     B          EPSC  (2:JDIM,2:KDIM),
     C          PDIL  (2:JDIM,2:KDIM),
     D          SIJ   (2:JDIM,2:KDIM,6),
     E          WIJ   (2:JDIM,2:KDIM,6),
     F          TAU   (2:JDIM,2:KDIM,6),
     G          SRC   (2:JDIM,2:KDIM,NF)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IMODEL(NMDL),
     1          IFDELQ(NQ)
C
      LOGICAL   ISAVPR
C
      ITURB  = IMODEL( 1)
      IFDL2K = IMODEL( 5)
      IASM   = IMODEL(10)
      IEPSC  = IMODEL(11)
      IPDIL  = IMODEL(12)
      ISRCW  = IMODEL(13)
C
C     Number of points in a cross-section to calculate source terms at
C
      NPTS   = (JDIM - 1) * (KDIM - 1)
C
C 1.  Calculate the Laplacian of the turbulent kinetic energy
C     Temporary variables used internally in DEL2KC:
C         TAU(1-3) for QJ  
C         SRC(1-3) for QK  
C         TAU(5->) for QPI 
C         SRC(4-6) for QMI 
C     QJ is dimensioned as JDIM,2:KDIM therefore need more room for it 
C     than 1 plane of TAU. Similarly for QK.  QPI and QMI require only 
C     one plane as they are dimensioned 2:JDIM,2:KDIM.
C     Use SRC for some of the temporary space to fix a memory allocation
C     problem.  Use DQDX,DQDY,DQDZ as temporary space here as it isn't
C     calculated until after this call.
C
      IF (IFDL2K .EQ. 1) THEN
         CALL DEL2KC (IDIM, JDIM, KDIM, I, ITURB,
     1                Q, S, VOL, IBCDAT, NBCS,
     2                TAU(2,2,1), SRC(2,2,1), TAU(2,2,5), SRC(2,2,4),
     3                DQDX, DQDY, DQDZ, DEL2K)
      ENDIF
C
C 2.  Calculate derivatives of Q at cell centers.
C     Temporary variables used internally in DELQC:
C         TAU(1-3) for QJ  
C         SRC(1-3) for QK  
C         TAU(5->) for QPI 
C         SRC(4-6) for QMI 
C     QJ is dimensioned as JDIM,2:KDIM therefore need more room for it 
C     than 1 plane of TAU. Similarly for QK.  QPI and QMI require only 
C     one plane as they are dimensioned 2:JDIM,2:KDIM.
C     Use SRC for some of the temporary space to fix a memory allocation
C     problem.
C
C     Calculate derivatives at interior cell centers only.
C
      NBOUND = 0
      CALL DELQC (IDIM, JDIM, KDIM, I, IFDELQ, NBOUND,
     1            Q, PROPS, S, VOL, IBCDAT, NBCS,
     2            TAU(2,2,1), SRC(2,2,1), TAU(2,2,5), SRC(2,2,4),
     3            QC, PROPSC, DQDX, DQDY, DQDZ)
C
C 3.  Calculate turbulent kinetic energy in 1D array
C
      CALL GETTKE (ITURB, NPTS, QC, TKE)
C
C 4.  Calculate the Reynold's stresses at the cell centers.  
C                                      ~
C     The array TAU contains {- rho (u" u")}.
C                                     i  j
C     For laminar flow:            TAU_ij => zero.
C         algebraic models:        TAU_ij => modeled with Boussinesq approx.
C         two-equation models:     TAU_ij => modeled with Boussinesq approx.
C         Reynold's stress models: TAU_ij => computed from PDE.
C     Pass in RMUT and DUDX,DUDY,DUDZ for Boussinesq approximation.
C     Pass in QC for Q data at same location as DQDX,... derivative data.
C     The Reynold's stresses are stored as:
C         TAU(1) = TAUXX, TAU(2) = TAUYY, TAU(3) = TAUZZ,
C         TAU(4) = TAUXY, TAU(5) = TAUXZ, TAU(6) = TAUYZ
C
      CALL TAUIJ (IMODEL, NPTS, QC, PROPSC, DQDX, DQDY, DQDZ, TAU)
C
C 5.  Store the value of CMU* from the local array back into the global
C     array if using Algebraic Reynolds stress model and ISAVPR is true
C
      IF (IASM .EQ. IASMGS  .AND.  ISAVPR) THEN
         CALL SAVEPR (IDIM, JDIM, KDIM, I, PROPSC, PROPS)
      ENDIF
C
C 6.  Calculate S_ij, W_ij
C
      CALL GETSW  (NPTS, QC, DQDX, DQDY, DQDZ, SIJ, WIJ)
C
C 7.  Turbulence model source terms
C
C 7-a. k-epsilon turbulence model
C
      IF (ITURB .EQ. ITKE) THEN
         CALL SRCKE  (IMODEL, NPTS, QC, PROPSC, DQDX, DQDY, DQDZ, TAU, 
     1                DEL2K, TKE, ETILDE, XSIW, FDAMP, 
     2                EPSC, PDIL, SIJ, WIJ, SRC, NPRLIM, NPRNEG)
C
C 7-b. k-omega turbulence model
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         CALL SRCKW  (NPTS, QC, PROPSC, DQDX, DQDY, DQDZ, TAU, SRC, 
     1                NPRLIM, NPRNEG)
C
C 7-c. Reynolds Stress turbulence model
C      SRCRS  calculates the source terms
C      SRCRSW calculates additional wall normal dependent source terms
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         CALL SRCRS  (IMODEL, NPTS, QC, PROPSC, DQDX, DQDY, DQDZ, 
     1                TAU, DEL2K, TKE, ETILDE, XSIW, FDAMP, 
     2                EPSC, PDIL, SIJ, WIJ, SRC, NPRLIM, NPRNEG)
C
C     Calculate wall normal dependent source terms for models that include them
C     Temporary storage: TAU    for SRCW
C                        ETILDE for YNI
C                        XSIW   for YNITOT
C
         IF (ISRCW .EQ. IYES) THEN
            CALL SRCRSW (IMODEL, QC, PROPSC, SN, DQDX, DQDY, DQDZ, 
     1                   TKE, EPSC, SRC, TAU, 
     2                   IDIM, JDIM, KDIM, I, NBCS, IBCDAT,
     3                   S, RC, ETILDE, XSIW)
         ENDIF
      ENDIF
C
C     Finished with the source terms at I
C
      RETURN
      END
