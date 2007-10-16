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
      SUBROUTINE SRCRES (IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, NPPTS, S, VOL, RC, RES,
     2                   QC, PROPSC, SN, DQDX, DQDY, DQDZ, DEL2K, 
     3                   TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ,
     4                   TAU, SRC, IBCDAT, NBCS, 
     5                   IFLUX, IMODEL, IFDELQ, NPRLIM, NPRNEG)
C
C     Routine to calculate the contriubtion to the residual from the
C     source terms.  This routine is used for either the 
C     time dependent calculation or the marching calculation.  
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     IPSTRT         : I Plane (cell center) to start calculations at 
C                      (Minimum Value: 2)
C     IPEND          : I Plane (cell center) to end caclulations at
C                      (Maximum Value: IDIM)
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
C     RC             : Position vector at cell centers
C     RES            : Residual
C     QC             : Storage for Q in one-dimensional array
C     PROPSC         : Storage for PROPS array in one-dimensional array
C     SN             : Storage for wall normals in one-dimensional array
C     DQDX,DQDY,DQDZ : Storage for derivatives of Q at cell centers
C     DEL2K          : Laplacian of turbulent kinetic energy
C     TKE            :
C     ETILDE         :
C     XSIW           :
C     FDAMP          :
C     EPSC           :
C     PDIL           :
C     SIJ            :
C     WIJ            :
C     TAU            : Storage for Reynolds stress terms
C     SRC            : Source terms for cell (J,K)
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IFLUX          : Flux evaluation scheme
C     IMODEL         : Turbulence model array
C     IFDELQ         : Array controlling calculation of delq
C     NPRLIM         : Number of times reset production
C NPRNEG         : Number of times invoke positivity preservation for production
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          QC    (2:JDIM,2:KDIM,NQ), 
     2          PROPSC(2:JDIM,2:KDIM,NP),
     3          SN    (2:JDIM,2:KDIM,3),
     4          DQDX  (2:JDIM,2:KDIM,NQ),
     5          DQDY  (2:JDIM,2:KDIM,NQ), 
     6          DQDZ  (2:JDIM,2:KDIM,NQ),
     7          DEL2K (2:JDIM,2:KDIM),
     8          TKE   (2:JDIM,2:KDIM),
     9          ETILDE(2:JDIM,2:KDIM),
     A          XSIW  (2:JDIM,2:KDIM),
     B          FDAMP (2:JDIM,2:KDIM,2),
     C          EPSC  (2:JDIM,2:KDIM),
     D          PDIL  (2:JDIM,2:KDIM),
     E          SIJ   (2:JDIM,2:KDIM,6),
     F          WIJ   (2:JDIM,2:KDIM,6),
     G          TAU   (2:JDIM,2:KDIM,6),
     H          SRC   (2:JDIM,2:KDIM,NF)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IMODEL(NMDL),
     1          IFDELQ(NQ)
C
C     Local variable:
C
      LOGICAL   ISAVPR
C
      ITURB  = IMODEL( 1)
      IFDL2K = IMODEL( 5)
      IASM   = IMODEL(10)
C
C     Set up control to store CMUSTAR back into global array if 
C     Algebraic Reynolds stress model
C
c-    Current EASM model for GS does not include variable C_mu for diffusion
c-    IF (IASM .EQ. IASMGS) THEN
c-       ISAVPR = .TRUE.
c-    ELSE
         ISAVPR = .FALSE.
c-    ENDIF
C
C     Number of points in a cross-section to calculate source terms at
C
      NPTS   = (JDIM - 1) * (KDIM - 1)
C
C     Loop over I calculating the source terms
C
      DO 300 I = IPSTRT, IPEND
C
         CALL SOURCE (IDIM, JDIM, KDIM, I, NPRLIM, NPRNEG,
     1                Q, PROPS, NPPTS, S, VOL, RC,
     2                QC, PROPSC, SN, DQDX, DQDY, DQDZ, DEL2K, 
     3                TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ, 
     4                TAU, SRC, IBCDAT, NBCS, 
     5                IFLUX, IMODEL, IFDELQ, ISAVPR)
C
C     Add the source term to the residual - include volume
C
         DO 250 L = 1, NF
            DO 240 K = 2, KDIM
               DO 230 J = 2, JDIM
                  RES(J,K,I,L) = RES(J,K,I,L) - SRC(J,K,L)*VOL(J,K,I)
  230          CONTINUE
  240       CONTINUE
  250    CONTINUE
  300 CONTINUE
C
C     Finished adding source terms to residual
C
      RETURN
      END
