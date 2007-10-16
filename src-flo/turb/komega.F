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






      SUBROUTINE RMUTKW (IDIM, JDIM, KDIM, I, Q, RMUT)
C
C     Routine to calculate the turbulent eddy viscosity  for the k-omega
C     turbulence model.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : Plane to do calculation at
C     Q              : Primitive variables at cell centers
C     RMUT           : Turbulent eddy viscosity at cell centers
C                      of the I-plane
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          RMUT  (0:JDIM+2,0:KDIM+2,0:IDIM+2) 
C
C     Calculate the turbulent eddy viscosity at cell centers from values
C     of turbulent kinetic energy (k) and specific dissipation rate (omega)
C
      DO 120 K = 0, KDIM + 2
         DO 110 J = 0, JDIM + 2
            RHO         = Q(J,K,I,1)
            TKE         = Q(J,K,I,6)
            OMEGA       = Q(J,K,I,7)
            RMUT(J,K,I) = CMU * RHO * TKE / (OMEGA + RSMALL) * 
     1                    RE / FSMACH
  110    CONTINUE
  120 CONTINUE
C
C
C     Finished calculating RMUT from k and omega
C
      RETURN
      END






      SUBROUTINE SRCKW (NPTS, Q, PROPS, DQDX, DQDY, DQDZ, TAU, SRC, 
     1                  NPRLIM, NPRNEG)
C
C Routine to calculate the source terms for the k-omega turbulence model.
C
C This routine is organized as follows:
C      1. Form the source terms for the k-epsilon model
C
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C PROPS          : Properties stored at cell centers
C                  PROPS(1) = RMU   molecular viscosity
C                  PROPS(2) = RMUT  turbulent eddy viscosity
C                  PROPS(3) = YPLUS Y+
C DQDX,DQDY,DQDZ : Derivatives of Q at cell centers
C TAU            : Reynolds stresses
C SRC            : Source terms for the k-omega model at the cell centers
C                  of the I-plane
C NPRLIM         : Counter for number of times production limited
C NPRNEG         : Number of times invoke positivity preservation for production
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          PROPS (NPTS,NP)
C
      DIMENSION DQDX  (NPTS,NQ),
     1          DQDY  (NPTS,NQ),
     2          DQDZ  (NPTS,NQ),
     3          TAU   (NPTS,6),
     4          SRC   (NPTS,NF)
C
C 1.  Form the source terms for the k-omega model
C
      DO 100 I = 1, NPTS
C
         RHO    = Q(I,1)
         TKE    = Q(I,6)
         OMEGA  = Q(I,7)
C
C Production term:
C    PROD = TAUXX*DUDX + TAUYY*DVDY + TAUZZ*DWDZ +
C           TAUXY*(DUDY+DVDX) + TAUXZ*(DUDZ+DWDX) + TAUYZ*(DVDZ+DWDY)
C
         PROD   =   TAU(I,1) * DQDX(I,2)
     1            + TAU(I,2) * DQDY(I,3)
     2            + TAU(I,3) * DQDZ(I,4)
     3            + TAU(I,4) * (DQDY(I,2) + DQDX(I,3))
     4            + TAU(I,5) * (DQDZ(I,2) + DQDX(I,4))
     5            + TAU(I,6) * (DQDZ(I,3) + DQDY(I,4))
C
C Alternate Form of Production Term: PROD = mu_t * Vorticity^2
C
C        WXY     = 0.5E0 * (DQDY(I,2) - DQDX(I,3))
C        WXZ     = 0.5E0 * (DQDZ(I,2) - DQDX(I,4))
C        WYZ     = 0.5E0 * (DQDZ(I,3) - DQDY(I,4))
C        WVORT   = 4.E0 * (WXY*WXY + WXZ*WXZ + WYZ*WYZ)
C        RMUT    = PROPS(I,2)
C        PROD    = RMUT * WVORT * FSMACH / RE
C
C Transition specification
C
         IF (PROPS(I,4) .LT. 0.E0) PROD = 0.E0
C
C Positivity preservation for production
C
         IF (POSPRD) THEN
            IF (PROD .LT. 0.E0) THEN
               NPRNEG = NPRNEG + 1
               PROD   = ABS (PROD)
            ENDIF
         ENDIF
C
         PRODUL = PROD
C
C Limit production term for robustness:
C
         TSTDIS = PRDLIM * BSTRKW * RHO * TKE * OMEGA 
         IF (PROD .GT. TSTDIS) THEN
            NPRLIM = NPRLIM + 1
            PROD   = TSTDIS
         ENDIF
C
C Production for omega equation: use either limited or unlimited form
C
         PRODW  = PRDE * PROD + PRDEM1 * PRODUL
C
C Calculate Source Term
C
         SRC(I,1) = 0.E0
         SRC(I,2) = 0.E0
         SRC(I,3) = 0.E0
         SRC(I,4) = 0.E0
         SRC(I,5) = 0.E0
         SRC(I,6) = PROD - BSTRKW * RHO * TKE * OMEGA 
         SRC(I,7) = GKW * OMEGA / TKE * PRODW - BKW * RHO * OMEGA*OMEGA
  100 CONTINUE
C
C Finished with k-omega source terms
C
      RETURN
      END
