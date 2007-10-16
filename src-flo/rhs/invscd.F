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
      SUBROUTINE FLUX (ITURB, NPTS, Q, S, TKE, H, F)
C
C     Routine to calculate interface fluxes as the full advection flux.
C     This routine calculates the fluxes for the Favre averaged 
C     Navier-Stokes equations for inviscid, laminar, algebraic turbulence 
C     model, two-equation turbulence model or Reynolds stress turbulence 
C     model.  It computes all of the inviscid fluxes for the above models.
C
C     The two-equation and Reynolds stress models need the turbulent kinetic 
C     energy (K) in the mean flow equations.  Therefore, this routine is
C     organized into the following sections:
C          1. Calculate Turbulent Kinetic Energy (TKE)
C             a. Inviscid, laminar, 
C                algebraic model       => K = 0
C             b. Two-Equation model    => K = K (Q(6))
C             c. Reynolds stress model => K = (Tau_xx+Tau_yy+Tau_zz)/2
C          2. Calculate Total Enthalpy (H)
C             a. Perfect gas    => H = gamma P / (RHO (GAMMA-1)) + (q^2)/2 + K
C          3. Calculate Advection Fluxes for Mean Flow Equations
C          4. Calculate additional Advection Fluxes
C             a. Inviscid, laminar, 
C                algebraic model       => No additional fluxes
C             b. Two-Equation model    => Fluxes for K and T2 equations
C             c. Reynolds Stress model => Fluxes for Tau_ij and T2 equations
C
C
C     ITURB          : Turbulence model
C     NPTS           : Size of array to calculate fluxes on.
C     Q              : Q at a cell interface to calculate fluxes at
C     S              : Metrics
C     TKE            : Turbulent kinetic energy 
C     H              : Total enthalpy 
C     F              : Flux at the interface
C
      include '../header/common.h'
C
      DIMENSION Q   (NPTS,NQ),
     1          S   (NPTS,3),
     2          TKE (NPTS),
     3          H   (NPTS),
     4          F   (NPTS,NF)
C
C     1. Calculate the Turbulence Kinetic Energy (TKE)
C
C     1-b. Turbulence Kinetic Energy for Two Equation Turbulence models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 100 I = 1, NPTS
            TKE(I) = Q(I,6)
  100    CONTINUE
C
C     1-c. Turbulence Kinetic Energy for Reynolds Stress Turbulence model
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 110 I = 1, NPTS
C
            TAUXX  = Q(I,6)
            TAUYY  = Q(I,7)
            TAUZZ  = Q(I,8)
            TKE(I) = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
  110    CONTINUE
C
C     1-a. Turbulence Kinetic Energy = 0 for Inviscid, Laminar and
C          Algebraic Turbulence models
C
      ELSE
         DO 120 I = 1, NPTS
            TKE(I) = 0.E0
  120    CONTINUE
      ENDIF
C
C     2. Calculate Total Enthalpy
C
C     2-a. Perfect Gas
C
      DO 200 I = 1, NPTS
C
         RHO   = Q(I,1)
         U     = Q(I,2) 
         V     = Q(I,3)
         W     = Q(I,4)
         P     = Q(I,5)
         H(I)  = GAMMA * P / (RHO * GAMM1) 
     1            + 0.5E0 * (U*U + V*V + W*W) + TKE(I)
  200 CONTINUE
C
C     3. Calculate F for Mean Flow Equations at the cell interfaces 
C
      DO 300 I = 1, NPTS
C
         RHO   = Q(I,1)
         U     = Q(I,2) 
         V     = Q(I,3)
         W     = Q(I,4)
         P     = Q(I,5)
C
         UBAR  = S(I,1)*U + S(I,2)*V + S(I,3)*W
C
C     Calculate F 
C     Area term is included in the RESID routine.
C
         F(I,1)  = RHO * UBAR 
         F(I,2)  = RHO * UBAR * U + S(I,1) * P
         F(I,3)  = RHO * UBAR * V + S(I,2) * P
         F(I,4)  = RHO * UBAR * W + S(I,3) * P
         F(I,5)  = RHO * UBAR * H(I)
 300  CONTINUE
C
C     4. Calculate Additional Advection Fluxes
C
C     4-a. Inviscid, Laminar, Algebraic Turbulence models => no additional flux
C
C     4-b. Two-Equation Turbulence models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 400 I = 1, NPTS
C
            RHO    = Q(I,1)
            U      = Q(I,2) 
            V      = Q(I,3)
            W      = Q(I,4)
C
            T2     = Q(I,7)
C
            UBAR   = S(I,1)*U + S(I,2)*V + S(I,3)*W
C
            F(I,6)  = RHO * UBAR * TKE(I)
            F(I,7)  = RHO * UBAR * T2 
  400    CONTINUE
C
C     4-c. Reynolds Stress Turbulence models
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 410 I = 1, NPTS
C
            RHO    = Q(I,1)
            U      = Q(I,2) 
            V      = Q(I,3)
            W      = Q(I,4)
C
            TAUXX  = Q(I,6)
            TAUYY  = Q(I,7)
            TAUZZ  = Q(I,8)
            TAUXY  = Q(I,9)
            TAUXZ  = Q(I,10)
            TAUYZ  = Q(I,11)
            T2     = Q(I,12)
C
            UBAR   = S(I,1)*U + S(I,2)*V + S(I,3)*W
C
            F(I,6)  = RHO * UBAR * TAUXX
            F(I,7)  = RHO * UBAR * TAUYY
            F(I,8)  = RHO * UBAR * TAUZZ
            F(I,9)  = RHO * UBAR * TAUXY
            F(I,10) = RHO * UBAR * TAUXZ
            F(I,11) = RHO * UBAR * TAUYZ
            F(I,12) = RHO * UBAR * T2
  410    CONTINUE
      ENDIF
C
C     Finished calculating full flux at the interface 
C
      RETURN
      END




      SUBROUTINE FDISS (ITURB, NPTS, Q, DQC, S, ALAM, F)
C
C Routine to calculate interface fluxes as the full advection flux.
C This routine calculates the fluxes for the Favre averaged 
C Navier-Stokes equations for inviscid, laminar, algebraic turbulence 
C model, two-equation turbulence model or Reynolds stress turbulence 
C model.  It computes all of the inviscid fluxes for the above models.
C
C The two-equation and Reynolds stress models need the turbulent kinetic 
C energy (K) in the mean flow equations.  Therefore, this routine is
C organized into the following sections:
C      1. Calculate artificial dissipation for 2nd or 4th order
C         a. Scalar dissipation
C            i.   Mean flow equations: Inviscid, laminar, algebraic model
C            ii.  Two-equation model
C            iii. Reynolds stress model
C
C ITURB          : Turbulence model
C NPTS           : Size of array to calculate fluxes on.
C Q              : Q at a cell interface to calculate fluxes at
C DQC            : Delta Q conserved at interfaces for dissipation
C                  = 2nd order -> third undivided difference
C                    4th order -> fifth undivided difference
C S              : Metrics
C ALAM           : Eigenvalues (Lambda)
C F              : Flux at the interface
C
      include '../header/common.h'
C
      DIMENSION Q   (NPTS,NQ),
     1          DQC (NPTS,NQ),
     2          S   (NPTS,3),
     3          ALAM(NPTS),
     4          F   (NPTS,NF)
C
C 1.   Calculate artificial dissipation for 2nd or 4th order
C 1-a. Scalar dissipation
C 1-a-i.   Mean flow equations: Inviscid, laminar, algebraic model
C
c     IF () THEN
         DO 100 I = 1, NPTS
C
         RHO     = Q(I,1)
         U       = Q(I,2) 
         V       = Q(I,3)
         W       = Q(I,4)
         P       = Q(I,5)
C
         UBAR    = S(I,1)*U + S(I,2)*V + S(I,3)*W
         A       = SQRT (GAMMA * P / RHO)
C
         ALAM(I) = ABS(UBAR) + A
C
C      Calculate F = F - F_diss
C      Area term is included in the RESID routine.
C
            F(I,1)  = F(I,1) - ADKP46 * ALAM(I) * DQC(I,1)
            F(I,2)  = F(I,2) - ADKP46 * ALAM(I) * DQC(I,2)
            F(I,3)  = F(I,3) - ADKP46 * ALAM(I) * DQC(I,3)
            F(I,4)  = F(I,4) - ADKP46 * ALAM(I) * DQC(I,4)
            F(I,5)  = F(I,5) - ADKP46 * ALAM(I) * DQC(I,5)
 100     CONTINUE
C
C 1-a-ii.  Two-equation model
C
         IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
            DO 110 I = 1, NPTS
               F(I,6)  = F(I,6) - ADKP46 * ALAM(I) * DQC(I,6)
               F(I,7)  = F(I,7) - ADKP46 * ALAM(I) * DQC(I,7)
 110        CONTINUE
C
C 1-a-iii. Reynolds stress model
C
         ELSE IF (ITURB .EQ. ITRS) THEN
            DO 120 I = 1, NPTS
               F(I, 6)  = F(I, 6) - ADKP46 * ALAM(I) * DQC(I, 6)
               F(I, 7)  = F(I, 7) - ADKP46 * ALAM(I) * DQC(I, 7)
               F(I, 8)  = F(I, 8) - ADKP46 * ALAM(I) * DQC(I, 8)
               F(I, 9)  = F(I, 9) - ADKP46 * ALAM(I) * DQC(I, 9)
               F(I,10)  = F(I,10) - ADKP46 * ALAM(I) * DQC(I,10)
               F(I,11)  = F(I,11) - ADKP46 * ALAM(I) * DQC(I,11)
               F(I,12)  = F(I,12) - ADKP46 * ALAM(I) * DQC(I,12)
 120        CONTINUE
         ENDIF
c     ELSE
c     ENDIF
C
C     Finished calculating artificial dissipation flux at interface
C
      RETURN
      END
