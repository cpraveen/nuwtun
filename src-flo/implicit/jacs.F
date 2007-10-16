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

      SUBROUTINE DFDQ (ITURB, NPTS, Q, S, TKE, H, A)
C
C     Routine to calculate Jacobians at the interface using Q
C     as DF/DQ for the mean equations and a either a two-equation 
C     turbulence model of the k-psi form or a Reynolds stress closure.
C     DF/DQ formulas match those given by Pulliam and Steger in 
C     AIAA 85-0360 in the Appendix for the special case of k = 0.
C
C
C     ITURB  : Turbulence model
C     NPTS   : Number of points at which to calculate DF/DQ
C     Q      : The primitive variables from which to calculate DF/DQ
C     S      : Metrics
C     TKE    : Turbulent kinetic energy
C     H      : Total enthalpy
C     A      : The Jacobian matrix DF/DQ
C
C     The two-equation and Reynolds stress models require the turbulent
C     kinetic energy in the mean flow equations.  Therefore, this routine
C     is organized into the following sections:
C
C          1. Calculate the Turbulent Kinetic Energy (TKE)
C             a. Inviscid, laminar, algebraic model => TKE = 0
C             b. Two-Equation model                 => TKE = Q(6)
C             c. Reynolds stress model              => TKE = (Q(7)+Q(8)+Q(9))/2
C          2. Calculate Total Enthalpy
C             a. Perfect Gas   => H = gamma P / (RHO(gamma-1)) + (q^2)/2 + TKE
C          3. Calculate Flux Jacobians for the Mean Flow Equations
C          4. Calculate Additional Flux Jacobians
C             a. Inviscid, laminar, algebraic model => No additional Jacobians
C             b. Two-Equation model                 => Jacobians for TKE, phi
C             c. Reynolds stress model              => Jacobians for Tau_ij,phi
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          S     (NPTS,4),
     2          TKE   (NPTS),
     3          H     (NPTS),
     4          A     (NPTS,NRANK,NRANK)
C
C     1. Calculate Turbulent Kinetic Energy (TKE)
C
C     1-b. Turbulent Kinetic Energy for Two-Equation models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 100 I = 1, NPTS
            TKE(I) = Q(I,6)
  100    CONTINUE
C
C     1-c. Turbulent Kinetic Energy for Reynolds Stress models
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 110 I = 1, NPTS
            TKE(I) = 0.5E0 * (Q(I,6) + Q(I,7) + Q(I,8))
  110    CONTINUE
C
C     1-a. Turbulent Kinetic Energy = 0 for inviscid, laminar, algebraic models
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
         RHO  = Q(I,1)
         U    = Q(I,2)
         V    = Q(I,3)
         W    = Q(I,4)
         P    = Q(I,5)
         H(I) = GAMMA * P / (RHO * GAMM1)
     1          + 0.5E0 * (U*U + V*V + W*W) + TKE(I)
  200 CONTINUE
C
C     3. Calculate DF/DQ at the cell interfaces for Mean Flow Equations
C
      DO 300 I = 1, NPTS
         RHO   = Q(I,1)
         U     = Q(I,2)
         V     = Q(I,3)
         W     = Q(I,4)
         P     = Q(I,5)
C
         UBAR  = S(I,1)*U + S(I,2)*V + S(I,3)*W
         PHISQ = 0.5E0 * GAMM1 * (U*U + V*V + W*W)
C
         A(I,1,1) =   0.0E0
         A(I,1,2) =   S(I,4) * S(I,1)
         A(I,1,3) =   S(I,4) * S(I,2)
         A(I,1,4) =   S(I,4) * S(I,3)
         A(I,1,5) =   0.0E0
C
         A(I,2,1) =   S(I,4)*( S(I,1)*PHISQ - U*UBAR )
         A(I,2,2) =   S(I,4)*( UBAR     - S(I,1)*(GAMMA-2.E0)*U )
         A(I,2,3) =   S(I,4)*( S(I,2)*U - S(I,1)*GAMM1*V )
         A(I,2,4) =   S(I,4)*( S(I,3)*U - S(I,1)*GAMM1*W )
         A(I,2,5) =   S(I,4)*S(I,1)*GAMM1
C
         A(I,3,1) =   S(I,4)*( S(I,2)*PHISQ - V*UBAR )
         A(I,3,2) =   S(I,4)*( S(I,1)*V - S(I,2)*GAMM1*U )
         A(I,3,3) =   S(I,4)*( UBAR     - S(I,2)*(GAMMA-2.E0)*V )
         A(I,3,4) =   S(I,4)*( S(I,3)*V - S(I,2)*GAMM1*W )
         A(I,3,5) =   S(I,4)*S(I,2)*GAMM1
C
         A(I,4,1) =   S(I,4)*( S(I,3)*PHISQ - W*UBAR )
         A(I,4,2) =   S(I,4)*( S(I,1)*W - S(I,3)*GAMM1*U )
         A(I,4,3) =   S(I,4)*( S(I,2)*W - S(I,3)*GAMM1*V )
         A(I,4,4) =   S(I,4)*( UBAR     - S(I,3)*(GAMMA-2.E0)*W )
         A(I,4,5) =   S(I,4)*S(I,3)*GAMM1
C
         A(I,5,1) =   S(I,4)*UBAR*( PHISQ - H(I) )
         A(I,5,2) =   S(I,4)*( S(I,1)*H(I) - GAMM1*UBAR*U )
         A(I,5,3) =   S(I,4)*( S(I,2)*H(I) - GAMM1*UBAR*V )
         A(I,5,4) =   S(I,4)*( S(I,3)*H(I) - GAMM1*UBAR*W )
         A(I,5,5) =   S(I,4)*GAMMA*UBAR 
C
  300 CONTINUE
C
C     4. Additional Jacobians
C
C     4-b. Two-Equation Models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 410 I = 1, NPTS
            RHO      = Q(I,1)
            U        = Q(I,2)
            V        = Q(I,3)
            W        = Q(I,4)
            P        = Q(I,5)
            T2       = Q(I,7)
            UBAR     = S(I,1)*U + S(I,2)*V + S(I,3)*W
C
            A(I,1,6) =   0.0E0
            A(I,1,7) =   0.0E0
C
            A(I,2,6) = - S(I,4)*S(I,1)*GAMM1
            A(I,2,7) =   0.0E0
C
            A(I,3,6) = - S(I,4)*S(I,2)*GAMM1
            A(I,3,7) =   0.0E0
C
            A(I,4,6) = - S(I,4)*S(I,3)*GAMM1
            A(I,4,7) =   0.0E0
C
            A(I,5,6) = - S(I,4)*GAMM1*UBAR 
            A(I,5,7) =   0.0E0
C
            A(I,6,1) = - S(I,4) * UBAR * TKE(I)
            A(I,6,2) =   S(I,4) * S(I,1) * TKE(I)
            A(I,6,3) =   S(I,4) * S(I,2) * TKE(I)
            A(I,6,4) =   S(I,4) * S(I,3) * TKE(I)
            A(I,6,5) =   0.0E0
            A(I,6,6) =   S(I,4) * UBAR
            A(I,6,7) =   0.0E0
C
            A(I,7,1) = - S(I,4) * UBAR * T2
            A(I,7,2) =   S(I,4) * S(I,1) * T2
            A(I,7,3) =   S(I,4) * S(I,2) * T2
            A(I,7,4) =   S(I,4) * S(I,3) * T2
            A(I,7,5) =   0.0E0
            A(I,7,6) =   0.0E0
            A(I,7,7) =   S(I,4) * UBAR
  410    CONTINUE
C
C     4-c. Additional Jacobians for Reynolds Stress models
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 420 I = 1, NPTS
C
            U          = Q(I,2)
            V          = Q(I,3)
            W          = Q(I,4)
            UBAR       = S(I,1)*U + S(I,2)*V + S(I,3)*W
C
            TAUXX      = Q(I, 6)
            TAUYY      = Q(I, 7)
            TAUZZ      = Q(I, 8)
            TAUXY      = Q(I, 9)
            TAUXZ      = Q(I,10)
            TAUYZ      = Q(I,11)
            T2         = Q(I,12)
C
            A(I, 1, 6) =   0.0E0
            A(I, 1, 7) =   0.0E0
            A(I, 1, 8) =   0.0E0
            A(I, 1, 9) =   0.0E0
            A(I, 1,10) =   0.0E0
            A(I, 1,11) =   0.0E0
            A(I, 1,12) =   0.0E0
C
            A(I, 2, 6) = - 0.5E0*S(I,4)*S(I,1)*GAMM1
            A(I, 2, 7) = - 0.5E0*S(I,4)*S(I,1)*GAMM1
            A(I, 2, 8) = - 0.5E0*S(I,4)*S(I,1)*GAMM1
            A(I, 2, 9) =   0.0E0
            A(I, 2,10) =   0.0E0
            A(I, 2,11) =   0.0E0
            A(I, 2,12) =   0.0E0
C
            A(I, 3, 6) = - 0.5E0*S(I,4)*S(I,2)*GAMM1
            A(I, 3, 7) = - 0.5E0*S(I,4)*S(I,2)*GAMM1
            A(I, 3, 8) = - 0.5E0*S(I,4)*S(I,2)*GAMM1
            A(I, 3, 9) =   0.0E0
            A(I, 3,10) =   0.0E0
            A(I, 3,11) =   0.0E0
            A(I, 3,12) =   0.0E0
C
            A(I, 4, 6) = - 0.5E0*S(I,4)*S(I,3)*GAMM1
            A(I, 4, 7) = - 0.5E0*S(I,4)*S(I,3)*GAMM1
            A(I, 4, 8) = - 0.5E0*S(I,4)*S(I,3)*GAMM1
            A(I, 4, 9) =   0.0E0
            A(I, 4,10) =   0.0E0
            A(I, 4,11) =   0.0E0
            A(I, 4,12) =   0.0E0
C
            A(I, 5, 6) = - 0.5E0*S(I,4)*UBAR*GAMM1
            A(I, 5, 7) = - 0.5E0*S(I,4)*UBAR*GAMM1
            A(I, 5, 8) = - 0.5E0*S(I,4)*UBAR*GAMM1
            A(I, 5, 9) =   0.0E0
            A(I, 5,10) =   0.0E0
            A(I, 5,11) =   0.0E0
            A(I, 5,12) =   0.0E0
C
            A(I, 6, 1) = - S(I,4)*UBAR*TAUXX
            A(I, 6, 2) =   S(I,4)*S(I,1)*TAUXX
            A(I, 6, 3) =   S(I,4)*S(I,2)*TAUXX
            A(I, 6, 4) =   S(I,4)*S(I,3)*TAUXX
            A(I, 6, 5) =   0.0E0
            A(I, 6, 6) =   S(I,4)*UBAR
            A(I, 6, 7) =   0.0E0
            A(I, 6, 8) =   0.0E0
            A(I, 6, 9) =   0.0E0
            A(I, 6,10) =   0.0E0
            A(I, 6,11) =   0.0E0
            A(I, 6,12) =   0.0E0
C
            A(I, 7, 1) = - S(I,4)*UBAR*TAUYY
            A(I, 7, 2) =   S(I,4)*S(I,1)*TAUYY
            A(I, 7, 3) =   S(I,4)*S(I,2)*TAUYY
            A(I, 7, 4) =   S(I,4)*S(I,3)*TAUYY
            A(I, 7, 5) =   0.0E0
            A(I, 7, 6) =   0.0E0
            A(I, 7, 7) =   S(I,4)*UBAR
            A(I, 7, 8) =   0.0E0
            A(I, 7, 9) =   0.0E0
            A(I, 7,10) =   0.0E0
            A(I, 7,11) =   0.0E0
            A(I, 7,12) =   0.0E0
C
            A(I, 8, 1) = - S(I,4)*UBAR*TAUZZ
            A(I, 8, 2) =   S(I,4)*S(I,1)*TAUZZ
            A(I, 8, 3) =   S(I,4)*S(I,2)*TAUZZ
            A(I, 8, 4) =   S(I,4)*S(I,3)*TAUZZ
            A(I, 8, 5) =   0.0E0
            A(I, 8, 6) =   0.0E0
            A(I, 8, 7) =   0.0E0
            A(I, 8, 8) =   S(I,4)*UBAR
            A(I, 8, 9) =   0.0E0
            A(I, 8,10) =   0.0E0
            A(I, 8,11) =   0.0E0
            A(I, 8,12) =   0.0E0
C
            A(I, 9, 1) = - S(I,4)*UBAR*TAUXY
            A(I, 9, 2) =   S(I,4)*S(I,1)*TAUXY
            A(I, 9, 3) =   S(I,4)*S(I,2)*TAUXY
            A(I, 9, 4) =   S(I,4)*S(I,3)*TAUXY
            A(I, 9, 5) =   0.0E0
            A(I, 9, 6) =   0.0E0
            A(I, 9, 7) =   0.0E0
            A(I, 9, 8) =   0.0E0
            A(I, 9, 9) =   S(I,4)*UBAR
            A(I, 9,10) =   0.0E0
            A(I, 9,11) =   0.0E0
            A(I, 9,12) =   0.0E0
C
            A(I,10, 1) = - S(I,4)*UBAR*TAUXZ
            A(I,10, 2) =   S(I,4)*S(I,1)*TAUXZ
            A(I,10, 3) =   S(I,4)*S(I,2)*TAUXZ
            A(I,10, 4) =   S(I,4)*S(I,3)*TAUXZ
            A(I,10, 5) =   0.0E0
            A(I,10, 6) =   0.0E0
            A(I,10, 7) =   0.0E0
            A(I,10, 8) =   0.0E0
            A(I,10, 9) =   0.0E0
            A(I,10,10) =   S(I,4)*UBAR
            A(I,10,11) =   0.0E0
            A(I,10,12) =   0.0E0
C
            A(I,11, 1) = - S(I,4)*UBAR*TAUYZ
            A(I,11, 2) =   S(I,4)*S(I,1)*TAUYZ
            A(I,11, 3) =   S(I,4)*S(I,2)*TAUYZ
            A(I,11, 4) =   S(I,4)*S(I,3)*TAUYZ
            A(I,11, 5) =   0.0E0
            A(I,11, 6) =   0.0E0
            A(I,11, 7) =   0.0E0
            A(I,11, 8) =   0.0E0
            A(I,11, 9) =   0.0E0
            A(I,11,10) =   0.0E0
            A(I,11,11) =   S(I,4)*UBAR
            A(I,11,12) =   0.0E0
C
            A(I,12, 1) = - S(I,4)*UBAR*T2
            A(I,12, 2) =   S(I,4)*S(I,1)*T2
            A(I,12, 3) =   S(I,4)*S(I,2)*T2
            A(I,12, 4) =   S(I,4)*S(I,3)*T2
            A(I,12, 5) =   0.0E0
            A(I,12, 6) =   0.0E0
            A(I,12, 7) =   0.0E0
            A(I,12, 8) =   0.0E0
            A(I,12, 9) =   0.0E0
            A(I,12,10) =   0.0E0
            A(I,12,11) =   0.0E0
            A(I,12,12) =   S(I,4)*UBAR
  420    CONTINUE
      ENDIF
C
C     Finished with DF/DQ calculation.
C
      RETURN
      END





      SUBROUTINE DFROE (ITURB, NPTS, QP, QM, S, AP, AM, AR, WORK)
C
C     Routine to calculate Jacobians by calling the 
C     appropriate Roe's implicit routine.
C
C
C     ITURB          : Turbulence model
C     NPTS           : Size of array to calculate fluxes on.
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     S              : Metrics
C     AP,AM          : The Jacobian matrices DF(Q+)/DQ and DF(Q-)/DQ
C     AR             : The Jacobian matrix T |Lambda| T(-1) using Roe's avg.
C     WORK           : Temporary storage for TKE,H in Jacobian routines
C
      include '../header/common.h'
C
      DIMENSION QP    (NPTS,NQ),
     1          QM    (NPTS,NQ),
     2          S     (NPTS,4),
     3          AP    (NPTS,NRANK,NRANK),
     4          AM    (NPTS,NRANK,NRANK),
     5          AR    (NPTS,NRANK,NRANK),
     6          WORK  (NPTS,4)
C
C     Calculate dF/dQ(QP), dF/dQ(QM) and T |Lambda| T[-1]
C
      CALL DFDQ   (ITURB, NPTS, QP, S, WORK(1,1), WORK(1,2), AP)
      CALL DFDQ   (ITURB, NPTS, QM, S, WORK(1,1), WORK(1,2), AM)
      CALL TLTROE (ITURB, NPTS, QP, QM, S, WORK(1,1), WORK(1,2), 
     1             WORK(1,3), WORK(1,4), AR)
C
C     Finished calculating Jacobians for Roe's scheme
C
      RETURN
      END






      SUBROUTINE DFVISJ (IDIM, JDIM, KDIM, I, ITURB,
     1                   Q, RMU, RMUT, S, VOL, IBCDAT, NBCS,
     2                   T, U, V, W, RMUJ, RMUTJ, AP, AM)
C
C     Routine to calculate the Jacobians of the viscous
C     terms for the J direction.  This routine calculates the terms
C     from the mean equations and adds the terms due to the 
C     turbulence model.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : Current I Plane
C     ITURB          : Turbulence model
C     Q              : Primitive variables at cell centers
C     RMU            : Molecular viscosity at cell centers
C     RMUT           : Turbulent eddy viscosity at cell centers
C     S              : Metrics
C     VOL            : Cell volumes
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     T              : Temporary storage of Temperature
C     U,V,W          : Temporary storage of U,V,W velocities
C     RMUJ           : Temporary storage of Viscosity at cell interfaces
C     RMUTJ          : Temporary storage of turbulent eddy viscosity 
C                      at cell interfaces
C     AP,AM          : Implicit terms
C                      AP = dF(i+1/2)/dQ(i+1)
C                      AM = dF(i+1/2)/dQ(i)
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          RMU   (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     2          RMUT  (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     3          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     4          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION T     (JDIM+1,2:KDIM),  
     1          RMUJ  (JDIM,2:KDIM), 
     2          RMUTJ (JDIM,2:KDIM),
     3          U     (JDIM,2:KDIM), 
     4          V     (JDIM,2:KDIM), 
     5          W     (JDIM,2:KDIM)
C
      DIMENSION AP    (JDIM,2:KDIM,NRANK,NRANK),
     1          AM    (JDIM,2:KDIM,NRANK,NRANK)
C
      DIMENSION IBCDAT(10,NBCS)
C
C     Zero the matrices
C
      DO 114 L = 1, NRANK
         DO 113 M = 1, NRANK
            DO 112 K = 2, KDIM
               DO 111 J = 1, JDIM
                  AP(J,K,L,M) = 0.E0
                  AM(J,K,L,M) = 0.E0
  111          CONTINUE
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
C
C     Calculate the thin-layer Navier-Stokes terms for the mean equations
C
      CALL DFVTLJ (IDIM, JDIM, KDIM, I, ITURB,
     1             Q, RMU, RMUT, S, VOL, IBCDAT, NBCS,
     2             T, U, V, W, RMUJ, RMUTJ, AP, AM)
C
C     Finished with implicit Thin Layer Viscous Terms in the J direction
C
      RETURN
      END






      SUBROUTINE DFVISK (IDIM, JDIM, KDIM, I, ITURB,
     1                   Q, RMU, RMUT, S, VOL, IBCDAT, NBCS,
     2                   T, U, V, W, RMUK, RMUTK, AP, AM)
C
C     Routine to calculate the thin layer Navier-Stokes terms for
C     the K direction.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : Current I Plane
C     ITURB          : Turbulence model
C     Q              : Primitive variables at cell centers
C     RMU            : Molecular viscosity at cell centers
C     RMUT           : Turbulent eddy viscosity at cell centers
C     S              : Metrics
C     VOL            : Cell volumes
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     T              : Temporary storage of Temperature
C     U,V,W          : Temporary storage of U,V,W velocities
C     RMUK           : Temporary storage of Viscosity at cell interfaces
C     RMUTK          : Temporary storage of turbulent eddy viscosity
C                      at cell interfaces
C     AP,AM          : Implicit terms
C                      AP = dF(i+1/2)/dQ(i+1)
C                      AM = dF(i+1/2)/dQ(i)
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          RMU   (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     2          RMUT  (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     3          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     4          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION T     (KDIM+1,2:JDIM),  
     1          RMUK  (KDIM,2:JDIM), 
     2          RMUTK (KDIM,2:JDIM),
     3          U     (KDIM,2:JDIM), 
     4          V     (KDIM,2:JDIM), 
     5          W     (KDIM,2:JDIM)
C
      DIMENSION AP    (KDIM,2:JDIM,NRANK,NRANK),
     1          AM    (KDIM,2:JDIM,NRANK,NRANK)
C
      DIMENSION IBCDAT(10,NBCS)
C
C     Zero the matrices
C
      DO 114 L = 1, NRANK
         DO 113 M = 1, NRANK
            DO 112 J = 2, JDIM
               DO 111 K = 1, KDIM
                  AP(K,J,L,M) = 0.E0
                  AM(K,J,L,M) = 0.E0
  111          CONTINUE
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
C
C     Calculate the implicit thin-layer Navier-Stokes terms for the mean eqns
C
      CALL DFVTLK (IDIM, JDIM, KDIM, I, 
     1             Q, RMU, RMUT, S, VOL, IBCDAT, NBCS,
     2             T, U, V, W, RMUK, RMUTK, AP, AM)
C
C     Finished 
C
      RETURN
      END






      SUBROUTINE DFVTLI (IDIM, JDIM, KDIM, K, IPSTRT, IPEND,
     1                   Q, RMUT, S, VOL, IBCDAT, NBCS,
     2                   T, U, V, W, RMU, RMUI, RMUTI, AP, AM)
C
C     Routine to calculate the implicit thin layer Navier-Stokes terms 
C     for the I direction.
C
C     This routine is not presently suited for the multiblock
C     approach.  It is close to being ready but is not.  The problem
C     lies that the metrics for one cell into the adjacent block
C     are needed to be completely general.  These are currently not available.
C     The routine metbc needs to be completed.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     K              : Current K Plane
C     IPSTRT         : I Plane (cell center) to start calculations at 
C                      (Minimum Value: 2)
C     IPEND          : I Plane (cell center) to end caclulations at
C                      (Maximum Value: IDIM)
C     Q              : Primitive variables at cell centers
C     RMUT           : Turbulent eddy viscosity at cell centers
C     S              : Metrics
C     VOL            : Cell volumes
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     T              : Temporary storage of Temperature
C     U,V,W          : Temporary storage of U,V,W velocities
C     RMU            : Temporary storage of Viscosity at cell centers
C     RMUI           : Temporary storage of Viscosity at cell interfaces
C     RMUTI          : Temporary storage of turbulent eddy viscosity 
C                      at cell interfaces
C     AP,AM          : Implicit terms
C                      AP = dF(i+1/2)/dQ(i+1)
C                      AM = dF(i+1/2)/dQ(i)
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          RMUT  (0:JDIM+2,0:KDIM+2,0:IDIM+2), 
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION T     (IPSTRT-1:IPEND+1,2:JDIM), 
     1          RMU   (IPSTRT-1:IPEND+1,2:JDIM), 
     2          RMUI  (IPSTRT-1:IPEND,2:JDIM), 
     3          RMUTI (IPSTRT-1:IPEND,2:JDIM),
     4          U     (IPSTRT-1:IPEND,2:JDIM), 
     5          V     (IPSTRT-1:IPEND,2:JDIM), 
     6          W     (IPSTRT-1:IPEND,2:JDIM)
C
      DIMENSION AP    (IPSTRT-1:IPEND,2:JDIM,NRANK,NRANK),
     1          AM    (IPSTRT-1:IPEND,2:JDIM,NRANK,NRANK)
C
      DIMENSION IBCDAT(10,NBCS)
C
C     Two-Thirds
C     Four-Thirds
C
      T23 = 2.E0 / 3.E0
      T43 = 4.E0 / 3.E0
C
C     Calculate T = T / TINF  at the cell centers
C
      DO 120 J = 2, JDIM
         DO 110 I = IPSTRT - 1, IPEND + 1
            T(I,J) = GAMMA * Q(J,K,I,5) / Q(J,K,I,1)
  110    CONTINUE
  120 CONTINUE
C
C     Calculate molecular viscosity at cell centers
C
      NPTS = (IPEND - IPSTRT + 1) * (JDIM - 1)
      CALL VSCSTY (NPTS, T, RMU)
C
C     Calculate U,V,W,RMUI,RMUTI at cell interfaces as averages of cell centers
C
      DO 220 J = 2, JDIM
         DO 210 I = IPSTRT - 1, IPEND
            U(I,J)     = 0.5E0 * ( Q(J,K,I,2)  + Q(J,K,I+1,2) )
            V(I,J)     = 0.5E0 * ( Q(J,K,I,3)  + Q(J,K,I+1,3) )
            W(I,J)     = 0.5E0 * ( Q(J,K,I,4)  + Q(J,K,I+1,4) )
            RMUI(I,J)  = 0.5E0 * ( RMU(I,J)    + RMU(I+1,J) )
            RMUTI(I,J) = 0.5E0 * ( RMUT(J,K,I) + RMUT(J,K,I+1) )
  210    CONTINUE
  220 CONTINUE
C
C     Loop over boundary conditions and fix the wall data
C
      DO 350 IBC = 1, NBCS
         IF (IBCDAT(1,IBC) .EQ. IWALL .AND. 
     1       IBCDAT(2,IBC) .EQ. 1) THEN
            INROUT = IBCDAT(3,IBC)
            JSTRT  = IBCDAT(6,IBC)
            JEND   = IBCDAT(7,IBC)
            KSTRT  = IBCDAT(8,IBC)
            KEND   = IBCDAT(9,IBC)
            IF (INROUT .EQ. 1) THEN
               IBCPT = 1
               INTR  = 1
            ELSE
               IBCPT = IDIM + 1
               INTR  = IDIM
            ENDIF
            IF (IBCPT .GE. IPSTRT-1 .AND. IBCPT .LE. IPEND .AND.
     1              K .GE. KSTRT    .AND.     K .LE. KEND      ) THEN
               DO 320 J = JSTRT, JEND
                  U(INTR,J)     = 0.E0
                  V(INTR,J)     = 0.E0
                  W(INTR,J)     = 0.E0
                  RMUI(INTR,J)  = RMU(IBCPT,J)
                  RMUTI(INTR,J) = RMUT(J,K,IBCPT)
  320          CONTINUE
            ENDIF
         ENDIF
  350 CONTINUE
C
C     Loop over J forming viscous terms along the I direction
C     on the I = IPSTRT - 1, IPEND faces 
C
C
      DO 520 J = 2, JDIM
         DO 510 I = IPSTRT - 1, IPEND
C
            AREAMR = S(J,K,I,4,1) * FSMACH / RE
C
C     Average metrics and cell volume
C
C     VOLAUX is the inverse of the average auxilliary cell volume
C
            VOLAUX = 2.0E0 / ( VOL(J,K,I) + VOL(J,K,I+1) )
            SPI1   = 0.5E0 * ( S(J,K,I,1,1)*S(J,K,I,4,1)
     1                       + S(J,K,I+1,1,1)*S(J,K,I+1,4,1) )
            SPI2   = 0.5E0 * ( S(J,K,I,2,1)*S(J,K,I,4,1)
     1                       + S(J,K,I+1,2,1)*S(J,K,I+1,4,1) )
            SPI3   = 0.5E0 * ( S(J,K,I,3,1)*S(J,K,I,4,1)
     1                       + S(J,K,I+1,3,1)*S(J,K,I+1,4,1) )
            SMI1   = 0.5E0 * ( S(J,K,I,1,1)*S(J,K,I,4,1)
     1                       + S(J,K,I-1,1,1)*S(J,K,I-1,4,1) )
            SMI2   = 0.5E0 * ( S(J,K,I,2,1)*S(J,K,I,4,1)
     1                       + S(J,K,I-1,2,1)*S(J,K,I-1,4,1) )
            SMI3   = 0.5E0 * ( S(J,K,I,3,1)*S(J,K,I,4,1)
     1                       + S(J,K,I-1,3,1)*S(J,K,I-1,4,1) )
C
            DUDX = VOLAUX * ( Q(J,K,I+1,2)*SPI1 - Q(J,K,I,2)*SMI1 )
            DUDY = VOLAUX * ( Q(J,K,I+1,2)*SPI2 - Q(J,K,I,2)*SMI2 )
            DUDZ = VOLAUX * ( Q(J,K,I+1,2)*SPI3 - Q(J,K,I,2)*SMI3 )
C
            DVDX = VOLAUX * ( Q(J,K,I+1,3)*SPI1 - Q(J,K,I,3)*SMI1 )
            DVDY = VOLAUX * ( Q(J,K,I+1,3)*SPI2 - Q(J,K,I,3)*SMI2 )
            DVDZ = VOLAUX * ( Q(J,K,I+1,3)*SPI3 - Q(J,K,I,3)*SMI3 )
C
            DWDX = VOLAUX * ( Q(J,K,I+1,4)*SPI1 - Q(J,K,I,4)*SMI1 )
            DWDY = VOLAUX * ( Q(J,K,I+1,4)*SPI2 - Q(J,K,I,4)*SMI2 )
            DWDZ = VOLAUX * ( Q(J,K,I+1,4)*SPI3 - Q(J,K,I,4)*SMI3 )
C
            DTDX = VOLAUX * ( T(I+1,J)*SPI1 - T(I,J)*SMI1 )
            DTDY = VOLAUX * ( T(I+1,J)*SPI2 - T(I,J)*SMI2 )
            DTDZ = VOLAUX * ( T(I+1,J)*SPI3 - T(I,J)*SMI3 )
C
C     Calculate Shear Terms
C     Use molecular viscosity (RMUI) and turbulent eddy viscosity (RMUTI)
C     averaged from cell centers to get viscosity at cell interface.
C     T23 is two-thirds
C
            VISC  = RMUI(I,J) + RMUTI(I,J)
            TAUXX = T23 * VISC * ( 2.E0*DUDX - DVDY - DWDZ )
            TAUYY = T23 * VISC * ( 2.E0*DVDY - DUDX - DWDZ )
            TAUZZ = T23 * VISC * ( 2.E0*DWDZ - DUDX - DVDY )
            TAUXY = VISC * ( DUDY + DVDX )
            TAUXZ = VISC * ( DWDX + DUDZ )
            TAUYZ = VISC * ( DWDY + DVDZ )
C
C     Heat Flux terms are not used in the current implementation of the
C     implicit viscous terms.  This is kept for compatibility reasons 
C     and in case they are later included in d(RMU) terms.
C     Calculate Heat Flux Terms
C
C           CNDCT = (RMUI(I,J)/PR + RMUTI(I,J)/PRT) / GAMM1
C           QDOTX = - CNDCT * DTDX
C           QDOTY = - CNDCT * DTDY
C           QDOTZ = - CNDCT * DTDZ
C
C     Region to calculate AM = dF(i+1/2)/dQ(i)  terms
C
C     Calculate DTau(xx,xy,xz,...)/DQ(1,2,3,...) terms
C
            RMRV   = VISC * VOLAUX / Q(J,K,I,1)
C
C     DTau(xx)
C
            DTXXD1 =   T23 * RMRV * ( 2.E0*SMI1*Q(J,K,I,2) 
     1                       - SMI2*Q(J,K,I,3) - SMI3*Q(J,K,I,4) )
            DTXXD2 = - T43 * RMRV * SMI1
            DTXXD3 =   T23 * RMRV * SMI2
            DTXXD4 =   T23 * RMRV * SMI3
C
C     DTau(yy)
C
            DTYYD1 =   T23 * RMRV * ( 2.E0*SMI2*Q(J,K,I,3) 
     1                       - SMI1*Q(J,K,I,2) - SMI3*Q(J,K,I,4) )
            DTYYD2 =   T23 * RMRV * SMI1
            DTYYD3 = - T43 * RMRV * SMI2
            DTYYD4 =   T23 * RMRV * SMI3
C
C     DTau(zz)
C
            DTZZD1 =   T23 * RMRV * ( 2.E0*SMI3*Q(J,K,I,4) 
     1                       - SMI1*Q(J,K,I,2) - SMI2*Q(J,K,I,3) )
            DTZZD2 =   T23 * RMRV * SMI1
            DTZZD3 =   T23 * RMRV * SMI2
            DTZZD4 = - T43 * RMRV * SMI3
C
C     DTau(xy)
C
            DTXYD1 =   RMRV * ( SMI1*Q(J,K,I,3) + SMI2*Q(J,K,I,2) )
            DTXYD2 = - RMRV * SMI2
            DTXYD3 = - RMRV * SMI1
C
C     DTau(xz)
C
            DTXZD1 =   RMRV * ( SMI1*Q(J,K,I,4) + SMI3*Q(J,K,I,2) )
            DTXZD2 = - RMRV * SMI3
            DTXZD4 = - RMRV * SMI1
C
C     DTau(yz)
C
            DTYZD1 =   RMRV * ( SMI2*Q(J,K,I,4) + SMI3*Q(J,K,I,3) )
            DTYZD3 = - RMRV * SMI3
            DTYZD4 = - RMRV * SMI2
C
C     DQdot(x)
C
            TEMP   = (RMUI(I,J)/PR + RMUTI(I,J)/PRT) * GAMMA * VOLAUX 
     1               / Q(J,K,I,1)
C
            DQDXD1 = - TEMP * SMI1 * ( Q(J,K,I,5)/GAMM1 
     1                   - 0.5E0*Q(J,K,I,1)*( Q(J,K,I,2)*Q(J,K,I,2)
     2                 +Q(J,K,I,3)*Q(J,K,I,3)+Q(J,K,I,4)*Q(J,K,I,4) ) )
     3                    / Q(J,K,I,1)
            DQDXD2 = - TEMP * SMI1 * Q(J,K,I,2) 
            DQDXD3 = - TEMP * SMI1 * Q(J,K,I,3) 
            DQDXD4 = - TEMP * SMI1 * Q(J,K,I,4) 
            DQDXD5 =   TEMP * SMI1 
C
C     DQdot(y)
C
            DQDYD1 = - TEMP * SMI2 * ( Q(J,K,I,5)/GAMM1
     1                   - 0.5E0*Q(J,K,I,1)*( Q(J,K,I,2)*Q(J,K,I,2)
     2                 +Q(J,K,I,3)*Q(J,K,I,3)+Q(J,K,I,4)*Q(J,K,I,4) ) )
     3                    / Q(J,K,I,1)
            DQDYD2 = - TEMP * SMI2 * Q(J,K,I,2)
            DQDYD3 = - TEMP * SMI2 * Q(J,K,I,3) 
            DQDYD4 = - TEMP * SMI2 * Q(J,K,I,4) 
            DQDYD5 =   TEMP * SMI2 
C
C     DQdot(z)
C
            DQDZD1 = - TEMP * SMI3 * ( Q(J,K,I,5)/GAMM1
     1                   - 0.5E0*Q(J,K,I,1)*( Q(J,K,I,2)*Q(J,K,I,2)
     2                 +Q(J,K,I,3)*Q(J,K,I,3)+Q(J,K,I,4)*Q(J,K,I,4) ) )
     3                    / Q(J,K,I,1)
            DQDZD2 = - TEMP * SMI3 * Q(J,K,I,2)
            DQDZD3 = - TEMP * SMI3 * Q(J,K,I,3)
            DQDZD4 = - TEMP * SMI3 * Q(J,K,I,4)
            DQDZD5 =   TEMP * SMI3
C
C     Calculate AM = dF(i+1/2)/dQ(i)
C
            AM(I,J,1,1) = 0.E0
            AM(I,J,1,2) = 0.E0
            AM(I,J,1,3) = 0.E0
            AM(I,J,1,4) = 0.E0
            AM(I,J,1,5) = 0.E0
C
            AM(I,J,2,1) = ( S(J,K,I,1,1)*DTXXD1 + S(J,K,I,2,1)*DTXYD1
     1                    + S(J,K,I,3,1)*DTXZD1 ) * AREAMR
            AM(I,J,2,2) = ( S(J,K,I,1,1)*DTXXD2 + S(J,K,I,2,1)*DTXYD2
     1                    + S(J,K,I,3,1)*DTXZD2 ) * AREAMR
            AM(I,J,2,3) = ( S(J,K,I,1,1)*DTXXD3 + S(J,K,I,2,1)*DTXYD3 )
     1                                            * AREAMR
            AM(I,J,2,4) = ( S(J,K,I,1,1)*DTXXD4 + S(J,K,I,3,1)*DTXZD4 ) 
     1                                            * AREAMR
            AM(I,J,2,5) = 0.E0
C
            AM(I,J,3,1) = ( S(J,K,I,1,1)*DTXYD1 + S(J,K,I,2,1)*DTYYD1
     1                    + S(J,K,I,3,1)*DTYZD1 ) * AREAMR
            AM(I,J,3,2) = ( S(J,K,I,1,1)*DTXYD2 + S(J,K,I,2,1)*DTYYD2 )
     1                                            * AREAMR
            AM(I,J,3,3) = ( S(J,K,I,1,1)*DTXYD3 + S(J,K,I,2,1)*DTYYD3
     1                    + S(J,K,I,3,1)*DTYZD3 ) * AREAMR
            AM(I,J,3,4) = ( S(J,K,I,2,1)*DTYYD4 + S(J,K,I,3,1)*DTYZD4 ) 
     1                                            * AREAMR
            AM(I,J,3,5) = 0.E0
C
            AM(I,J,4,1) = ( S(J,K,I,1,1)*DTXZD1 + S(J,K,I,2,1)*DTYZD1
     1                    + S(J,K,I,3,1)*DTZZD1 ) * AREAMR
            AM(I,J,4,2) = ( S(J,K,I,1,1)*DTXZD2 + S(J,K,I,3,1)*DTZZD2 ) 
     1                                            * AREAMR
            AM(I,J,4,3) = ( S(J,K,I,2,1)*DTYZD3 + S(J,K,I,3,1)*DTZZD3 ) 
     1                                            * AREAMR
            AM(I,J,4,4) = ( S(J,K,I,1,1)*DTXZD4 + S(J,K,I,2,1)*DTYZD4
     1                    + S(J,K,I,3,1)*DTZZD4 ) * AREAMR
            AM(I,J,4,5) = 0.E0
C
            AM(I,J,5,1) = AREAMR *
     1 ( S(J,K,I,1,1)*(U(I,J)*DTXXD1+V(I,J)*DTXYD1+W(I,J)*DTXZD1-DQDXD1)
     2 + S(J,K,I,2,1)*(U(I,J)*DTXYD1+V(I,J)*DTYYD1+W(I,J)*DTYZD1-DQDYD1)
     3 + S(J,K,I,3,1)*(U(I,J)*DTXZD1+V(I,J)*DTYZD1+W(I,J)*DTZZD1-DQDZD1)
     4 - 0.5E0*Q(J,K,I,2)/Q(J,K,I,1)*(S(J,K,I,1,1)*TAUXX
     5                    +S(J,K,I,2,1)*TAUXY+S(J,K,I,3,1)*TAUXZ)
     6 - 0.5E0*Q(J,K,I,3)/Q(J,K,I,1)*(S(J,K,I,1,1)*TAUXY
     7                    +S(J,K,I,2,1)*TAUYY+S(J,K,I,3,1)*TAUYZ)
     8 - 0.5E0*Q(J,K,I,4)/Q(J,K,I,1)*(S(J,K,I,1,1)*TAUXZ
     9                    +S(J,K,I,2,1)*TAUYZ+S(J,K,I,3,1)*TAUZZ) )
            AM(I,J,5,2) = AREAMR *
     1 ( S(J,K,I,1,1)*(U(I,J)*DTXXD2+V(I,J)*DTXYD2+W(I,J)*DTXZD2-DQDXD2)
     2 + S(J,K,I,2,1)*(U(I,J)*DTXYD2+V(I,J)*DTYYD2              -DQDYD2)
     3 + S(J,K,I,3,1)*(U(I,J)*DTXZD2              +W(I,J)*DTZZD2-DQDZD2)
     4 + 0.5E0/Q(J,K,I,1)*(S(J,K,I,1,1)*TAUXX
     5                    +S(J,K,I,2,1)*TAUXY+S(J,K,I,3,1)*TAUXZ) )
            AM(I,J,5,3) = AREAMR *
     1 ( S(J,K,I,1,1)*(U(I,J)*DTXXD3+V(I,J)*DTXYD3              -DQDXD3)
     2 + S(J,K,I,2,1)*(U(I,J)*DTXYD3+V(I,J)*DTYYD3+W(I,J)*DTYZD3-DQDYD3)
     3 + S(J,K,I,3,1)*(              V(I,J)*DTYZD3+W(I,J)*DTZZD3-DQDZD3)
     4 + 0.5E0/Q(J,K,I,1)*(S(J,K,I,1,1)*TAUXY
     5                    +S(J,K,I,2,1)*TAUYY+S(J,K,I,3,1)*TAUYZ) )
            AM(I,J,5,4) = AREAMR *
     1 ( S(J,K,I,1,1)*(U(I,J)*DTXXD4              +W(I,J)*DTXZD4-DQDXD4)
     2 + S(J,K,I,2,1)*(              V(I,J)*DTYYD4+W(I,J)*DTYZD4-DQDYD4)
     3 + S(J,K,I,3,1)*(U(I,J)*DTXZD4+V(I,J)*DTYZD4+W(I,J)*DTZZD4-DQDZD4)
     4 + 0.5E0/Q(J,K,I,1)*(S(J,K,I,1,1)*TAUXZ
     5                    +S(J,K,I,2,1)*TAUYZ+S(J,K,I,3,1)*TAUZZ) )
            AM(I,J,5,5) = - AREAMR *
     1                      ( S(J,K,I,1,1)*DQDXD5 + S(J,K,I,2,1)*DQDYD5 
     2                      + S(J,K,I,3,1)*DQDZD5 )
C
C     Region to calculate AP = dF(i+1/2)/dQ(i)  terms
C
C     Calculate DTau(xx,xy,xz,...)/DQ(1,2,3,...) terms
C
            RMRV   = VISC * VOLAUX / Q(J,K,I+1,1)
C
C     DTau(xx)
C
            DTXXD1 = - T23 * RMRV * ( 2.E0*SPI1*Q(J,K,I+1,2) 
     1                     - SPI2*Q(J,K,I+1,3) - SPI3*Q(J,K,I+1,4) )
            DTXXD2 =   T43 * RMRV * SPI1
            DTXXD3 = - T23 * RMRV * SPI2
            DTXXD4 = - T23 * RMRV * SPI3
C
C     DTau(yy)
C
            DTYYD1 = - T23 * RMRV * ( 2.E0*SPI2*Q(J,K,I+1,3) 
     1                     - SPI1*Q(J,K,I+1,2) - SPI3*Q(J,K,I+1,4) )
            DTYYD2 = - T23 * RMRV * SPI1
            DTYYD3 =   T43 * RMRV * SPI2
            DTYYD4 = - T23 * RMRV * SPI3
C
C     DTau(zz)
C
            DTZZD1 = - T23 * RMRV * ( 2.E0*SPI3*Q(J,K,I+1,4) 
     1                     - SPI1*Q(J,K,I+1,2) - SPI2*Q(J,K,I+1,3) )
            DTZZD2 = - T23 * RMRV * SPI1
            DTZZD3 = - T23 * RMRV * SPI2
            DTZZD4 =   T43 * RMRV * SPI3
C
C     DTau(xy)
C
            DTXYD1 = - RMRV * ( SPI1*Q(J,K,I+1,3) + SPI2*Q(J,K,I+1,2) )
            DTXYD2 =   RMRV * SPI2
            DTXYD3 =   RMRV * SPI1
C
C     DTau(xz)
C
            DTXZD1 = - RMRV * ( SPI1*Q(J,K,I+1,4) + SPI3*Q(J,K,I+1,2) )
            DTXZD2 =   RMRV * SPI3
            DTXZD4 =   RMRV * SPI1
C
C     DTau(yz)
C
            DTYZD1 = - RMRV * ( SPI2*Q(J,K,I+1,4) + SPI3*Q(J,K,I+1,3) )
            DTYZD3 =   RMRV * SPI3
            DTYZD4 =   RMRV * SPI2
C
C     DQdot(x)
C
            TEMP   = (RMUI(I,J)/PR + RMUTI(I,J)/PRT) * GAMMA * VOLAUX 
     1               / Q(J,K,I+1,1)
C
            DQDXD1 =   TEMP * SPI1 * ( Q(J,K,I+1,5)/GAMM1 
     1             - 0.5E0*Q(J,K,I+1,1)*( Q(J,K,I+1,2)*Q(J,K,I+1,2)
     2         +Q(J,K,I+1,3)*Q(J,K,I+1,3)+Q(J,K,I+1,4)*Q(J,K,I+1,4) ) )
     3                    / Q(J,K,I+1,1)
            DQDXD2 =   TEMP * SPI1 * Q(J,K,I+1,2) 
            DQDXD3 =   TEMP * SPI1 * Q(J,K,I+1,3) 
            DQDXD4 =   TEMP * SPI1 * Q(J,K,I+1,4) 
            DQDXD5 = - TEMP * SPI1 
C
C     DQdot(y)
C
            DQDYD1 =   TEMP * SPI2 * ( Q(J,K,I+1,5)/GAMM1
     1             - 0.5E0*Q(J,K,I+1,1)*( Q(J,K,I+1,2)*Q(J,K,I+1,2)
     2         +Q(J,K,I+1,3)*Q(J,K,I+1,3)+Q(J,K,I+1,4)*Q(J,K,I+1,4) ) )
     3                    / Q(J,K,I+1,1)
            DQDYD2 =   TEMP * SPI2 * Q(J,K,I+1,2)
            DQDYD3 =   TEMP * SPI2 * Q(J,K,I+1,3) 
            DQDYD4 =   TEMP * SPI2 * Q(J,K,I+1,4) 
            DQDYD5 = - TEMP * SPI2 
C
C     DQdot(z)
C
            DQDZD1 =   TEMP * SPI3 * ( Q(J,K,I+1,5)/GAMM1
     1             - 0.5E0*Q(J,K,I+1,1)*( Q(J,K,I+1,2)*Q(J,K,I+1,2)
     2         +Q(J,K,I+1,3)*Q(J,K,I+1,3)+Q(J,K,I+1,4)*Q(J,K,I+1,4) ) )
     3                    / Q(J,K,I+1,1)
            DQDZD2 =   TEMP * SPI3 * Q(J,K,I+1,2)
            DQDZD3 =   TEMP * SPI3 * Q(J,K,I+1,3)
            DQDZD4 =   TEMP * SPI3 * Q(J,K,I+1,4)
            DQDZD5 = - TEMP * SPI3
C
C     Calculate AP = dF(i+1/2)/dQ(i+1)
C
            AP(I,J,1,1) = 0.E0
            AP(I,J,1,2) = 0.E0
            AP(I,J,1,3) = 0.E0
            AP(I,J,1,4) = 0.E0
            AP(I,J,1,5) = 0.E0
C
            AP(I,J,2,1) = ( S(J,K,I,1,1)*DTXXD1 + S(J,K,I,2,1)*DTXYD1
     1                    + S(J,K,I,3,1)*DTXZD1 ) * AREAMR
            AP(I,J,2,2) = ( S(J,K,I,1,1)*DTXXD2 + S(J,K,I,2,1)*DTXYD2
     1                    + S(J,K,I,3,1)*DTXZD2 ) * AREAMR
            AP(I,J,2,3) = ( S(J,K,I,1,1)*DTXXD3 + S(J,K,I,2,1)*DTXYD3 )
     1                                            * AREAMR
            AP(I,J,2,4) = ( S(J,K,I,1,1)*DTXXD4 + S(J,K,I,3,1)*DTXZD4 ) 
     1                                            * AREAMR
            AP(I,J,2,5) = 0.E0
C
            AP(I,J,3,1) = ( S(J,K,I,1,1)*DTXYD1 + S(J,K,I,2,1)*DTYYD1
     1                    + S(J,K,I,3,1)*DTYZD1 ) * AREAMR
            AP(I,J,3,2) = ( S(J,K,I,1,1)*DTXYD2 + S(J,K,I,2,1)*DTYYD2 )
     1                                            * AREAMR
            AP(I,J,3,3) = ( S(J,K,I,1,1)*DTXYD3 + S(J,K,I,2,1)*DTYYD3
     1                    + S(J,K,I,3,1)*DTYZD3 ) * AREAMR
            AP(I,J,3,4) = ( S(J,K,I,2,1)*DTYYD4 + S(J,K,I,3,1)*DTYZD4 ) 
     1                                            * AREAMR
            AP(I,J,3,5) = 0.E0
C
            AP(I,J,4,1) = ( S(J,K,I,1,1)*DTXZD1 + S(J,K,I,2,1)*DTYZD1
     1                    + S(J,K,I,3,1)*DTZZD1 ) * AREAMR
            AP(I,J,4,2) = ( S(J,K,I,1,1)*DTXZD2 + S(J,K,I,3,1)*DTZZD2 ) 
     1                                            * AREAMR
            AP(I,J,4,3) = ( S(J,K,I,2,1)*DTYZD3 + S(J,K,I,3,1)*DTZZD3 ) 
     1                                            * AREAMR
            AP(I,J,4,4) = ( S(J,K,I,1,1)*DTXZD4 + S(J,K,I,2,1)*DTYZD4
     1                    + S(J,K,I,3,1)*DTZZD4 ) * AREAMR
            AP(I,J,4,5) = 0.E0
C
            AP(I,J,5,1) = AREAMR *
     1 ( S(J,K,I,1,1)*(U(I,J)*DTXXD1+V(I,J)*DTXYD1+W(I,J)*DTXZD1-DQDXD1)
     2 + S(J,K,I,2,1)*(U(I,J)*DTXYD1+V(I,J)*DTYYD1+W(I,J)*DTYZD1-DQDYD1)
     3 + S(J,K,I,3,1)*(U(I,J)*DTXZD1+V(I,J)*DTYZD1+W(I,J)*DTZZD1-DQDZD1)
     4 - 0.5E0*Q(J,K,I+1,2)/Q(J,K,I+1,1)*(S(J,K,I,1,1)*TAUXX
     5                    +S(J,K,I,2,1)*TAUXY+S(J,K,I,3,1)*TAUXZ)
     6 - 0.5E0*Q(J,K,I+1,3)/Q(J,K,I+1,1)*(S(J,K,I,1,1)*TAUXY
     7                    +S(J,K,I,2,1)*TAUYY+S(J,K,I,3,1)*TAUYZ)
     8 - 0.5E0*Q(J,K,I+1,4)/Q(J,K,I+1,1)*(S(J,K,I,1,1)*TAUXZ
     9                    +S(J,K,I,2,1)*TAUYZ+S(J,K,I,3,1)*TAUZZ) )
            AP(I,J,5,2) = AREAMR *
     1 ( S(J,K,I,1,1)*(U(I,J)*DTXXD2+V(I,J)*DTXYD2+W(I,J)*DTXZD2-DQDXD2)
     2 + S(J,K,I,2,1)*(U(I,J)*DTXYD2+V(I,J)*DTYYD2              -DQDYD2)
     3 + S(J,K,I,3,1)*(U(I,J)*DTXZD2              +W(I,J)*DTZZD2-DQDZD2)
     4 + 0.5E0/Q(J,K,I+1,1)*(S(J,K,I,1,1)*TAUXX
     5                      +S(J,K,I,2,1)*TAUXY+S(J,K,I,3,1)*TAUXZ) )
            AP(I,J,5,3) = AREAMR *
     1 ( S(J,K,I,1,1)*(U(I,J)*DTXXD3+V(I,J)*DTXYD3              -DQDXD3)
     2 + S(J,K,I,2,1)*(U(I,J)*DTXYD3+V(I,J)*DTYYD3+W(I,J)*DTYZD3-DQDYD3)
     3 + S(J,K,I,3,1)*(              V(I,J)*DTYZD3+W(I,J)*DTZZD3-DQDZD3)
     4 + 0.5E0/Q(J,K,I+1,1)*(S(J,K,I,1,1)*TAUXY
     5                      +S(J,K,I,2,1)*TAUYY+S(J,K,I,3,1)*TAUYZ) )
            AP(I,J,5,4) = AREAMR *
     1 ( S(J,K,I,1,1)*(U(I,J)*DTXXD4              +W(I,J)*DTXZD4-DQDXD4)
     2 + S(J,K,I,2,1)*(              V(I,J)*DTYYD4+W(I,J)*DTYZD4-DQDYD4)
     3 + S(J,K,I,3,1)*(U(I,J)*DTXZD4+V(I,J)*DTYZD4+W(I,J)*DTZZD4-DQDZD4)
     4 + 0.5E0/Q(J,K,I+1,1)*(S(J,K,I,1,1)*TAUXZ
     5                      +S(J,K,I,2,1)*TAUYZ+S(J,K,I,3,1)*TAUZZ) )
            AP(I,J,5,5) = - AREAMR *
     1                      ( S(J,K,I,1,1)*DQDXD5 + S(J,K,I,2,1)*DQDYD5 
     2                      + S(J,K,I,3,1)*DQDZD5 )
  510    CONTINUE
  520 CONTINUE
C
C     Finished with implicit Thin Layer Viscous Terms in the I direction
C
      RETURN
      END






      SUBROUTINE DFVTLJ (IDIM, JDIM, KDIM, I, ITURB,
     1                   Q, RMU, RMUT, S, VOL, IBCDAT, NBCS,
     2                   T, U, V, W, RMUJ, RMUTJ, AP, AM)
C
C     Routine to calculate the Jacobians of the thin layer Navier-Stokes 
C     terms for the J direction.
C
C     This routine is not presently suited for the multiblock
C     approach.  It is close to being ready but is not.  The problem
C     lies that the metrics for one cell into the adjacent block
C     are needed to be completely general.  These are currently not available.
C     The routine metbc needs to be completed.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : Current I Plane
C     ITURB          : Turbulence model
C     Q              : Primitive variables at cell centers
C     RMU            : Molecular viscosity at cell centers
C     RMUT           : Turbulent eddy viscosity at cell centers
C     S              : Metrics
C     VOL            : Cell volumes
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     T              : Temporary storage of Temperature
C     U,V,W          : Temporary storage of U,V,W velocities
C     RMUJ           : Temporary storage of Viscosity at cell interfaces
C     RMUTJ          : Temporary storage of turbulent eddy viscosity 
C                      at cell interfaces
C     AP,AM          : Implicit terms
C                      AP = dF(i+1/2)/dQ(i+1)
C                      AM = dF(i+1/2)/dQ(i)
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          RMU   (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     2          RMUT  (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     3          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     4          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION T     (JDIM+1,2:KDIM),  
     1          RMUJ  (JDIM,2:KDIM), 
     2          RMUTJ (JDIM,2:KDIM),
     3          U     (JDIM,2:KDIM), 
     4          V     (JDIM,2:KDIM), 
     5          W     (JDIM,2:KDIM)
C
      DIMENSION AP    (JDIM,2:KDIM,NRANK,NRANK),
     1          AM    (JDIM,2:KDIM,NRANK,NRANK)
C
      DIMENSION IBCDAT(10,NBCS)
C
C     Two-Thirds
C     Four-Thirds
C
      T23 = 2.E0 / 3.E0
      T43 = 4.E0 / 3.E0
C
C     Calculate T = T / TINF  at the cell centers
C
      DO 120 K = 2, KDIM
         DO 110 J = 1, JDIM + 1
            T(J,K) = GAMMA * Q(J,K,I,5) / Q(J,K,I,1)
  110    CONTINUE
  120 CONTINUE
C
C     Calculate U,V,W,RMUJ,RMUTJ at cell interfaces as averages of cell centers
C
      DO 220 K = 2, KDIM
         DO 210 J = 1, JDIM
            U(J,K)     = 0.5E0 * ( Q(J,K,I,2)  + Q(J+1,K,I,2) )
            V(J,K)     = 0.5E0 * ( Q(J,K,I,3)  + Q(J+1,K,I,3) )
            W(J,K)     = 0.5E0 * ( Q(J,K,I,4)  + Q(J+1,K,I,4) )
            RMUJ(J,K)  = 0.5E0 * ( RMU(J,K,I)  + RMU(J+1,K,I) )
            RMUTJ(J,K) = 0.5E0 * ( RMUT(J,K,I) + RMUT(J+1,K,I) )
  210    CONTINUE
  220 CONTINUE
C
C     Loop over boundary conditions and fix the wall data
C
      DO 350 IBC = 1, NBCS
         IF (IBCDAT(1,IBC) .EQ. IWALL .AND. 
     1       IBCDAT(2,IBC) .EQ. 2) THEN
            INROUT = IBCDAT(3,IBC)
            ISTRT  = IBCDAT(4,IBC)
            IEND   = IBCDAT(5,IBC)
            KSTRT  = IBCDAT(8,IBC)
            KEND   = IBCDAT(9,IBC)
            IF (INROUT .EQ. 1) THEN
               JBCPT = 1
               JNTR  = 1
            ELSE
               JBCPT = JDIM + 1
               JNTR  = JDIM
            ENDIF
            IF (I .GE. ISTRT .AND. I .LE. IEND) THEN
               DO 330 K = KSTRT, KEND
                  U(JNTR,K)     = 0.E0
                  V(JNTR,K)     = 0.E0
                  W(JNTR,K)     = 0.E0
                  RMUJ(JNTR,K)  = RMU(JBCPT,K,I)
                  RMUTJ(JNTR,K) = RMUT(JBCPT,K,I)
  330          CONTINUE
            ENDIF
         ENDIF
  350 CONTINUE
C
C     Loop over K forming viscous terms along the J direction
C     on the J = 1, JDIM faces 
C
      DO 520 K = 2, KDIM
         DO 510 J = 1, JDIM
C
            AREAMR = S(J,K,I,4,2) * FSMACH / RE
C
C     Average metrics and cell volume
C
C     VOLAUX is the inverse of the average auxilliary cell volume
C
            VOLAUX = 2.0E0 / ( VOL(J,K,I) + VOL(J+1,K,I) )
            SPJ1   = 0.5E0 * ( S(J+1,K,I,1,2)*S(J+1,K,I,4,2)
     1                       + S(J,K,I,1,2)*S(J,K,I,4,2) )
            SPJ2   = 0.5E0 * ( S(J+1,K,I,2,2)*S(J+1,K,I,4,2)
     1                       + S(J,K,I,2,2)*S(J,K,I,4,2) )
            SPJ3   = 0.5E0 * ( S(J+1,K,I,3,2)*S(J+1,K,I,4,2)
     1                       + S(J,K,I,3,2)*S(J,K,I,4,2) )
            SMJ1   = 0.5E0 * ( S(J-1,K,I,1,2)*S(J-1,K,I,4,2)
     1                       + S(J,K,I,1,2)*S(J,K,I,4,2) )
            SMJ2   = 0.5E0 * ( S(J-1,K,I,2,2)*S(J-1,K,I,4,2)
     1                       + S(J,K,I,2,2)*S(J,K,I,4,2) )
            SMJ3   = 0.5E0 * ( S(J-1,K,I,3,2)*S(J-1,K,I,4,2)
     1                       + S(J,K,I,3,2)*S(J,K,I,4,2) )
C
C     Calculate first derivative terms using Divergence Theorem
C
            DUDX = VOLAUX * (Q(J+1,K,I,2)*SPJ1 - Q(J,K,I,2)*SMJ1)
            DUDY = VOLAUX * (Q(J+1,K,I,2)*SPJ2 - Q(J,K,I,2)*SMJ2)
            DUDZ = VOLAUX * (Q(J+1,K,I,2)*SPJ3 - Q(J,K,I,2)*SMJ3)
C
            DVDX = VOLAUX * (Q(J+1,K,I,3)*SPJ1 - Q(J,K,I,3)*SMJ1)
            DVDY = VOLAUX * (Q(J+1,K,I,3)*SPJ2 - Q(J,K,I,3)*SMJ2)
            DVDZ = VOLAUX * (Q(J+1,K,I,3)*SPJ3 - Q(J,K,I,3)*SMJ3)
C
            DWDX = VOLAUX * (Q(J+1,K,I,4)*SPJ1 - Q(J,K,I,4)*SMJ1)
            DWDY = VOLAUX * (Q(J+1,K,I,4)*SPJ2 - Q(J,K,I,4)*SMJ2)
            DWDZ = VOLAUX * (Q(J+1,K,I,4)*SPJ3 - Q(J,K,I,4)*SMJ3)
C
            DTDX = VOLAUX * (T(J+1,K)*SPJ1 - T(J,K)*SMJ1)
            DTDY = VOLAUX * (T(J+1,K)*SPJ2 - T(J,K)*SMJ2)
            DTDZ = VOLAUX * (T(J+1,K)*SPJ3 - T(J,K)*SMJ3)
C
C     Calculate Shear Terms
C     Use molecular viscosity (RMUJ) and turbulent eddy viscosity (RMUTJ)
C     averaged from cell centers to get viscosity at cell interface.
C     T23 is two-thirds.
C
            VISC  = RMUJ(J,K) + RMUTJ(J,K)
            TAUXX = T23 * VISC * ( 2.E0*DUDX - DVDY - DWDZ )
            TAUYY = T23 * VISC * ( 2.E0*DVDY - DUDX - DWDZ )
            TAUZZ = T23 * VISC * ( 2.E0*DWDZ - DUDX - DVDY )
            TAUXY = VISC * ( DUDY + DVDX )
            TAUXZ = VISC * ( DWDX + DUDZ )
            TAUYZ = VISC * ( DWDY + DVDZ )
C
C     Heat Flux terms are not used in the current implementation of the
C     implicit viscous terms.  This is kept for compatibility reasons
C     and in case they are later included in d(RMU) terms.
C     Calculate Heat Flux Terms
C
C           CNDCT = (RMUJ(J,K)/PR + RMUTJ(J,K)/PRT) / GAMM1
C           QDOTX = - CNDCT * DTDX
C           QDOTY = - CNDCT * DTDY
C           QDOTZ = - CNDCT * DTDZ
C
C     Region to calculate AM = dF(i+1/2)/dQ(i)  terms
C
C     Calculate DTau(xx,xy,xz,...)/DQ(1,2,3,...) terms
C
            RMRV    = VISC * VOLAUX / Q(J,K,I,1)
C
C     DTau(xx)
C
            DTXXD1 =   T23 * RMRV * ( 2.E0*SMJ1*Q(J,K,I,2) 
     1                       - SMJ2*Q(J,K,I,3) - SMJ3*Q(J,K,I,4) )
            DTXXD2 = - T43 * RMRV * SMJ1
            DTXXD3 =   T23 * RMRV * SMJ2
            DTXXD4 =   T23 * RMRV * SMJ3
C
C     DTau(yy)
C
            DTYYD1 =   T23 * RMRV * ( 2.E0*SMJ2*Q(J,K,I,3) 
     1                       - SMJ1*Q(J,K,I,2) - SMJ3*Q(J,K,I,4) )
            DTYYD2 =   T23 * RMRV * SMJ1
            DTYYD3 = - T43 * RMRV * SMJ2
            DTYYD4 =   T23 * RMRV * SMJ3
C
C     DTau(zz)
C
            DTZZD1 =   T23 * RMRV * ( 2.E0*SMJ3*Q(J,K,I,4) 
     1                       - SMJ1*Q(J,K,I,2) - SMJ2*Q(J,K,I,3) )
            DTZZD2 =   T23 * RMRV * SMJ1
            DTZZD3 =   T23 * RMRV * SMJ2
            DTZZD4 = - T43 * RMRV * SMJ3
C
C     DTau(xy)
C
            DTXYD1 =   RMRV * ( SMJ1*Q(J,K,I,3) + SMJ2*Q(J,K,I,2) )
            DTXYD2 = - RMRV * SMJ2
            DTXYD3 = - RMRV * SMJ1
C
C     DTau(xz)
C
            DTXZD1 =   RMRV * ( SMJ1*Q(J,K,I,4) + SMJ3*Q(J,K,I,2) )
            DTXZD2 = - RMRV * SMJ3
            DTXZD4 = - RMRV * SMJ1
C
C     DTau(yz)
C
            DTYZD1 =   RMRV * ( SMJ2*Q(J,K,I,4) + SMJ3*Q(J,K,I,3) )
            DTYZD3 = - RMRV * SMJ3
            DTYZD4 = - RMRV * SMJ2
C
C     DQdot(x)
C
            TEMP   =   (RMUJ(J,K)/PR + RMUTJ(J,K)/PRT) * GAMMA * VOLAUX
     1                  / Q(J,K,I,1)
C
            DQDXD1 = - TEMP * SMJ1 * ( Q(J,K,I,5)/GAMM1 
     1                   - 0.5E0*Q(J,K,I,1)*( Q(J,K,I,2)*Q(J,K,I,2)
     2                 +Q(J,K,I,3)*Q(J,K,I,3)+Q(J,K,I,4)*Q(J,K,I,4) ) )
     3                    / Q(J,K,I,1)
            DQDXD2 = - TEMP * SMJ1 * Q(J,K,I,2) 
            DQDXD3 = - TEMP * SMJ1 * Q(J,K,I,3) 
            DQDXD4 = - TEMP * SMJ1 * Q(J,K,I,4) 
            DQDXD5 =   TEMP * SMJ1 
C
C     DQdot(y)
C
            DQDYD1 = - TEMP * SMJ2 * ( Q(J,K,I,5)/GAMM1
     1                   - 0.5E0*Q(J,K,I,1)*( Q(J,K,I,2)*Q(J,K,I,2)
     2                 +Q(J,K,I,3)*Q(J,K,I,3)+Q(J,K,I,4)*Q(J,K,I,4) ) )
     3                    / Q(J,K,I,1)
            DQDYD2 = - TEMP * SMJ2 * Q(J,K,I,2)
            DQDYD3 = - TEMP * SMJ2 * Q(J,K,I,3) 
            DQDYD4 = - TEMP * SMJ2 * Q(J,K,I,4) 
            DQDYD5 =   TEMP * SMJ2 
C
C     DQdot(z)
C
            DQDZD1 = - TEMP * SMJ3 * ( Q(J,K,I,5)/GAMM1
     1                   - 0.5E0*Q(J,K,I,1)*( Q(J,K,I,2)*Q(J,K,I,2)
     2                 +Q(J,K,I,3)*Q(J,K,I,3)+Q(J,K,I,4)*Q(J,K,I,4) ) )
     3                    / Q(J,K,I,1)
            DQDZD2 = - TEMP * SMJ3 * Q(J,K,I,2)
            DQDZD3 = - TEMP * SMJ3 * Q(J,K,I,3)
            DQDZD4 = - TEMP * SMJ3 * Q(J,K,I,4)
            DQDZD5 =   TEMP * SMJ3
C
C     Calculate AM = dF(i+1/2)/dQ(i)
C
            AM(J,K,1,1) = 0.E0
            AM(J,K,1,2) = 0.E0
            AM(J,K,1,3) = 0.E0
            AM(J,K,1,4) = 0.E0
            AM(J,K,1,5) = 0.E0
C
            AM(J,K,2,1) = ( S(J,K,I,1,2)*DTXXD1 + S(J,K,I,2,2)*DTXYD1
     1                    + S(J,K,I,3,2)*DTXZD1 ) * AREAMR
            AM(J,K,2,2) = ( S(J,K,I,1,2)*DTXXD2 + S(J,K,I,2,2)*DTXYD2
     1                    + S(J,K,I,3,2)*DTXZD2 ) * AREAMR
            AM(J,K,2,3) = ( S(J,K,I,1,2)*DTXXD3 + S(J,K,I,2,2)*DTXYD3 )
     1                                            * AREAMR
            AM(J,K,2,4) = ( S(J,K,I,1,2)*DTXXD4 + S(J,K,I,3,2)*DTXZD4 ) 
     1                                            * AREAMR
            AM(J,K,2,5) = 0.E0
C
            AM(J,K,3,1) = ( S(J,K,I,1,2)*DTXYD1 + S(J,K,I,2,2)*DTYYD1
     1                    + S(J,K,I,3,2)*DTYZD1 ) * AREAMR
            AM(J,K,3,2) = ( S(J,K,I,1,2)*DTXYD2 + S(J,K,I,2,2)*DTYYD2 )
     1                                            * AREAMR
            AM(J,K,3,3) = ( S(J,K,I,1,2)*DTXYD3 + S(J,K,I,2,2)*DTYYD3
     1                    + S(J,K,I,3,2)*DTYZD3 ) * AREAMR
            AM(J,K,3,4) = ( S(J,K,I,2,2)*DTYYD4 + S(J,K,I,3,2)*DTYZD4 ) 
     1                                            * AREAMR
            AM(J,K,3,5) = 0.E0
C
            AM(J,K,4,1) = ( S(J,K,I,1,2)*DTXZD1 + S(J,K,I,2,2)*DTYZD1
     1                    + S(J,K,I,3,2)*DTZZD1 ) * AREAMR
            AM(J,K,4,2) = ( S(J,K,I,1,2)*DTXZD2 + S(J,K,I,3,2)*DTZZD2 ) 
     1                                            * AREAMR
            AM(J,K,4,3) = ( S(J,K,I,2,2)*DTYZD3 + S(J,K,I,3,2)*DTZZD3 ) 
     1                                            * AREAMR
            AM(J,K,4,4) = ( S(J,K,I,1,2)*DTXZD4 + S(J,K,I,2,2)*DTYZD4
     1                    + S(J,K,I,3,2)*DTZZD4 ) * AREAMR
            AM(J,K,4,5) = 0.E0
C
            AM(J,K,5,1) = AREAMR *
     1 ( S(J,K,I,1,2)*(U(J,K)*DTXXD1+V(J,K)*DTXYD1+W(J,K)*DTXZD1-DQDXD1)
     2 + S(J,K,I,2,2)*(U(J,K)*DTXYD1+V(J,K)*DTYYD1+W(J,K)*DTYZD1-DQDYD1)
     3 + S(J,K,I,3,2)*(U(J,K)*DTXZD1+V(J,K)*DTYZD1+W(J,K)*DTZZD1-DQDZD1)
     4 - 0.5E0*Q(J,K,I,2)/Q(J,K,I,1)*(S(J,K,I,1,2)*TAUXX
     5                    +S(J,K,I,2,2)*TAUXY+S(J,K,I,3,2)*TAUXZ)
     6 - 0.5E0*Q(J,K,I,3)/Q(J,K,I,1)*(S(J,K,I,1,2)*TAUXY
     7                    +S(J,K,I,2,2)*TAUYY+S(J,K,I,3,2)*TAUYZ)
     8 - 0.5E0*Q(J,K,I,4)/Q(J,K,I,1)*(S(J,K,I,1,2)*TAUXZ
     9                    +S(J,K,I,2,2)*TAUYZ+S(J,K,I,3,2)*TAUZZ) )
            AM(J,K,5,2) = AREAMR *
     1 ( S(J,K,I,1,2)*(U(J,K)*DTXXD2+V(J,K)*DTXYD2+W(J,K)*DTXZD2-DQDXD2)
     2 + S(J,K,I,2,2)*(U(J,K)*DTXYD2+V(J,K)*DTYYD2              -DQDYD2)
     3 + S(J,K,I,3,2)*(U(J,K)*DTXZD2              +W(J,K)*DTZZD2-DQDZD2)
     4 + 0.5E0/Q(J,K,I,1)*(S(J,K,I,1,2)*TAUXX
     5                    +S(J,K,I,2,2)*TAUXY+S(J,K,I,3,2)*TAUXZ) )
            AM(J,K,5,3) = AREAMR *
     1 ( S(J,K,I,1,2)*(U(J,K)*DTXXD3+V(J,K)*DTXYD3              -DQDXD3)
     2 + S(J,K,I,2,2)*(U(J,K)*DTXYD3+V(J,K)*DTYYD3+W(J,K)*DTYZD3-DQDYD3)
     3 + S(J,K,I,3,2)*(              V(J,K)*DTYZD3+W(J,K)*DTZZD3-DQDZD3)
     4 + 0.5E0/Q(J,K,I,1)*(S(J,K,I,1,2)*TAUXY
     5                    +S(J,K,I,2,2)*TAUYY+S(J,K,I,3,2)*TAUYZ) )
            AM(J,K,5,4) = AREAMR *
     1 ( S(J,K,I,1,2)*(U(J,K)*DTXXD4              +W(J,K)*DTXZD4-DQDXD4)
     2 + S(J,K,I,2,2)*(              V(J,K)*DTYYD4+W(J,K)*DTYZD4-DQDYD4)
     3 + S(J,K,I,3,2)*(U(J,K)*DTXZD4+V(J,K)*DTYZD4+W(J,K)*DTZZD4-DQDZD4)
     4 + 0.5E0/Q(J,K,I,1)*(S(J,K,I,1,2)*TAUXZ
     5                    +S(J,K,I,2,2)*TAUYZ+S(J,K,I,3,2)*TAUZZ) )
            AM(J,K,5,5) = - AREAMR *
     1                      ( S(J,K,I,1,2)*DQDXD5 + S(J,K,I,2,2)*DQDYD5 
     2                      + S(J,K,I,3,2)*DQDZD5 )
  510    CONTINUE
         DO 515 J = 1, JDIM
C
            AREAMR = S(J,K,I,4,2) * FSMACH / RE
C
C     Average metrics and cell volume
C
C     VOLAUX is the inverse of the average auxilliary cell volume
C
            VOLAUX = 2.0E0 / ( VOL(J,K,I) + VOL(J+1,K,I) )
            SPJ1   = 0.5E0 * ( S(J+1,K,I,1,2)*S(J+1,K,I,4,2)
     1                       + S(J,K,I,1,2)*S(J,K,I,4,2) )
            SPJ2   = 0.5E0 * ( S(J+1,K,I,2,2)*S(J+1,K,I,4,2)
     1                       + S(J,K,I,2,2)*S(J,K,I,4,2) )
            SPJ3   = 0.5E0 * ( S(J+1,K,I,3,2)*S(J+1,K,I,4,2)
     1                       + S(J,K,I,3,2)*S(J,K,I,4,2) )
            SMJ1   = 0.5E0 * ( S(J-1,K,I,1,2)*S(J-1,K,I,4,2)
     1                       + S(J,K,I,1,2)*S(J,K,I,4,2) )
            SMJ2   = 0.5E0 * ( S(J-1,K,I,2,2)*S(J-1,K,I,4,2)
     1                       + S(J,K,I,2,2)*S(J,K,I,4,2) )
            SMJ3   = 0.5E0 * ( S(J-1,K,I,3,2)*S(J-1,K,I,4,2)
     1                       + S(J,K,I,3,2)*S(J,K,I,4,2) )
C
C     Calculate first derivative terms using Divergence Theorem
C
            DUDX = VOLAUX * (Q(J+1,K,I,2)*SPJ1 - Q(J,K,I,2)*SMJ1)
            DUDY = VOLAUX * (Q(J+1,K,I,2)*SPJ2 - Q(J,K,I,2)*SMJ2)
            DUDZ = VOLAUX * (Q(J+1,K,I,2)*SPJ3 - Q(J,K,I,2)*SMJ3)
C
            DVDX = VOLAUX * (Q(J+1,K,I,3)*SPJ1 - Q(J,K,I,3)*SMJ1)
            DVDY = VOLAUX * (Q(J+1,K,I,3)*SPJ2 - Q(J,K,I,3)*SMJ2)
            DVDZ = VOLAUX * (Q(J+1,K,I,3)*SPJ3 - Q(J,K,I,3)*SMJ3)
C
            DWDX = VOLAUX * (Q(J+1,K,I,4)*SPJ1 - Q(J,K,I,4)*SMJ1)
            DWDY = VOLAUX * (Q(J+1,K,I,4)*SPJ2 - Q(J,K,I,4)*SMJ2)
            DWDZ = VOLAUX * (Q(J+1,K,I,4)*SPJ3 - Q(J,K,I,4)*SMJ3)
C
            DTDX = VOLAUX * (T(J+1,K)*SPJ1 - T(J,K)*SMJ1)
            DTDY = VOLAUX * (T(J+1,K)*SPJ2 - T(J,K)*SMJ2)
            DTDZ = VOLAUX * (T(J+1,K)*SPJ3 - T(J,K)*SMJ3)
C
C     Calculate Shear Terms
C     Use molecular viscosity (RMUJ) and turbulent eddy viscosity (RMUTJ)
C     averaged from cell centers to get viscosity at cell interface.
C     T23 is two-thirds.
C
            VISC  = RMUJ(J,K) + RMUTJ(J,K)
            TAUXX = T23 * VISC * ( 2.E0*DUDX - DVDY - DWDZ )
            TAUYY = T23 * VISC * ( 2.E0*DVDY - DUDX - DWDZ )
            TAUZZ = T23 * VISC * ( 2.E0*DWDZ - DUDX - DVDY )
            TAUXY = VISC * ( DUDY + DVDX )
            TAUXZ = VISC * ( DWDX + DUDZ )
            TAUYZ = VISC * ( DWDY + DVDZ )
C
C     Heat Flux terms are not used in the current implementation of the
C     implicit viscous terms.  This is kept for compatibility reasons
C     and in case they are later included in d(RMU) terms.
C     Calculate Heat Flux Terms
C
C           CNDCT = (RMUJ(J,K)/PR + RMUTJ(J,K)/PRT) / GAMM1
C           QDOTX = - CNDCT * DTDX
C           QDOTY = - CNDCT * DTDY
C           QDOTZ = - CNDCT * DTDZ
C
C     Region to calculate AP = dF(i+1/2)/dQ(i)  terms
C
C     Calculate DTau(xx,xy,xz,...)/DQ(1,2,3,...) terms
C
            RMRV    = VISC * VOLAUX / Q(J+1,K,I,1)
C
C     DTau(xx)
C
            DTXXD1 = - T23 * RMRV * ( 2.E0*SPJ1*Q(J+1,K,I,2) 
     1                     - SPJ2*Q(J+1,K,I,3) - SPJ3*Q(J+1,K,I,4) )
            DTXXD2 =   T43 * RMRV * SPJ1
            DTXXD3 = - T23 * RMRV * SPJ2
            DTXXD4 = - T23 * RMRV * SPJ3
C
C     DTau(yy)
C
            DTYYD1 = - T23 * RMRV * ( 2.E0*SPJ2*Q(J+1,K,I,3) 
     1                     - SPJ1*Q(J+1,K,I,2) - SPJ3*Q(J+1,K,I,4) )
            DTYYD2 = - T23 * RMRV * SPJ1
            DTYYD3 =   T43 * RMRV * SPJ2
            DTYYD4 = - T23 * RMRV * SPJ3
C
C     DTau(zz)
C
            DTZZD1 = - T23 * RMRV * ( 2.E0*SPJ3*Q(J+1,K,I,4) 
     1                     - SPJ1*Q(J+1,K,I,2) - SPJ2*Q(J+1,K,I,3) )
            DTZZD2 = - T23 * RMRV * SPJ1
            DTZZD3 = - T23 * RMRV * SPJ2
            DTZZD4 =   T43 * RMRV * SPJ3
C
C     DTau(xy)
C
            DTXYD1 = - RMRV * ( SPJ1*Q(J+1,K,I,3) + SPJ2*Q(J+1,K,I,2) )
            DTXYD2 =   RMRV * SPJ2
            DTXYD3 =   RMRV * SPJ1
C
C     DTau(xz)
C
            DTXZD1 = - RMRV * ( SPJ1*Q(J+1,K,I,4) + SPJ3*Q(J+1,K,I,2) )
            DTXZD2 =   RMRV * SPJ3
            DTXZD4 =   RMRV * SPJ1
C
C     DTau(yz)
C
            DTYZD1 = - RMRV * ( SPJ2*Q(J+1,K,I,4) + SPJ3*Q(J+1,K,I,3) )
            DTYZD3 =   RMRV * SPJ3
            DTYZD4 =   RMRV * SPJ2
C
C     DQdot(x)
C
            TEMP   =   (RMUJ(J,K)/PR + RMUTJ(J,K)/PRT) * GAMMA * VOLAUX
     1                  / Q(J+1,K,I,1)
C
            DQDXD1 =   TEMP * SPJ1 * ( Q(J+1,K,I,5)/GAMM1 
     1             - 0.5E0*Q(J+1,K,I,1)*( Q(J+1,K,I,2)*Q(J+1,K,I,2)
     2         +Q(J+1,K,I,3)*Q(J+1,K,I,3)+Q(J+1,K,I,4)*Q(J+1,K,I,4) ) )
     3                    / Q(J+1,K,I,1)
            DQDXD2 =   TEMP * SPJ1 * Q(J+1,K,I,2) 
            DQDXD3 =   TEMP * SPJ1 * Q(J+1,K,I,3) 
            DQDXD4 =   TEMP * SPJ1 * Q(J+1,K,I,4) 
            DQDXD5 = - TEMP * SPJ1 
C
C     DQdot(y)
C
            DQDYD1 =   TEMP * SPJ2 * ( Q(J+1,K,I,5)/GAMM1
     1             - 0.5E0*Q(J+1,K,I,1)*( Q(J+1,K,I,2)*Q(J+1,K,I,2)
     2         +Q(J+1,K,I,3)*Q(J+1,K,I,3)+Q(J+1,K,I,4)*Q(J+1,K,I,4) ) )
     3                    / Q(J+1,K,I,1)
            DQDYD2 =   TEMP * SPJ2 * Q(J+1,K,I,2)
            DQDYD3 =   TEMP * SPJ2 * Q(J+1,K,I,3) 
            DQDYD4 =   TEMP * SPJ2 * Q(J+1,K,I,4) 
            DQDYD5 = - TEMP * SPJ2 
C
C     DQdot(z)
C
            DQDZD1 =   TEMP * SPJ3 * ( Q(J+1,K,I,5)/GAMM1
     1             - 0.5E0*Q(J+1,K,I,1)*( Q(J+1,K,I,2)*Q(J+1,K,I,2)
     2         +Q(J+1,K,I,3)*Q(J+1,K,I,3)+Q(J+1,K,I,4)*Q(J+1,K,I,4) ) )
     3                    / Q(J+1,K,I,1)
            DQDZD2 =   TEMP * SPJ3 * Q(J+1,K,I,2)
            DQDZD3 =   TEMP * SPJ3 * Q(J+1,K,I,3)
            DQDZD4 =   TEMP * SPJ3 * Q(J+1,K,I,4)
            DQDZD5 = - TEMP * SPJ3
C
C     Calculate AP = dF(i+1/2)/dQ(i+1)
C
            AP(J,K,1,1) = 0.E0
            AP(J,K,1,2) = 0.E0
            AP(J,K,1,3) = 0.E0
            AP(J,K,1,4) = 0.E0
            AP(J,K,1,5) = 0.E0
C
            AP(J,K,2,1) = ( S(J,K,I,1,2)*DTXXD1 + S(J,K,I,2,2)*DTXYD1
     1                    + S(J,K,I,3,2)*DTXZD1 ) * AREAMR
            AP(J,K,2,2) = ( S(J,K,I,1,2)*DTXXD2 + S(J,K,I,2,2)*DTXYD2
     1                    + S(J,K,I,3,2)*DTXZD2 ) * AREAMR
            AP(J,K,2,3) = ( S(J,K,I,1,2)*DTXXD3 + S(J,K,I,2,2)*DTXYD3 )
     1                                            * AREAMR
            AP(J,K,2,4) = ( S(J,K,I,1,2)*DTXXD4 + S(J,K,I,3,2)*DTXZD4 ) 
     1                                            * AREAMR
            AP(J,K,2,5) = 0.E0
C
            AP(J,K,3,1) = ( S(J,K,I,1,2)*DTXYD1 + S(J,K,I,2,2)*DTYYD1
     1                    + S(J,K,I,3,2)*DTYZD1 ) * AREAMR
            AP(J,K,3,2) = ( S(J,K,I,1,2)*DTXYD2 + S(J,K,I,2,2)*DTYYD2 )
     1                                            * AREAMR
            AP(J,K,3,3) = ( S(J,K,I,1,2)*DTXYD3 + S(J,K,I,2,2)*DTYYD3
     1                    + S(J,K,I,3,2)*DTYZD3 ) * AREAMR
            AP(J,K,3,4) = ( S(J,K,I,2,2)*DTYYD4 + S(J,K,I,3,2)*DTYZD4 ) 
     1                                            * AREAMR
            AP(J,K,3,5) = 0.E0
C
            AP(J,K,4,1) = ( S(J,K,I,1,2)*DTXZD1 + S(J,K,I,2,2)*DTYZD1
     1                    + S(J,K,I,3,2)*DTZZD1 ) * AREAMR
            AP(J,K,4,2) = ( S(J,K,I,1,2)*DTXZD2 + S(J,K,I,3,2)*DTZZD2 ) 
     1                                            * AREAMR
            AP(J,K,4,3) = ( S(J,K,I,2,2)*DTYZD3 + S(J,K,I,3,2)*DTZZD3 ) 
     1                                            * AREAMR
            AP(J,K,4,4) = ( S(J,K,I,1,2)*DTXZD4 + S(J,K,I,2,2)*DTYZD4
     1                    + S(J,K,I,3,2)*DTZZD4 ) * AREAMR
            AP(J,K,4,5) = 0.E0
C
            AP(J,K,5,1) = AREAMR *
     1 ( S(J,K,I,1,2)*(U(J,K)*DTXXD1+V(J,K)*DTXYD1+W(J,K)*DTXZD1-DQDXD1)
     2 + S(J,K,I,2,2)*(U(J,K)*DTXYD1+V(J,K)*DTYYD1+W(J,K)*DTYZD1-DQDYD1)
     3 + S(J,K,I,3,2)*(U(J,K)*DTXZD1+V(J,K)*DTYZD1+W(J,K)*DTZZD1-DQDZD1)
     4 - 0.5E0*Q(J+1,K,I,2)/Q(J+1,K,I,1)*(S(J,K,I,1,2)*TAUXX
     5                    +S(J,K,I,2,2)*TAUXY+S(J,K,I,3,2)*TAUXZ)
     6 - 0.5E0*Q(J+1,K,I,3)/Q(J+1,K,I,1)*(S(J,K,I,1,2)*TAUXY
     7                    +S(J,K,I,2,2)*TAUYY+S(J,K,I,3,2)*TAUYZ)
     8 - 0.5E0*Q(J+1,K,I,4)/Q(J+1,K,I,1)*(S(J,K,I,1,2)*TAUXZ
     9                    +S(J,K,I,2,2)*TAUYZ+S(J,K,I,3,2)*TAUZZ) )
            AP(J,K,5,2) = AREAMR *
     1 ( S(J,K,I,1,2)*(U(J,K)*DTXXD2+V(J,K)*DTXYD2+W(J,K)*DTXZD2-DQDXD2)
     2 + S(J,K,I,2,2)*(U(J,K)*DTXYD2+V(J,K)*DTYYD2              -DQDYD2)
     3 + S(J,K,I,3,2)*(U(J,K)*DTXZD2              +W(J,K)*DTZZD2-DQDZD2)
     4 + 0.5E0/Q(J+1,K,I,1)*(S(J,K,I,1,2)*TAUXX
     5                      +S(J,K,I,2,2)*TAUXY+S(J,K,I,3,2)*TAUXZ) )
            AP(J,K,5,3) = AREAMR *
     1 ( S(J,K,I,1,2)*(U(J,K)*DTXXD3+V(J,K)*DTXYD3              -DQDXD3)
     2 + S(J,K,I,2,2)*(U(J,K)*DTXYD3+V(J,K)*DTYYD3+W(J,K)*DTYZD3-DQDYD3)
     3 + S(J,K,I,3,2)*(              V(J,K)*DTYZD3+W(J,K)*DTZZD3-DQDZD3)
     4 + 0.5E0/Q(J+1,K,I,1)*(S(J,K,I,1,2)*TAUXY
     5                      +S(J,K,I,2,2)*TAUYY+S(J,K,I,3,2)*TAUYZ) )
            AP(J,K,5,4) = AREAMR *
     1 ( S(J,K,I,1,2)*(U(J,K)*DTXXD4              +W(J,K)*DTXZD4-DQDXD4)
     2 + S(J,K,I,2,2)*(              V(J,K)*DTYYD4+W(J,K)*DTYZD4-DQDYD4)
     3 + S(J,K,I,3,2)*(U(J,K)*DTXZD4+V(J,K)*DTYZD4+W(J,K)*DTZZD4-DQDZD4)
     4 + 0.5E0/Q(J+1,K,I,1)*(S(J,K,I,1,2)*TAUXZ
     5                      +S(J,K,I,2,2)*TAUYZ+S(J,K,I,3,2)*TAUZZ) )
            AP(J,K,5,5) = - AREAMR *
     1                      ( S(J,K,I,1,2)*DQDXD5 + S(J,K,I,2,2)*DQDYD5 
     2                      + S(J,K,I,3,2)*DQDZD5 )
  515    CONTINUE
  520 CONTINUE
C
C     Add additional terms for turbulence model
C
C     2-i.  Two-Equation models
C           N.B. This is an approximation for the Jacobian including only
C                the terms on the diagonal!
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 820 K = 2, KDIM
            DO 810 J = 1, JDIM
C
               VOLAUX = 2.0E0 / ( VOL(J,K,I) + VOL(J+1,K,I) )
               TERM   = S(J,K,I,4,2) * S(J,K,I,4,2) * VOLAUX 
     1                  * FSMACH / RE
C
               AM(J,K,6,6) = - (RMUJ(J,K) + SIGK  * RMUTJ(J,K)) * TERM 
     1                         / Q(J  ,K,I,1)
               AP(J,K,6,6) =   (RMUJ(J,K) + SIGK  * RMUTJ(J,K)) * TERM 
     1                         / Q(J+1,K,I,1)
               AM(J,K,7,7) = - (RMUJ(J,K) + SIGT2 * RMUTJ(J,K)) * TERM 
     1                         / Q(J  ,K,I,1)
               AP(J,K,7,7) =   (RMUJ(J,K) + SIGT2 * RMUTJ(J,K)) * TERM 
     1                         / Q(J+1,K,I,1)
  810       CONTINUE
  820    CONTINUE
C
C     2-ii. Reynolds Stress model
C           N.B. This is an approximation for the Jacobian including only
C                the terms on the diagonal!
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 920 K = 2, KDIM
            DO 910 J = 1, JDIM
C
               VOLAUX = 2.0E0 / ( VOL(J,K,I) + VOL(J+1,K,I) )
               TERM   = S(J,K,I,4,2) * S(J,K,I,4,2) * VOLAUX 
     1                  * FSMACH / RE
C
               AM(J,K, 6, 6) = - (RMUJ(J,K) + SIGK  * RMUTJ(J,K)) * TERM 
     1                           / Q(J  ,K,I,1)
               AP(J,K, 6, 6) =   (RMUJ(J,K) + SIGK  * RMUTJ(J,K)) * TERM 
     1                           / Q(J+1,K,I,1)
               AM(J,K, 7, 7) = - (RMUJ(J,K) + SIGK  * RMUTJ(J,K)) * TERM 
     1                           / Q(J  ,K,I,1)
               AP(J,K, 7, 7) =   (RMUJ(J,K) + SIGK  * RMUTJ(J,K)) * TERM 
     1                           / Q(J+1,K,I,1)
               AM(J,K, 8, 8) = - (RMUJ(J,K) + SIGK  * RMUTJ(J,K)) * TERM 
     1                           / Q(J  ,K,I,1)
               AP(J,K, 8, 8) =   (RMUJ(J,K) + SIGK  * RMUTJ(J,K)) * TERM 
     1                           / Q(J+1,K,I,1)
               AM(J,K, 9, 9) = - (RMUJ(J,K) + SIGK  * RMUTJ(J,K)) * TERM 
     1                           / Q(J  ,K,I,1)
               AP(J,K, 9, 9) =   (RMUJ(J,K) + SIGK  * RMUTJ(J,K)) * TERM 
     1                           / Q(J+1,K,I,1)
               AM(J,K,10,10) = - (RMUJ(J,K) + SIGK  * RMUTJ(J,K)) * TERM 
     1                           / Q(J  ,K,I,1)
               AP(J,K,10,10) =   (RMUJ(J,K) + SIGK  * RMUTJ(J,K)) * TERM 
     1                           / Q(J+1,K,I,1)
               AM(J,K,11,11) = - (RMUJ(J,K) + SIGK  * RMUTJ(J,K)) * TERM 
     1                           / Q(J  ,K,I,1)
               AP(J,K,11,11) =   (RMUJ(J,K) + SIGK  * RMUTJ(J,K)) * TERM 
     1                           / Q(J+1,K,I,1)
C
               AM(J,K,12,12) = - (RMUJ(J,K) + SIGT2 * RMUTJ(J,K)) * TERM 
     1                           / Q(J  ,K,I,1)
               AP(J,K,12,12) =   (RMUJ(J,K) + SIGT2 * RMUTJ(J,K)) * TERM 
     1                           / Q(J+1,K,I,1)
  910       CONTINUE
  920    CONTINUE
      ENDIF
C
C     Finished with implicit Thin Layer Viscous Terms in the J direction
C
      RETURN
      END






      SUBROUTINE DFVTLK (IDIM, JDIM, KDIM, I,
     1                   Q, RMU, RMUT, S, VOL, IBCDAT, NBCS,
     2                   T, U, V, W, RMUK, RMUTK, AP, AM)
C
C     Routine to calculate the thin layer Navier-Stokes terms for
C     the K direction.
C
C     This routine is not presently suited for the multiblock
C     approach.  It is close to being ready but is not.  The problem
C     lies that the metrics for one cell into the adjacent block
C     are needed to be completely general.  These are currently not available.
C     The routine metbc needs to be completed.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : Current I Plane
C     Q              : Primitive variables at cell centers
C     RMU            : Molecular viscosity at cell centers
C     RMUT           : Turbulent eddy viscosity at cell centers
C     S              : Metrics
C     VOL            : Cell volumes
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     T              : Temporary storage of Temperature
C     U,V,W          : Temporary storage of U,V,W velocities
C     RMUK           : Temporary storage of Viscosity at cell interfaces
C     RMUTK          : Temporary storage of turbulent eddy viscosity
C                      at cell interfaces
C     AP,AM          : Implicit terms
C                      AP = dF(i+1/2)/dQ(i+1)
C                      AM = dF(i+1/2)/dQ(i)
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          RMU   (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     2          RMUT  (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     3          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     4          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION T     (KDIM+1,2:JDIM),  
     1          RMUK  (KDIM,2:JDIM), 
     2          RMUTK (KDIM,2:JDIM),
     3          U     (KDIM,2:JDIM), 
     4          V     (KDIM,2:JDIM), 
     5          W     (KDIM,2:JDIM)
C
      DIMENSION AP    (KDIM,2:JDIM,NRANK,NRANK),
     1          AM    (KDIM,2:JDIM,NRANK,NRANK)
C
      DIMENSION IBCDAT(10,NBCS)
C
C     Two-Thirds
C     Four-Thirds
C
      T23 = 2.E0 / 3.E0
      T43 = 4.E0 / 3.E0
C
C     Calculate T = T / TINF  at the cell centers
C
      DO 120 J = 2, JDIM
         DO 110 K = 1, KDIM + 1
            T(K,J) = GAMMA * Q(J,K,I,5) / Q(J,K,I,1)
  110    CONTINUE
  120 CONTINUE
C
C     Calculate U,V,W,RMUK,RMUTK at cell interfaces as averages of cell centers
C
      DO 220 J = 2, JDIM
         DO 210 K = 1, KDIM
            U(K,J)     = 0.5E0 * ( Q(J,K,I,2)  + Q(J,K+1,I,2) )
            V(K,J)     = 0.5E0 * ( Q(J,K,I,3)  + Q(J,K+1,I,3) )
            W(K,J)     = 0.5E0 * ( Q(J,K,I,4)  + Q(J,K+1,I,4) )
            RMUK(K,J)  = 0.5E0 * ( RMU(J,K,I)  + RMU(J,K+1,I) )
            RMUTK(K,J) = 0.5E0 * ( RMUT(J,K,I) + RMUT(J,K+1,I) )
  210    CONTINUE
  220 CONTINUE
C
C     Loop over boundary conditions and fix the wall data
C
      DO 350 IBC = 1, NBCS
         IF (IBCDAT(1,IBC) .EQ. IWALL .AND. 
     1       IBCDAT(2,IBC) .EQ. 3) THEN
            INROUT = IBCDAT(3,IBC)
            ISTRT  = IBCDAT(4,IBC)
            IEND   = IBCDAT(5,IBC)
            JSTRT  = IBCDAT(6,IBC)
            JEND   = IBCDAT(7,IBC)
            IF (INROUT .EQ. 1) THEN
               KBCPT = 1
               KNTR  = 1
            ELSE
               KBCPT = KDIM + 1
               KNTR  = KDIM
            ENDIF
            IF (I .GE. ISTRT .AND. I .LE. IEND) THEN
               DO 330 J = JSTRT, JEND
                  U(KNTR,J)     = 0.E0
                  V(KNTR,J)     = 0.E0
                  W(KNTR,J)     = 0.E0
                  RMUK(KNTR,J)  = RMU(J,KBCPT,I)
                  RMUTK(KNTR,J) = RMUT(J,KBCPT,I)
  330          CONTINUE
            ENDIF
         ENDIF
  350 CONTINUE
C
C     Loop over J forming viscous terms along the K direction
C     on the K = 1, KDIM faces 
C
      DO 520 J = 2, JDIM
         DO 510 K = 1, KDIM
C
            AREAMR = S(J,K,I,4,3) * FSMACH / RE
C
C     Average metrics and cell volume
C
C     VOLAUX is the inverse of the average auxilliary cell volume
C
            VOLAUX = 2.0E0 / ( VOL(J,K,I) + VOL(J,K+1,I) )
            SPK1   = 0.5E0 * ( S(J,K+1,I,1,3)*S(J,K+1,I,4,3)
     1                       + S(J,K,I,1,3)*S(J,K,I,4,3) )
            SPK2   = 0.5E0 * ( S(J,K+1,I,2,3)*S(J,K+1,I,4,3)
     1                       + S(J,K,I,2,3)*S(J,K,I,4,3) )
            SPK3   = 0.5E0 * ( S(J,K+1,I,3,3)*S(J,K+1,I,4,3)
     1                       + S(J,K,I,3,3)*S(J,K,I,4,3) )
            SMK1   = 0.5E0 * ( S(J,K-1,I,1,3)*S(J,K-1,I,4,3)
     1                       + S(J,K,I,1,3)*S(J,K,I,4,3) )
            SMK2   = 0.5E0 * ( S(J,K-1,I,2,3)*S(J,K-1,I,4,3)
     1                       + S(J,K,I,2,3)*S(J,K,I,4,3) )
            SMK3   = 0.5E0 * ( S(J,K-1,I,3,3)*S(J,K-1,I,4,3)
     1                       + S(J,K,I,3,3)*S(J,K,I,4,3) )
C
C     Calculate first derivative terms using Divergence Theorem
C
            DUDX = VOLAUX * (Q(J,K+1,I,2)*SPK1 - Q(J,K,I,2)*SMK1)
            DUDY = VOLAUX * (Q(J,K+1,I,2)*SPK2 - Q(J,K,I,2)*SMK2)
            DUDZ = VOLAUX * (Q(J,K+1,I,2)*SPK3 - Q(J,K,I,2)*SMK3)
C
            DVDX = VOLAUX * (Q(J,K+1,I,3)*SPK1 - Q(J,K,I,3)*SMK1)
            DVDY = VOLAUX * (Q(J,K+1,I,3)*SPK2 - Q(J,K,I,3)*SMK2)
            DVDZ = VOLAUX * (Q(J,K+1,I,3)*SPK3 - Q(J,K,I,3)*SMK3)
C
            DWDX = VOLAUX * (Q(J,K+1,I,4)*SPK1 - Q(J,K,I,4)*SMK1)
            DWDY = VOLAUX * (Q(J,K+1,I,4)*SPK2 - Q(J,K,I,4)*SMK2)
            DWDZ = VOLAUX * (Q(J,K+1,I,4)*SPK3 - Q(J,K,I,4)*SMK3)
C
            DTDX = VOLAUX * (T(K+1,J)*SPK1 - T(K,J)*SMK1)
            DTDY = VOLAUX * (T(K+1,J)*SPK2 - T(K,J)*SMK2)
            DTDZ = VOLAUX * (T(K+1,J)*SPK3 - T(K,J)*SMK3)
C
C     Calculate Shear Terms
C     Use molecular viscosity (RMUK) and turbulent eddy viscosity (RMUTK)
C     averaged from cell centers to get viscosity at cell interface.
C     T23 is two-thirds.
C
            VISC  = RMUK(K,J) + RMUTK(K,J)
            TAUXX = T23 * VISC * ( 2.E0*DUDX - DVDY - DWDZ )
            TAUYY = T23 * VISC * ( 2.E0*DVDY - DUDX - DWDZ )
            TAUZZ = T23 * VISC * ( 2.E0*DWDZ - DUDX - DVDY )
            TAUXY = VISC * ( DUDY + DVDX )
            TAUXZ = VISC * ( DWDX + DUDZ )
            TAUYZ = VISC * ( DWDY + DVDZ )
C
C     Heat Flux terms are not used in the current implementation of the
C     implicit viscous terms.  This is kept for compatibility reasons
C     and in case they are later included in d(RMU) terms.
C     Calculate Heat Flux Terms
C
C           CNDCT = (RMUK(K,J)/PR + RMUTK(K,J)/PRT) / GAMM1
C           QDOTX = - CNDCT * DTDX
C           QDOTY = - CNDCT * DTDY
C           QDOTZ = - CNDCT * DTDZ
C
C     Region to calculate AM = dF(i+1/2)/dQ(i)  terms
C
C     Calculate DTau(xx,xy,xz,...)/DQ(1,2,3,...) terms
C
            RMRV   = VISC * VOLAUX / Q(J,K,I,1)
C
C     DTau(xx)
C
            DTXXD1 =   T23 * RMRV * ( 2.E0*SMK1*Q(J,K,I,2) 
     1                       - SMK2*Q(J,K,I,3) - SMK3*Q(J,K,I,4) )
            DTXXD2 = - T43 * RMRV * SMK1
            DTXXD3 =   T23 * RMRV * SMK2
            DTXXD4 =   T23 * RMRV * SMK3
C
C     DTau(yy)
C
            DTYYD1 =   T23 * RMRV * ( 2.E0*SMK2*Q(J,K,I,3) 
     1                       - SMK1*Q(J,K,I,2) - SMK3*Q(J,K,I,4) )
            DTYYD2 =   T23 * RMRV * SMK1
            DTYYD3 = - T43 * RMRV * SMK2
            DTYYD4 =   T23 * RMRV * SMK3
C
C     DTau(zz)
C
            DTZZD1 =   T23 * RMRV * ( 2.E0*SMK3*Q(J,K,I,4) 
     1                       - SMK1*Q(J,K,I,2) - SMK2*Q(J,K,I,3) )
            DTZZD2 =   T23 * RMRV * SMK1
            DTZZD3 =   T23 * RMRV * SMK2
            DTZZD4 = - T43 * RMRV * SMK3
C
C     DTau(xy)
C
            DTXYD1 =   RMRV * ( SMK1*Q(J,K,I,3) + SMK2*Q(J,K,I,2) )
            DTXYD2 = - RMRV * SMK2
            DTXYD3 = - RMRV * SMK1
C
C     DTau(xz)
C
            DTXZD1 =   RMRV * ( SMK1*Q(J,K,I,4) + SMK3*Q(J,K,I,2) )
            DTXZD2 = - RMRV * SMK3
            DTXZD4 = - RMRV * SMK1
C
C     DTau(yz)
C
            DTYZD1 =   RMRV * ( SMK2*Q(J,K,I,4) + SMK3*Q(J,K,I,3) )
            DTYZD3 = - RMRV * SMK3
            DTYZD4 = - RMRV * SMK2
C
C     DQdot(x)
C
            TEMP   =   (RMUK(K,J)/PR + RMUTK(K,J)/PRT) * GAMMA * VOLAUX
     1                  / Q(J,K,I,1)
C
            DQDXD1 = - TEMP * SMK1 * ( Q(J,K,I,5)/GAMM1 
     1                   - 0.5E0*Q(J,K,I,1)*( Q(J,K,I,2)*Q(J,K,I,2)
     2                 +Q(J,K,I,3)*Q(J,K,I,3)+Q(J,K,I,4)*Q(J,K,I,4) ) )
     3                    / Q(J,K,I,1)
            DQDXD2 = - TEMP * SMK1 * Q(J,K,I,2) 
            DQDXD3 = - TEMP * SMK1 * Q(J,K,I,3) 
            DQDXD4 = - TEMP * SMK1 * Q(J,K,I,4) 
            DQDXD5 =   TEMP * SMK1 
C
C     DQdot(y)
C
            DQDYD1 = - TEMP * SMK2 * ( Q(J,K,I,5)/GAMM1
     1                   - 0.5E0*Q(J,K,I,1)*( Q(J,K,I,2)*Q(J,K,I,2)
     2                 +Q(J,K,I,3)*Q(J,K,I,3)+Q(J,K,I,4)*Q(J,K,I,4) ) )
     3                    / Q(J,K,I,1)
            DQDYD2 = - TEMP * SMK2 * Q(J,K,I,2)
            DQDYD3 = - TEMP * SMK2 * Q(J,K,I,3) 
            DQDYD4 = - TEMP * SMK2 * Q(J,K,I,4) 
            DQDYD5 =   TEMP * SMK2 
C
C     DQdot(z)
C
            DQDZD1 = - TEMP * SMK3 * ( Q(J,K,I,5)/GAMM1
     1                   - 0.5E0*Q(J,K,I,1)*( Q(J,K,I,2)*Q(J,K,I,2)
     2                 +Q(J,K,I,3)*Q(J,K,I,3)+Q(J,K,I,4)*Q(J,K,I,4) ) )
     3                    / Q(J,K,I,1)
            DQDZD2 = - TEMP * SMK3 * Q(J,K,I,2)
            DQDZD3 = - TEMP * SMK3 * Q(J,K,I,3)
            DQDZD4 = - TEMP * SMK3 * Q(J,K,I,4)
            DQDZD5 =   TEMP * SMK3
C
C     Calculate AM = dF(i+1/2)/dQ(i)
C
            AM(K,J,1,1) = 0.E0
            AM(K,J,1,2) = 0.E0
            AM(K,J,1,3) = 0.E0
            AM(K,J,1,4) = 0.E0
            AM(K,J,1,5) = 0.E0
C
            AM(K,J,2,1) = ( S(J,K,I,1,3)*DTXXD1 + S(J,K,I,2,3)*DTXYD1
     1                    + S(J,K,I,3,3)*DTXZD1 ) * AREAMR
            AM(K,J,2,2) = ( S(J,K,I,1,3)*DTXXD2 + S(J,K,I,2,3)*DTXYD2
     1                    + S(J,K,I,3,3)*DTXZD2 ) * AREAMR
            AM(K,J,2,3) = ( S(J,K,I,1,3)*DTXXD3 + S(J,K,I,2,3)*DTXYD3 )
     1                                            * AREAMR
            AM(K,J,2,4) = ( S(J,K,I,1,3)*DTXXD4 + S(J,K,I,3,3)*DTXZD4 ) 
     1                                            * AREAMR
            AM(K,J,2,5) = 0.E0
C
            AM(K,J,3,1) = ( S(J,K,I,1,3)*DTXYD1 + S(J,K,I,2,3)*DTYYD1
     1                    + S(J,K,I,3,3)*DTYZD1 ) * AREAMR
            AM(K,J,3,2) = ( S(J,K,I,1,3)*DTXYD2 + S(J,K,I,2,3)*DTYYD2 )
     1                                            * AREAMR
            AM(K,J,3,3) = ( S(J,K,I,1,3)*DTXYD3 + S(J,K,I,2,3)*DTYYD3
     1                    + S(J,K,I,3,3)*DTYZD3 ) * AREAMR
            AM(K,J,3,4) = ( S(J,K,I,2,3)*DTYYD4 + S(J,K,I,3,3)*DTYZD4 ) 
     1                                            * AREAMR
            AM(K,J,3,5) = 0.E0
C
            AM(K,J,4,1) = ( S(J,K,I,1,3)*DTXZD1 + S(J,K,I,2,3)*DTYZD1
     1                    + S(J,K,I,3,3)*DTZZD1 ) * AREAMR
            AM(K,J,4,2) = ( S(J,K,I,1,3)*DTXZD2 + S(J,K,I,3,3)*DTZZD2 ) 
     1                                            * AREAMR
            AM(K,J,4,3) = ( S(J,K,I,2,3)*DTYZD3 + S(J,K,I,3,3)*DTZZD3 ) 
     1                                            * AREAMR
            AM(K,J,4,4) = ( S(J,K,I,1,3)*DTXZD4 + S(J,K,I,2,3)*DTYZD4
     1                    + S(J,K,I,3,3)*DTZZD4 ) * AREAMR
            AM(K,J,4,5) = 0.E0
C
            AM(K,J,5,1) = AREAMR *
     1 ( S(J,K,I,1,3)*(U(K,J)*DTXXD1+V(K,J)*DTXYD1+W(K,J)*DTXZD1-DQDXD1)
     2 + S(J,K,I,2,3)*(U(K,J)*DTXYD1+V(K,J)*DTYYD1+W(K,J)*DTYZD1-DQDYD1)
     3 + S(J,K,I,3,3)*(U(K,J)*DTXZD1+V(K,J)*DTYZD1+W(K,J)*DTZZD1-DQDZD1)
     4 - 0.5E0*Q(J,K,I,2)/Q(J,K,I,1)*(S(J,K,I,1,3)*TAUXX
     5                    +S(J,K,I,2,3)*TAUXY+S(J,K,I,3,3)*TAUXZ)
     6 - 0.5E0*Q(J,K,I,3)/Q(J,K,I,1)*(S(J,K,I,1,3)*TAUXY
     7                    +S(J,K,I,2,3)*TAUYY+S(J,K,I,3,3)*TAUYZ)
     8 - 0.5E0*Q(J,K,I,4)/Q(J,K,I,1)*(S(J,K,I,1,3)*TAUXZ
     9                    +S(J,K,I,2,3)*TAUYZ+S(J,K,I,3,3)*TAUZZ) )
            AM(K,J,5,2) = AREAMR *
     1 ( S(J,K,I,1,3)*(U(K,J)*DTXXD2+V(K,J)*DTXYD2+W(K,J)*DTXZD2-DQDXD2)
     2 + S(J,K,I,2,3)*(U(K,J)*DTXYD2+V(K,J)*DTYYD2              -DQDYD2)
     3 + S(J,K,I,3,3)*(U(K,J)*DTXZD2              +W(K,J)*DTZZD2-DQDZD2)
     4 + 0.5E0/Q(J,K,I,1)*(S(J,K,I,1,3)*TAUXX
     5                    +S(J,K,I,2,3)*TAUXY+S(J,K,I,3,3)*TAUXZ) )
            AM(K,J,5,3) = AREAMR *
     1 ( S(J,K,I,1,3)*(U(K,J)*DTXXD3+V(K,J)*DTXYD3              -DQDXD3)
     2 + S(J,K,I,2,3)*(U(K,J)*DTXYD3+V(K,J)*DTYYD3+W(K,J)*DTYZD3-DQDYD3)
     3 + S(J,K,I,3,3)*(              V(K,J)*DTYZD3+W(K,J)*DTZZD3-DQDZD3)
     4 + 0.5E0/Q(J,K,I,1)*(S(J,K,I,1,3)*TAUXY
     5                    +S(J,K,I,2,3)*TAUYY+S(J,K,I,3,3)*TAUYZ) )
            AM(K,J,5,4) = AREAMR *
     1 ( S(J,K,I,1,3)*(U(K,J)*DTXXD4              +W(K,J)*DTXZD4-DQDXD4)
     2 + S(J,K,I,2,3)*(              V(K,J)*DTYYD4+W(K,J)*DTYZD4-DQDYD4)
     3 + S(J,K,I,3,3)*(U(K,J)*DTXZD4+V(K,J)*DTYZD4+W(K,J)*DTZZD4-DQDZD4)
     4 + 0.5E0/Q(J,K,I,1)*(S(J,K,I,1,3)*TAUXZ
     5                    +S(J,K,I,2,3)*TAUYZ+S(J,K,I,3,3)*TAUZZ) )
            AM(K,J,5,5) = - AREAMR *
     1                      ( S(J,K,I,1,3)*DQDXD5 + S(J,K,I,2,3)*DQDYD5 
     2                      + S(J,K,I,3,3)*DQDZD5 )
  510    CONTINUE
         DO 515 K = 1, KDIM
C
            AREAMR = S(J,K,I,4,3) * FSMACH / RE
C
C     Average metrics and cell volume
C
C     VOLAUX is the inverse of the average auxilliary cell volume
C
            VOLAUX = 2.0E0 / ( VOL(J,K,I) + VOL(J,K+1,I) )
            SPK1   = 0.5E0 * ( S(J,K+1,I,1,3)*S(J,K+1,I,4,3)
     1                       + S(J,K,I,1,3)*S(J,K,I,4,3) )
            SPK2   = 0.5E0 * ( S(J,K+1,I,2,3)*S(J,K+1,I,4,3)
     1                       + S(J,K,I,2,3)*S(J,K,I,4,3) )
            SPK3   = 0.5E0 * ( S(J,K+1,I,3,3)*S(J,K+1,I,4,3)
     1                       + S(J,K,I,3,3)*S(J,K,I,4,3) )
            SMK1   = 0.5E0 * ( S(J,K-1,I,1,3)*S(J,K-1,I,4,3)
     1                       + S(J,K,I,1,3)*S(J,K,I,4,3) )
            SMK2   = 0.5E0 * ( S(J,K-1,I,2,3)*S(J,K-1,I,4,3)
     1                       + S(J,K,I,2,3)*S(J,K,I,4,3) )
            SMK3   = 0.5E0 * ( S(J,K-1,I,3,3)*S(J,K-1,I,4,3)
     1                       + S(J,K,I,3,3)*S(J,K,I,4,3) )
C
C     Calculate first derivative terms using Divergence Theorem
C
            DUDX = VOLAUX * (Q(J,K+1,I,2)*SPK1 - Q(J,K,I,2)*SMK1)
            DUDY = VOLAUX * (Q(J,K+1,I,2)*SPK2 - Q(J,K,I,2)*SMK2)
            DUDZ = VOLAUX * (Q(J,K+1,I,2)*SPK3 - Q(J,K,I,2)*SMK3)
C
            DVDX = VOLAUX * (Q(J,K+1,I,3)*SPK1 - Q(J,K,I,3)*SMK1)
            DVDY = VOLAUX * (Q(J,K+1,I,3)*SPK2 - Q(J,K,I,3)*SMK2)
            DVDZ = VOLAUX * (Q(J,K+1,I,3)*SPK3 - Q(J,K,I,3)*SMK3)
C
            DWDX = VOLAUX * (Q(J,K+1,I,4)*SPK1 - Q(J,K,I,4)*SMK1)
            DWDY = VOLAUX * (Q(J,K+1,I,4)*SPK2 - Q(J,K,I,4)*SMK2)
            DWDZ = VOLAUX * (Q(J,K+1,I,4)*SPK3 - Q(J,K,I,4)*SMK3)
C
            DTDX = VOLAUX * (T(K+1,J)*SPK1 - T(K,J)*SMK1)
            DTDY = VOLAUX * (T(K+1,J)*SPK2 - T(K,J)*SMK2)
            DTDZ = VOLAUX * (T(K+1,J)*SPK3 - T(K,J)*SMK3)
C
C     Calculate Shear Terms
C     Use molecular viscosity (RMUK) and turbulent eddy viscosity (RMUTK)
C     averaged from cell centers to get viscosity at cell interface.
C     T23 is two-thirds.
C
            VISC  = RMUK(K,J) + RMUTK(K,J)
            TAUXX = T23 * VISC * ( 2.E0*DUDX - DVDY - DWDZ )
            TAUYY = T23 * VISC * ( 2.E0*DVDY - DUDX - DWDZ )
            TAUZZ = T23 * VISC * ( 2.E0*DWDZ - DUDX - DVDY )
            TAUXY = VISC * ( DUDY + DVDX )
            TAUXZ = VISC * ( DWDX + DUDZ )
            TAUYZ = VISC * ( DWDY + DVDZ )
C
C     Heat Flux terms are not used in the current implementation of the
C     implicit viscous terms.  This is kept for compatibility reasons
C     and in case they are later included in d(RMU) terms.
C     Calculate Heat Flux Terms
C
C           CNDCT = (RMUK(K,J)/PR + RMUTK(K,J)/PRT) / GAMM1
C           QDOTX = - CNDCT * DTDX
C           QDOTY = - CNDCT * DTDY
C           QDOTZ = - CNDCT * DTDZ
C
C     Region to calculate AP = dF(i+1/2)/dQ(i)  terms
C
C     Calculate DTau(xx,xy,xz,...)/DQ(1,2,3,...) terms
C
            RMRV   = VISC * VOLAUX / Q(J,K+1,I,1)
C
C     DTau(xx)
C
            DTXXD1 = - T23 * RMRV * ( 2.E0*SPK1*Q(J,K+1,I,2) 
     1                     - SPK2*Q(J,K+1,I,3) - SPK3*Q(J,K+1,I,4) )
            DTXXD2 =   T43 * RMRV * SPK1
            DTXXD3 = - T23 * RMRV * SPK2
            DTXXD4 = - T23 * RMRV * SPK3
C
C     DTau(yy)
C
            DTYYD1 = - T23 * RMRV * ( 2.E0*SPK2*Q(J,K+1,I,3) 
     1                     - SPK1*Q(J,K+1,I,2) - SPK3*Q(J,K+1,I,4) )
            DTYYD2 = - T23 * RMRV * SPK1
            DTYYD3 =   T43 * RMRV * SPK2
            DTYYD4 = - T23 * RMRV * SPK3
C
C     DTau(zz)
C
            DTZZD1 = - T23 * RMRV * ( 2.E0*SPK3*Q(J,K+1,I,4) 
     1                     - SPK1*Q(J,K+1,I,2) - SPK2*Q(J,K+1,I,3) )
            DTZZD2 = - T23 * RMRV * SPK1
            DTZZD3 = - T23 * RMRV * SPK2
            DTZZD4 =   T43 * RMRV * SPK3
C
C     DTau(xy)
C
            DTXYD1 = - RMRV * ( SPK1*Q(J,K+1,I,3) + SPK2*Q(J,K+1,I,2) )
            DTXYD2 =   RMRV * SPK2
            DTXYD3 =   RMRV * SPK1
C
C     DTau(xz)
C
            DTXZD1 = - RMRV * ( SPK1*Q(J,K+1,I,4) + SPK3*Q(J,K+1,I,2) )
            DTXZD2 =   RMRV * SPK3
            DTXZD4 =   RMRV * SPK1
C
C     DTau(yz)
C
            DTYZD1 = - RMRV * ( SPK2*Q(J,K+1,I,4) + SPK3*Q(J,K+1,I,3) )
            DTYZD3 =   RMRV * SPK3
            DTYZD4 =   RMRV * SPK2
C
C     DQdot(x)
C
            TEMP   =   (RMUK(K,J)/PR + RMUTK(K,J)/PRT) * GAMMA * VOLAUX
     1                  / Q(J,K+1,I,1)
C
            DQDXD1 =   TEMP * SPK1 * ( Q(J,K+1,I,5)/GAMM1 
     1             - 0.5E0*Q(J,K+1,I,1)*( Q(J,K+1,I,2)*Q(J,K+1,I,2)
     2         +Q(J,K+1,I,3)*Q(J,K+1,I,3)+Q(J,K+1,I,4)*Q(J,K+1,I,4) ) )
     3                    / Q(J,K+1,I,1)
            DQDXD2 =   TEMP * SPK1 * Q(J,K+1,I,2) 
            DQDXD3 =   TEMP * SPK1 * Q(J,K+1,I,3) 
            DQDXD4 =   TEMP * SPK1 * Q(J,K+1,I,4) 
            DQDXD5 = - TEMP * SPK1 
C
C     DQdot(y)
C
            DQDYD1 =   TEMP * SPK2 * ( Q(J,K+1,I,5)/GAMM1
     1             - 0.5E0*Q(J,K+1,I,1)*( Q(J,K+1,I,2)*Q(J,K+1,I,2)
     2         +Q(J,K+1,I,3)*Q(J,K+1,I,3)+Q(J,K+1,I,4)*Q(J,K+1,I,4) ) )
     3                    / Q(J,K+1,I,1)
            DQDYD2 =   TEMP * SPK2 * Q(J,K+1,I,2)
            DQDYD3 =   TEMP * SPK2 * Q(J,K+1,I,3) 
            DQDYD4 =   TEMP * SPK2 * Q(J,K+1,I,4) 
            DQDYD5 = - TEMP * SPK2 
C
C     DQdot(z)
C
            DQDZD1 =   TEMP * SPK3 * ( Q(J,K+1,I,5)/GAMM1
     1             - 0.5E0*Q(J,K+1,I,1)*( Q(J,K+1,I,2)*Q(J,K+1,I,2)
     2         +Q(J,K+1,I,3)*Q(J,K+1,I,3)+Q(J,K+1,I,4)*Q(J,K+1,I,4) ) )
     3                    / Q(J,K+1,I,1)
            DQDZD2 =   TEMP * SPK3 * Q(J,K+1,I,2)
            DQDZD3 =   TEMP * SPK3 * Q(J,K+1,I,3)
            DQDZD4 =   TEMP * SPK3 * Q(J,K+1,I,4)
            DQDZD5 = - TEMP * SPK3
C
C     Calculate AP = dF(i+1/2)/dQ(i+1)
C
            AP(K,J,1,1) = 0.E0
            AP(K,J,1,2) = 0.E0
            AP(K,J,1,3) = 0.E0
            AP(K,J,1,4) = 0.E0
            AP(K,J,1,5) = 0.E0
C
            AP(K,J,2,1) = ( S(J,K,I,1,3)*DTXXD1 + S(J,K,I,2,3)*DTXYD1
     1                    + S(J,K,I,3,3)*DTXZD1 ) * AREAMR
            AP(K,J,2,2) = ( S(J,K,I,1,3)*DTXXD2 + S(J,K,I,2,3)*DTXYD2
     1                    + S(J,K,I,3,3)*DTXZD2 ) * AREAMR
            AP(K,J,2,3) = ( S(J,K,I,1,3)*DTXXD3 + S(J,K,I,2,3)*DTXYD3 )
     1                                            * AREAMR
            AP(K,J,2,4) = ( S(J,K,I,1,3)*DTXXD4 + S(J,K,I,3,3)*DTXZD4 ) 
     1                                            * AREAMR
            AP(K,J,2,5) = 0.E0
C
            AP(K,J,3,1) = ( S(J,K,I,1,3)*DTXYD1 + S(J,K,I,2,3)*DTYYD1
     1                    + S(J,K,I,3,3)*DTYZD1 ) * AREAMR
            AP(K,J,3,2) = ( S(J,K,I,1,3)*DTXYD2 + S(J,K,I,2,3)*DTYYD2 )
     1                                            * AREAMR
            AP(K,J,3,3) = ( S(J,K,I,1,3)*DTXYD3 + S(J,K,I,2,3)*DTYYD3
     1                    + S(J,K,I,3,3)*DTYZD3 ) * AREAMR
            AP(K,J,3,4) = ( S(J,K,I,2,3)*DTYYD4 + S(J,K,I,3,3)*DTYZD4 ) 
     1                                            * AREAMR
            AP(K,J,3,5) = 0.E0
C
            AP(K,J,4,1) = ( S(J,K,I,1,3)*DTXZD1 + S(J,K,I,2,3)*DTYZD1
     1                    + S(J,K,I,3,3)*DTZZD1 ) * AREAMR
            AP(K,J,4,2) = ( S(J,K,I,1,3)*DTXZD2 + S(J,K,I,3,3)*DTZZD2 ) 
     1                                            * AREAMR
            AP(K,J,4,3) = ( S(J,K,I,2,3)*DTYZD3 + S(J,K,I,3,3)*DTZZD3 ) 
     1                                            * AREAMR
            AP(K,J,4,4) = ( S(J,K,I,1,3)*DTXZD4 + S(J,K,I,2,3)*DTYZD4
     1                    + S(J,K,I,3,3)*DTZZD4 ) * AREAMR
            AP(K,J,4,5) = 0.E0
C
            AP(K,J,5,1) = AREAMR *
     1 ( S(J,K,I,1,3)*(U(K,J)*DTXXD1+V(K,J)*DTXYD1+W(K,J)*DTXZD1-DQDXD1)
     2 + S(J,K,I,2,3)*(U(K,J)*DTXYD1+V(K,J)*DTYYD1+W(K,J)*DTYZD1-DQDYD1)
     3 + S(J,K,I,3,3)*(U(K,J)*DTXZD1+V(K,J)*DTYZD1+W(K,J)*DTZZD1-DQDZD1)
     4 - 0.5E0*Q(J,K+1,I,2)/Q(J,K+1,I,1)*(S(J,K,I,1,3)*TAUXX
     5                    +S(J,K,I,2,3)*TAUXY+S(J,K,I,3,3)*TAUXZ)
     6 - 0.5E0*Q(J,K+1,I,3)/Q(J,K+1,I,1)*(S(J,K,I,1,3)*TAUXY
     7                    +S(J,K,I,2,3)*TAUYY+S(J,K,I,3,3)*TAUYZ)
     8 - 0.5E0*Q(J,K+1,I,4)/Q(J,K+1,I,1)*(S(J,K,I,1,3)*TAUXZ
     9                    +S(J,K,I,2,3)*TAUYZ+S(J,K,I,3,3)*TAUZZ) )
            AP(K,J,5,2) = AREAMR *
     1 ( S(J,K,I,1,3)*(U(K,J)*DTXXD2+V(K,J)*DTXYD2+W(K,J)*DTXZD2-DQDXD2)
     2 + S(J,K,I,2,3)*(U(K,J)*DTXYD2+V(K,J)*DTYYD2              -DQDYD2)
     3 + S(J,K,I,3,3)*(U(K,J)*DTXZD2              +W(K,J)*DTZZD2-DQDZD2)
     4 + 0.5E0/Q(J,K+1,I,1)*(S(J,K,I,1,3)*TAUXX
     5                      +S(J,K,I,2,3)*TAUXY+S(J,K,I,3,3)*TAUXZ) )
            AP(K,J,5,3) = AREAMR *
     1 ( S(J,K,I,1,3)*(U(K,J)*DTXXD3+V(K,J)*DTXYD3              -DQDXD3)
     2 + S(J,K,I,2,3)*(U(K,J)*DTXYD3+V(K,J)*DTYYD3+W(K,J)*DTYZD3-DQDYD3)
     3 + S(J,K,I,3,3)*(              V(K,J)*DTYZD3+W(K,J)*DTZZD3-DQDZD3)
     4 + 0.5E0/Q(J,K+1,I,1)*(S(J,K,I,1,3)*TAUXY
     5                      +S(J,K,I,2,3)*TAUYY+S(J,K,I,3,3)*TAUYZ) )
            AP(K,J,5,4) = AREAMR *
     1 ( S(J,K,I,1,3)*(U(K,J)*DTXXD4              +W(K,J)*DTXZD4-DQDXD4)
     2 + S(J,K,I,2,3)*(              V(K,J)*DTYYD4+W(K,J)*DTYZD4-DQDYD4)
     3 + S(J,K,I,3,3)*(U(K,J)*DTXZD4+V(K,J)*DTYZD4+W(K,J)*DTZZD4-DQDZD4)
     4 + 0.5E0/Q(J,K+1,I,1)*(S(J,K,I,1,3)*TAUXZ
     5                      +S(J,K,I,2,3)*TAUYZ+S(J,K,I,3,3)*TAUZZ) )
            AP(K,J,5,5) = - AREAMR *
     1                      ( S(J,K,I,1,3)*DQDXD5 + S(J,K,I,2,3)*DQDYD5 
     2                      + S(J,K,I,3,3)*DQDZD5 )
  515    CONTINUE
  520 CONTINUE
C
C     Finished with implicit Thin Layer Viscous Terms in the K direction
C
      RETURN
      END






      SUBROUTINE DSRCDQ (IDIM, JDIM, KDIM, I, IMODEL, 
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   QC, PROPSC, DQDX, DQDY, DQDZ, DEL2K, 
     3                   TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ,
     4                   TAU, SRC, IBCDAT, NBCS, ASRC)
C
C Routine to calculate Jacobians for the source terms.  The source terms have
C previously been calculated and are passed into this routine.
C
C This routine is organized as follows:
C     1.   Initialize source Jacobians to zero and update non-zero elements.
C     2.   Call the appropriate source term Jacobian routine
C          a. k-epsilon Two-Equation turbulence model
C          b. k-omega   Two-Equation turbulence model
C          c. Reynolds Stress turbulence model
C          d. Error
C     3.   Scale Jacobian by the volume
C
C IDIM,JDIM,KDIM : Dimensions of current block
C I              : The current plane
C IMODEL         : Turbulence model array
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
C QC             : Storage for Q in one-dimensional array
C PROPSC         : Storage for PROPS array in one-dimensional array
C DQDX,DQDY,DQDZ : Storage for derivatives of Q at cell centers
C DEL2K          : Laplacian of turbulent kinetic energy
C TKE            : Turbulnet kinetic energy
C ETILDE         : Epsilon^tilde
C XSIW           :
C FDAMP          : Damping functions F1, F2
C EPSC           :
C PDIL           :
C SIJ            :
C WIJ            :
C TAU            : Storage for Reynolds stress terms
C SRC            : Source terms for cell (J,K)
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C ASRC           : The Jacobian matrix for the source terms
C
      include '../header/common.h'
C
      DIMENSION IMODEL(NMDL)
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION QC    (2:JDIM,2:KDIM,NQ),
     1          PROPSC(2:JDIM,2:KDIM,NP),
     2          DQDX  (2:JDIM,2:KDIM,NQ),
     3          DQDY  (2:JDIM,2:KDIM,NQ),
     4          DQDZ  (2:JDIM,2:KDIM,NQ),
     5          DEL2K (2:JDIM,2:KDIM),
     6          TKE   (2:JDIM,2:KDIM),
     7          ETILDE(2:JDIM,2:KDIM),
     8          XSIW  (2:JDIM,2:KDIM),
     9          FDAMP (2:JDIM,2:KDIM,2),
     A          EPSC  (2:JDIM,2:KDIM),
     B          PDIL  (2:JDIM,2:KDIM),
     C          SIJ   (2:JDIM,2:KDIM,6),
     D          WIJ   (2:JDIM,2:KDIM,6),
     E          TAU   (2:JDIM,2:KDIM,6),
     F          SRC   (2:JDIM,2:KDIM,NF)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION ASRC  (2:JDIM,2:KDIM,NRANK,NRANK)
C
      ITURB  = IMODEL(1)
C
C     Number of points in a cross-section to calculate source terms Jacobians
C
      NPTS   = (JDIM - 1) * (KDIM - 1)
C
C 1.   Initialize source Jacobians to zero and update only non-zero elements.
C
      DO 140 M = 1, NRANK
         DO 130 L = 1, NRANK
            DO 120 K = 2, KDIM
               DO 110 J = 2, JDIM
                  ASRC(J,K,L,M) = 0.E0
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
C
C 2.   Call the appropriate source term Jacobian routine
C 2-a. k-epsilon Two-Equation turbulence model
C
      IF (ITURB .EQ. ITKE) THEN
         CALL DSRCKE (NPTS, QC, PROPSC, DQDX, DQDY, DQDZ, TAU, 
     1                TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ, 
     2                SRC, ASRC)
C
C 2-b. k-omega Two-Equation turbulence model
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         CALL DSRCKW (NPTS, QC, PROPSC, DQDX, DQDY, DQDZ, TAU, ASRC)
C
C 2-c. Reynolds Stress turbulence model
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         CALL DSRCRS (IMODEL, NPTS, QC, PROPSC, DQDX, DQDY, DQDZ, TAU, 
     1                TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ, 
     2                SRC, ASRC)
C
C 2-d. Error
C
      ELSE
         WRITE (IOUT,1000) ITURB
         STOP
      ENDIF
C
C 3.   Scale Jacobian by the volume
C
      DO 340 M = 1, NRANK
         DO 330 L = 1, NRANK
            DO 320 K = 2, KDIM
               DO 310 J = 2, JDIM
                  ASRC(J,K,L,M) = ASRC(J,K,L,M) * VOL(J,K,I)
  310          CONTINUE
  320       CONTINUE
  330    CONTINUE
  340 CONTINUE
C
C     Finished calculating Jacobians
C
      RETURN
 1000 FORMAT (' ','DSRCDQ: ERROR-> Invalid turbulence model: ',I5,
     1       /' ',8X,'RUN ABORTING!')
      END






      SUBROUTINE DSRCKE (NPTS, Q, PROPS, DQDX, DQDY, DQDZ, TAU, 
     1                   TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ,
     2                   SRC, ASRC)
C
C Routine to calculate the implicit source term Jacobians for the 
C k-epsilon turbulence model.
C
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C PROPS          : Properties stored at cell centers
C                  PROPS(1) = RMU   molecular viscosity
C                  PROPS(2) = RMUT  turbulent eddy viscosity
C                  PROPS(3) = YPLUS Y+
C DQDX,DQDY,DQDZ : Derivatives of Q at cell centers
C TAU            : Reynolds stresses
C TKE            : Turbulent kinetic energy
C ETILDE         : Epsilon^tilde
C XSIW           :
C FDAMP          : Damping functions F1, F2
C EPSC           :
C PDIL           :
C SIJ            :
C WIJ            :
C SRC            :
C ASRC           : Source terms for the k-omega model at the cell centers
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
     4          TKE   (NPTS),
     5          ETILDE(NPTS),
     6          XSIW  (NPTS),
     7          FDAMP (NPTS,2),
     8          EPSC  (NPTS),
     9          PDIL  (NPTS),
     A          SIJ   (NPTS,6),
     B          WIJ   (NPTS,6),
     C          SRC   (NPTS,NF),
     D          ASRC  (NPTS,NRANK,NRANK)
C
C     One-Third
C     Two-Thirds
C
      T13 = 1.E0 / 3.E0
      T23 = 2.E0 / 3.E0
C
C     Form the source Jacobians for the k-epsilon model.
C     Update only the non-zero elements.
C
      DO 100 I = 1, NPTS
C
         RHO   = Q    (I,1)
         P     = Q    (I,5)
         EPSLN = Q    (I,7)
C
         RMU   = PROPS(I,1)
         RMUT  = PROPS(I,2)
         YPLUS = PROPS(I,3)
C
         SXX   = SIJ  (I,1)
         SYY   = SIJ  (I,2)
         SZZ   = SIJ  (I,3)
         SXY   = SIJ  (I,4)
         SXZ   = SIJ  (I,5)
         SYZ   = SIJ  (I,6)
C
         SKK   = SXX + SYY + SZZ
         TERM  =         SXX*SXX + SYY*SYY + SZZ*SZZ
     1           + 2.E0*(SXY*SXY + SXZ*SXZ + SYZ*SYZ) - T13*SKK*SKK
C
         DPDK  =   4.E0*RMUT/(RHO*TKE(I))*TERM*FSMACH/RE - T23*SKK
         DPDE  = - 2.E0*RMUT/(RHO*EPSLN )*TERM*FSMACH/RE
C
C     Calculate Implicit Source Term
C
c        write (iout,'(a,i5,2(1pe14.4))') ' SRC : i,', i, 
c    1               src(i,6), src(i,7)
C-----Begin DS/DQ (complete)------------------------------
c        ASRC(I,6,6) = DPDK
c        ASRC(I,6,7) = DPDE - 1.E0
C
c        ASRC(I,7,6) = CEPS1*FDAMP(I,1)*EPSLN/TKE(I)*DPDK
c    1                 - SRC(I,7) / (RHO*TKE(I))
c        ASRC(I,7,7) = EPSLN/TKE(I)*( CEPS1*FDAMP(I,1)*DPDE
c    1                              - CEPS2*FDAMP(I,2) )
c    2                 + SRC(I,7) / (RHO*EPSLN)
c        write (iout,'(a,i5,4(1pe14.4))') ' Full: i,', i, 
c    1               asrc(i,6,6), asrc(i,6,7), asrc(i,7,6), asrc(i,7,7)
C-----End   DS/DQ (complete)------------------------------
C
C-----Begin DS/DQ (van Drommme)---------------------------
         TIMSCL      =   EPSLN / TKE(I)
         ASRC(I,6,6) = - 2.E0 * TIMSCL
         ASRC(I,6,7) =   0.E0
C
         ASRC(I,7,6) =          CEPS2 * FDAMP(I,2) * TIMSCL * TIMSCL
         ASRC(I,7,7) = - 2.E0 * CEPS2 * FDAMP(I,2) * TIMSCL
c        write (iout,'(a,i5,4(1pe14.4))') ' VanD: i,', i, 
c    1               asrc(i,6,6), asrc(i,6,7), asrc(i,7,6), asrc(i,7,7)
C-----End   DS/DQ (van Drommme)---------------------------
C-----Begin DS/DQ (Blended DS(-)/DQ always + DS(+)/DQ if SRC < 0)-----
c        S6          = 0.5E0 * (1.E0 - SIGN (1.E0, SRC(I,6)))
c        S7          = 0.5E0 * (1.E0 - SIGN (1.E0, SRC(I,7)))
C
c        TIMSCL      =   EPSLN / TKE(I)
c        ASRC(I,6,6) = S6 * DPDK - 2.E0 * TIMSCL
c        ASRC(I,6,7) = S6 * DPDE
C
c        ASRC(I,7,6) = S7 * ( CEPS1*FDAMP(I,1)*EPSLN/TKE(I)*DPDK
c    1                      - SRC(I,7) / (RHO*TKE(I)) )
c    2               + (1.E0-S7) * CEPS2 * FDAMP(I,2) * TIMSCL * TIMSCL
c        ASRC(I,7,7) = S7 * ( EPSLN/TKE(I)*( CEPS1*FDAMP(I,1)*DPDE
c    1                                   - CEPS2*FDAMP(I,2) )
c    2                      + SRC(I,7) / (RHO*EPSLN) )
c    3               + (1.E0-S7) * (-2.E0 * CEPS2 * FDAMP(I,2) * TIMSCL)
c        write (iout,'(a,i5,4(1pe14.4))') ' Blen: i,', i, 
c    1               asrc(i,6,6), asrc(i,6,7), asrc(i,7,6), asrc(i,7,7)
C-----End   DS/DQ (Blended DS(-)/DQ always + DS(+)/DQ if SRC < 0)-----
C
  100 CONTINUE
C
C     Finished with implicit k-epsilon source terms
C
      RETURN
      END






      SUBROUTINE DSRCKW (NPTS, Q, PROPS, DQDX, DQDY, DQDZ, TAU, ASRC)
C
C     Routine to calculate the implicit source Jacobians for the k-omega
C     turbulence model.
C
C     NPTS           : Number of points to calculate source terms at
C     Q              : Primitive variables at cell centers
C     PROPS          : Properties stored at cell centers
C                      PROPS(1) = RMU   molecular viscosity
C                      PROPS(2) = RMUT  turbulent eddy viscosity
C                      PROPS(3) = YPLUS Y+
C     DQDX,DQDY,DQDZ : Derivatives of Q at cell centers
C     TAU            : Reynolds stresses
C     ASRC           : Source Jacobians for the k-omega model at cell centers
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
     4          ASRC  (NPTS,NRANK,NRANK)
C
C
C     Two-Thirds
C
      T23 = 2.E0 / 3.E0
C
C     Form the implicit source Jacobians for the k-omega model.
C     Modify only the non-zero elements.
C
      DO 100 I = 1, NPTS
C
         RHO   = Q(I,1)
         TKE   = Q(I,6)
         OMEGA = Q(I,7)
C
C     Production term:
C        PROD = TAUXX*DUDX + TAUYY*DVDY + TAUZZ*DWDZ +
C               TAUXY*(DUDY+DVDX) + TAUXZ*(DUDZ+DWDX) + TAUYZ*(DVDZ+DWDY)
C
         PROD  =   TAU(I,1) * DQDX(I,2)
     1           + TAU(I,2) * DQDY(I,3)
     2           + TAU(I,3) * DQDZ(I,4)
     3           + TAU(I,4) * (DQDY(I,2) + DQDX(I,3))
     4           + TAU(I,5) * (DQDZ(I,2) + DQDX(I,4))
     5           + TAU(I,6) * (DQDZ(I,3) + DQDY(I,4))
C
C     Calculate Implicit Source Term
C
         ASRC(I,6,1) =   BSTRKW * TKE * OMEGA
         ASRC(I,6,6) = - BSTRKW * OMEGA
         ASRC(I,6,7) = - BSTRKW * TKE
C
         ASRC(I,7,1) =   BKW * OMEGA * OMEGA
         ASRC(I,7,7) = - 2.E0 * BKW * OMEGA
  100 CONTINUE
C
C     Finished with implicit k-omega source terms
C
      RETURN
      END





      SUBROUTINE DSRCRS (IMODEL, NPTS, Q, PROPS, DQDX, DQDY, DQDZ, TAU, 
     1                   TKE, ETILDE, XSIW, FDAMP, EPSC, PDIL, SIJ, WIJ,
     2                   SRC, ASRC)
C
C Routine to calculate the source term Jacobians for the Reynolds stress
C turbulence model.
C
C IMODEL         : Turbulence model array
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C PROPS          : Properties stored at cell centers
C                  PROPS(1) = RMU   molecular viscosity
C                  PROPS(2) = RMUT  turbulent eddy viscosity
C                  PROPS(3) = YPLUS Y+
C DQDX,DQDY,DQDZ : Derivatives of Q at cell centers
C TAU            : Reynolds stresses
C TKE            : Turbulent kinetic energy
C ETILDE         : Epsilon^tilde
C XSIW           :
C FDAMP          : Damping functions F1, F2
C EPSC           :
C PDIL           :
C SIJ            :
C WIJ            :
C SRC            :
C ASRC           : Source terms for the k-omega model at the cell centers
C
      include '../header/common.h'
C
      DIMENSION IMODEL(NMDL)
C
      DIMENSION Q     (NPTS,NQ),
     1          PROPS (NPTS,NP)
C
      DIMENSION DQDX  (NPTS,NQ),
     1          DQDY  (NPTS,NQ),
     2          DQDZ  (NPTS,NQ),
     3          TAU   (NPTS,6),
     4          TKE   (NPTS),
     5          ETILDE(NPTS),
     6          XSIW  (NPTS),
     7          FDAMP (NPTS,2),
     8          EPSC  (NPTS),
     9          PDIL  (NPTS),
     A          SIJ   (NPTS,6),
     B          WIJ   (NPTS,6),
     C          SRC   (NPTS,NF),
     D          ASRC  (NPTS,NRANK,NRANK)
C
C     Two-Thirds
C
      T23 = 2.E0 / 3.E0
C
      ITURB  = IMODEL(1)
      IPSMDL = IMODEL(6)
      IEPSIJ = IMODEL(7)
C
C     Call routine to calculate Pressure-Strain Correlation Jacobians
C
      CALL DPSTRN (IPSMDL, NPTS, Q, PROPS, DQDX, DQDY, DQDZ, ASRC)
C
C     Form the source terms for the Reynolds stress model 
C
      DO 100 I = 1, NPTS
C
         RHO    = Q(I,1)
         P      = Q(I,5)
         TAUXX  = Q(I,6)
         TAUYY  = Q(I,7)
         TAUZZ  = Q(I,8)
         TAUXY  = Q(I,9)
         TAUXZ  = Q(I,10)
         TAUYZ  = Q(I,11)
         EPSLN  = Q(I,12)
C
         RMU    = PROPS(I,1)
         RMUT   = PROPS(I,2)
         YPLUS  = PROPS(I,3)
         YNORML = ABS (PROPS(I,4))
C
C     DU_k/DX_k = (DU/DX + DV/DY + DW/DZ) / 3
C
         DUKDXK = (DQDX(I,2) + DQDY(I,3) + DQDZ(I,4)) / 3.E0
C
C     Production term for epsilon:
C        PRODE = - CEPS1*epsilon/k*rho*tau_ij*(Du_i/Dx_j-Du_k/Dx_k*delta_ij/3)
C                - 4/3*rho*epsilon*Du_k/Dx_k
C        from: Speziale and Sarkar, ICASE 91-9
C        wall damping from: Shima, F. Fluids Eng., Vol 111
C
         ARG    = RHO * SQRT(TKE(I)) * YNORML / RMU * RE / FSMACH
         FW     = EXP( -(0.015E0 * ARG)**4 )
         CEPSTR = CEPS1 * (1.E0 + CEPS4*FW)
         CTAUDU = CEPSTR * RHO *
     1             ( TAUXX * (DQDX(I,2) - DUKDXK)
     2             + TAUYY * (DQDY(I,3) - DUKDXK)
     3             + TAUZZ * (DQDZ(I,4) - DUKDXK)
     4             + TAUXY * (DQDY(I,2) + DQDX(I,3)) 
     5             + TAUXZ * (DQDZ(I,2) + DQDX(I,4))
     6             + TAUYZ * (DQDZ(I,3) + DQDY(I,4)) )
C
C     Turbulent Mach Number (from Sarkar et al.) for compressibility
C     correction to dissipation rate
C
         ASQ    = GAMMA * P / RHO
         RMTSQ  = 2.E0 * TKE(I) / ASQ
C
C     Dissipation
C
         CDISS  = T23 * (1.E0 + ALF1 * RMTSQ) 
         DISS   = T23 * (1.E0 + ALF1 * RMTSQ) * RHO * EPSLN
C
C     Dissipation for epsilon from Shima
C
         DISSE  = CEPS2 * RHO * EPSLN * EPSLN / TKE(I) - RHO * XSIW(I)
C
C     Calculate the Pressure-Dilatation from Sarkar, ICASE 91-42
C
c        PDIL  = T23 * (-ALF2*PROD + ALF3*RHO*EPSLN) * RMTSQ
C
C     Calculate Source Jacobians due to Production (P_ij)
C     and Dissipation terms
C     Calculate only non-zero terms.  Add these to terms calculated in
C     the pressure-strain Jacobian routine.
C
         ASRC(I, 6, 6) = ASRC(I, 6, 6) - 2.E0 * DQDX(I,2)
         ASRC(I, 6, 9) = ASRC(I, 6, 9) - 2.E0 * DQDY(I,2)
         ASRC(I, 6,10) = ASRC(I, 6,10) - 2.E0 * DQDZ(I,2)
         ASRC(I, 6,12) = ASRC(I, 6,12) - CDISS
C
         ASRC(I, 7, 7) = ASRC(I, 7, 7) - 2.E0 * DQDY(I,3)
         ASRC(I, 7, 9) = ASRC(I, 7, 9) - 2.E0 * DQDX(I,3)
         ASRC(I, 7,11) = ASRC(I, 7,11) - 2.E0 * DQDZ(I,3)
         ASRC(I, 7,12) = ASRC(I, 7,12) - CDISS
C
         ASRC(I, 8, 8) = ASRC(I, 8, 8) - 2.E0 * DQDZ(I,4)
         ASRC(I, 8,10) = ASRC(I, 8,10) - 2.E0 * DQDX(I,4)
         ASRC(I, 8,11) = ASRC(I, 8,11) - 2.E0 * DQDY(I,4)
         ASRC(I, 8,12) = ASRC(I, 8,12) - CDISS
C
         ASRC(I, 9, 6) = ASRC(I, 9, 6) -  DQDX(I,3)
         ASRC(I, 9, 7) = ASRC(I, 9, 7) -  DQDY(I,2)
         ASRC(I, 9, 9) = ASRC(I, 9, 9) - (DQDY(I,3) + DQDX(I,2))
         ASRC(I, 9,10) = ASRC(I, 9,10) -  DQDZ(I,3)
         ASRC(I, 9,11) = ASRC(I, 9,11) -  DQDZ(I,2)
C
         ASRC(I,10, 6) = ASRC(I,10, 6) -  DQDX(I,4)
         ASRC(I,10, 8) = ASRC(I,10, 8) -  DQDZ(I,2)
         ASRC(I,10, 9) = ASRC(I,10, 9) -  DQDY(I,4)
         ASRC(I,10,10) = ASRC(I,10,10) - (DQDZ(I,4) + DQDX(I,2))
         ASRC(I,10,11) = ASRC(I,10,11) -  DQDY(I,2)
C
         ASRC(I,11, 7) = ASRC(I,11, 7) -  DQDY(I,4)
         ASRC(I,11, 8) = ASRC(I,11, 8) -  DQDZ(I,3)
         ASRC(I,11, 9) = ASRC(I,11, 9) -  DQDX(I,4)
         ASRC(I,11,10) = ASRC(I,11,10) -  DQDX(I,3)
         ASRC(I,11,11) = ASRC(I,11,11) - (DQDZ(I,4) + DQDY(I,3))
C
         ASRC(I,12, 6) = ASRC(I,12, 6) 
     1              - CEPSTR * EPSLN / TKE(I) * (DQDX(I,2) - DUKDXK)
     2              + 0.5E0 * CTAUDU * EPSLN / (RHO * TKE(I) * TKE(I))
     3              + 0.5E0 * CEPS2 * EPSLN * EPSLN / (TKE(I) * TKE(I))
         ASRC(I,12, 7) = ASRC(I,12, 7) 
     1              - CEPSTR * EPSLN / TKE(I) * (DQDY(I,3) - DUKDXK)
     2              + 0.5E0 * CTAUDU * EPSLN / (RHO * TKE(I) * TKE(I))
     3              + 0.5E0 * CEPS2 * EPSLN * EPSLN / (TKE(I) * TKE(I))
         ASRC(I,12, 8) = ASRC(I,12, 8) 
     1              - CEPSTR * EPSLN / TKE(I) * (DQDZ(I,4) - DUKDXK)
     2              + 0.5E0 * CTAUDU * EPSLN / (RHO * TKE(I) * TKE(I))
     3              + 0.5E0 * CEPS2 * EPSLN * EPSLN / (TKE(I) * TKE(I))
         ASRC(I,12, 9) = ASRC(I,12, 9) 
     1              - CEPSTR * EPSLN / TKE(I) * (DQDY(I,2)+DQDX(I,3)) 
         ASRC(I,12,10) = ASRC(I,12,10) 
     1              - CEPSTR * EPSLN / TKE(I) * (DQDZ(I,2)+DQDX(I,4))
         ASRC(I,12,11) = ASRC(I,12,11) 
     1              - CEPSTR * EPSLN / TKE(I) * (DQDZ(I,3)+DQDY(I,4))
         ASRC(I,12,12) = ASRC(I,12,12) 
     1              - CTAUDU / (RHO * TKE(I))
     2              - 4.E0 * DUKDXK
     3              - 2.E0 * CEPS2 * EPSLN / TKE(I)
  100 CONTINUE
C
C     Finished with Reynolds Stress source terms
C
      RETURN
      END






      SUBROUTINE DPSTRN (IPSMDL, NPTS, Q, PROPS, DQDX, DQDY, DQDZ, ASRC)
C
C     Routine to calculate the Pressure-Strain Correlation Jacobians for a 
C     Reynolds stress turbulence model.
C
C     IPSMDL         : Pressure-strain model
C     NPTS           : Number of points to calculate source terms at
C     Q              : Primitive variables at cell centers
C     PROPS          : Properties stored at cell centers
C                      PROPS(1) = RMU   molecular viscosity
C                      PROPS(2) = RMUT  turbulent eddy viscosity
C                      PROPS(3) = YPLUS Y+
C     DQDX,DQDY,DQDZ : Derivatives of Q at cell centers
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          PROPS (NPTS,NP)
C
      DIMENSION DQDX  (NPTS,NQ),
     1          DQDY  (NPTS,NQ),
     2          DQDZ  (NPTS,NQ),
     3          ASRC  (NPTS,NRANK,NRANK)
C
C     One-Third
C     Two-Thirds
C     Four-Thirds
C
      T13 = 1.E0 / 3.E0
      T23 = 2.E0 / 3.E0
      T43 = 4.E0 / 3.E0
C
C     Calculate Pressure-Strain Correlation
C
      IF (IPSMDL .EQ. IPGL) THEN
C
C     Gison-Launder Pressure-Strain Correlation
C     Launder-Reece-Rodi "Basic Model" Pressure-Strain Correlation.
C     The form for this model is from Sarkar&Balakrishnan, ICASE 90-18
C
         write (iout,'(a)') ' Implicit Gibson-Launder not supported '
         stop
C
C     Shima Pressure Strain Correlation
C     J. Fluids Engineering, Vol. 110
C*NOTE ALSO USE FOR LAI-SO, SSG, FLT FOR NOW
C
      ELSE IF (ipsmdl .eq. ipls .or.
     1         ipsmdl .eq. iplrr1 .or. ipsmdl .eq. ipssg .or.
     2         ipsmdl .eq. ipflt ) THEN
         CPHI2A = (PSC2 + 8.E0) / 11.E0
         CPHI2B = (30.E0*PSC2 - 2.E0) / 55.E0
         CPHI2C = (8.E0 *PSC2 - 2.E0) / 11.E0
         DO 200 I = 1, NPTS
C
            RHO    = Q(I,1)
            P      = Q(I,5)
            TAUXX  = Q(I,6)
            TAUYY  = Q(I,7)
            TAUZZ  = Q(I,8)
            TAUXY  = Q(I,9)
            TAUXZ  = Q(I,10)
            TAUYZ  = Q(I,11)
            EPSLN  = Q(I,12)
C
            TKE    = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
            RMU    = PROPS(I,1)
            YNORML = ABS (PROPS(I,4))
C
C     Anisotropy tensor
C
            BXX    = 0.5E0 * TAUXX / TKE - T13
            BYY    = 0.5E0 * TAUYY / TKE - T13
            BZZ    = 0.5E0 * TAUZZ / TKE - T13
            BXY    = 0.5E0 * TAUXY / TKE
            BXZ    = 0.5E0 * TAUXZ / TKE
            BYZ    = 0.5E0 * TAUYZ / TKE
C
C     DU_k/DX_k = DU/DX + DV/DY + DW/DZ
C
            DUKDXK = DQDX(I,2) + DQDY(I,3) + DQDZ(I,4)
C
C      Wall damping
C
            ARG    = RHO * SQRT(TKE) * YNORML / RMU * RE / FSMACH
            FW     = EXP( -(0.015E0 * ARG)**4 )
            C1STAR = PSC1 * (1.E0 - (1.E0-1.E0/PSC1)*FW)
            CROTTA = 0.25E0 * C1STAR * EPSLN / TKE / TKE
C
C      Add implicit Jacobian for pressure-strain to source Jacobian.
C      DPI_ij = DPhi*_ij,1 + DPhi_ij,2 + DPhi*_ij,w
C
            ASRC(I, 6, 6) = ASRC(I, 6, 6) - CROTTA * (2.E0*TKE - TAUXX)
     1                      - (PSALFA*FW - CPHI2A) * T43 * DQDX(I,2)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (2.E0 * DQDX(I,2) - T23 * DUKDXK)
     4                      - (PSBETA*FW - CPHI2C) * T43 * DQDX(I,2)
            ASRC(I, 6, 7) = ASRC(I, 6, 7) + CROTTA * TAUXX
     1                      + (PSALFA*FW - CPHI2A) * T23 * DQDY(I,3)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (2.E0 * DQDX(I,2) - T23 * DUKDXK)
     4                      + (PSBETA*FW - CPHI2C) * T23 * DQDY(I,3)
            ASRC(I, 6, 8) = ASRC(I, 6, 8) + CROTTA * TAUXX
     1                      + (PSALFA*FW - CPHI2A) * T23 * DQDZ(I,4)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (2.E0 * DQDX(I,2) - T23 * DUKDXK)
     4                      + (PSBETA*FW - CPHI2C) * T23 * DQDZ(I,4)
            ASRC(I, 6, 9) = ASRC(I, 6, 9)
     1                      + (PSALFA*FW - CPHI2A) * T23 * 
     2                                 (DQDX(I,3) - 2.E0 * DQDY(I,2))
     3                      + (PSBETA*FW - CPHI2C) * T23 * 
     4                                 (DQDY(I,2) - 2.E0 * DQDX(I,3))
            ASRC(I, 6,10) = ASRC(I, 6,10)
     1                      + (PSALFA*FW - CPHI2A) * T23 * 
     2                                 (DQDX(I,4) - 2.E0 * DQDZ(I,2))
     3                      + (PSBETA*FW - CPHI2C) * T23 * 
     4                                 (DQDZ(I,2) - 2.E0 * DQDX(I,4))
            ASRC(I, 6,11) = ASRC(I, 6,11)
     1                      + (PSALFA*FW - CPHI2A) * T23 * 
     2                                 (DQDZ(I,3) + DQDY(I,4))
     3                      + (PSBETA*FW - CPHI2C) * T23 * 
     4                                 (DQDZ(I,3) + DQDY(I,4))
            ASRC(I, 6,12) = ASRC(I, 6,12) - C1STAR * BXX
C
            ASRC(I, 7, 6) = ASRC(I, 7, 6) + CROTTA * TAUYY
     1                      + (PSALFA*FW - CPHI2A) * T23 * DQDX(I,2)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (2.E0 * DQDY(I,3) - T23 * DUKDXK)
     4                      + (PSBETA*FW - CPHI2C) * T23 * DQDX(I,2)
            ASRC(I, 7, 7) = ASRC(I, 7, 7) - CROTTA * (2.E0*TKE - TAUYY)
     1                      - (PSALFA*FW - CPHI2A) * T43 * DQDY(I,3)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (2.E0 * DQDY(I,3) - T23 * DUKDXK)
     4                      - (PSBETA*FW - CPHI2C) * T43 * DQDY(I,3)
            ASRC(I, 7, 8) = ASRC(I, 7, 8) + CROTTA * TAUYY
     1                      + (PSALFA*FW - CPHI2A) * T23 * DQDZ(I,4)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (2.E0 * DQDY(I,3) - T23 * DUKDXK)
     4                      + (PSBETA*FW - CPHI2C) * T23 * DQDZ(I,4)
            ASRC(I, 7, 9) = ASRC(I, 7, 9)
     1                      + (PSALFA*FW - CPHI2A) * T23 * 
     2                                 (DQDY(I,2) - 2.E0 * DQDX(I,3))
     3                      + (PSBETA*FW - CPHI2C) * T23 * 
     4                                 (DQDX(I,3) - 2.E0 * DQDY(I,2))
            ASRC(I, 7,10) = ASRC(I, 7,10)
     1                      + (PSALFA*FW - CPHI2A) * T23 * 
     2                                 (DQDZ(I,2) + DQDX(I,4))
     3                      + (PSBETA*FW - CPHI2C) * T23 * 
     4                                 (DQDZ(I,2) + DQDX(I,4))
            ASRC(I, 7,11) = ASRC(I, 7,11)
     1                      + (PSALFA*FW - CPHI2A) * T23 * 
     2                                 (DQDY(I,4) - 2.E0 * DQDZ(I,3))
     3                      + (PSBETA*FW - CPHI2C) * T23 * 
     4                                 (DQDZ(I,3) - 2.E0 * DQDY(I,4))
            ASRC(I, 7,12) = ASRC(I, 7,12) - C1STAR * BYY
C
            ASRC(I, 8, 6) = ASRC(I, 8, 6) + CROTTA * TAUZZ
     1                      + (PSALFA*FW - CPHI2A) * T23 * DQDX(I,2)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (2.E0 * DQDZ(I,4) - T23 * DUKDXK)
     4                      + (PSBETA*FW - CPHI2C) * T23 * DQDX(I,2)
            ASRC(I, 8, 7) = ASRC(I, 8, 7) + CROTTA * TAUZZ
     1                      + (PSALFA*FW - CPHI2A) * T23 * DQDY(I,3)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (2.E0 * DQDZ(I,4) - T23 * DUKDXK)
     4                      + (PSBETA*FW - CPHI2C) * T23 * DQDY(I,3)
            ASRC(I, 8, 8) = ASRC(I, 8, 8) - CROTTA * (2.E0*TKE - TAUZZ)
     1                      - (PSALFA*FW - CPHI2A) * T43 * DQDZ(I,4)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (2.E0 * DQDZ(I,4) - T23 * DUKDXK)
     4                      - (PSBETA*FW - CPHI2C) * T43 * DQDZ(I,4)
            ASRC(I, 8, 9) = ASRC(I, 8, 9)
     1                      + (PSALFA*FW - CPHI2A) * T23 * 
     2                                 (DQDY(I,2) + DQDX(I,3))
     3                      + (PSBETA*FW - CPHI2C) * T23 * 
     4                                 (DQDX(I,3) + DQDY(I,2))
            ASRC(I, 8,10) = ASRC(I, 8,10)
     1                      + (PSALFA*FW - CPHI2A) * T23 * 
     2                                 (DQDZ(I,2) - 2.E0 * DQDX(I,4))
     3                      + (PSBETA*FW - CPHI2C) * T23 * 
     4                                 (DQDX(I,4) - 2.E0 * DQDZ(I,2))
            ASRC(I, 8,11) = ASRC(I, 8,11)
     1                      + (PSALFA*FW - CPHI2A) * T23 * 
     2                                 (DQDZ(I,3) - 2.E0 * DQDY(I,4))
     3                      + (PSBETA*FW - CPHI2C) * T23 * 
     4                                 (DQDY(I,4) - 2.E0 * DQDZ(I,3))
            ASRC(I, 8,12) = ASRC(I, 8,12) - C1STAR * BZZ
C
            ASRC(I, 9, 6) = ASRC(I, 9, 6) + CROTTA * TAUXY
     1                      - (PSALFA*FW - CPHI2A) * DQDX(I,3)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (DQDY(I,2) + DQDX(I,3))
     4                      - (PSBETA*FW - CPHI2C) * DQDY(I,2)
            ASRC(I, 9, 7) = ASRC(I, 9, 7) + CROTTA * TAUXY
     1                      - (PSALFA*FW - CPHI2A) * DQDY(I,2)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (DQDY(I,2) + DQDX(I,3))
     4                      - (PSBETA*FW - CPHI2C) * DQDX(I,3)
            ASRC(I, 9, 8) = ASRC(I, 9, 8) + CROTTA * TAUXY
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (DQDY(I,2) + DQDX(I,3))
            ASRC(I, 9, 9) = ASRC(I, 9, 9) - 0.5E0 * C1STAR / TKE / TKE
     1                      - (PSALFA*FW - CPHI2A) * 
     2                                 (DQDY(I,3) + DQDX(I,2))
     3                      - (PSBETA*FW - CPHI2C) * 
     4                                 (DQDY(I,3) + DQDX(I,2))
            ASRC(I, 9,10) = ASRC(I, 9,10)
     1                      + (PSALFA*FW - CPHI2A) * DQDZ(I,3)
     2                      + (PSBETA*FW - CPHI2C) * DQDY(I,4)
            ASRC(I, 9,11) = ASRC(I, 9,11)
     1                      + (PSALFA*FW - CPHI2A) * DQDZ(I,2)
     2                      + (PSBETA*FW - CPHI2C) * DQDX(I,4)
            ASRC(I, 9,12) = ASRC(I, 9,12) - C1STAR * BXY
C
            ASRC(I,10, 6) = ASRC(I,10, 6) + CROTTA * TAUXZ
     1                      - (PSALFA*FW - CPHI2A) * DQDX(I,4)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (DQDZ(I,2) + DQDX(I,4))
     4                      - (PSBETA*FW - CPHI2C) * DQDZ(I,2)
            ASRC(I,10, 7) = ASRC(I,10, 7) + CROTTA * TAUXZ
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (DQDZ(I,2) + DQDX(I,4))
            ASRC(I,10, 8) = ASRC(I,10, 8) + CROTTA * TAUXZ
     1                      - (PSALFA*FW - CPHI2A) * DQDZ(I,2)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (DQDZ(I,2) + DQDX(I,4))
     4                      - (PSBETA*FW - CPHI2C) * DQDX(I,4)
            ASRC(I,10, 9) = ASRC(I,10, 9)
     1                      + (PSALFA*FW - CPHI2A) * DQDY(I,4)
     2                      + (PSBETA*FW - CPHI2C) * DQDZ(I,3)
            ASRC(I,10,10) = ASRC(I,10,10) - 0.5E0 * C1STAR / TKE / TKE
     1                      - (PSALFA*FW - CPHI2A) * 
     2                                 (DQDZ(I,4) + DQDX(I,2))
     3                      - (PSBETA*FW - CPHI2C) * 
     4                                 (DQDZ(I,4) + DQDX(I,2))
            ASRC(I,10,11) = ASRC(I,10,11)
     1                      + (PSALFA*FW - CPHI2A) * DQDY(I,2)
     2                      + (PSBETA*FW - CPHI2C) * DQDX(I,3)
            ASRC(I,10,12) = ASRC(I,10,12) - C1STAR * BXZ
C
            ASRC(I,11, 6) = ASRC(I,11, 6) + CROTTA * TAUYZ
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (DQDZ(I,3) + DQDY(I,4))
            ASRC(I,11, 7) = ASRC(I,11, 7) + CROTTA * TAUYZ
     1                      - (PSALFA*FW - CPHI2A) * DQDY(I,4)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (DQDZ(I,3) + DQDY(I,4))
     4                      - (PSBETA*FW - CPHI2C) * DQDZ(I,3)
            ASRC(I,11, 8) = ASRC(I,11, 8) + CROTTA * TAUYZ
     1                      - (PSALFA*FW - CPHI2A) * DQDZ(I,3)
     2                      + 0.5E0*(PSGAMA*FW - CPHI2B) * 
     3                                 (DQDZ(I,3) + DQDY(I,4))
     4                      - (PSBETA*FW - CPHI2C) * DQDY(I,4)
            ASRC(I,11, 9) = ASRC(I,11, 9)
     1                      + (PSALFA*FW - CPHI2A) * DQDX(I,4)
     2                      + (PSBETA*FW - CPHI2C) * DQDZ(I,2)
            ASRC(I,11,10) = ASRC(I,11,10)
     1                      + (PSALFA*FW - CPHI2A) * DQDX(I,3)
     2                      + (PSBETA*FW - CPHI2C) * DQDY(I,2)
            ASRC(I,11,11) = ASRC(I,11,11) - 0.5E0 * C1STAR / TKE / TKE
     1                      - (PSALFA*FW - CPHI2A) * 
     2                                 (DQDZ(I,4) + DQDY(I,3))
     3                      - (PSBETA*FW - CPHI2C) * 
     4                                 (DQDZ(I,4) + DQDY(I,3))
            ASRC(I,11,12) = ASRC(I,11,12) - C1STAR * BYZ
C
  200    CONTINUE
C
C     Error
C
      ELSE
         WRITE (IOUT,1000) IPSMDL
         STOP
      ENDIF
C
C     Finished calculating pressure-strain correlation Jacobians
C
      RETURN
 1000 FORMAT (' ','ERROR - RUN ABORTING',
     1       /' ','Incorrect Pressure-Strain Correlation ',
     2            'Specified in routine DPSTRN = ',I10)
      END







      SUBROUTINE TLTROE (ITURB, NPTS, QP, QM, S, TKEP, TKEM, 
     1                   HP, HM, A)
C
C     Routine to calculate T |Lambda| T(inverse) using
C     Roe averaged variables at the interface for the mean equations
C     with either a two equation turbulence model of the form k-psi
C     or a Reynolds stress closure.
C
C
C     ITURB     : Turbulence model
C     NPTS      : Number of points to calculate T*LAMBDA*T(-1)
C     QP,QM     : Q(+),Q(-) at a cell interface for the MUSCL scheme
C     S         : Metrics
C     TKEP,TKEM : Turbulent kinetic energy for QP,QM
C     HP,HM     : Total enthalpy for QP,QM
C     A         : The T*LAMBDA*T(-1) matrix
C
C     The two-equation and Reynolds stress models require the turbulent
C     kinetic energy in the mean flow equations.  Therefore, this routine
C     is organized into the following sections:
C
C          1. Calculate the Turbulent Kinetic Energy (TKE)
C             a. Inviscid, laminar, algebraic model => TKE = 0
C             b. Two-Equation model                 => TKE = Q(6)
C             c. Reynolds stress model              => TKE = (Q(7)+Q(8)+Q(9))/2
C          2. Calculate Total Enthalpy
C             a. Perfect Gas   => H = gamma P / (RHO(gamma-1)) + (q^2)/2 + TKE
C          3. Calculate Flux Jacobians for the Mean Flow Equations
C          4. Calculate Additional Flux Jacobians
C             a. Inviscid, laminar, algebraic model => No additional Jacobians
C             b. Two-Equation model                 => Jacobians for TKE, phi
C             c. Reynolds stress model              => Jacobians for Tau_ij,phi
C
      include '../header/common.h'
C
      DIMENSION QP    (NPTS,NQ),
     1          QM    (NPTS,NQ),
     2          S     (NPTS,4),
     3          TKEP  (NPTS),
     4          TKEM  (NPTS),
     5          HP    (NPTS),
     6          HM    (NPTS),
     7          A     (NPTS,NRANK,NRANK)
C
C     1. Calculate Turbulent Kinetic Energy (TKEP,TKEM)
C
C     1-b. Turbulent Kinetic Energy for Two-Equation models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 100 I = 1, NPTS
            TKEP(I) = QP(I,6)
            TKEM(I) = QM(I,6)
  100    CONTINUE
C
C     1-c. Turbulent Kinetic Energy for Reynolds Stress models
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 110 I = 1, NPTS
            TKEP(I) = 0.5E0 * (QP(I,6) + QP(I,7) + QP(I,8))
            TKEM(I) = 0.5E0 * (QM(I,6) + QM(I,7) + QM(I,8))
  110    CONTINUE
C
C     1-a. Turbulent Kinetic Energy = 0 for inviscid, laminar, algebraic models
C
      ELSE
         DO 120 I = 1, NPTS
            TKEP(I) = 0.E0
            TKEM(I) = 0.E0
  120    CONTINUE
      ENDIF
C
C     2. Calculate Total Enthalpy
C
C     2-a. Perfect Gas
C
      DO 200 I = 1, NPTS
         RHOP  = QP(I,1)
         UP    = QP(I,2)
         VP    = QP(I,3)
         WP    = QP(I,4)
         PP    = QP(I,5)
         HP(I) = GAMMA * PP / (RHOP * GAMM1)
     1          + 0.5E0 * (UP*UP + VP*VP + WP*WP) + TKEP(I)
C
         RHOM  = QM(I,1)
         UM    = QM(I,2)
         VM    = QM(I,3)
         WM    = QM(I,4)
         PM    = QM(I,5)
         HM(I) = GAMMA * PM / (RHOM * GAMM1)
     1          + 0.5E0 * (UM*UM + VM*VM + WM*WM) + TKEM(I)
  200 CONTINUE
C
C     3. Calculate T |Lambda| T(inverse) at the cell interfaces
C        for Mean Flow Equations
C
      DO 300 I = 1, NPTS
C
C     Calculate +
C
         RHOP   = QP(I,1)
         UP     = QP(I,2) 
         VP     = QP(I,3)
         WP     = QP(I,4)
         PP     = QP(I,5)
C
         UBARP  = S(I,1)*UP + S(I,2)*VP + S(I,3)*WP
         ASQP   = GAMMA * PP / RHOP
C
C     Calculate -
C
         RHOM   = QM(I,1)
         UM     = QM(I,2) 
         VM     = QM(I,3)
         WM     = QM(I,4)
         PM     = QM(I,5)
C
         UBARM  = S(I,1)*UM + S(I,2)*VM + S(I,3)*WM
         ASQM   = GAMMA * PM / RHOM
C
C     Roe averaging
C
         RHOMSR = SQRT (RHOM)
         RHOPSR = SQRT (RHOP)
         RHOAV  = RHOMSR + RHOPSR
         UR     = ( RHOMSR*UM + RHOPSR*UP ) / RHOAV
         VR     = ( RHOMSR*VM + RHOPSR*VP ) / RHOAV
         WR     = ( RHOMSR*WM + RHOPSR*WP ) / RHOAV
         UBARR  = S(I,1)*UR + S(I,2)*VR + S(I,3)*WR
         HR     = ( RHOMSR*HM(I) + RHOPSR*HP(I) ) / RHOAV
C
         TKER   = ( RHOMSR*TKEM(I) + RHOPSR*TKEP(I) ) / RHOAV
C
         QSQR   = UR*UR + VR*VR + WR*WR
         ASQR   = GAMM1 * (HR - 0.5E0*QSQR - TKER)
         AR     = SQRT (ASQR)
         RHOR   = RHOMSR*RHOPSR
C
C     Calculate the eigenvalues.
C
         ALAM1  = ABS (UBARR)
         ALAM4  = ABS (UBARR + AR)
         ALAM5  = ABS (UBARR - AR)
C
C     Add the contribution due to the fourth eigenvalue.
C
C     TLT41  - contribution to T |Lambda| T(-1) from 4th eigenvalue to
C              the first column
C     TEMP42 - factor to multiply TLT4 by on the 2nd row
C
         TLT41  =   ALAM4 * (0.25E0*GAMM1*QSQR/ASQR - 0.5E0*UBARR/AR)
         TLT42  =   ALAM4 * 0.5E0 * (S(I,1)/AR - GAMM1*UR/ASQR)
         TLT43  =   ALAM4 * 0.5E0 * (S(I,2)/AR - GAMM1*VR/ASQR)
         TLT44  =   ALAM4 * 0.5E0 * (S(I,3)/AR - GAMM1*WR/ASQR)
         TLT45  =   ALAM4 * 0.5E0 * GAMM1 / ASQR
C
         TEMP42 =   UR + AR*S(I,1)
         TEMP43 =   VR + AR*S(I,2)
         TEMP44 =   WR + AR*S(I,3)
         TEMP45 =   HR + AR*UBARR
C
C     Add the contribution due to the fifth eigenvalue.
C
         TLT51  =   ALAM5 * (0.25E0*GAMM1*QSQR/ASQR + 0.5E0*UBARR/AR)
         TLT52  = - ALAM5 * 0.5E0 * (S(I,1)/AR + GAMM1*UR/ASQR)
         TLT53  = - ALAM5 * 0.5E0 * (S(I,2)/AR + GAMM1*VR/ASQR)
         TLT54  = - ALAM5 * 0.5E0 * (S(I,3)/AR + GAMM1*WR/ASQR)
         TLT55  =   ALAM5 * 0.5E0 * GAMM1 / ASQR
C
         TEMP52 =   UR - AR*S(I,1)
         TEMP53 =   VR - AR*S(I,2)
         TEMP54 =   WR - AR*S(I,3)
         TEMP55 =   HR - AR*UBARR
C
         A(I,1,1) = ( (1.E0 - 0.5E0 * GAMM1 * QSQR / ASQR) * ALAM1
     1                + TLT41 + TLT51 ) * S(I,4)
         A(I,1,2) = ( GAMM1 * UR / ASQR * ALAM1
     1                + TLT42 + TLT52 ) * S(I,4)
         A(I,1,3) = ( GAMM1 * VR / ASQR * ALAM1
     1                + TLT43 + TLT53 ) * S(I,4)
         A(I,1,4) = ( GAMM1 * WR / ASQR * ALAM1
     1                + TLT44 + TLT54 ) * S(I,4)
         A(I,1,5) = ( - GAMM1 / ASQR * ALAM1
     1                + TLT45 + TLT55 ) * S(I,4)
C
         A(I,2,1) = ( (S(I,1)*UBARR - 0.5E0*GAMM1*UR*QSQR/ASQR) * ALAM1
     1                + TLT41*TEMP42 + TLT51*TEMP52 ) * S(I,4)
         A(I,2,2) = ( (GAMM1*UR*UR/ASQR - S(I,1)*S(I,1) + 1.E0) * ALAM1
     1                + TLT42*TEMP42 + TLT52*TEMP52 ) * S(I,4)
         A(I,2,3) = ( (GAMM1*UR*VR/ASQR - S(I,1)*S(I,2)) * ALAM1
     1                + TLT43*TEMP42 + TLT53*TEMP52 ) * S(I,4)
         A(I,2,4) = ( (GAMM1*UR*WR/ASQR - S(I,1)*S(I,3)) * ALAM1
     1                + TLT44*TEMP42 + TLT54*TEMP52 ) * S(I,4)
         A(I,2,5) = ( - GAMM1 * UR / ASQR * ALAM1
     1                + TLT45*TEMP42 + TLT55*TEMP52 ) * S(I,4)
C
         A(I,3,1) = ( (S(I,2)*UBARR - 0.5E0*GAMM1*VR*QSQR/ASQR) * ALAM1
     1                + TLT41*TEMP43 + TLT51*TEMP53 ) * S(I,4)
         A(I,3,2) = ( (GAMM1*VR*UR/ASQR - S(I,2)*S(I,1)) * ALAM1
     1                + TLT42*TEMP43 + TLT52*TEMP53 ) * S(I,4)
         A(I,3,3) = ( (GAMM1*VR*VR/ASQR - S(I,2)*S(I,2) + 1.E0) * ALAM1
     1                + TLT43*TEMP43 + TLT53*TEMP53 ) * S(I,4)
         A(I,3,4) = ( (GAMM1*VR*WR/ASQR - S(I,2)*S(I,3)) * ALAM1
     1                + TLT44*TEMP43 + TLT54*TEMP53 ) * S(I,4)
         A(I,3,5) = ( - GAMM1 * VR / ASQR * ALAM1
     1                + TLT45*TEMP43 + TLT55*TEMP53 ) * S(I,4)
C
         A(I,4,1) = ( (S(I,3)*UBARR - 0.5E0*GAMM1*WR*QSQR/ASQR) * ALAM1
     1                + TLT41*TEMP44 + TLT51*TEMP54 ) * S(I,4)
         A(I,4,2) = ( (GAMM1*WR*UR/ASQR - S(I,3)*S(I,1)) * ALAM1
     1                + TLT42*TEMP44 + TLT52*TEMP54 ) * S(I,4)
         A(I,4,3) = ( (GAMM1*WR*VR/ASQR - S(I,3)*S(I,2)) * ALAM1
     1                + TLT43*TEMP44 + TLT53*TEMP54 ) * S(I,4)
         A(I,4,4) = ( (GAMM1*WR*WR/ASQR - S(I,3)*S(I,3) + 1.E0) * ALAM1
     1                + TLT44*TEMP44 + TLT54*TEMP54 ) * S(I,4)
         A(I,4,5) = ( - GAMM1 * WR / ASQR * ALAM1
     1                + TLT45*TEMP44 + TLT55*TEMP54 ) * S(I,4)
C
         A(I,5,1) = ( (UBARR*UBARR - 0.50E0*QSQR 
     1                             - 0.5E0*GAMM1*QSQR*TKER/ASQR
     2                            - 0.25E0*GAMM1*QSQR*QSQR/ASQR)*ALAM1
     3                + TLT41*TEMP45 + TLT51*TEMP55 ) * S(I,4)
         A(I,5,2) = ( (0.5E0*GAMM1*QSQR*UR/ASQR + UR 
     1                 - S(I,1)*UBARR + GAMM1*UR*TKER/ASQR) * ALAM1
     2                + TLT42*TEMP45 + TLT52*TEMP55 ) * S(I,4)
         A(I,5,3) = ( (0.5E0*GAMM1*QSQR*VR/ASQR + VR 
     1                 - S(I,2)*UBARR + GAMM1*VR*TKER/ASQR) * ALAM1
     2                + TLT43*TEMP45 + TLT53*TEMP55 ) * S(I,4)
         A(I,5,4) = ( (0.5E0*GAMM1*QSQR*WR/ASQR + WR
     1                 - S(I,3)*UBARR + GAMM1*WR*TKER/ASQR) * ALAM1
     2                + TLT44*TEMP45 + TLT54*TEMP55 ) * S(I,4)
         A(I,5,5) = ( - (0.5E0*GAMM1*QSQR/ASQR + GAMM1*TKER/ASQR)*ALAM1
     1                + TLT45*TEMP45 + TLT55*TEMP55 ) * S(I,4)
  300 CONTINUE
C
C     4. Additional Jacobians
C
C     4-b. Two-Equation Models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 410 I = 1, NPTS
C
C     Calculate +
C
            RHOP   = QP(I,1)
            UP     = QP(I,2) 
            VP     = QP(I,3)
            WP     = QP(I,4)
            PP     = QP(I,5)
C
            T2P    = QP(I,7)
C
            UBARP  = S(I,1)*UP + S(I,2)*VP + S(I,3)*WP
            ASQP   = GAMMA * PP / RHOP
C
C     Calculate -
C
            RHOM   = QM(I,1)
            UM     = QM(I,2) 
            VM     = QM(I,3)
            WM     = QM(I,4)
            PM     = QM(I,5)
C
            T2M    = QM(I,7)
C
            UBARM  = S(I,1)*UM + S(I,2)*VM + S(I,3)*WM
            ASQM   = GAMMA * PM / RHOM
C
C     Roe averaging
C
            RHOMSR = SQRT (RHOM)
            RHOPSR = SQRT (RHOP)
            RHOAV  = RHOMSR + RHOPSR
            UR     = ( RHOMSR*UM + RHOPSR*UP ) / RHOAV
            VR     = ( RHOMSR*VM + RHOPSR*VP ) / RHOAV
            WR     = ( RHOMSR*WM + RHOPSR*WP ) / RHOAV
            UBARR  = S(I,1)*UR + S(I,2)*VR + S(I,3)*WR
            HR     = ( RHOMSR*HM(I) + RHOPSR*HP(I) ) / RHOAV
C
            TKER   = ( RHOMSR*TKEM(I) + RHOPSR*TKEP(I) ) / RHOAV
            T2R    = ( RHOMSR*T2M + RHOPSR*T2P ) / RHOAV
C
            QSQR   = UR*UR + VR*VR + WR*WR
            ASQR   = GAMM1 * (HR - 0.5E0*QSQR - TKER)
            AR     = SQRT (ASQR)
            RHOR   = RHOMSR*RHOPSR
C
C     Calculate the eigenvalues.
C
            ALAM1  = ABS (UBARR)
            ALAM4  = ABS (UBARR + AR)
            ALAM5  = ABS (UBARR - AR)
C
C     Add the contribution due to the fourth eigenvalue.
C
C     TLT41  - contribution to T |Lambda| T(-1) from 4th eigenvalue to
C              the first column
C     TEMP42 - factor to multiply TLT4 by on the 2nd row
C
            TLT41  =   ALAM4 * (0.25E0*GAMM1*QSQR/ASQR - 0.5E0*UBARR/AR)
            TLT42  =   ALAM4 * 0.5E0 * (S(I,1)/AR - GAMM1*UR/ASQR)
            TLT43  =   ALAM4 * 0.5E0 * (S(I,2)/AR - GAMM1*VR/ASQR)
            TLT44  =   ALAM4 * 0.5E0 * (S(I,3)/AR - GAMM1*WR/ASQR)
            TLT45  =   ALAM4 * 0.5E0 * GAMM1 / ASQR
            TLT46  = - ALAM4 * 0.5E0 * GAMM1 / ASQR
C
            TEMP42 =   UR + AR*S(I,1)
            TEMP43 =   VR + AR*S(I,2)
            TEMP44 =   WR + AR*S(I,3)
            TEMP45 =   HR + AR*UBARR
            TEMP46 =   TKER
            TEMP47 =   T2R
C
C     Add the contribution due to the fifth eigenvalue.
C
            TLT51  =   ALAM5 * (0.25E0*GAMM1*QSQR/ASQR + 0.5E0*UBARR/AR)
            TLT52  = - ALAM5 * 0.5E0 * (S(I,1)/AR + GAMM1*UR/ASQR)
            TLT53  = - ALAM5 * 0.5E0 * (S(I,2)/AR + GAMM1*VR/ASQR)
            TLT54  = - ALAM5 * 0.5E0 * (S(I,3)/AR + GAMM1*WR/ASQR)
            TLT55  =   ALAM5 * 0.5E0 * GAMM1 / ASQR
            TLT56  = - ALAM5 * 0.5E0 * GAMM1 / ASQR
C
            TEMP52 =   UR - AR*S(I,1)
            TEMP53 =   VR - AR*S(I,2)
            TEMP54 =   WR - AR*S(I,3)
            TEMP55 =   HR - AR*UBARR
            TEMP56 =   TKER
            TEMP57 =   T2R
C
            A(I,1,6) = ( GAMM1 / ASQR * ALAM1
     1                   + TLT46 + TLT56 ) * S(I,4)
            A(I,1,7) = 0.0E0
C
            A(I,2,6) = ( GAMM1 * UR / ASQR * ALAM1
     1                   + TLT46*TEMP42 + TLT56*TEMP52 ) * S(I,4)
            A(I,2,7) = 0.0E0
C
            A(I,3,6) = ( GAMM1 * VR / ASQR * ALAM1
     1                   + TLT46*TEMP43 + TLT56*TEMP53 ) * S(I,4)
            A(I,3,7) = 0.0E0
C
            A(I,4,6) = ( GAMM1 * WR / ASQR * ALAM1
     1                   + TLT46*TEMP44 + TLT56*TEMP54 ) * S(I,4)
            A(I,4,7) = 0.0E0
C
            A(I,5,6) = ( (1.E0+0.5E0*GAMM1*QSQR/ASQR+GAMM1*TKER/ASQR)* 
     1                   ALAM1 + TLT45*TEMP45 + TLT55*TEMP55 ) * S(I,4)
            A(I,5,7) = 0.0E0
C
            A(I,6,1) = ( ( - 0.5E0 * GAMM1 * QSQR * TKER / ASQR) * ALAM1
     1                   + TLT41*TEMP46 + TLT51*TEMP56 ) * S(I,4)
            A(I,6,2) = ( GAMM1 * UR * TKER / ASQR * ALAM1
     1                   + TLT42*TEMP46 + TLT52*TEMP56 ) * S(I,4)
            A(I,6,3) = ( GAMM1 * VR * TKER / ASQR * ALAM1
     1                   + TLT43*TEMP46 + TLT53*TEMP56 ) * S(I,4)
            A(I,6,4) = ( GAMM1 * WR * TKER / ASQR * ALAM1
     1                   + TLT44*TEMP46 + TLT54*TEMP56 ) * S(I,4)
            A(I,6,5) = ( - GAMM1 * TKER / ASQR * ALAM1
     1                   + TLT45*TEMP46 + TLT55*TEMP56 ) * S(I,4)
            A(I,6,6) = ( (1.E0 + GAMM1 * TKER / ASQR) * ALAM1
     1                   + TLT46*TEMP46 + TLT56*TEMP56 ) * S(I,4)
            A(I,6,7) = 0.0E0
C
            A(I,7,1) = ( ( - 0.5E0 * GAMM1 * QSQR * T2R / ASQR) * ALAM1
     1                   + TLT41*TEMP47 + TLT51*TEMP57 ) * S(I,4)
            A(I,7,2) = ( GAMM1 * UR * T2R / ASQR * ALAM1
     1                   + TLT42*TEMP47 + TLT52*TEMP57 ) * S(I,4)
            A(I,7,3) = ( GAMM1 * VR * T2R / ASQR * ALAM1
     1                   + TLT43*TEMP47 + TLT53*TEMP57 ) * S(I,4)
            A(I,7,4) = ( GAMM1 * WR * T2R / ASQR * ALAM1
     1                   + TLT44*TEMP47 + TLT54*TEMP57 ) * S(I,4)
            A(I,7,5) = ( - GAMM1 * T2R / ASQR * ALAM1
     1                   + TLT45*TEMP47 + TLT55*TEMP57 ) * S(I,4)
            A(I,7,6) = ( GAMM1 * T2R / ASQR * ALAM1
     1                   + TLT46*TEMP47 + TLT56*TEMP57 ) * S(I,4)
            A(I,7,7) = ALAM1 * S(I,4)
  410    CONTINUE
C
C     4-c. Additional Jacobians for Reynolds Stress models
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 420 I = 1, NPTS
C
C     Calculate +
C
            RHOP   = QP(I, 1)
            UP     = QP(I, 2) 
            VP     = QP(I, 3)
            WP     = QP(I, 4)
            PP     = QP(I, 5)
C
            TAUXXP = QP(I, 6)
            TAUYYP = QP(I, 7)
            TAUZZP = QP(I, 8)
            TAUXYP = QP(I, 9)
            TAUXZP = QP(I,10)
            TAUYZP = QP(I,11)
            T2P    = QP(I,12)
C
            UBARP  = S(I,1)*UP + S(I,2)*VP + S(I,3)*WP
            ASQP   = GAMMA * PP / RHOP
C
C     Calculate -
C
            RHOM   = QM(I, 1)
            UM     = QM(I, 2) 
            VM     = QM(I, 3)
            WM     = QM(I, 4)
            PM     = QM(I, 5)
C
            TAUXXM = QM(I, 6)
            TAUYYM = QM(I, 7)
            TAUZZM = QM(I, 8)
            TAUXYM = QM(I, 9)
            TAUXZM = QM(I,10)
            TAUYZM = QM(I,11)
            T2M    = QM(I,12)
C
            UBARM  = S(I,1)*UM + S(I,2)*VM + S(I,3)*WM
            ASQM   = GAMMA * PM / RHOM
C
C     Roe averaging
C
            RHOMSR = SQRT (RHOM)
            RHOPSR = SQRT (RHOP)
            RHOAV  = RHOMSR + RHOPSR
            UR     = ( RHOMSR*UM + RHOPSR*UP ) / RHOAV
            VR     = ( RHOMSR*VM + RHOPSR*VP ) / RHOAV
            WR     = ( RHOMSR*WM + RHOPSR*WP ) / RHOAV
            UBARR  = S(I,1)*UR + S(I,2)*VR + S(I,3)*WR
            HR     = ( RHOMSR*HM(I) + RHOPSR*HP(I) ) / RHOAV
C
            TKER   = ( RHOMSR*TKEM(I) + RHOPSR*TKEP(I) ) / RHOAV
            TAUXXR = ( RHOMSR*TAUXXM  + RHOPSR*TAUXXP  ) / RHOAV
            TAUYYR = ( RHOMSR*TAUYYM  + RHOPSR*TAUYYP  ) / RHOAV
            TAUZZR = ( RHOMSR*TAUZZM  + RHOPSR*TAUZZP  ) / RHOAV
            TAUXYR = ( RHOMSR*TAUXYM  + RHOPSR*TAUXYP  ) / RHOAV
            TAUXZR = ( RHOMSR*TAUXZM  + RHOPSR*TAUXZP  ) / RHOAV
            TAUYZR = ( RHOMSR*TAUYZM  + RHOPSR*TAUYZP  ) / RHOAV
            T2R    = ( RHOMSR*T2M     + RHOPSR*T2P     ) / RHOAV
C
            QSQR   = UR*UR + VR*VR + WR*WR
            ASQR   = GAMM1 * (HR - 0.5E0*QSQR - TKER)
            AR     = SQRT (ASQR)
            RHOR   = RHOMSR*RHOPSR
C
C     Calculate the eigenvalues.
C
            ALAM1  = ABS (UBARR)
            ALAM4  = ABS (UBARR + AR)
            ALAM5  = ABS (UBARR - AR)
C
C     Add the contribution due to the fourth eigenvalue.
C
C     TLT41  - contribution to T |Lambda| T(-1) from 4th eigenvalue to
C              the first column
C     TEMP42 - factor to multiply TLT4 by on the 2nd row
C
            TLT41  =   ALAM4 * (0.25E0*GAMM1*QSQR/ASQR - 0.5E0*UBARR/AR)
            TLT42  =   ALAM4 * 0.50E0 * (S(I,1)/AR - GAMM1*UR/ASQR)
            TLT43  =   ALAM4 * 0.50E0 * (S(I,2)/AR - GAMM1*VR/ASQR)
            TLT44  =   ALAM4 * 0.50E0 * (S(I,3)/AR - GAMM1*WR/ASQR)
            TLT45  =   ALAM4 * 0.50E0 * GAMM1 / ASQR
            TLT46  = - ALAM4 * 0.25E0 * GAMM1 / ASQR
C
            TEMP42 =   UR + AR*S(I,1)
            TEMP43 =   VR + AR*S(I,2)
            TEMP44 =   WR + AR*S(I,3)
            TEMP45 =   HR + AR*UBARR
C
C     Add the contribution due to the fifth eigenvalue.
C
            TLT51  =   ALAM5 * (0.25E0*GAMM1*QSQR/ASQR + 0.5E0*UBARR/AR)
            TLT52  = - ALAM5 * 0.50E0 * (S(I,1)/AR + GAMM1*UR/ASQR)
            TLT53  = - ALAM5 * 0.50E0 * (S(I,2)/AR + GAMM1*VR/ASQR)
            TLT54  = - ALAM5 * 0.50E0 * (S(I,3)/AR + GAMM1*WR/ASQR)
            TLT55  =   ALAM5 * 0.50E0 * GAMM1 / ASQR
            TLT56  = - ALAM5 * 0.25E0 * GAMM1 / ASQR
C
            TEMP52 =   UR - AR*S(I,1)
            TEMP53 =   VR - AR*S(I,2)
            TEMP54 =   WR - AR*S(I,3)
            TEMP55 =   HR - AR*UBARR
C
            A(I, 1, 6) = ( 0.5E0 * GAMM1 / ASQR * ALAM1
     1                   + TLT46 + TLT56 ) * S(I,4)
            A(I, 1, 7) = A(I, 1, 6)
            A(I, 1, 8) = A(I, 1, 6)
            A(I, 1, 9) = 0.0E0
            A(I, 1,10) = 0.0E0
            A(I, 1,11) = 0.0E0
            A(I, 1,12) = 0.0E0
C
            A(I, 2, 6) = ( 0.5E0 * GAMM1 * UR / ASQR * ALAM1
     1                   + TLT46*TEMP42 + TLT56*TEMP52 ) * S(I,4)
            A(I, 2, 7) = A(I, 2, 6)
            A(I, 2, 8) = A(I, 2, 6)
            A(I, 2, 9) = 0.0E0
            A(I, 2,10) = 0.0E0
            A(I, 2,11) = 0.0E0
            A(I, 2,12) = 0.0E0
C
            A(I, 3, 6) = ( 0.5E0 * GAMM1 * VR / ASQR * ALAM1
     1                   + TLT46*TEMP43 + TLT56*TEMP53 ) * S(I,4)
            A(I, 3, 7) = A(I, 3, 6)
            A(I, 3, 8) = A(I, 3, 6)
            A(I, 3, 9) = 0.0E0
            A(I, 3,10) = 0.0E0
            A(I, 3,11) = 0.0E0
            A(I, 3,12) = 0.0E0
C
            A(I, 4, 6) = ( 0.5E0 * GAMM1 * WR / ASQR * ALAM1
     1                   + TLT46*TEMP44 + TLT56*TEMP54 ) * S(I,4)
            A(I, 4, 7) = A(I, 4, 6)
            A(I, 4, 8) = A(I, 4, 6)
            A(I, 4, 9) = 0.0E0
            A(I, 4,10) = 0.0E0
            A(I, 4,11) = 0.0E0
            A(I, 4,12) = 0.0E0
C
            A(I, 5, 6) = ( (0.5E0+0.25E0*GAMM1*QSQR/ASQR
     1                     +0.5E0*GAMM1*TKER/ASQR)* ALAM1
     2                   + TLT46*TEMP45 + TLT56*TEMP55 ) * S(I,4)
            A(I, 5, 7) = A(I, 5, 6)
            A(I, 5, 8) = A(I, 5, 6)
            A(I, 5, 9) = 0.0E0
            A(I, 5,10) = 0.0E0
            A(I, 5,11) = 0.0E0
            A(I, 5,12) = 0.0E0
C
            A(I, 6, 1) = ( ( - 0.5E0 * GAMM1 * QSQR * TAUXXR/ASQR)*ALAM1
     1                   + TLT41*TAUXXR + TLT51*TAUXXR ) * S(I,4)
            A(I, 6, 2) = ( GAMM1 * UR * TAUXXR / ASQR * ALAM1
     1                   + TLT42*TAUXXR + TLT52*TAUXXR ) * S(I,4)
            A(I, 6, 3) = ( GAMM1 * VR * TAUXXR / ASQR * ALAM1
     1                   + TLT43*TAUXXR + TLT53*TAUXXR ) * S(I,4)
            A(I, 6, 4) = ( GAMM1 * WR * TAUXXR / ASQR * ALAM1
     1                   + TLT44*TAUXXR + TLT54*TAUXXR ) * S(I,4)
            A(I, 6, 5) = ( - GAMM1 * TAUXXR / ASQR * ALAM1
     1                   + TLT45*TAUXXR + TLT55*TAUXXR ) * S(I,4)
            A(I, 6, 6) = ( (1.E0 + 0.5E0 * GAMM1 * TAUXXR / ASQR)*ALAM1
     1                   + TLT46*TAUXXR + TLT56*TAUXXR ) * S(I,4)
            A(I, 6, 7) = ( 0.5E0 * GAMM1 * TAUXXR / ASQR * ALAM1
     1                   + TLT46*TAUXXR + TLT56*TAUXXR ) * S(I,4)
            A(I, 6, 8) = ( 0.5E0 * GAMM1 * TAUXXR / ASQR * ALAM1
     1                   + TLT46*TAUXXR + TLT56*TAUXXR ) * S(I,4)
            A(I, 6, 9) = 0.0E0
            A(I, 6,10) = 0.0E0
            A(I, 6,11) = 0.0E0
            A(I, 6,12) = 0.0E0
C
            A(I, 7, 1) = ( ( - 0.5E0 * GAMM1 * QSQR * TAUYYR/ASQR)*ALAM1
     1                   + TLT41*TAUYYR + TLT51*TAUYYR ) * S(I,4)
            A(I, 7, 2) = ( GAMM1 * UR * TAUYYR / ASQR * ALAM1
     1                   + TLT42*TAUYYR + TLT52*TAUYYR ) * S(I,4)
            A(I, 7, 3) = ( GAMM1 * VR * TAUYYR / ASQR * ALAM1
     1                   + TLT43*TAUYYR + TLT53*TAUYYR ) * S(I,4)
            A(I, 7, 4) = ( GAMM1 * WR * TAUYYR / ASQR * ALAM1
     1                   + TLT44*TAUYYR + TLT54*TAUYYR ) * S(I,4)
            A(I, 7, 5) = ( - GAMM1 * TAUYYR / ASQR * ALAM1
     1                   + TLT45*TAUYYR + TLT55*TAUYYR ) * S(I,4)
            A(I, 7, 6) = ( 0.5E0 * GAMM1 * TAUYYR / ASQR * ALAM1
     1                   + TLT46*TAUYYR + TLT56*TAUYYR ) * S(I,4)
            A(I, 7, 7) = ( (1.E0 + 0.5E0 * GAMM1 * TAUYYR / ASQR)*ALAM1
     1                   + TLT46*TAUYYR + TLT56*TAUYYR ) * S(I,4)
            A(I, 7, 8) = ( 0.5E0 * GAMM1 * TAUYYR / ASQR * ALAM1
     1                   + TLT46*TAUYYR + TLT56*TAUYYR ) * S(I,4)
            A(I, 7, 9) = 0.0E0
            A(I, 7,10) = 0.0E0
            A(I, 7,11) = 0.0E0
            A(I, 7,12) = 0.0E0
C
            A(I, 8, 1) = ( ( - 0.5E0 * GAMM1 * QSQR * TAUZZR/ASQR)*ALAM1
     1                   + TLT41*TAUZZR + TLT51*TAUZZR ) * S(I,4)
            A(I, 8, 2) = ( GAMM1 * UR * TAUZZR / ASQR * ALAM1
     1                   + TLT42*TAUZZR + TLT52*TAUZZR ) * S(I,4)
            A(I, 8, 3) = ( GAMM1 * VR * TAUZZR / ASQR * ALAM1
     1                   + TLT43*TAUZZR + TLT53*TAUZZR ) * S(I,4)
            A(I, 8, 4) = ( GAMM1 * WR * TAUZZR / ASQR * ALAM1
     1                   + TLT44*TAUZZR + TLT54*TAUZZR ) * S(I,4)
            A(I, 8, 5) = ( - GAMM1 * TAUZZR / ASQR * ALAM1
     1                   + TLT45*TAUZZR + TLT55*TAUZZR ) * S(I,4)
            A(I, 8, 6) = ( 0.5E0 * GAMM1 * TAUZZR / ASQR * ALAM1
     1                   + TLT46*TAUZZR + TLT56*TAUZZR ) * S(I,4)
            A(I, 8, 7) = ( 0.5E0 * GAMM1 * TAUZZR / ASQR * ALAM1
     1                   + TLT46*TAUZZR + TLT56*TAUZZR ) * S(I,4)
            A(I, 8, 8) = ( (1.E0 + 0.5E0 * GAMM1 * TAUZZR / ASQR)*ALAM1
     1                   + TLT46*TAUZZR + TLT56*TAUZZR ) * S(I,4)
            A(I, 8, 9) = 0.0E0
            A(I, 8,10) = 0.0E0
            A(I, 8,11) = 0.0E0
            A(I, 8,12) = 0.0E0
C
            A(I, 9, 1) = ( ( - 0.5E0 * GAMM1 * QSQR * TAUXYR/ASQR)*ALAM1
     1                   + TLT41*TAUXYR + TLT51*TAUXYR ) * S(I,4)
            A(I, 9, 2) = ( GAMM1 * UR * TAUXYR / ASQR * ALAM1
     1                   + TLT42*TAUXYR + TLT52*TAUXYR ) * S(I,4)
            A(I, 9, 3) = ( GAMM1 * VR * TAUXYR / ASQR * ALAM1
     1                   + TLT43*TAUXYR + TLT53*TAUXYR ) * S(I,4)
            A(I, 9, 4) = ( GAMM1 * WR * TAUXYR / ASQR * ALAM1
     1                   + TLT44*TAUXYR + TLT54*TAUXYR ) * S(I,4)
            A(I, 9, 5) = ( - GAMM1 * TAUXYR / ASQR * ALAM1
     1                   + TLT45*TAUXYR + TLT55*TAUXYR ) * S(I,4)
            A(I, 9, 6) = ( 0.5E0 * GAMM1 * TAUXYR / ASQR * ALAM1
     1                   + TLT46*TAUXYR + TLT56*TAUXYR ) * S(I,4)
            A(I, 9, 7) = ( 0.5E0 * GAMM1 * TAUXYR / ASQR * ALAM1
     1                   + TLT46*TAUXYR + TLT56*TAUXYR ) * S(I,4)
            A(I, 9, 8) = ( 0.5E0 * GAMM1 * TAUXYR / ASQR * ALAM1
     1                   + TLT46*TAUXYR + TLT56*TAUXYR ) * S(I,4)
            A(I, 9, 9) = ALAM1 * S(I,4)
            A(I, 9,10) = 0.0E0
            A(I, 9,11) = 0.0E0
            A(I, 9,12) = 0.0E0
C
            A(I,10, 1) = ( ( - 0.5E0 * GAMM1 * QSQR * TAUXZR/ASQR)*ALAM1
     1                   + TLT41*TAUXZR + TLT51*TAUXZR ) * S(I,4)
            A(I,10, 2) = ( GAMM1 * UR * TAUXZR / ASQR * ALAM1
     1                   + TLT42*TAUXZR + TLT52*TAUXZR ) * S(I,4)
            A(I,10, 3) = ( GAMM1 * VR * TAUXZR / ASQR * ALAM1
     1                   + TLT43*TAUXZR + TLT53*TAUXZR ) * S(I,4)
            A(I,10, 4) = ( GAMM1 * WR * TAUXZR / ASQR * ALAM1
     1                   + TLT44*TAUXZR + TLT54*TAUXZR ) * S(I,4)
            A(I,10, 5) = ( - GAMM1 * TAUXZR / ASQR * ALAM1
     1                   + TLT45*TAUXZR + TLT55*TAUXZR ) * S(I,4)
            A(I,10, 6) = ( 0.5E0 * GAMM1 * TAUXZR / ASQR * ALAM1
     1                   + TLT46*TAUXZR + TLT56*TAUXZR ) * S(I,4)
            A(I,10, 7) = ( 0.5E0 * GAMM1 * TAUXZR / ASQR * ALAM1
     1                   + TLT46*TAUXZR + TLT56*TAUXZR ) * S(I,4)
            A(I,10, 8) = ( 0.5E0 * GAMM1 * TAUXZR / ASQR * ALAM1
     1                   + TLT46*TAUXZR + TLT56*TAUXZR ) * S(I,4)
            A(I,10, 9) = 0.0E0
            A(I,10,10) = ALAM1 * S(I,4)
            A(I,10,11) = 0.0E0
            A(I,10,12) = 0.0E0
C
            A(I,11, 1) = ( ( - 0.5E0 * GAMM1 * QSQR * TAUYZR/ASQR)*ALAM1
     1                   + TLT41*TAUYZR + TLT51*TAUYZR ) * S(I,4)
            A(I,11, 2) = ( GAMM1 * UR * TAUYZR / ASQR * ALAM1
     1                   + TLT42*TAUYZR + TLT52*TAUYZR ) * S(I,4)
            A(I,11, 3) = ( GAMM1 * VR * TAUYZR / ASQR * ALAM1
     1                   + TLT43*TAUYZR + TLT53*TAUYZR ) * S(I,4)
            A(I,11, 4) = ( GAMM1 * WR * TAUYZR / ASQR * ALAM1
     1                   + TLT44*TAUYZR + TLT54*TAUYZR ) * S(I,4)
            A(I,11, 5) = ( - GAMM1 * TAUYZR / ASQR * ALAM1
     1                   + TLT45*TAUYZR + TLT55*TAUYZR ) * S(I,4)
            A(I,11, 6) = ( 0.5E0 * GAMM1 * TAUYZR / ASQR * ALAM1
     1                   + TLT46*TAUYZR + TLT56*TAUYZR ) * S(I,4)
            A(I,11, 7) = ( 0.5E0 * GAMM1 * TAUYZR / ASQR * ALAM1
     1                   + TLT46*TAUYZR + TLT56*TAUYZR ) * S(I,4)
            A(I,11, 8) = ( 0.5E0 * GAMM1 * TAUYZR / ASQR * ALAM1
     1                   + TLT46*TAUYZR + TLT56*TAUYZR ) * S(I,4)
            A(I,11, 9) = 0.0E0
            A(I,11,10) = 0.0E0
            A(I,11,11) = ALAM1 * S(I,4)
            A(I,11,12) = 0.0E0
C
            A(I,12, 1) = ( ( - 0.5E0 * GAMM1 * QSQR * T2R/ASQR)*ALAM1
     1                   + TLT41*T2R + TLT51*T2R ) * S(I,4)
            A(I,12, 2) = ( GAMM1 * UR * T2R / ASQR * ALAM1
     1                   + TLT42*T2R + TLT52*T2R ) * S(I,4)
            A(I,12, 3) = ( GAMM1 * VR * T2R / ASQR * ALAM1
     1                   + TLT43*T2R + TLT53*T2R ) * S(I,4)
            A(I,12, 4) = ( GAMM1 * WR * T2R / ASQR * ALAM1
     1                   + TLT44*T2R + TLT54*T2R ) * S(I,4)
            A(I,12, 5) = ( - GAMM1 * T2R / ASQR * ALAM1
     1                   + TLT45*T2R + TLT55*T2R ) * S(I,4)
            A(I,12, 6) = ( 0.5E0 * GAMM1 * T2R / ASQR * ALAM1
     1                   + TLT46*T2R + TLT56*T2R ) * S(I,4)
            A(I,12, 7) = ( 0.5E0 * GAMM1 * T2R / ASQR * ALAM1
     1                   + TLT46*T2R + TLT56*T2R ) * S(I,4)
            A(I,12, 8) = ( 0.5E0 * GAMM1 * T2R / ASQR * ALAM1
     1                   + TLT46*T2R + TLT56*T2R ) * S(I,4)
            A(I,12, 9) = 0.0E0
            A(I,12,10) = 0.0E0
            A(I,12,11) = 0.0E0
            A(I,12,12) = ALAM1 * S(I,4)
C
  420    CONTINUE
      ENDIF
C
C     Finished with T*LAMBDA*T(-1) calculation.
C
      RETURN
      END
