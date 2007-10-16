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
      SUBROUTINE FLUROE(ITURB, QP, QM, S, F)
C
C     Routine to calculate interface fluxes using Roe's
C     approximate Riemann solver.  This routine calculates the
C     fluxes for the Favre averaged Navier-Stokes equations for inviscid,
C     laminar, algebraic turbulence model, two-equation turbulence model
C     or Reynolds stress turbulence model.  It computes all of the inviscid
C     fluxes for the above models.
C
C     The two-equation and Reynolds stress models need the turbulent kinetic 
C     energy (K) in the mean flow equations.  Therefore, this routine is
C     organized into the following sections:
C          1. Calculate Turbulent Kinetic Energy (TKEP,TKEM)
C             a. Inviscid, laminar, 
C                algebraic model       => K = 0
C             b. Two-Equation model    => K = K (Q(6))
C             c. Reynolds stress model => K = (Tau_xx+Tau_yy+Tau_zz)/2
C          2. Calculate Total Enthalpy (HP,HM)
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
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     S              : Metrics
C     TKEP,TKEM      : Turbulent kinetic energy for QP,QM
C     HP,HM          : Total enthalpy for QP,QM
C     WORK           : Temporary storage for the terms common to the 
C                      T |Lambda| T[-1] DQ term
C     F              : Flux at the interface
C
C---------------------------------------------------------------------- 
      include '../header/common.h'
      DIMENSION QP(*), QM(*), S(*), F(*)

      DIMENSION WORK(4)
C
C     1. Calculate the Turbulence Kinetic Energy (TKEP,TKEM)
C
C     1-b. Turbulence Kinetic Energy for Two Equation Turbulence models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         TKEP = QP(6)
         TKEM = QM(6)
C
C     1-c. Turbulence Kinetic Energy for Reynolds Stress Turbulence model
C
      ELSE IF (ITURB .EQ. ITRS) THEN
C
         TAUXXP = QP(6)
         TAUYYP = QP(7)
         TAUZZP = QP(8)
         TKEP   = 0.5E0 * (TAUXXP + TAUYYP + TAUZZP)
C
         TAUXXM = QM(6)
         TAUYYM = QM(7)
         TAUZZM = QM(8)
         TKEM   = 0.5E0 * (TAUXXM + TAUYYM + TAUZZM)
C
C     1-a. Turbulence Kinetic Energy = 0 for Inviscid, Laminar and
C          Algebraic Turbulence models
C
      ELSE
         TKEP = 0.E0
         TKEM = 0.E0
      ENDIF
C
C     2. Calculate Total Enthalpy
C
C     2-a. Perfect Gas
C
      RHOP   = QP(1)
      UP     = QP(2) 
      VP     = QP(3)
      WP     = QP(4)
      PP     = QP(5)
      HP     = GAMMA * PP / (RHOP * GAMM1) 
     1         + 0.5E0 * (UP*UP + VP*VP + WP*WP) + TKEP
      UBARP  = S(1)*UP + S(2)*VP + S(3)*WP
C
      RHOM   = QM(1)
      UM     = QM(2) 
      VM     = QM(3)
      WM     = QM(4)
      PM     = QM(5)
      HM     = GAMMA * PM / (RHOM * GAMM1)
     1         + 0.5E0 * (UM*UM + VM*VM + WM*WM) + TKEM
      UBARM  = S(1)*UM + S(2)*VM + S(3)*WM
C
C     3. Calculate F for Mean Flow Equations at the cell interfaces 
C        using Roe's approximate Riemann Solver
C
C
C     Calculate T * ABS(LAMBDA) * T(-1) * DQ
C
      DRHO   = RHOP   - RHOM
      DU     = UP     - UM
      DV     = VP     - VM
      DW     = WP     - WM
      DUBAR  = UBARP  - UBARM
      DP     = PP     - PM
C
      DK     = TKEP - TKEM
C
      RHOMSR = SQRT (RHOM) 
      RHOPSR = SQRT (RHOP)
      RHOAV  = RHOMSR + RHOPSR
      RHOAVI = 1.0 / RHOAV
      UR     = ( RHOMSR*UM + RHOPSR*UP ) * RHOAVI
      VR     = ( RHOMSR*VM + RHOPSR*VP ) * RHOAVI
      WR     = ( RHOMSR*WM + RHOPSR*WP ) * RHOAVI
      UBARR  = S(1)*UR + S(2)*VR + S(3)*WR
      HR     = ( RHOMSR*HM   + RHOPSR*HP )   * RHOAVI
C
      TKER   = ( RHOMSR*TKEM + RHOPSR*TKEP ) * RHOAVI
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
C     Entropy Fix: Based on Yee, Warming, and Harten (1985) and Yee ()
C       Apply fix only to acoustic eigenvalues to preserve boundary layer
C       calculations. Avoidance of the carbuncle phenomena requires fix
C       applied to convective eigenvalue also. Details based on conversations
C       with Jeff White.
C       Begin {N.B. Comment out this section to remove the entropy fix
C
      TSTLAM = CNTRPY * (SQRT(QSQR) + AR)
      TSTLAM = MAX (TSTLAM, RSMALL)
      IF (ALAM4 .LT. TSTLAM) ALAM4=0.5E0*(ALAM4*ALAM4/TSTLAM+TSTLAM)
      IF (ALAM5 .LT. TSTLAM) ALAM5=0.5E0*(ALAM5*ALAM5/TSTLAM+TSTLAM)
C
C       End} Entropy Fix
C
C     Temporary constants used in T * ABS(LAMBDA) * T(-1) * DQ term
C
      WORK(1) = ALAM1 * (DRHO - DP/ASQR)
      WORK(2) = ALAM1 * RHOR
      WORK(3) = ALAM4 * 0.5E0 * (DP/ASQR + RHOR*DUBAR/AR)
      WORK(4) = ALAM5 * 0.5E0 * (DP/ASQR - RHOR*DUBAR/AR)
C
C     Calculate F = 0.5 * (F+ + F-) - 0.5 * T * ABS(LAMBDA) * T(-1) * DQ
C     Area term is included in the RESID routine.
C
      F(1)  = 0.5E0 * 
     1          ( RHOP*UBARP
     2          + RHOM*UBARM
     3          - ( WORK(1)
     4            + WORK(3)
     5            + WORK(4) ) )
      F(2)  = 0.5E0 *
     1          ( (RHOP*UBARP*UP + S(1)*PP)
     2          + (RHOM*UBARM*UM + S(1)*PM)
     3          - ( UR * WORK(1)
     4            + (DU - S(1)*DUBAR) * WORK(2)
     5            + (UR + AR*S(1))    * WORK(3)
     6            + (UR - AR*S(1))    * WORK(4) ) )
      F(3)  = 0.5E0 * 
     1          ( (RHOP*UBARP*VP + S(2)*PP)
     2          + (RHOM*UBARM*VM + S(2)*PM)
     3          - ( VR * WORK(1)
     4            + (DV - S(2)*DUBAR) * WORK(2)
     5            + (VR + AR*S(2))    * WORK(3)
     6            + (VR - AR*S(2))    * WORK(4) ) )
      F(4)  = 0.5E0 * 
     1          ( (RHOP*UBARP*WP + S(3)*PP)
     2          + (RHOM*UBARM*WM + S(3)*PM)
     3          - ( WR * WORK(1)
     4            + (DW - S(3)*DUBAR) * WORK(2)
     5            + (WR + AR*S(3))    * WORK(3)
     6            + (WR - AR*S(3))    * WORK(4) ) )
      F(5)  = 0.5E0 *
     1          ( RHOP*UBARP*HP
     2          + RHOM*UBARM*HM
     3          - ( (0.5E0*QSQR + TKER) * WORK(1)
     4            + (UR*DU+VR*DV+WR*DW-UBARR*DUBAR+DK)*WORK(2)
     5            + (HR + AR*UBARR)    * WORK(3)
     6            + (HR - AR*UBARR)    * WORK(4) ) )
C
C     4. Calculate Additional Advection Fluxes
C
C     4-a. Inviscid, Laminar, Algebraic Turbulence models => no additional flux
C
C     4-b. Two-Equation Turbulence models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
C
         T2P     = QP(7)
         T2M     = QM(7)
         DT2     = T2P - T2M
         T2R     = ( RHOMSR*T2M + RHOPSR*T2P ) * RHOAVI
C
         F(6)  = 0.5E0 *
     1             ( RHOP*UBARP*TKEP
     2             + RHOM*UBARM*TKEM
     3             - ( TKER * WORK(1)
     4               + DK   * WORK(2)
     5               + TKER * WORK(3)
     6               + TKER * WORK(4) ) )
         F(7)  = 0.5E0 *
     1             ( RHOP*UBARP*T2P
     2             + RHOM*UBARM*T2M
     3             - ( T2R * WORK(1)
     4               + DT2 * WORK(2)
     5               + T2R * WORK(3)
     6               + T2R * WORK(4) ) )
C
C     4-c. Reynolds Stress Turbulence models
C
      ELSE IF (ITURB .EQ. ITRS) THEN
C
         TAUXXP  = QP(6)
         TAUYYP  = QP(7)
         TAUZZP  = QP(8)
         TAUXYP  = QP(9)
         TAUXZP  = QP(10)
         TAUYZP  = QP(11)
         T2P     = QP(12)
C
         TAUXXM  = QM(6)
         TAUYYM  = QM(7)
         TAUZZM  = QM(8)
         TAUXYM  = QM(9)
         TAUXZM  = QM(10)
         TAUYZM  = QM(11)
         T2M     = QM(12)
C
         DTAUXX  = TAUXXP - TAUXXM
         DTAUYY  = TAUYYP - TAUYYM
         DTAUZZ  = TAUZZP - TAUZZM
         DTAUXY  = TAUXYP - TAUXYM
         DTAUXZ  = TAUXZP - TAUXZM
         DTAUYZ  = TAUYZP - TAUYZM
         DT2     = T2P    - T2M
C
         TAUXXR  = ( RHOMSR*TAUXXM + RHOPSR*TAUXXP ) * RHOAVI
         TAUYYR  = ( RHOMSR*TAUYYM + RHOPSR*TAUYYP ) * RHOAVI
         TAUZZR  = ( RHOMSR*TAUZZM + RHOPSR*TAUZZP ) * RHOAVI
         TAUXYR  = ( RHOMSR*TAUXYM + RHOPSR*TAUXYP ) * RHOAVI
         TAUXZR  = ( RHOMSR*TAUXZM + RHOPSR*TAUXZP ) * RHOAVI
         TAUYZR  = ( RHOMSR*TAUYZM + RHOPSR*TAUYZP ) * RHOAVI
         T2R     = ( RHOMSR*T2M    + RHOPSR*T2P )    * RHOAVI
C
         F(6)  = 0.5E0 *
     1             ( RHOP*UBARP*TAUXXP
     2             + RHOM*UBARM*TAUXXM
     3             - ( TAUXXR * WORK(1)
     4               + DTAUXX * WORK(2)
     5               + TAUXXR * WORK(3)
     6               + TAUXXR * WORK(4) ) )
         F(7)  = 0.5E0 *
     1             ( RHOP*UBARP*TAUYYP
     2             + RHOM*UBARM*TAUYYM
     3             - ( TAUYYR * WORK(1)
     4               + DTAUYY * WORK(2)
     5               + TAUYYR * WORK(3)
     6               + TAUYYR * WORK(4) ) )
         F(8)  = 0.5E0 *
     1             ( RHOP*UBARP*TAUZZP
     2                + RHOM*UBARM*TAUZZM
     3                - ( TAUZZR * WORK(1)
     4                  + DTAUZZ * WORK(2)
     5                  + TAUZZR * WORK(3)
     6                  + TAUZZR * WORK(4) ) )
         F(9)  = 0.5E0 *
     1             ( RHOP*UBARP*TAUXYP
     2             + RHOM*UBARM*TAUXYM
     3             - ( TAUXYR * WORK(1)
     4               + DTAUXY * WORK(2)
     5               + TAUXYR * WORK(3)
     6               + TAUXYR * WORK(4) ) )
         F(10) = 0.5E0 *
     1             ( RHOP*UBARP*TAUXZP
     2             + RHOM*UBARM*TAUXZM
     3             - ( TAUXZR * WORK(1)
     4               + DTAUXZ * WORK(2)
     5               + TAUXZR * WORK(3)
     6               + TAUXZR * WORK(4) ) )
         F(11) = 0.5E0 *
     1             ( RHOP*UBARP*TAUYZP
     2             + RHOM*UBARM*TAUYZM
     3             - ( TAUYZR * WORK(1)
     4               + DTAUYZ * WORK(2)
     5               + TAUYZR * WORK(3)
     6               + TAUYZR * WORK(4) ) )
         F(12) = 0.5E0 *
     1             ( RHOP*UBARP*T2P
     2             + RHOM*UBARM*T2M
     3             - ( T2R * WORK(1)
     4               + DT2 * WORK(2)
     5               + T2R * WORK(3)
     6               + T2R * WORK(4) ) )

      ENDIF

      RETURN
      END
