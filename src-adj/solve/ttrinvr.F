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
      SUBROUTINE TTrINVR  (ITURB, NPTS, Q, S, EL, EM, TKE, R, TIR)
C
C     Routine to calculate T[-1]^T R.  This routine is done in two steps:
C          I.  Calculate T[-1]^T R for mean-flow equations
C          II. Calculate T[-1]^T R for additional equations
C              a. Inviscid, Laminar, Algebraic turb. => None
C              b. Two-Equation turb. model           => k, phi equations
C              c. Reynolds stress model              => tau_ij, epsilon
C
C     ITURB          : Turbulence model
C     NPTS           : Size of array to calculate on
C     Q              : Flow variables
C     S              : Metrics - including area term
C     EL, EM         : Unit tangency vectors
C     TKE            : Tubulent kinetic energy
C     R              : Residual to be multiplied by T[-1]
C     TIR            : T[-1] R
C
      include '../../src-flo/header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          S     (NPTS,4),
     2          EL    (NPTS,3),
     3          EM    (NPTS,3),
     4          TKE   (NPTS),
     5          R     (NPTS,NF),
     6          TIR   (NPTS,NF)
C
C     I. Calculate the T[-1] R for the mean-flow equations
C
      DO 100 I = 1, NPTS
         RHO   = Q(I,1)
         U     = Q(I,2)
         V     = Q(I,3)
         W     = Q(I,4)
         P     = Q(I,5)
C
         QSQ   = U*U + V*V + W*W
         ASQ   = GAMMA * P / RHO
         A     = SQRT ( ASQ )
C
         UBAR  = S (I,1) * U + S (I,2) * V + S (I,3) * W
         VBAR  = EL(I,1) * U + EL(I,2) * V + EL(I,3) * W
         WBAR  = EM(I,1) * U + EM(I,2) * V + EM(I,3) * W
C
         CAPG  = GAMM1 / ASQ
C
         TIR(I, 1) = (1.E0 - 0.5E0*CAPG*QSQ)       * R(I,1)
     1             - VBAR                          * R(I,2) 
     2             - WBAR                          * R(I,3) 
     3             + 0.5E0*(0.5E0*CAPG*QSQ-UBAR/A) * R(I,4) 
     4             + 0.5E0*(0.5E0*CAPG*QSQ+UBAR/A) * R(I,5)
         TIR(I, 2) = CAPG*U                        * R(I,1) 
     1             + EL(I,1)                       * R(I,2) 
     2             + EM(I,1)                       * R(I,3) 
     3             + 0.5E0*(S(I,1)/A-CAPG*U)       * R(I,4)
     4             - 0.5E0*(S(I,1)/A+CAPG*U)       * R(I,5)
         TIR(I, 3) = CAPG*V                        * R(I,1) 
     1             + EL(I,2)                       * R(I,2) 
     2             + EM(I,2)                       * R(I,3) 
     3             + 0.5E0*(S(I,2)/A-CAPG*V)       * R(I,4)
     4             - 0.5E0*(S(I,2)/A+CAPG*V)       * R(I,5)
         TIR(I, 4) = CAPG*W                        * R(I,1)
     1             + EL(I,3)                       * R(I,2)
     2             + EM(I,3)                       * R(I,3)
     3             + 0.5E0*(S(I,3)/A-CAPG*W)       * R(I,4)
     4             - 0.5E0*(S(I,3)/A+CAPG*W)       * R(I,5)
         TIR(I, 5) = -CAPG                         * R(I,1)
     1             + 0.5E0*CAPG                    * R(I,4)
     2             + 0.5E0*CAPG                    * R(I,5)
  100 CONTINUE
C
C     II.   Additional equations for T[-1] R
C
C     II-b. Two-Equation model
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         STOP "TTrINVR: viscous not finished"
         DO 210 I = 1, NPTS
            RHO   = Q(I, 1)
            P     = Q(I, 5)
            EPSLN = Q(I, 7)
C
            ASQ   = GAMMA * P / RHO
            A     = SQRT ( ASQ )
C
            CAPG  = GAMM1 / ASQ
C
            TIR(I, 1) =  TIR(I, 1) +         CAPG * R(I, 6)
            TIR(I, 4) =  TIR(I, 4) - 0.5E0 * CAPG * R(I, 6)
            TIR(I, 5) =  TIR(I, 5) - 0.5E0 * CAPG * R(I, 6)
            TIR(I, 6) = -TKE(I)*R(I, 1) + R(I, 6)
            TIR(I, 7) = -EPSLN *R(I, 1) + R(I, 7)
  210    CONTINUE
C
C     II-c. Reynolds Stress model
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         STOP "TTrINVR: viscous not finished"
         DO 220 I = 1, NPTS
            RHO   = Q(I, 1)
            P     = Q(I, 5)
            TAUXX = Q(I, 6)
            TAUYY = Q(I, 7)
            TAUZZ = Q(I, 8)
            TAUXY = Q(I, 9)
            TAUXZ = Q(I,10)
            TAUYZ = Q(I,11)
            EPSLN = Q(I,12)
C
            ASQ   = GAMMA * P / RHO
            A     = SQRT ( ASQ )
C
            CAPG  = GAMM1 / ASQ
C
            TRACE     =  R(I,6) + R(I,7) + R(I,8)
            TIR(I, 1) =  TIR(I, 1) + 0.50E0 * CAPG * TRACE
            TIR(I, 4) =  TIR(I, 4) - 0.25E0 * CAPG * TRACE
            TIR(I, 5) =  TIR(I, 5) - 0.25E0 * CAPG * TRACE
            TIR(I, 6) = -TAUXX*R(I, 1) + R(I, 6)
            TIR(I, 7) = -TAUYY*R(I, 1) + R(I, 7)
            TIR(I, 8) = -TAUZZ*R(I, 1) + R(I, 8)
            TIR(I, 9) = -TAUXY*R(I, 1) + R(I, 9)
            TIR(I,10) = -TAUXZ*R(I, 1) + R(I,10)
            TIR(I,11) = -TAUYZ*R(I, 1) + R(I,11)
            TIR(I,12) = -EPSLN*R(I, 1) + R(I,12)
  220    CONTINUE
      ENDIF
C
      RETURN
      END
