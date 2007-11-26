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
      SUBROUTINE TTrR     (ITURB, NPTS, Q, S, EL, EM, TKE, R, TTR)
C
C     Routine to calculate T^T R.  This routine is done in two steps:
C          I.  Calculate T^T R for mean-flow equations
C          II. Calculate T^T R for additional equations
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
C     TTR            : T R (T times R)
C
      include '../../src-flo/header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          S     (NPTS,4),
     2          EL    (NPTS,3),
     3          EM    (NPTS,3),
     4          TKE   (NPTS),
     5          R     (NPTS,NF),
     6          TTR   (NPTS,NF)
C
C     I. Calculate the T R for the mean-flow equations
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
         H     = ASQ / GAMM1 + 0.5E0 * QSQ + TKE(I)
C
C
         UBAR  = S (I,1) * U + S (I,2) * V + S (I,3) * W
         VBAR  = EL(I,1) * U + EL(I,2) * V + EL(I,3) * W
         WBAR  = EM(I,1) * U + EM(I,2) * V + EM(I,3) * W
C
         CAPG  = GAMM1 / ASQ
C
         TTR(I, 1) =                      R(I,1) 
     1             + U                  * R(I,2)
     2             + V                  * R(I,3)
     3             + W                  * R(I,4)
     4             + (0.5E0*QSQ+TKE(I)) * R(I,5)
         TTR(I, 2) =
     1             + EL(I,1)            * R(I,2)
     2             + EL(I,2)            * R(I,3)
     3             + EL(I,3)            * R(I,4)
     4             + VBAR               * R(I,5)
         TTR(I, 3) =
     1             + EM(I,1)            * R(I,2)
     2             + EM(I,2)            * R(I,3)
     3             + EM(I,3)            * R(I,4)
     4             + WBAR               * R(I,5)
         TTR(I, 4) =                      R(I,1) 
     1             + (U+S(I,1)*A)       * R(I,2)
     2             + (V+S(I,2)*A)       * R(I,3)
     3             + (W+S(I,3)*A)       * R(I,4)
     4             + (H+UBAR*A)         * R(I,5)
         TTR(I, 5) =                      R(I,1) 
     1             + (U-S(I,1)*A)       * R(I,2)
     2             + (V-S(I,2)*A)       * R(I,3)
     3             + (W-S(I,3)*A)       * R(I,4)
     4             + (H-UBAR*A)         * R(I,5)
  100 CONTINUE
C
C     II.   Additional equations for T R
C
C     II-b. Two-Equation model
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         STOP "TTrR: viscous not finished"
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
            TTR(I, 5) = TTR(I, 5) + R(I, 6)
            TTR(I, 6) = TKE(I)*(R(I, 1) + R(I, 4) + R(I, 5)) + R(I, 6)
            TTR(I, 7) = EPSLN *(R(I, 1) + R(I, 4) + R(I, 5)) + R(I, 7)
  210    CONTINUE
C
C     II-c. Reynolds Stress model
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         STOP "TTrR: viscous not finished"
         DO 220 I = 1, NPTS
            TAUXX = Q(I, 6)
            TAUYY = Q(I, 7)
            TAUZZ = Q(I, 8)
            TAUXY = Q(I, 9)
            TAUXZ = Q(I,10)
            TAUYZ = Q(I,11)
            EPSLN = Q(I,12)
C
            TERM      = R(I, 1) + R(I, 4) + R(I, 5)
            TRACE     = R(I, 6) + R(I, 7) + R(I, 8)
            TTR(I, 5) = TTR(I, 5) + 0.5E0 * TRACE
            TTR(I, 6) = TAUXX*TERM + R(I, 6)
            TTR(I, 7) = TAUYY*TERM + R(I, 7)
            TTR(I, 8) = TAUZZ*TERM + R(I, 8)
            TTR(I, 9) = TAUXY*TERM + R(I, 9)
            TTR(I,10) = TAUXZ*TERM + R(I,10)
            TTR(I,11) = TAUYZ*TERM + R(I,11)
            TTR(I,12) = EPSLN*TERM + R(I,12)
  220    CONTINUE
      ENDIF
C
      RETURN
      END
