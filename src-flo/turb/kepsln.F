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





      SUBROUTINE DAMPMU (IMODEL, NPTS, Q, PROPS, TKE, FMU)
C
C Routine to calculate damping functions for turbulence models
C
C This routine is organized as follows:
C     1. Calculate FMU from appropriate damping function
C        a. High Re No.   (FMU = 1)
C        b. Speziale, Abid, Anderson
C        c. Zhang, So, Gatski, Speziale
C     2. For Algebraic Stress Model, modify FMU = FMU * CMUSTR / CMU
C        to account for variable CMUSTR in eddy-viscosity
C
C IMODEL         : Turbulence model control array
C NPTS           : Number of points to calculate damping terms at
C Q              : Primitive variables at cell centers
C PROPS          : Properties stored at cell centers
C                  PROPS(1) = RMU   molecular viscosity
C                  PROPS(2) = RMUT  turbulent eddy viscosity
C                  PROPS(3) = YPLUS Y+
C TKE            : Turbulent kinetic energy
C FMU            :
C
      include '../header/common.h'
C
      DIMENSION IMODEL(NMDL)
C
      DIMENSION Q     (NPTS,NQ),
     1          PROPS (NPTS,NP),
     2          TKE   (NPTS),
     3          FMU   (NPTS)
C
C Constant to prevent division by zero
C
      DATA SMALL / 1.0E-10 /
C
      ITURB = IMODEL( 1)
      IDAMP = IMODEL( 3)
      IASM  = IMODEL(10)
C
C 1. Calculate FMU
C
C 1-a. High Reynolds Number form => FDAMP = 1
C
      IF (IDAMP .EQ. IFMHR) THEN
         DO 100 I = 1, NPTS
            FMU(I) = 1.E0
  100    CONTINUE
C
C 1-b. Speziale et al. (AIAA 90-1481)
C      RET is the turbulent Reynolds number (appropriately non-dimensional)
C
      ELSE IF (IDAMP .EQ. IFMSAA) THEN
         DO 200 I = 1, NPTS
            RHO    = Q    (I,1)
            EPSLN  = Q    (I,IEPSLN)
            RMU    = PROPS(I,1)
            YPLUS  = PROPS(I,3)
C
C Turbulent Reynolds number.
C Check for positivity.
C
            RET    = RHO * TKE(I) * TKE(I) / (EPSLN * RMU + RSMALL) * 
     1               RE / FSMACH
            RET    = MAX (RET, SMALL)
C
C Damping function = O(1/y) near wall.  Fmu -> 1 far from wall.
C The form of Fmu in Speziale et al. needs to be limited to 1 as get
C away from wall as RET approaches a small value outside of bl (Abid)
C
            FMU(I) = (1.E0 + 3.45E0 / SQRT(RET)) * TANH( YPLUS / 70.E0 )
            FMU(I) = MIN (FMU(I), 1.E0)
  200    CONTINUE
C
C 1-c. Zhang, So, Gatski, Speziale (Near-Wall Turbulent Flows)
C      RET is the turbulent Reynolds number (appropriately non-dimensional)
C
      ELSE IF (IDAMP .EQ. IFMZSG) THEN
         DO 300 I = 1, NPTS
            RHO    = Q    (I,1)
            EPSLN  = Q    (I,IEPSLN)
            RMU    = PROPS(I,1)
            YPLUS  = PROPS(I,3)
C
C Turbulent Reynolds number.
C Check for positivity.
C
            RET    = RHO * TKE(I) * TKE(I) / (EPSLN*RMU + RSMALL) * 
     1               RE / FSMACH
            RET    = MAX (RET, SMALL)
C
C Damping function = O(1/y) near wall.  Fmu -> 1 far from wall.
C The form of Fmu in Speziale et al. needs to be limited to 1 as get
C away from wall as RET approaches a small value outside of bl (Abid)
C
            FMU(I) = (1.E0 + 3.45E0 / SQRT(RET)) * TANH(YPLUS / 115.E0)
            FMU(I) = MIN (FMU(I), 1.E0)
  300    CONTINUE
      ELSE
         WRITE (IOUT,1010) IDAMP
         STOP
      ENDIF
C
C 2. Algebraic Stress Model  
C    Set FMU = FMU * CMUSTR / CMU to account for variable CMUSTR
C
c-N.B. Currently not including cmu* in heat transfer and diffusion for ASM
c-    IF (IASM .EQ. IASMGS) THEN
c-       DO 900 I = 1, NPTS
c-          CMUSTR = PROPS(I,ICMUST)
c-          FMU(I) = FMU(I) * CMUSTR / CMU
c-900    CONTINUE
c-    ENDIF
C
C Finished with damping functions
C
      RETURN
 1000 FORMAT (' ','DAMPMU  : Invalid turbulence model specified ',
     1            'ITURB = ',I5,
     2       /' ','Run Aborting.')
 1010 FORMAT (' ','DAMPMU  : Invalid damping functions specified ',
     1            'IDAMP = ',I5,
     2       /' ','Run Aborting.')
      END





      SUBROUTINE DAMPEE (IMODEL, NPTS, Q, PROPS, DEL2K, TKE,
     1                   TAU, DQDX, DQDY, DQDZ, ETILDE, XSIW, FDAMP)
C
C Routine to calculate the damping functions F1, F2 and epsilon^tilde for
C the source terms for the k-epsilon turbulence model.
C This routine also implements the Anisotropic Dissipation Rate Model.
C This routine is organized as:
C      I. Calculate damping functions
C         a. High Reynolds number form
C         b. Speziale-Abid-Anderson model
C         c. Zhang-So-Gatski-Speziale model
C         d. Zhang-So-Speziale-Lai model
C         e. RNG (Yakhot-Orszag-Thangam-Gatski-Speziale) model
C         f. So et al. 1995 model
C         g. Abid's algebraic stress(Abid,Morrison,Gatski,Speziale AIAA 96-0565)
C      II. Anisotropic dissipation rate model
C
C IMODEL         : Turbulence model control array
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C PROPS          : Properties stored at cell centers
C                  PROPS(1) = RMU   molecular viscosity
C                  PROPS(2) = RMUT  turbulent eddy viscosity
C                  PROPS(3) = YPLUS Y+
C DEL2K          :
C TKE            : Turbulent kinetic energy
C ETILDE         :
C XSIW           :
C FDAMP          :
C
      include '../header/common.h'
C
      DIMENSION IMODEL(NMDL)
C
      DIMENSION Q     (NPTS,NQ),
     1          PROPS (NPTS,NP)
C
      DIMENSION DEL2K (NPTS),
     1          TKE   (NPTS),
     2          TAU   (NPTS,6),
     3          DQDX  (NPTS,NQ),
     4          DQDY  (NPTS,NQ),
     5          DQDZ  (NPTS,NQ),
     6          ETILDE(NPTS),
     7          XSIW  (NPTS),
     8          FDAMP (NPTS,2)
C
      T23    = 2.E0 / 3.E0
C
      ITURB  = IMODEL(1)
      IEPSEQ = IMODEL(2)
      IADRM  = IMODEL(4)
C
C I. Form the damping functions for the source terms epsilon equation 
C
C I-a. High Reynolds number form
C
      IF (IEPSEQ .EQ. IEEHR) THEN
         DO 100 I = 1, NPTS
            EPSLN      = Q(I,IEPSLN)
            FDAMP(I,1) = 1.E0
            FDAMP(I,2) = 1.E0
            ETILDE(I)  = EPSLN
            XSIW(I)    = 0.E0
  100    CONTINUE
C
C I-b. Speziale-Abid-Anderson. AIAA 90-1481 for wall bounded flows
C
      ELSE IF (IEPSEQ .EQ. IEESAA) THEN
         DO 200 I = 1, NPTS
C
            RHO   = Q(I,1)
            P     = Q(I,5)
            EPSLN = Q(I,IEPSLN)
C
            RMU   = PROPS(I,1)
            RMUT  = PROPS(I,2)
            YPLUS = PROPS(I,3)
C
C RET = turbulent Reynolds number
C F1  = 1
C F2  = [(1 - exp(-Y+/A2))**2] * (1 - 2/9 exp(-(RET/6)**2) )
C
            FDAMP(I,1) = 1.E0
            RET        = RHO * TKE(I) * TKE(I) / (EPSLN * RMU)*RE/FSMACH
            F2         = ( 1.E0 - EXP(-YPLUS / A2KE) )**2
            CE2        = ( 1.E0 - 2.E0/9.E0 * EXP(-RET*RET/36.E0) )
            IF (YPLUS .GT. 500.E0) CE2 = 1.E0
            FDAMP(I,2) = CE2 * F2 
            ETILDE(I)  = EPSLN
            XSIW(I)    = 0.E0
  200    CONTINUE
C
C I-c. Zhang-So-Gatski-Speziale. Near-Wall Turbulent Flows
C
      ELSE IF (IEPSEQ .EQ. IEEZSG) THEN
         DO 300 I = 1, NPTS
C
            RHO    = Q(I,1)
            EPSLN  = Q(I,IEPSLN)
C
            RMU    = PROPS(I,1)
            RMUT   = PROPS(I,2)
            YPLUS  = PROPS(I,3)
            YNORML = ABS (PROPS(I,4))
C
C RET = turbulent Reynolds number
C
            RET        = RHO * TKE(I) * TKE(I) / (EPSLN * RMU)*RE/FSMACH
            FW2        = EXP( -(RET /  40.E0)**2 )
            IF (YPLUS .GT. 250.E0) FW2 = 0.E0
            FDAMP(I,1) = 1.E0
            FDAMP(I,2) = 1.E0
            ETILDE(I)  = EPSLN - RMU/RHO*DEL2K(I)*FSMACH/RE
            ETILDE(I)  = MAX (ETILDE(I), 0.E0)
C
            ESTAR      = EPSLN - 2.E0*RMU/RHO*TKE(I)/(YNORML*YNORML)
     1                               *FSMACH/RE
            ESTAR      = MAX (ESTAR,  0.E0)
C
C Production term: (TAUij = - rho tau_ij)
C    PROD = TAUXX*DUDX + TAUYY*DVDY + TAUZZ*DWDZ +
C           TAUXY*(DUDY+DVDX) + TAUXZ*(DUDZ+DWDX) + TAUYZ*(DVDZ+DWDY)
C
            PROD       =   TAU(I,1) * DQDX(I,2)
     1                   + TAU(I,2) * DQDY(I,3)
     2                   + TAU(I,3) * DQDZ(I,4)
     3                   + TAU(I,4) * (DQDY(I,2) + DQDX(I,3))
     4                   + TAU(I,5) * (DQDZ(I,2) + DQDX(I,4))
     5                   + TAU(I,6) * (DQDZ(I,3) + DQDY(I,4))
C
            IF (POSPRD) PROD = ABS (PROD)
C
            PRODUL     = PROD
C
C Limit production term for robustness:
C
            TSTDIS     = PRDLIM * RHO * EPSLN
            PROD       = MIN (PROD, TSTDIS)
C
            PRODE      = PRDE * PROD + PRDEM1 * PRODUL
C
            XSIW(I)    = FW2 * ( -2.E0 * RHO * EPSLN * ETILDE(I)
     1                         + 1.5E0 * RHO * ESTAR * ESTAR 
     2                         - 1.5E0 * CEPS1 * EPSLN * PRODE) / TKE(I)
  300    CONTINUE
C
C I-d. Zhang-So-Speziale-Lai. (AIAA J. Vol. 31, No. 1, Jan. 1993)
C
      ELSE IF (IEPSEQ .EQ. IEEZSL) THEN
         DO 400 I = 1, NPTS
C
            RHO    = Q(I,1)
            EPSLN  = Q(I,IEPSLN)
C
            RMU    = PROPS(I,1)
            RMUT   = PROPS(I,2)
            YPLUS  = PROPS(I,3)
            YNORML = ABS (PROPS(I,4))
C
C RET = turbulent Reynolds number
C
            RET        = RHO * TKE(I) * TKE(I) / (EPSLN * RMU)*RE/FSMACH
            FW2        = EXP( -(RET /  64.E0)**2 )
            IF (YPLUS .GT. 250.E0) FW2 = 0.E0
            FDAMP(I,1) = 1.E0
            FDAMP(I,2) = 1.E0
            ETILDE(I)  = EPSLN - RMU/RHO*DEL2K(I)*FSMACH/RE
            ETILDE(I)  = MAX (ETILDE(I), 0.E0)
C
            ESTAR      = EPSLN - 2.E0*RMU/RHO*TKE(I)/(YNORML*YNORML)
     1                               *FSMACH/RE
            ESTAR      = MAX (ESTAR,  0.E0)
C
            XSIW(I)    = FW2 * ( -2.E0 * RHO * EPSLN * ETILDE(I)
     1                         + 1.5E0 * RHO * ESTAR * ESTAR  ) / TKE(I)
  400    CONTINUE
C
C I-e. RNG Yakhot-Orszag-Thangam-Gatski-Speziale Phys. Fluids A Vol. 4
C      No. 7, July 1992
C      Coding based on original code of Prof. J. M. Cimbala.
C
C     eta_0 = sqrt [(C_eps2 - 1) / (C_mu * (C_eps1 - 1) )]
C
      ELSE IF (IEPSEQ .EQ. IEERNG) THEN
         ETA0   = 4.38E0
         BETRNG = 0.012E0
         DO 500 I = 1, NPTS
C
            RHO    = Q(I,1)
            EPSLN  = Q(I,IEPSLN)
C
            RMU    = PROPS(I,1)
            RMUT   = PROPS(I,2)
C
C Production term: (TAUij = - rho tau_ij)
C    PROD = TAUXX*DUDX + TAUYY*DVDY + TAUZZ*DWDZ +
C           TAUXY*(DUDY+DVDX) + TAUXZ*(DUDZ+DWDX) + TAUYZ*(DVDZ+DWDY)
C
            PROD       =   TAU(I,1) * DQDX(I,2)
     1                   + TAU(I,2) * DQDY(I,3)
     2                   + TAU(I,3) * DQDZ(I,4)
     3                   + TAU(I,4) * (DQDY(I,2) + DQDX(I,3))
     4                   + TAU(I,5) * (DQDZ(I,2) + DQDX(I,4))
     5                   + TAU(I,6) * (DQDZ(I,3) + DQDY(I,4))
C
            IF (POSPRD) PROD = ABS (PROD)
C
            PRODUL     = PROD
C
C Limit production term for robustness:
C
            TSTDIS     = PRDLIM * RHO * EPSLN
            PROD       = MIN (PROD, TSTDIS)
C
            PRODE      = PRDE * PROD + PRDEM1 * PRODUL
C
C S_ij
C
            SSXX       = DQDX(I,2)
            SSYY       = DQDY(I,3)
            SSZZ       = DQDZ(I,4)
            SSXY       = 0.5E0 * (DQDY(I,2) + DQDX(I,3))
            SSXZ       = 0.5E0 * (DQDZ(I,2) + DQDX(I,4))
            SSYZ       = 0.5E0 * (DQDZ(I,3) + DQDY(I,4))
C
            SIJSIJ     =     ( SSXX*SSXX + SSYY*SSYY + SSZZ*SSZZ
     1                 + 2.E0*(SSXY*SSXY + SSXZ*SSXZ + SSYZ*SSYZ) )
C
            ETA        = SQRT (2.E0 * SIJSIJ) * TKE(I) / EPSLN
C
C Polynomial in eta for RNG
C
            FACT       = ETA * (1.E0 - ETA / ETA0) 
     1                   / (1.E0 + BETRNG * ETA * ETA * ETA)
C
C XSIW = - R
C
            FDAMP(I,1) = 1.E0
            FDAMP(I,2) = 1.E0
            ETILDE(I)  = EPSLN
            XSIW(I)    = - FACT * PRODE * EPSLN / TKE(I)
  500    CONTINUE
C
C I-f. So et al. 1995
C
      ELSE IF (IEPSEQ .EQ. IEES95) THEN
         DO 600 I = 1, NPTS
C
            RHO    = Q(I,1)
            EPSLN  = Q(I,IEPSLN)
C
            RMU    = PROPS(I,1)
            RMUT   = PROPS(I,2)
            YPLUS  = PROPS(I,3)
            YNORML = ABS (PROPS(I,4))
C
C RET = turbulent Reynolds number
C
            RET         = RHO * TKE(I) * TKE(I) / (EPSLN*RMU)*RE/FSMACH
            FW2         = EXP( -(RET /  40.E0)**2 )
            FDAMP (I,1) = 1.E0
            FDAMP (I,2) = 1.E0 - 2.E0/9.E0 * EXP( -(RET / 6.E0)**2 )
            ETILDE(I)   = EPSLN
C
C Production term: (TAUij = - rho tau_ij)
C    PROD = TAUXX*DUDX + TAUYY*DVDY + TAUZZ*DWDZ +
C           TAUXY*(DUDY+DVDX) + TAUXZ*(DUDZ+DWDX) + TAUYZ*(DVDZ+DWDY)
C
            PROD        =   TAU(I,1) * DQDX(I,2)
     1                    + TAU(I,2) * DQDY(I,3)
     2                    + TAU(I,3) * DQDZ(I,4)
     3                    + TAU(I,4) * (DQDY(I,2) + DQDX(I,3))
     4                    + TAU(I,5) * (DQDZ(I,2) + DQDX(I,4))
     5                    + TAU(I,6) * (DQDZ(I,3) + DQDY(I,4))
C
            IF (POSPRD) PROD = ABS (PROD)
C
            PRODUL      = PROD
C
C Limit production term for robustness:
C
            TSTDIS      = PRDLIM * RHO * EPSLN
            PROD        = MIN (PROD, TSTDIS)
C 
            PRODE       = PRDE * PROD + PRDEM1 * PRODUL
C
C N.B. DEL2K = 2 {|del[sqrt(k)]|}^2  so that a 1/2 is needed on the ceps4 term
C
            XSIW  (I)   = EPSLN / TKE(I) * 
     1                  ( CEPS3 * SQRT(RET) * RHO * EPSLN 
     2                  + CEPS4 * 0.5E0 * DEL2K(I) * RMU * FSMACH / RE
     3                  + CEPS5 * FW2 * PRODE )
  600    CONTINUE
C
C I-g. Abid's algebraic stress (Abid,Morrison,Gatski,Speziale AIAA 96-0565)
C
      ELSE IF (IEPSEQ .EQ. IEEABD) THEN
         DO 700 I = 1, NPTS
C
            RHO    = Q(I,1)
            EPSLN  = Q(I,IEPSLN)
C
            RMU    = PROPS(I,1)
            YPLUS  = PROPS(I,3)
            YNORML = ABS (PROPS(I,4))
C
C F1  = 1
C F2  = [(1 - exp(-Y+/5.5))**2]
C
            FDAMP(I,1) = 1.E0
            FDAMP(I,2) = ( 1.E0 - EXP(-YPLUS / 5.5E0) )**2
c-----      REY        = RHO * SQRT(TKE(I)) * YNORML / RMU * RE / FSMACH
c-----      FDAMP(I,2) = ( 1.E0 - EXP(-REY / 12.E0) )
            ETILDE(I)  = EPSLN
            XSIW(I)    = 0.E0
  700    CONTINUE
C
C Error
C
      ELSE
         WRITE (IOUT,1000) IEPSEQ
         STOP
      ENDIF
C
C II. Anisotropic Dissipation Rate Model (ADRM)
C
      	IF (IADRM .EQ. IADRGS) THEN
C
C     Constants based on SSG pressure-strain
C
c----       C1     = 6.80E0
            C2     = 0.36E0
            C3     = 1.25E0
            C4     = 0.40E0
C
            CE5    = 5.80E0
            BET3   = 0.60E0
C
            BET    = 0.75E0 * (14.E0*BET3/11.E0 - 16.E0/33.E0)
C
            TBP1   = 2.E0 / 15.E0 * (1.E0 + BET)
C
            DO 900 I = 1, NPTS
               RHO        = Q(I,1)
               EPSLN      = Q(I,IEPSLN)
C
               TIMSCL     = TKE(I) / EPSLN
C
C Production term: (TAUij = - rho tau_ij)
C    PROD = TAUXX*DUDX + TAUYY*DVDY + TAUZZ*DWDZ +
C           TAUXY*(DUDY+DVDX) + TAUXZ*(DUDZ+DWDX) + TAUYZ*(DVDZ+DWDY)
C
               PROD       =   TAU(I,1) * DQDX(I,2)
     1                      + TAU(I,2) * DQDY(I,3)
     2                      + TAU(I,3) * DQDZ(I,4)
     3                      + TAU(I,4) * (DQDY(I,2) + DQDX(I,3))
     4                      + TAU(I,5) * (DQDZ(I,2) + DQDX(I,4))
     5                      + TAU(I,6) * (DQDZ(I,3) + DQDY(I,4))
C
               IF (POSPRD) PROD = ABS (PROD)
C
C Limit production term for robustness:
C
               TSTDIS     = PRDLIM * RHO * EPSLN
               PROD       = MIN (PROD, TSTDIS)
C
C     Constants based on SSG pressure-strain
C
               POVERE     = PROD / (RHO*EPSLN)
C
               C1         = 3.40E0 + 1.8E0 * POVERE
               CGEPS      = POVERE
               CG         = POVERE
C
               GSMALL     = 1.E0 / (0.5E0 * C1 + CG - 1.E0)
C
               AL1        = 0.5E0 * (4.E0 / 3.E0 - C2) * GSMALL
               AL2        = 0.5E0 * (2.E0        - C4) * GSMALL
               AL3        = 0.5E0 * (2.E0        - C3) * GSMALL
C
               GEPS       = 1.E0 / (CE5 + CGEPS - 1.E0)
C
               BET1       = ( 7.E0 * BET3 + 1.E0) / 11.E0 * GEPS
               BET2       = (15.E0 * BET3 - 1.E0) / 11.E0 * GEPS
C
               TB2A3      = 2.E0 / 3.E0  * (BET2 / AL3)**2
               TB1A2      = 2.E0         * (BET1 / AL2)**2
C
C     S_ij 
C
               SXX        = DQDX(I,2)
               SYY        = DQDY(I,3)
               SZZ        = DQDZ(I,4)
               SXY        = 0.5E0 * (DQDY(I,2) + DQDX(I,3))
               SXZ        = 0.5E0 * (DQDZ(I,2) + DQDX(I,4))
               SYZ        = 0.5E0 * (DQDZ(I,3) + DQDY(I,4))
C
C     W_ij
C
               WYX        = 0.5E0 * (DQDX(I,3) - DQDY(I,2))
               WZX        = 0.5E0 * (DQDX(I,4) - DQDZ(I,2))
               WZY        = 0.5E0 * (DQDY(I,4) - DQDZ(I,3))
               WXY        = 0.5E0 * (DQDY(I,2) - DQDX(I,3))
               WXZ        = 0.5E0 * (DQDZ(I,2) - DQDX(I,4))
               WYZ        = 0.5E0 * (DQDZ(I,3) - DQDY(I,4))
C
               ETASQ      =      AL3 * AL3 * TIMSCL * TIMSCL *
     1                          ( SXX*SXX + SYY*SYY + SZZ*SZZ
     2                    + 2.E0*(SXY*SXY + SXZ*SXZ + SYZ*SYZ) )
               ZETSQ      =      AL2 * AL2 * TIMSCL * TIMSCL *
     1                      2.E0*(WXY*WXY + WXZ*WXZ + WYZ*WYZ)
C
               CMUSTR     = 1.E0 / (1.E0 - TB2A3*ETASQ + TB1A2*ZETSQ)
               CEPSTR     = CEPS1 + TBP1 * GEPS * CMUSTR / CMU
C
               FDAMP(I,1) = FDAMP(I,1) * CEPSTR / CEPS1
  900      CONTINUE
        ENDIF
C
C Finished with k-epsilon damping terms
C
      RETURN
 1000 FORMAT (' ','DAMPEE: Invalid epsilon model specified = ',I5,
     1       /' ','Run Aborting.')
      END





      SUBROUTINE EPSCMP (IMODEL, NPTS, Q, PROPS, TKE, SIJ, WIJ, EPSC)
C
C Routine to calculate compressible dissipation rate where the dissipation
C rate has been decomposed into a solenoidal and compressible portion as
C      eps = eps_s + eps_c 
C where eps_s is the solenoidal (incompressible) portion and eps_c is the
C dilatational (compressible) portion.  This routine returns eps_c.
C This routine is used for two-equation (k-epsilon) and Reynolds stress models.
C
C This routine is organized as follows:
C     1. Calculate EPSC from appropriate model
C        a. No compressibility correction => epsc = 0
C        b. Sarkar's correction           => epsc = const * M_T^2 * eps
C        c. Ristorcelli's correction
C
C IMODEL         : Turbulence model control array
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C PROPS          : Properties stored at cell centers
C                  PROPS(1) = RMU   molecular viscosity
C                  PROPS(2) = RMUT  turbulent eddy viscosity
C                  PROPS(3) = YPLUS Y+
C TKE            : Turbulent Kinetic Energy
C SIJ            : 
C WIJ            :
C EPSC           : Compressible dissipation rate
C
      include '../header/common.h'
C
      DIMENSION IMODEL(NMDL)
C
      DIMENSION Q     (NPTS,NQ),
     1          PROPS (NPTS,NP),
     2          TKE   (NPTS),
     3          SIJ   (NPTS,6),
     4          WIJ   (NPTS,6),
     5          EPSC  (NPTS)
C
      IEPSC = IMODEL(11)
C
C 1. Calculate EPSC
C
C 1-a. No compressibility correction
C
      IF (IEPSC .EQ. IECNO) THEN
         DO 100 I = 1, NPTS
            EPSC(I) = 0.E0
  100    CONTINUE
C
C 1-b. Sarkar's compressible dissipation
C
      ELSE IF (IEPSC .EQ. IECSAR) THEN
         DO 200 I = 1, NPTS
            RHO     = Q    (I,1)
            P       = Q    (I,5)
            EPSLN   = Q    (I,IEPSLN)
C
C Turbulent Mach Number (from Sarkar et al.) 
C
            ASQ     = GAMMA * P / RHO
            RMTSQ   = 2.E0 * TKE(I) / ASQ
C
            EPSC(I) = ALF1 * RMTSQ * EPSLN
  200    CONTINUE
C
C 1-c. Ristorcelli's compressible dissipation
C
      ELSE IF (IEPSC .EQ. IECRIS) THEN
         T13   = 1.E0 / 3.E0
         T23   = 2.E0 / 3.E0
         T49   = 4.E0 / 9.E0
C
C Ristorcelli's integral constants - from experimental data
C     Slow Integral 1,2,3
C     Rapid Integral 1,2,3
C
         SI1   =  0.300E0
         SI2   = 13.768E0
         SI3   =  2.623E0
C
         RI1   =  1.392E0
         RI2   =  0.E0
         RI3   =  3.E0
C
         RALFA = 1.E0
C
         RCON1 = 16.E0 / (3.E0 * RALFA)
         RCON2 = T23 ** 5
         RCON3 = (1.E0 / 15.E0) ** 2
C
         DO 300 I = 1, NPTS
            RHO     = Q    (I,1)
            P       = Q    (I,5)
            EPSLN   = Q    (I,IEPSLN)
            RMU     = PROPS(I,1)
C
C Turbulent Mach Number, Reynolds Number (from Ristorcelli's definition)
C
            ASQ     = GAMMA * P / RHO
            RMTSQ   = T23 * TKE(I) / ASQ
            RMT4    = RMTSQ * RMTSQ
            RET     = T49 * RHO * TKE(I) * TKE(I) / (EPSLN * RMU) * 
     1                RE / FSMACH
C
            SXX     = SIJ(I,1)
            SYY     = SIJ(I,2)
            SZZ     = SIJ(I,3)
            SXY     = SIJ(I,4)
            SXZ     = SIJ(I,5)
            SYZ     = SIJ(I,6)
C
            WYX     = WIJ(I,1)
            WZX     = WIJ(I,2)
            WZY     = WIJ(I,3)
            WXY     = WIJ(I,4)
            WXZ     = WIJ(I,5)
            WYZ     = WIJ(I,6)
C
C S_ij * S_ij (symmetric tensor), W_ij * W_ij (anti-symmetric tensor)
C Note: Ristorcelli's notes define S_ij = .5 (U_i,j + U_j,i - 1/3 U_k,k)
C       the S_ij used here do NOT subtract off the trace.  The additional
C       term 1/3 trace^2 subtracted off of SIJSIJ accounts for this difference.
C
            TRACES  = SXX + SYY + SZZ
            SIJSIJ  =         SXX*SXX + SYY*SYY + SZZ*SZZ
     1                + 2.E0*(SXY*SXY + SXZ*SXZ + SYZ*SYZ)
     2                - T13 * TRACES * TRACES
            WIJWIJ  = - 2.E0*(WXY*WYX + WXZ*WZX + WYZ*WZY)
C
            SHATSQ  = TKE(I) * TKE(I) / (EPSLN * EPSLN) * SIJSIJ
            WHATSQ  = TKE(I) * TKE(I) / (EPSLN * EPSLN) * WIJWIJ
C
            EPSC(I) = EPSLN * 
     1            ( RCON1 * RMT4 / RET * (SI2 + 6.E0 * SI1 * SI3)
     2            + RCON2 * RMT4 / RET * (3.E0*SHATSQ+5.E0*WHATSQ)
     3              * (0.6E0*RI3 + RCON3*(13.E0*SHATSQ+15.E0*WHATSQ)
     4                                  *RALFA*RI1) )
  300    CONTINUE
C
C Error
C
      ELSE
         WRITE (IOUT,1000) IEPSC
         STOP
      ENDIF
C
C Finished with compressible dissipation
C
      RETURN
 1000 FORMAT (' ','EPSCMP  : Invalid compressible dissipation model ',
     1            'specified IEPSC = ',I5,
     2       /' ','Run Aborting.')
      END







      SUBROUTINE PRDIL  (IMODEL, NPTS, Q, PROPS, DQDX, DQDY, DQDZ, TAU, 
     1                   TKE, SIJ, WIJ, EPSC, PDIL)
C
C Routine to calculate pressure-dilatation model.  This routine returns p_dil.
C This routine is used for two-equation (k-epsilon) and Reynolds stress models.
C
C This routine is organized as follows:
C     1. Calculate PDIL from appropriate model
C        a. No pressure-dilatation            => pdil = 0
C        b. Sarkar's pressure-dilatation      => pdil = const * M_T^2 * eps
C        c. Ristorcelli's pressure-dilatation => pdil = chi   * M_T^2 * (P-e)
C
C IMODEL         : Turbulence model control array
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C PROPS          : Properties stored at cell centers
C                  PROPS(1) = RMU   molecular viscosity
C                  PROPS(2) = RMUT  turbulent eddy viscosity
C                  PROPS(3) = YPLUS Y+
C DQDX,DQDY,DQDZ : Derivatives of Q at cell centers
C TAU            : Reynolds stresses
C TKE            : Turbulent Kinetic Energy
C SIJ            : 
C WIJ            :
C EPSC           : Compressible dissipation rate
C PDIL           : Pressure-Dilatation
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
     5          SIJ   (NPTS,6),
     6          WIJ   (NPTS,6),
     7          EPSC  (NPTS),
     8          PDIL  (NPTS)
C
      IPDIL = IMODEL(12)
C
C 1. Calculate PDIL
C
C 1-a. No Pressure-Dilatation
C
      IF (IPDIL .EQ. IPDNO) THEN
         DO 100 I = 1, NPTS
            PDIL(I) = 0.E0
  100    CONTINUE
C
C 1-b. Sarkar's Pressure-Dilatation
C
      ELSE IF (IPDIL .EQ. IPDSAR) THEN
         DO 200 I = 1, NPTS
            RHO     = Q    (I,1)
            P       = Q    (I,5)
            EPSLN   = Q    (I,IEPSLN)
C
C Production term:
C    PROD = TAUXX*DUDX + TAUYY*DVDY + TAUZZ*DWDZ +
C           TAUXY*(DUDY+DVDX) + TAUXZ*(DUDZ+DWDX) + TAUYZ*(DVDZ+DWDY)
C
            PROD    =   TAU(I,1) * DQDX(I,2)
     1                + TAU(I,2) * DQDY(I,3)
     2                + TAU(I,3) * DQDZ(I,4)
     3                + TAU(I,4) * (DQDY(I,2) + DQDX(I,3))
     4                + TAU(I,5) * (DQDZ(I,2) + DQDX(I,4))
     5                + TAU(I,6) * (DQDZ(I,3) + DQDY(I,4))
C
            IF (POSPRD) PROD = ABS (PROD)
C
C Limit production term for robustness:
C
            TSTDIS  = PRDLIM * RHO * EPSLN
            PROD    = MIN (PROD, TSTDIS)
C
C Turbulent Mach Number (from Sarkar et al.) 
C
            ASQ     = GAMMA * P / RHO
            RMTSQ   = 2.E0 * TKE(I) / ASQ
            RMT     = SQRT (RMTSQ)
C
C Pressure-Dilatation
C This form is from Sarkar, ICASE 91-42
C
            PDIL(I) = - ALF2 * PROD * RMTSQ + ALF3 * RHO * EPSLN * RMTSQ
C
  200    CONTINUE
C
C 1-c. Ristorcelli's Pressure-Dilatation
C
      ELSE IF (IPDIL .EQ. IPDRIS) THEN
         T13   = 1.E0 / 3.E0
         T23   = 2.E0 / 3.E0
         T49   = 4.E0 / 9.E0
C
C Ristorcelli's integral constants - from experimental data
C     Slow Integral 1,2,3
C     Rapid Integral 1,2,3
C
         SI1   =  0.300E0
         SI2   = 13.768E0
         SI3   =  2.623E0
C
         RI1   =  1.392E0
         RI2   =  0.E0
         RI3   =  3.E0
C
         RALFA = 1.E0
C
         DO 300 I = 1, NPTS
            RHO     = Q    (I,1)
            P       = Q    (I,5)
            EPSLN   = Q    (I,IEPSLN)
            RMU     = PROPS(I,1)
C
C Production term:
C    PROD = TAUXX*DUDX + TAUYY*DVDY + TAUZZ*DWDZ +
C           TAUXY*(DUDY+DVDX) + TAUXZ*(DUDZ+DWDX) + TAUYZ*(DVDZ+DWDY)
C
            PROD    =   TAU(I,1) * DQDX(I,2)
     1                + TAU(I,2) * DQDY(I,3)
     2                + TAU(I,3) * DQDZ(I,4)
     3                + TAU(I,4) * (DQDY(I,2) + DQDX(I,3))
     4                + TAU(I,5) * (DQDZ(I,2) + DQDX(I,4))
     5                + TAU(I,6) * (DQDZ(I,3) + DQDY(I,4))
C
            IF (POSPRD) PROD = ABS (PROD)
C
C Limit production term for robustness:
C
            TSTDIS  = PRDLIM * RHO * EPSLN
            PROD    = MIN (PROD, TSTDIS)
C
C Turbulent Mach Number, Reynolds Number (from Ristorcelli's definition)
C
            ASQ     = GAMMA * P / RHO
            RMTSQ   = T23 * TKE(I) / ASQ
            RMT4    = RMTSQ * RMTSQ
C
            SXX     = SIJ(I,1)
            SYY     = SIJ(I,2)
            SZZ     = SIJ(I,3)
            SXY     = SIJ(I,4)
            SXZ     = SIJ(I,5)
            SYZ     = SIJ(I,6)
C
            WYX     = WIJ(I,1)
            WZX     = WIJ(I,2)
            WZY     = WIJ(I,3)
            WXY     = WIJ(I,4)
            WXZ     = WIJ(I,5)
            WYZ     = WIJ(I,6)
C
C S_ij * S_ij (symmetric tensor), W_ij * W_ij (anti-symmetric tensor)
C Note: Ristorcelli's notes define S_ij = .5 (U_i,j + U_j,i - 1/3 U_k,k)
C       the S_ij used here do NOT subtract off the trace.  The additional
C       term 1/3 trace^2 subtracted off of SIJSIJ accounts for this difference.
C
            TRACES  = SXX + SYY + SZZ
            SIJSIJ  =         SXX*SXX + SYY*SYY + SZZ*SZZ
     1                + 2.E0*(SXY*SXY + SXZ*SXZ + SYZ*SYZ)
     2                - T13 * TRACES * TRACES
            WIJWIJ  = - 2.E0*(WXY*WYX + WXZ*WZX + WYZ*WZY)
C
            SHATSQ  = TKE(I) * TKE(I) / (EPSLN * EPSLN) * SIJSIJ
            WHATSQ  = TKE(I) * TKE(I) / (EPSLN * EPSLN) * WIJWIJ
C
            PDIR    = T23 * T23 * T23 / 30.E0 * RALFA * RALFA * RI1
            PDI     = T23 * SI1 + PDIR * (3.E0*SHATSQ+5.E0*WHATSQ)
            CHIPD   = 2.E0 * PDI / 
     1                (1.E0+2.E0*PDI*RMTSQ+1.5E0*PDIR*RMT4*GAMMA*GAMM1)
C
C Pressure-Dilatation
C
            PDIL(I) = - CHIPD * RMTSQ * (PROD - RHO * (EPSLN + EPSC(I)))
  300    CONTINUE
C
C Error
C
      ELSE
         WRITE (IOUT,1000) IPDIL
         STOP
      ENDIF
C
C Finished with pressure-dilatation
C
      RETURN
 1000 FORMAT (' ','PRDIL   : Invalid Pressure-Dilatation model ',
     1            'specified IPDIL = ',I5,
     2       /' ','Run Aborting.')
      END







      SUBROUTINE RMUTKE (IMODEL, NPTS, Q, PROPS, TKE, FMU)
C
C     Routine to calculate the turbulent eddy viscosity  for the k-epsilon
C     turbulence model.
C
C     IMODEL         : Tur
C     NPTS           : Number of points 
C     Q              : Primitive variables at cell centers
C     PROPS          : Properties at cell centers
C     TKE            : Turbulent kinetic energy
C     NDAMP          :
C     FMU            : Damping functions for mut
C
      include '../header/common.h'
C
      DIMENSION IMODEL(NMDL)
C
      DIMENSION Q     (NPTS,NQ),
     1          PROPS (NPTS,NP),
     2          TKE   (NPTS),
     3          FMU   (NPTS)
C
C     Constant to prevent division by zero
C
      DATA SMALL / 1.0E-10 /
C
      ITURB = IMODEL(1)
C
C     Calculate damping functions
C
      CALL DAMPMU (IMODEL, NPTS, Q, PROPS, TKE, FMU)
C
C     Calculate the turbulent eddy viscosity at cell centers from values
C     of turbulent kinetic energy (k) and dissipation (epsilon).
C
      DO 100 I = 1, NPTS
         RHO        = Q(I,1)
         EPSLN      = Q(I,IEPSLN)
C
         PROPS(I,2) = CMU * FMU(I) * RHO * TKE(I) * TKE(I) / 
     1                (EPSLN + RSMALL) * RE / FSMACH
  100 CONTINUE
C
C     Finished calculating RMUT from k and epsilon
C
      RETURN
 1000 FORMAT (' ','RMUTKE: Invalid turbulence model specified ',
     1            'ITURB = ',I5,
     2       /' ','Run Aborting.')
      END






      SUBROUTINE SRCKE (IMODEL, NPTS, Q, PROPS, DQDX, DQDY, DQDZ, TAU, 
     1                  DEL2K, TKE, ETILDE, XSIW, FDAMP, 
     2                  EPSC, PDIL, SIJ, WIJ, SRC, NPRLIM, NPRNEG)
C
C Routine to calculate the source terms for the k-epsilon turbulence model.
C
C This routine is organized as follows:
C      1. Calculate damping functions and epsilon^tilde
C      2. Calculate compressible dissipation
C      3. Calculate pressure-dilatation
C      4. Form the source terms for the k-epsilon model 
C
C IMODEL         : Turbulence model control array
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C PROPS          : Properties stored at cell centers
C                  PROPS(1) = RMU   molecular viscosity
C                  PROPS(2) = RMUT  turbulent eddy viscosity
C                  PROPS(3) = YPLUS Y+
C DQDX,DQDY,DQDZ : Derivatives of Q at cell centers
C TAU            : Reynolds stresses
C DEL2K          :
C TKE            :
C ETILDE         :
C XSIW           :
C FDAMP          :
C EPSC           :
C PDIL           :
C SIJ            :
C WIJ            :
C SRC            : Source terms for the k-omega model at the cell centers
C                  of the I-plane
C NPRLIM         : Number of limitations to the production term
C NPRNEG         : Number of times invoke positivity preservation for production
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
     4          DEL2K (NPTS),
     5          TKE   (NPTS),
     6          ETILDE(NPTS),
     7          XSIW  (NPTS),
     8          FDAMP (NPTS,2),
     9          EPSC  (NPTS),
     A          PDIL  (NPTS),
     B          SIJ   (NPTS,6),
     C          WIJ   (NPTS,6),
     D          SRC   (NPTS,NF)
C
C 1.  Calculate damping functions and epsilon^tilde
C
      CALL DAMPEE (IMODEL, NPTS, Q, PROPS, DEL2K, TKE,
     1             TAU, DQDX, DQDY, DQDZ, ETILDE, XSIW, FDAMP)
C
C 2.  Calculate compressible dissipation
C
      CALL EPSCMP (IMODEL, NPTS, Q, PROPS, TKE, SIJ, WIJ, EPSC)
C
C 3.  Calculate pressure-dilatation
C
      CALL PRDIL  (IMODEL, NPTS, Q, PROPS, DQDX, DQDY, DQDZ, TAU, 
     1             TKE, SIJ, WIJ, EPSC, PDIL)
C
C 4.  Form the source terms for the k-epsilon model 
C
      DO 100 I = 1, NPTS
C
         RHO    = Q(I,1)
         P      = Q(I,5)
         EPSLN  = Q(I,7)
C
         RMU    = PROPS(I,1)
         RMUT   = PROPS(I,2)
         YPLUS  = PROPS(I,3)
C
C DU_k/DX_k = (DU/DX + DV/DY + DW/DZ) / 3
C
         DUKDXK = (DQDX(I,2) + DQDY(I,3) + DQDZ(I,4)) / 3.E0
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
         TSTDIS = PRDLIM * RHO * EPSLN
         IF (PROD .GT. TSTDIS) THEN
            NPRLIM = NPRLIM + 1
            PROD   = TSTDIS
         ENDIF
C
C Production term for epsilon: use either limited or unlimited form
C->  PRODE = - CEPS1*epsilon/k*rho*tau_ij*(Du_i/Dx_j-Du_k/Dx_k*delta_ij/3)
C->          - 4/3*rho*epsilon*Du_k/Dx_k
C->  from: Speziale and Sarkar, ICASE 91-9
C
         PRODE = CEPS1 * FDAMP(I,1) * EPSLN / TKE(I) 
     1           * (PRDE * PROD + PRDEM1 * PRODUL)
C
c->  1           ( TAU(I,1) * (DQDX(I,2) - DUKDXK)
c->  2           + TAU(I,2) * (DQDY(I,3) - DUKDXK)
c->  3           + TAU(I,3) * (DQDZ(I,4) - DUKDXK)
c->  4           + TAU(I,4) * (DQDY(I,2) + DQDX(I,3)) 
c->  5           + TAU(I,5) * (DQDZ(I,2) + DQDX(I,4))
c->  6           + TAU(I,6) * (DQDZ(I,3) + DQDY(I,4)) )
c->  7         - 4.E0 * RHO * EPSLN * DUKDXK
C
C Dissipation
C
         DISSK = RHO * (EPSLN + EPSC(I))
C
         DISSE = CEPS2 * FDAMP(I,2) * RHO * EPSLN * ETILDE(I) / TKE(I)
C
C Calculate the Pressure-Velocity correlation term:
C    PVEL  = CPVELK * mu_T * sigma_rho / rho^2 * (DrhoDx_i DPDx_i)
C
         PVEL  = CPVELK * RMUT * SIGRHO / (RHO * RHO) * FSMACH / RE *
     1           ( DQDX(I,1) * DQDX(I,5) + DQDY(I,1) * DQDY(I,5) 
     2           + DQDZ(I,1) * DQDZ(I,5) ) 
C
C Calculate Source Term
C
         SRC(I,1) = 0.E0
         SRC(I,2) = 0.E0
         SRC(I,3) = 0.E0
         SRC(I,4) = 0.E0
         SRC(I,5) = 0.E0
         SRC(I,6) = PROD  - DISSK - PVEL + PDIL(I)
         SRC(I,7) = PRODE - DISSE + XSIW(I)
  100 CONTINUE
C
C Finished with k-epsilon source terms
C
      RETURN
      END
