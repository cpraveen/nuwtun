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


      SUBROUTINE CHAPND (STRING1, STRING2)
C
C Routine to append character STRING2 to STRING1.
C
C N.B. This routine assumes that STRING1 is right padded with spaces. 
C      Unpredictable results will be returned if STRING1 is not right padded.
C
C STRING1 : Input base character string (with at least one terminating space)
C STRING2 : Input character string to be appended to STRING1
C
      CHARACTER*(*) STRING1,
     1              STRING2
C
      LEN1 = ICHLNG (STRING1)
C
      STRING1 = STRING1(1:LEN1) // STRING2
C
C Finished with CHAPND
C
      RETURN
      END



      INTEGER FUNCTION ICHLNG (STRING)
C
C Function to return the length of a string defined as the location where
C the first white space is located.
C
C STRING : Input string
C
      CHARACTER*80 STRING
C
      LOC = 1
  100 CONTINUE
      IF (STRING(LOC:LOC) .EQ. ' ') THEN
         GO TO 200
      ENDIF
      LOC = LOC + 1
      IF (LOC .GT. 80) THEN
         GO TO 200
      ENDIF
      GO TO 100
C
C Found space or end of string
C
  200 CONTINUE
      ICHLNG = LOC - 1
C
      RETURN
      END



      SUBROUTINE MONITR (IDIM, JDIM, KDIM, Q, NPRMON, IPRMON)
C
C     Routine to monitor output during the iteration cycle
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     Q              : Primitive variables at cell centers
C     NPRMON         : Number of monitor print statements for current block
C     IPRMON         : Monitor print control vector for current block
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION IPRMON(10,NPRMON)
C
      DIMENSION ISTRT(3), IEND(3), INC(3)
C
      CHARACTER*10 TITLE(3)
      DATA TITLE / 'J, K =    ',
     1             'K, I =    ',
     2             'I, J =    ' /
C
C     Loop over all print statements
C
      do 130 i = 0, idim+2
         do 120 k = 0, kdim+2
            do 110 j = 0, jdim+2
               if (q(j,k,i,4) .ne. 0.e0) then
                  write (iout,'(a,3i5,e15.7)') ' i,j,k,w ',
     1                                           i,j,k,q(j,k,i,4)
               endif
  110       continue
  120    continue
  130 continue
C
C     Output a few blank lines to make output look better
C
      WRITE (IOUT,1015)
C
      RETURN
 1000 FORMAT (//' ',20X,A10,2I5/
     1        ' ','   I   J   K',
     2            6X,'RHO',11X,'U',12X,'V',
     3            12X,'W',9X,'P/PINF',8X,'T/TINF')
 1010 FORMAT (' ',3I4,7E13.5)
 1015 FORMAT (//' ')
      END



      SUBROUTINE SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW)
C
C This routine sets the following constants for the various pressure-strain 
C correlations:
C      IEPSIJ, IEPSEQ, IDAMP, ISRCW, PSC1, PSC2, PSCW, PSALFA, PSBETA, PSGAMA
C
C IPSMDL         : Input pressure strain model to set constants for
C IEPSIJ         : Epsilon_ij model
C IEPSEQ         : Epsilon equation model
C IDAMP          : Damping used for mu_turbulent
C ISRCW          : Control of additional wall normal dependent source routines
C
      include '../header/common.h'
C
      IF (IPSMDL .EQ. IPLRR1) THEN
         IEPSIJ = IEISO
         IEPSEQ = IEEHR 
         IDAMP  = IFMHR 
         ISRCW  = INO
         PSC1   = 3.0E0
         PSC2   = 0.4E0
         PSCW   = 1.0E0
         PSALFA = 0.E0
         PSBETA = 0.E0
         PSGAMA = 0.E0
      ELSE IF (IPSMDL .EQ. IPLRR2) THEN
         IEPSIJ = IEISO
         IEPSEQ = IEEHR 
         IDAMP  = IFMHR 
         ISRCW  = INO
         PSC1   = 3.0E0
         PSC2   = 0.4E0
         PSCW   = 0.E0
         PSALFA = 0.E0
         PSBETA = 0.E0
         PSGAMA = 0.E0
      ELSE IF (IPSMDL .EQ. IPSSG ) THEN
         IEPSIJ = IEISO
         IEPSEQ = IEEHR 
         IDAMP  = IFMHR 
         ISRCW  = INO
         PSALFA = 0.E0
         PSBETA = 0.E0
         PSGAMA = 0.E0
      ELSE IF (IPSMDL .EQ. IPFLT ) THEN
         IEPSIJ = IEFLT
         IEPSEQ = IEEHR 
         IDAMP  = IFMHR 
         ISRCW  = INO
         PSALFA = 0.E0
         PSBETA = 0.E0
         PSGAMA = 0.E0
      ELSE IF (IPSMDL .EQ. IPLS  ) THEN
         IEPSIJ = IELAIS
         IEPSEQ = IEEZSG
         IDAMP  = IFMZSG
         ISRCW  = IYES
         PSC1   = 3.0E0
         PSC2   = 0.4E0
         PSCW   = 0.01363540E0
         PSALFA = 0.45E0
         PSBETA = 0.E0
         PSGAMA = 0.E0
      ELSE IF (IPSMDL .EQ. IPSO95) THEN
         IEPSIJ = IESO95
         IEPSEQ = IEES95
         IDAMP  = IFMZSG
         ISRCW  = INO
         PSC1   = 3.0E0
         PSC2   = 0.4E0
         PSCW   = 0.0E0
         PSALFA = - 0.36E0
         PSBETA =   0.0005E0
         PSGAMA =   0.075E0
      ELSE IF (IPSMDL .EQ. IPGL  ) THEN
         IEPSIJ = IEISO
         IEPSEQ = IEEHR 
         IDAMP  = IFMHR 
         ISRCW  = INO
         PSC1   = 3.6E0
         PSC2   = 0.6E0
         PSCW   = 0.E0
         PSALFA = 0.E0
         PSBETA = 0.E0
         PSGAMA = 0.E0
      ELSE
         WRITE (IOUT,1000) IPSMDL
      ENDIF
C
C Finished
C
      RETURN
 1000 FORMAT (' ','SETPIJ: ERROR-> Invalid pressure-strain model: ',
     1            'IPSMDL = ',I5)
      END




      SUBROUTINE SFHT  (IDIM, JDIM, KDIM, IPSTRT, IPEND, ITURB,
     1                  Q, PROPS, S, VOL, IBCDAT, RC, IBLOCK)
C
C Routine to output the skin friction, heat transfer, and boundary layer
C quantities.
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C ITURB          : Turbulence model
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C S              : Metrics
C VOL            : Cell volumes
C IBCDAT         : Data controlling boundary conditions:
C                  IBCDAT(1,ibc) = Boundary Condition Type
C                         2      = Direction (1=i, 2=j, 3=k)
C                         3      = i(1 or idim)
C                         4      = istrt
C                         5      = iend
C                         6      = jstrt
C                         7      = jend
C                         8      = kstrt
C                         9      = kend
C                         10     = profile identifier
C RC             : Cell center grid
C IBLOCK         : Block number to identify output
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
C
C Set control for higher order evaluation of du/dn and dT/dn
C     IHIORD = 0 for first  order
C     IHIORD = 1 for second order
C
      IHIORD = 1
C
C Output headers
C
      IF (IHIORD .EQ. 0) THEN
         WRITE (IOUT,1000) IBLOCK, '1st'
      ELSE
         WRITE (IOUT,1000) IBLOCK, '2nd'
      ENDIF
      WRITE (IOUT,1005) 
C
C 1. Setup memory for looping.
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
      IJKDIM(1) = IDIM
      IJKDIM(2) = JDIM
      IJKDIM(3) = KDIM
C
C Store boundary condition data in mnemonic names
C
      IBCTYP   = IBCDAT(1)
      IDIR     = IBCDAT(2)
      INROUT   = IBCDAT(3)
      ISTRT(1) = MAX (IBCDAT(4), IPSTRT)
      IEND (1) = MIN (IBCDAT(5), IPEND)
      ISTRT(2) = IBCDAT(6)
      IEND (2) = IBCDAT(7)
      ISTRT(3) = IBCDAT(8)
      IEND (3) = IBCDAT(9)
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCPT2 is the boundary condition location outside the block from IBCPT
C      IBCN   is the point to the inside of the block from IBCPT
C      IBCM   is the location of the metrics
C
C
      IF (INROUT .EQ. 1) THEN
         IBCPT  = 1
         IBCPT2 = 0
         IBCN   = 2
         IBCN2  = 3
         IBCM   = 1
         IBCEND = IJKDIM(IDIR) + 1
         IBCINC = 1
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCN2  = IJKDIM(IDIR) - 1
         IBCM   = IJKDIM(IDIR)
         IBCEND = 1
         IBCINC = - 1
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2. Loop over all points on the boundary surface and calculate skin friction,
C    heat transfer, and boundary layer quantities for output.
C
      DO 220 I2 = ISTRT(IC2), IEND(IC2)
         DO 210 I3 = ISTRT(IC3), IEND(IC3)
C
            IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            IN   = IBCN  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JN   = IBCN  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KN   = IBCN  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            IN2  = IBCN2 *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JN2  = IBCN2 *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KN2  = IBCN2 *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            RHOW = Q(JBC,KBC,IBC,1) 
            RHON = Q(JN ,KN ,IN ,1) 
            RHON2= Q(JN2,KN2,IN2,1) 
C
            PW   = Q(JBC,KBC,IBC,5) 
            PN   = Q(JN ,KN ,IN ,5) 
            PN2  = Q(JN2,KN2,IN2,5) 
C
C Calculate T = T / TINF  at the wall and cell center off the wall
C
            TBC  = GAMMA * PW  / RHOW
            TN   = GAMMA * PN  / RHON
            TN2  = GAMMA * PN2 / RHON2
C
C Calculate the velocity parallel to the wall at the two cell centers off 
C the wall (to allow for second order accurate expressions).
C UBAR = U(parallel) = U - UDOTN where U is the total velocity vector
C and UDOTN is the velocity at the cell center normal to the wall.
C (Velocity on the wall should be identically zero)
C
            UDOTN  = S(JM,KM,IM,1,IDIR) * Q(JN ,KN ,IN ,2)
     1             + S(JM,KM,IM,2,IDIR) * Q(JN ,KN ,IN ,3)
     2             + S(JM,KM,IM,3,IDIR) * Q(JN ,KN ,IN ,4)
            UBARN  = SQRT (
     1               (Q(JN ,KN ,IN ,2) - UDOTN *S(JM,KM,IM,1,IDIR))**2
     2             + (Q(JN ,KN ,IN ,3) - UDOTN *S(JM,KM,IM,2,IDIR))**2
     3             + (Q(JN ,KN ,IN ,4) - UDOTN *S(JM,KM,IM,3,IDIR))**2 )
C
            UDOTN2 = S(JM,KM,IM,1,IDIR) * Q(JN2,KN2,IN2,2)
     1             + S(JM,KM,IM,2,IDIR) * Q(JN2,KN2,IN2,3)
     2             + S(JM,KM,IM,3,IDIR) * Q(JN2,KN2,IN2,4)
            UBARN2 = SQRT (
     1               (Q(JN2,KN2,IN2,2) - UDOTN2*S(JM,KM,IM,1,IDIR))**2
     2             + (Q(JN2,KN2,IN2,3) - UDOTN2*S(JM,KM,IM,2,IDIR))**2
     3             + (Q(JN2,KN2,IN2,4) - UDOTN2*S(JM,KM,IM,3,IDIR))**2 )
C
C Viscosity (RMUW) at the wall and the cell center
C
            RMUW   = PROPS(JBC,KBC,IBC,1)
            RMUN   = PROPS(JN, KN ,IN ,1)
C
C Distances
C
            DN     = SQRT ( (RC(JN ,KN ,IN ,1) - RC(JBC,KBC,IBC,1))**2
     1                    + (RC(JN ,KN ,IN ,2) - RC(JBC,KBC,IBC,2))**2
     2                    + (RC(JN ,KN ,IN ,3) - RC(JBC,KBC,IBC,3))**2 )
C
            DN2    = SQRT ( (RC(JN2,KN2,IN2,1) - RC(JN ,KN ,IN ,1))**2
     1                    + (RC(JN2,KN2,IN2,2) - RC(JN ,KN ,IN ,2))**2
     2                    + (RC(JN2,KN2,IN2,3) - RC(JN ,KN ,IN ,3))**2 )
C
C Calculate first derivative terms using first or second order differencing
C
            DUDN   = (1-IHIORD) * (UBARN      / DN)
     1               +  IHIORD  * (- DN*DN*UBARN2 + (DN+DN2)**2*UBARN)
     2                            / (DN*DN2*(DN+DN2))
            DTDN   = (1-IHIORD) * ((TN - TBC) / DN)
     1               +  IHIORD  * (- DN*DN*TN2 + (DN+DN2)**2*TN 
     2                             - DN2*(2.E0*DN+DN2)*TBC)
     3                            / (DN*DN2*(DN+DN2))
C
C Pressure Coefficients Pw/Qinf =  P_w          / ( 1/2 * rhoinf * uinf ** 2 )
C                       Cp      = (P_w - P_inf) / ( 1/2 * rhoinf * uinf ** 2 )
C
            PWQINF =  PW         * 2.E0 / ( FSMACH * FSMACH )
            CP     = (PW - PREF) * 2.E0 / ( FSMACH * FSMACH )
C
C Find the maximum velocity magnitude in the profile
C
            UMAX   = 0.E0
CDIR$ IVDEP
            DO 100 I1 = IBCPT, IBCEND, IBCINC
C
               II   = I1    *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JJ   = I1    *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KK   = I1    *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               UMAG = SQRT ( Q(JJ,KK,II,2) * Q(JJ,KK,II,2) +
     1                       Q(JJ,KK,II,3) * Q(JJ,KK,II,3) +
     2                       Q(JJ,KK,II,4) * Q(JJ,KK,II,4) )
               UMAX = MAX  ( UMAX, UMAG )
  100       CONTINUE
C
C Calculate the edge velocity.  Based on Chris Rumsey's notes with the
C following assumptions: 
C      1. H   = H_inf 
C      2. P   = P_wall 
C      3. Isentropic from edge to wall to get density
C      4. TKE = TKE_inf
C If these assumptions result in imaginary UE, then set UE to UMAX
C
            HEDGE  = EINF + PREF / RHOINF
            PEDGE  = PW
            RHOEDG = RHOINF * (PEDGE / PREF)**(1.E0/GAMMA)
            TKEDGE = TKEINF
C
            ARGUE  = 2.E0*(HEDGE - GAMMA/GAMM1*PEDGE/RHOEDG - TKEDGE)
            IF (ARGUE .LE. 0.E0) THEN
               UE  = UMAX
            ELSE
               UE  = SQRT (ARGUE)
            ENDIF
C
C Make sure that UE is bounded by UMAX
C
            UE     = MIN (UE, UMAX)
C
C Skin Friction Cf      = Tau_wall / ( 1/2 * rho_inf  * u_inf  ** 2 )
C               Cf_edge = Tau_wall / ( 1/2 * rho_edge * u_edge ** 2 )
C
            IF ( IFWALF ) THEN
               TAUW   = PROPS (JBC2,KBC2,IBC2,ITQTAU)
            ELSE
               TAUW   = RMUW * DUDN
            ENDIF
            CF     = 2.E0 * TAUW / ( RE * FSMACH )
            CFEDGE = CF * (RHOINF/RHOEDG) * (FSMACH/UE)**2
C
C Heat Transfer Ch      = Qdot_wall / ( rho_inf  * u_inf  * (H_inf  - H_wall) )
C               Ch_edge = Qdot_wall / ( rho_edge * u_edge * (H_edge - H_wall) )
C
            IF ( IFWALF ) THEN
               QW     = PROPS (JBC2,KBC2,IBC2,ITQDOT)
            ELSE
               QW     = RMUW * DTDN / PR
            ENDIF
            DENOM = RE*(1.E0-TWALL/TINF+0.5E0*GAMM1*FSMACH*FSMACH)
C
C Check to see if can be non-dimensionalized this way
C
            IF (DENOM .EQ. 0.E0) THEN
               CH     = RMAX
               CHEDGE = RMAX
            ELSE
               CH     = QW / DENOM
               CHEDGE = CH * (RHOINF/RHOEDG) * (FSMACH/UE)**2
            ENDIF
C
C Ustar / ainf = sqrt (Tauwall / rho) * sqrt (Minf / Re)
C
            USTAR = SQRT (ABS(TAUW) / RHOW) * SQRT (FSMACH / RE)
C
C Yplus = y ustar / nu : ustar = sqrt (Tauwall / rho) 
C
            YPLUS = PROPS(JN ,KN ,IN ,3)
C
C Calculate boundary layer data:
C 1. Find the boundary layer edge defined as:
C    i.   location where U = .995 U_edge
C 2. Integrate using trapezoidal rule to find:
C    i.   boundary layer thickness (DELTA)
C    ii.  displacement thickness   (DSTAR)
C    iii. momentum thickness       (THETA)
C Formulas for DELTA, DSTAR, THETA from White, 'Viscous Fluid Flow', p.606
C
C JEDGE is the last cell-center before the boundary layer edge
C
            UEDGE  = .995E0 * UE
            DELTA  = 0.E0
            DSTAR  = 0.E0
            THETA  = 0.E0
CDIR$ IVDEP
            DO 110 I1 = IBCPT, IBCEND, IBCINC
C
               II    = I1    *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JJ    = I1    *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KK    = I1    *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IP    = II + IBCINC*IKD(IC1,1)
               JP    = JJ + IBCINC*IKD(IC1,2)
               KP    = KK + IBCINC*IKD(IC1,3)
C
               U     = SQRT ( Q(JJ,KK,II,2) * Q(JJ,KK,II,2) + 
     1                        Q(JJ,KK,II,3) * Q(JJ,KK,II,3) + 
     2                        Q(JJ,KK,II,4) * Q(JJ,KK,II,4) ) 
               UP    = SQRT ( Q(JP,KP,IP,2) * Q(JP,KP,IP,2) + 
     1                        Q(JP,KP,IP,3) * Q(JP,KP,IP,3) + 
     2                        Q(JP,KP,IP,4) * Q(JP,KP,IP,4) ) 
               IF (U  .LT. UEDGE .AND. UP .GE. UEDGE) THEN
                  I1EDGE = I1
                  GO TO 119
               ENDIF
C
               DY    = SQRT ( (RC(JP,KP,IP,1)-RC(JJ,KK,II,1))**2
     1                      + (RC(JP,KP,IP,2)-RC(JJ,KK,II,2))**2
     2                      + (RC(JP,KP,IP,3)-RC(JJ,KK,II,3))**2 )
C
               U     = U  / FSMACH
               UP    = UP / FSMACH
               RHOU  = Q(JJ,KK,II,1) * U
               RHOUP = Q(JP,KP,IP,1) * UP
C
               DELTA = DELTA + DY
               DSTAR = DSTAR + 0.5E0 * ( (1.E0 - RHOU) 
     1                                 + (1.E0 - RHOUP) )*DY
               THETA = THETA + 0.5E0 * ( (1.E0 - U) * RHOU
     1                                 + (1.E0 - UP)* RHOUP )*DY
  110       CONTINUE
C
C Could not find bl edge.  Output error.
C
            DELTA  = - RMAX
            DSTAR  = - RMAX
            THETA  = - RMAX
C
  119       CONTINUE
C
C Output Cf, Ch, Cp, Delta, Delta*, Theta
C 
  209       CONTINUE
            WRITE (IOUT,1010) IBC, JBC, KBC, RC(JBC,KBC,IBC,1), 
     1                        RC(JBC,KBC,IBC,2), RC(JBC,KBC,IBC,3), 
     2                        PWQINF, CP, CF, CFEDGE, CH, CHEDGE, DN, 
     3                        YPLUS, USTAR, DELTA, DSTAR, THETA, UE
  210    CONTINUE
  220 CONTINUE
C
C Finished with Skin Friction and Heat Transfer
C
      RETURN
 1000 FORMAT (//' ',' Surface Data for BLOCK ',I5,10X,'(',A,
     1              ' Order accurate derivatives)')
 1005 FORMAT ( /' ',3X,'I',3X,'J',3X,'K',6X,'X',12X,'Y',12X,'Z',9X,
     1              'Pw/Qinf',8X,'Cp',11X,'Cf',9X,'Cf_edge',8X, 
     2              'Ch',9X,'Ch_edge',8X,'DN',10X,'Yplus',8X,'Ustar',8X,
     3              'Delta',8X,'Delta*',7X,'Theta',8X,'U_edge'/)
 1010 FORMAT (' ',3I4,16(1PE13.5))
      END




      SUBROUTINE xSFHTJ (IDIM, JDIM, KDIM, IPSTRT, IPEND, ITURB,
     1                  Q, PROPS, S, VOL, IBCDAT, RC, IBLOCK)
C
C Routine to calculate the skin friction and heat transfer using 
C the thin layer Navier-Stokes assumption for the J direction.
C
C This routine is not presently suited for the multiblock
C approach.  It is close to being ready but is not.  The problem
C lies that the metrics for one cell into the adjacent block
C are needed to be completely general.  These are currently not available.
C The routine metbc needs to be completed.
C
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C ITURB          : Turbulence model
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C S              : Metrics
C VOL            : Cell volumes
C IBCDAT         : Data controlling boundary conditions:
C                  IBCDAT(1,ibc) = Boundary Condition Type
C                         2      = Direction (1=i, 2=j, 3=k)
C                         3      = i(1 or idim)
C                         4      = istrt
C                         5      = iend
C                         6      = jstrt
C                         7      = jend
C                         8      = kstrt
C                         9      = kend
C RC             : Cell center grid
C IBLOCK         : Block number to identify output
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      DIMENSION IBCDAT(10)
C
C Store boundary condition data in mnemonic names
C
      IBCTYP = IBCDAT(1)
      IDIR   = IBCDAT(2)
      INROUT = IBCDAT(3)
      ISTRT  = MAX (IBCDAT(4), IPSTRT)
      IEND   = MIN (IBCDAT(5), IPEND)
      JSTRT  = IBCDAT(6)
      JEND   = IBCDAT(7)
      KSTRT  = IBCDAT(8)
      KEND   = IBCDAT(9)
c
      WRITE (IOUT,1000) IBLOCK
      WRITE (IOUT,1005) 
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBC  is the boundary condition location
C      IBC2 is the boundary condition location outside the block from IBC
C      IN   is the point to the inside of the block from IBC
C      IM   is the location of the metrics
C
C
C J Direction
C
      IF (IDIR .EQ. 2) THEN
         IF (INROUT .EQ. 1) THEN
            JBC  = 1
            JBC2 = 0
            JN   = 2
            JM   = 1
            JMP1 = 2
            JBEG = 1
            JEND = JDIM + 1
            JINC = 1
         ELSE
            JBC  = JDIM + 1
            JBC2 = JDIM + 2
            JN   = JDIM
            JM   = JDIM
            JMP1 = JDIM - 1
            JBEG = JDIM + 1
            JEND = 1
            JINC = -1
         ENDIF
C
C Loop over the J = JBC plane
C
         DO 220 I = ISTRT, IEND
            DO 210 K = KSTRT, KEND
C
C Calculate T = T / TINF  at the wall and cell center off the wall
C
               TBC = GAMMA * Q(JBC,K,I,5) / Q(JBC,K,I,1)
               TN  = GAMMA * Q(JN,K,I,5) / Q(JN,K,I,1)
C
C Calculate the velocity parallel to the wall at 
C the cell center off the wall.
C U(parallel) = U - UBAR where U is the total velocity vector
C and UBAR is the velocity normal to the wall.
C (Velocity on the wall should be identically zero)
C
               UBARN = S(JM,K,I,1,2) * Q(JN,K,I,2)
     1               + S(JM,K,I,2,2) * Q(JN,K,I,3)
     2               + S(JM,K,I,3,2) * Q(JN,K,I,4)
               VELN  = SQRT ( (Q(JN,K,I,2) - UBARN*S(JM,K,I,1,2))**2
     1                      + (Q(JN,K,I,3) - UBARN*S(JM,K,I,2,2))**2
     2                      + (Q(JN,K,I,4) - UBARN*S(JM,K,I,3,2))**2 )
C
C Viscosity (RMUW) at the wall and the cell center
C
               RMUW = PROPS(JBC,K,I,1)
               RMUN = PROPS(JN, K,I,1)
C
C Average metrics and cell volume
C
               VOLAUX = 0.5E0 * VOL(JMP1,K,I)
               SPJ1   = 0.5E0 * ( S(JMP1,K,I,1,2)*S(JMP1,K,I,4,2)
     1                          + S(JM,K,I,1,2)*S(JM,K,I,4,2) )
               SPJ2   = 0.5E0 * ( S(JMP1,K,I,2,2)*S(JMP1,K,I,4,2)
     1                          + S(JM,K,I,2,2)*S(JM,K,I,4,2) )
               SPJ3   = 0.5E0 * ( S(JMP1,K,I,3,2)*S(JMP1,K,I,4,2)
     1                          + S(JM,K,I,3,2)*S(JM,K,I,4,2) )
               SMJ1   = S(JM,K,I,1,2)*S(JM,K,I,4,2)
               SMJ2   = S(JM,K,I,2,2)*S(JM,K,I,4,2)
               SMJ3   = S(JM,K,I,3,2)*S(JM,K,I,4,2)
C
C Calculate first derivative terms using Divergence Theorem
C
               DVELDX = VELN*SPJ1 / VOLAUX
               DVELDY = VELN*SPJ2 / VOLAUX
               DVELDZ = VELN*SPJ3 / VOLAUX
               DTDX   = (TN*SPJ1 - TBC*SMJ1) / VOLAUX
               DTDY   = (TN*SPJ2 - TBC*SMJ2) / VOLAUX
               DTDZ   = (TN*SPJ3 - TBC*SMJ3) / VOLAUX
C
C Pressure Coefficients Pw/Qinf =  P_w          / ( 1/2 * rhoinf * uinf ** 2 )
C                       Cp      = (P_w - P_inf) / ( 1/2 * rhoinf * uinf ** 2 )
C
               PW     = Q(JBC,K,I,5) 
               PWQINF =  PW         * 2.E0 / ( FSMACH * FSMACH )
               CP     = (PW - PREF) * 2.E0 / ( FSMACH * FSMACH )
C
C Calculate the edge velocity.  Based on Chris Rumsey's notes with the
C following assumptions: 
C      1. H   = H_inf 
C      2. P   = P_wall 
C      3. Isentropic from edge to wall to get density
C      4. TKE = TKE_inf
C
               HEDGE  = EINF + PREF / RHOINF
               PEDGE  = PW
               RHOEDG = RHOINF * (PEDGE / PREF)**(1.E0/GAMMA)
               TKEDGE = TKEINF
C
               ARGUE  = 2.E0*(HEDGE - GAMMA/GAMM1*PEDGE/RHOEDG - TKEDGE)
               IF (ARGUE .LE. 0.E0) THEN
                  UE     = FSMACH
               ELSE
                  UE     = SQRT (ARGUE)
               ENDIF
C
C Skin Friction Cf      = Tau_wall / ( 1/2 * rho_inf  * u_inf  ** 2 )
C               Cf_edge = Tau_wall / ( 1/2 * rho_edge * u_edge ** 2 )
C
               DVELDN = DVELDX * S(JM,K,I,1,2) + DVELDY * S(JM,K,I,2,2)
     1                + DVELDZ * S(JM,K,I,3,2)
               IF ( IFWALF ) THEN
                  TAUW   = PROPS (JBC2,K,I,ITQTAU)
               ELSE
                  TAUW   = RMUW * DVELDN
               ENDIF
               CF     = 2.E0 * TAUW / ( RE * FSMACH )
               CFEDGE = CF * (RHOINF/RHOEDG) * (FSMACH/UE)**2
C
C Heat Transfer Ch      = Qdot_wall / ( rho_inf  * u_inf  * (H_inf  - H_wall) )
C               Ch_edge = Qdot_wall / ( rho_edge * u_edge * (H_edge - H_wall) )
C
               DTDN  = DTDX * S(JM,K,I,1,2) + DTDY * S(JM,K,I,2,2)
     1               + DTDZ * S(JM,K,I,3,2)
               IF ( IFWALF ) THEN
                  QW     = PROPS (JBC2,K,I,ITQDOT)
               ELSE
                  QW     = RMUW * DTDN / PR
               ENDIF
               DENOM = RE*(1.E0-TWALL/TINF+0.5E0*GAMM1*FSMACH*FSMACH)
C
C Check to see if can be non-dimensionalized this way
C
               IF (DENOM .EQ. 0.E0) THEN
                  CH     = RMAX
                  CHEDGE = RMAX
               ELSE
                  CH     = QW / DENOM
                  CHEDGE = CH * (RHOINF/RHOEDG) * (FSMACH/UE)**2
               ENDIF
C
C Ustar / ainf = sqrt (Tauwall / rho) * sqrt (Minf / Re)
C
               USTAR = SQRT (TAUW / Q(JBC,K,I,1)) *
     1                 SQRT (FSMACH / RE)
C
C Yplus = y ustar / nu : ustar = sqrt (Tauwall / rho) 
C
               DN    = VOLAUX / (0.5E0*(S(JMP1,K,I,4,2)+S(JM,K,I,4,2)))
               YPLUS = PROPS(JN,K,I,3)
C
C Find the boundary layer edge defined as:
C      1. Location where U = .995 U_edge
C JEDGE is the last cell-center before the boundary layer edge
C
               UE     = .995E0 * UE
               DO 110 J = JBEG, JEND-JINC, JINC
                  U  = SQRT ( Q(J     ,K,I,2)*Q(J     ,K,I,2) + 
     1                        Q(J     ,K,I,3)*Q(J     ,K,I,3) + 
     2                        Q(J     ,K,I,4)*Q(J     ,K,I,4) ) 
                  UP = SQRT ( Q(J+JINC,K,I,2)*Q(J+JINC,K,I,2) + 
     1                        Q(J+JINC,K,I,3)*Q(J+JINC,K,I,3) + 
     2                        Q(J+JINC,K,I,4)*Q(J+JINC,K,I,4) ) 
                  IF (U  .LT. UE .AND. UP .GE. UE) THEN
                     JEDGE = J
                     GO TO 119
                  ENDIF
  110          CONTINUE
C
C Could not find bl edge.  Output error.
C
               DELTA = - RMAX
               DSTAR = - RMAX
               THETA = - RMAX
               GO TO 209
C
C Integrate to find:
C    the boundary layer thickness (DELTA),
C    the displacement thickness (DSTAR),
C    the momentum thickness (THETA),
C using the trapezoidal rule.  Formulas for DELTA, DSTAR, THETA
C given on page 606 of White's 'Viscous Fluid Flow' for compressible fluids.
C
  119          CONTINUE
               DELTA = 0.E0
               DSTAR = 0.E0
               THETA = 0.E0
               DO 120 J = JBEG, JEDGE-JINC, JINC
                  DY    = SQRT ( (RC(J+JINC,K,I,1)-RC(J,K,I,1))**2
     1                         + (RC(J+JINC,K,I,2)-RC(J,K,I,2))**2
     2                         + (RC(J+JINC,K,I,3)-RC(J,K,I,3))**2 )
                  U     = SQRT (Q(J     ,K,I,2)*Q(J     ,K,I,2) + 
     1                          Q(J     ,K,I,3)*Q(J     ,K,I,3) + 
     2                          Q(J     ,K,I,4)*Q(J     ,K,I,4) )/FSMACH
                  RHOU  = Q(J,K,I,1) * U
                  UP    = SQRT (Q(J+JINC,K,I,2)*Q(J+JINC,K,I,2) + 
     1                          Q(J+JINC,K,I,3)*Q(J+JINC,K,I,3) + 
     2                          Q(J+JINC,K,I,4)*Q(J+JINC,K,I,4) )/FSMACH
                  RHOUP = Q(J+JINC,K,I,1) * UP
                  DELTA = DELTA + DY
                  DSTAR = DSTAR + 0.5E0 * ( (1.E0 - RHOU) 
     1                                    + (1.E0 - RHOUP) )*DY
                  THETA = THETA + 0.5E0 * ( (1.E0 - U) * RHOU
     1                                    + (1.E0 - UP)* RHOUP )*DY
  120          CONTINUE
C
C Add contribution from JEDGE to DELTA to DELTA, DSTAR and THETA
C
               UJE   = SQRT ( Q(JEDGE     ,K,I,2)*Q(JEDGE     ,K,I,2) + 
     1                        Q(JEDGE     ,K,I,3)*Q(JEDGE     ,K,I,3) + 
     2                        Q(JEDGE     ,K,I,4)*Q(JEDGE     ,K,I,4) )
               UJEP  = SQRT ( Q(JEDGE+JINC,K,I,2)*Q(JEDGE+JINC,K,I,2) + 
     1                        Q(JEDGE+JINC,K,I,3)*Q(JEDGE+JINC,K,I,3) + 
     2                        Q(JEDGE+JINC,K,I,4)*Q(JEDGE+JINC,K,I,4) )
               FACT  = (UE - UJE) / (UJEP - UJE)
               DY    = FACT *
     1              SQRT ( (RC(JEDGE+JINC,K,I,1)-RC(JEDGE,K,I,1))**2
     2                   + (RC(JEDGE+JINC,K,I,2)-RC(JEDGE,K,I,2))**2
     3                   + (RC(JEDGE+JINC,K,I,3)-RC(JEDGE,K,I,3))**2 )
               RHOU  = Q(JEDGE,K,I,1) * UJE / FSMACH
               RHOE  = (1.E0 - FACT)*Q(JEDGE,K,I,1) 
     1                 + FACT*Q(JEDGE+JINC,K,I,1)
               RHOUP = RHOE * UE / FSMACH
               DELTA = DELTA + DY
               DSTAR = DSTAR + 0.5E0 * ( (1.E0 - RHOU) 
     1                                 + (1.E0 - RHOUP) )*DY
               THETA = THETA + 0.5E0 * ( (1.E0 - RHOU) * RHOU
     1                                 + (1.E0 - RHOUP)* RHOUP )*DY
C
C Print data
C 
  209          CONTINUE
               WRITE (IOUT,1010) I, JBC, K, RC(JBC,K,I,1), 
     1                           RC(JBC,K,I,2), RC(JBC,K,I,3), 
     2                           PWQINF, CP, CF, CFEDGE, CH, CHEDGE, DN, 
     3                           YPLUS, USTAR, DELTA, DSTAR, THETA
  210       CONTINUE
  220    CONTINUE
      ENDIF
C
C Finished with Skin Friction and Heat Transfer
C
      RETURN
 1000 FORMAT (//' ',' Surface Data for BLOCK ',I5)
 1005 FORMAT ( /' ',3X,'I',3X,'J',3X,'K',6X,'X',12X,'Y',12X,'Z',9X,
     1              'Pw/Qinf',8X,'Cp',11X,'Cf',9X,'Cf_edge',8X,
     2              'Ch',9X,'Ch_edge',8X,'DN',10X,'Yplus',8X,
     3              'Ustar',8X,'Delta',8X,'Delta*',7X,'Theta'/)
 1010 FORMAT (' ',3I4,15(1PE13.5))
      END
