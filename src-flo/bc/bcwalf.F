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
      SUBROUTINE BCWALF (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, VOL, IBCDAT)
C
C Routine to do TURBULENT WALL FUNCTION boundary condition.
C
C This routine is organized as follows:
C     1. Setup constants for turbulent wall functions
C     2. Setup variables for looping
C     3. Set boundary conditions for the mean variables
C     4. Set additional boundary conditions for turbulence variables:
C        a. Two-Equation variables
C           i.  K-epsilon
C           ii. K-omega
C        b. Reynolds stress variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
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
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3),
     3          ISTRTI(3),  
     4          IENDI (3)
C
C 1. Control of iterative process
C
      MAXIT  = 50
C
C 32-Bit machine can't converge 1x10-7
C
      TOLER  = 1.E-06
C
C No underrelaxation results in more planes not converging,
C underrelaxation with RELAX = .9 converges more rapidly (11 to 20 cycles)
C than none but some planes still do not converge (get in limit cycle)
C RELAX = .5 requires 14 iterations instead of 11 but is more robust
C
      RELAX  = 0.5E0
C
C Constants for wall functions
C
      RKAPVK = 0.41E0
      ESTAR  = 9.0E0
      YV     = 11.6E0
      SCALE  = SQRT (RE / FSMACH)
C
C 2. Setup memory for looping.
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
C ISTRTI, IENDI correspond to the Interior of the boundary condition face
C i.e., this does not extend to the edge/corner 
C
      DO 5 L = 1, 3
         ISTRTI(L) = ISTRT(L)
         IENDI (L) = IEND (L)
    5 CONTINUE
C
C Modify beginning, ending indices to extend boundary condition to edge/corner
C Reset kstrt, kend to 2 for improved efficiency on two-dimensional flows.
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
      IF (.NOT. THREED .AND. .NOT. AXISYM) THEN
         ISTRT (3) = 2
         IEND  (3) = 2
         ISTRTI(3) = 2
         IENDI (3) = 2
      ENDIF
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
         IBCM2  = 2
         SGNM   = 1.E0
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCN2  = IJKDIM(IDIR) - 1
         IBCM   = IJKDIM(IDIR)
         IBCM2  = IJKDIM(IDIR) - 1
         SGNM   = - 1.E0
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2-b. Set the boundary conditions on the edges/corners 
C
      IF (ISTRT(IC2) .EQ. 1) THEN
         I2 = ISTRT(IC2)
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
            RHON = Q(JN ,KN ,IN ,1)
            PN   = Q(JN ,KN ,IN ,5)
C
            Q(JBC ,KBC ,IBC ,1) = RHON
            Q(JBC ,KBC ,IBC ,2) = 0.E0
            Q(JBC ,KBC ,IBC ,3) = 0.E0
            Q(JBC ,KBC ,IBC ,4) = 0.E0
            Q(JBC ,KBC ,IBC ,5) = PN
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
C
            PROPS (JBC2,KBC2,IBC2,ITQTAU) = 0.E0
            PROPS (JBC2,KBC2,IBC2,ITQDOT) = 0.E0
  210    CONTINUE
      ENDIF
C
      IF (IEND (IC2) .EQ. (IJKDIM(IC2)+1)) THEN
         I2 = IEND (IC2)
         DO 220 I3 = ISTRT(IC3), IEND(IC3)
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
            RHON = Q(JN ,KN ,IN ,1)
            PN   = Q(JN ,KN ,IN ,5)
C
            Q(JBC ,KBC ,IBC ,1) = RHON
            Q(JBC ,KBC ,IBC ,2) = 0.E0
            Q(JBC ,KBC ,IBC ,3) = 0.E0
            Q(JBC ,KBC ,IBC ,4) = 0.E0
            Q(JBC ,KBC ,IBC ,5) = PN
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
C
            PROPS (JBC2,KBC2,IBC2,ITQTAU) = 0.E0
            PROPS (JBC2,KBC2,IBC2,ITQDOT) = 0.E0
  220    CONTINUE
      ENDIF
C
      IF (ISTRT(IC3) .EQ. 1) THEN
         I3 = ISTRT(IC3)
         DO 230 I2 = ISTRT(IC2), IEND(IC2)
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
            RHON = Q(JN ,KN ,IN ,1)
            PN   = Q(JN ,KN ,IN ,5)
C
            Q(JBC ,KBC ,IBC ,1) = RHON
            Q(JBC ,KBC ,IBC ,2) = 0.E0
            Q(JBC ,KBC ,IBC ,3) = 0.E0
            Q(JBC ,KBC ,IBC ,4) = 0.E0
            Q(JBC ,KBC ,IBC ,5) = PN
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
C
            PROPS (JBC2,KBC2,IBC2,ITQTAU) = 0.E0
            PROPS (JBC2,KBC2,IBC2,ITQDOT) = 0.E0
  230    CONTINUE
      ENDIF
C
      IF (IEND (IC3) .EQ. (IJKDIM(IC3)+1)) THEN
         I3 = IEND (IC3)
         DO 240 I2 = ISTRT(IC2), IEND(IC2)
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
            RHON = Q(JN ,KN ,IN ,1)
            PN   = Q(JN ,KN ,IN ,5)
C
            Q(JBC ,KBC ,IBC ,1) = RHON
            Q(JBC ,KBC ,IBC ,2) = 0.E0
            Q(JBC ,KBC ,IBC ,3) = 0.E0
            Q(JBC ,KBC ,IBC ,4) = 0.E0
            Q(JBC ,KBC ,IBC ,5) = PN
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
C
            PROPS (JBC2,KBC2,IBC2,ITQTAU) = 0.E0
            PROPS (JBC2,KBC2,IBC2,ITQDOT) = 0.E0
  240    CONTINUE
      ENDIF
C
C 3. Boundary conditions for Mean Variables using Wall Functions
C    Compute using the wall functions only in the INTERIOR of the boundary 
C
      DO 320 I3 = ISTRTI(IC3), IENDI(IC3)
CDIR$ IVDEP
         DO 310 I2 = ISTRTI(IC2), IENDI(IC2)
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
            IM2  = IBCM2 *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JM2  = IBCM2 *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KM2  = IBCM2 *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
c---------- RHOW   = Q(JBC,KBC,IBC,1)
            RHON   = Q(JN ,KN ,IN ,1)
            PN     = Q(JN ,KN ,IN ,5)
C
C Calculate velocity magnitude parallel to the wall
C
            UDOTN  = ( Q(JN,KN,IN,2) * S(JM,KM,IM,1,IDIR) +
     1                 Q(JN,KN,IN,3) * S(JM,KM,IM,2,IDIR) +
     2                 Q(JN,KN,IN,4) * S(JM,KM,IM,3,IDIR) ) * SGNM
            UPN    = Q(JN,KN,IN,2) - UDOTN * SGNM * S(JM,KM,IM,1,IDIR)
            VPN    = Q(JN,KN,IN,3) - UDOTN * SGNM * S(JM,KM,IM,2,IDIR)
            WPN    = Q(JN,KN,IN,4) - UDOTN * SGNM * S(JM,KM,IM,3,IDIR)
            UPARN  = SQRT( UPN*UPN + VPN*VPN + WPN*WPN )
C
            UDOTN2 = ( Q(JN2,KN2,IN2,2) * S(JM2,KM2,IM2,1,IDIR) +
     1                 Q(JN2,KN2,IN2,3) * S(JM2,KM2,IM2,2,IDIR) +
     2                 Q(JN2,KN2,IN2,4) * S(JM2,KM2,IM2,3,IDIR) ) * SGNM
            UPN2   = Q(JN2,KN2,IN2,2)-UDOTN2*SGNM*S(JM2,KM2,IM2,1,IDIR)
            VPN2   = Q(JN2,KN2,IN2,3)-UDOTN2*SGNM*S(JM2,KM2,IM2,2,IDIR)
            WPN2   = Q(JN2,KN2,IN2,4)-UDOTN2*SGNM*S(JM2,KM2,IM2,3,IDIR)
            UPARN2 = SQRT( UPN2*UPN2 + VPN2*VPN2 + WPN2*WPN2 )
C
            TN     = GAMMA * PN / RHON
C
            RMUW   = PROPS(JBC,KBC,IBC,1)
            RMUN   = PROPS(JN ,KN ,IN ,1)
            YNORML = ABS (PROPS(JN ,KN ,IN ,4))
            YNORML = MAX (YNORML, RSMALL)
C
C Calculate normal Reynolds stress
C
            IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
               TKEN   = Q(JN,KN,IN,6)
               TAUNN  = 2.E0 / 3.E0 * TKEN
            ELSE IF (ITURB .EQ. ITRS) THEN
C
               RNX    = S(JM,KM,IM,1,IDIR) * SGNM
               RNY    = S(JM,KM,IM,2,IDIR) * SGNM
               RNZ    = S(JM,KM,IM,3,IDIR) * SGNM
C
               TAUXX  = Q(JN,KN,IN, 6)
               TAUYY  = Q(JN,KN,IN, 7)
               TAUZZ  = Q(JN,KN,IN, 8)
               TAUXY  = Q(JN,KN,IN, 9)
               TAUXZ  = Q(JN,KN,IN,10)
               TAUYZ  = Q(JN,KN,IN,11)
C
C Transform stresses to local streamwise, normal coordinates
C
               TAUNN  =   RNX*RNX*TAUXX + RNY*RNY*TAUYY + RNZ*RNZ*TAUZZ 
     1            + 2.E0*(RNX*RNY*TAUXY + RNX*RNZ*TAUXZ + RNY*RNZ*TAUYZ)
C
            ELSE
               TAUNN  = 0.E0
            ENDIF
c+++++++
c+++        IF (S(JM,KM,IM,4,IDIR) .LE. 0.E0 .AND.
c+++ 1          VOL(JM,KM,IM)        .LE. 0.E0) THEN
c+++           WRITE (iout,'(a,3i5)') ' i, j, k ', im, jm, km
c+++        END IF
c+++++++
C
C Wall pressure boundary condition and wall density
C
            PW     = PN + RHON * TAUNN
            IF ( TWSPE ) THEN
               TW    = TWALL / TINF
            ELSE
               TW    = TN + 0.5E0 * UPARN * UPARN * GAMM1 * PRT
            ENDIF
            RHOW = GAMMA * PW / TW
C
C Initialize Tau_w for iterative procedure to previous result
C
            TAUW   = PROPS (JBC2,KBC2,IBC2,ITQTAU)
C
C Iterate on tau_w
C
            DO 300 ITER = 1, MAXIT
               YPLUSN = SCALE * YNORML * SQRT (ABS(TAUW) * RHOW) / RMUW
C
C Case 1: Y+_N > YV ==> u+ = Ln(y+) / kappa + B
C
               IF (YPLUSN .GT. YV) THEN
                  RMUTN = SCALE * RKAPVK * YNORML * SQRT(ABS(TAUW)*RHON)
                  IF ( TWSPE ) THEN
                     TW    = TWALL / TINF
                     QW    = TAUW * ( (TW - TN) / (PRT * UPARN + RSMALL) 
     1                              - 0.5E0 * GAMM1 * UPARN ) 
                  ELSE
                     QW    = 0.E0
                     TW    = TN + 0.5E0 * UPARN * UPARN * GAMM1 * PRT
                  ENDIF
                  RHOW = GAMMA * PW / TW
C
C Calculate U_C_N
C
                  DEN1  = SQRT( QW*QW + 2.E0*TW/(GAMM1*PRT)*TAUW*TAUW )
                  DEN1  = MAX (DEN1, 1.e-08)
                  ARG11 = (QW + UPARN*TAUW) / DEN1
                  ARG12 = (QW) / DEN1
                  ARG11 = MIN (  1.E0, ARG11 )
                  ARG11 = MAX ( -1.E0, ARG11 )
                  ARG12 = MIN (  1.E0, ARG12 )
                  ARG12 = MAX ( -1.E0, ARG12 )
                  UCN   = SQRT( 2.E0*TW/(GAMM1*PRT) ) 
     1                    * ( ASIN(ARG11) - ASIN(ARG12) )
C
C Calculate Tau_w
C
                  TAUWN = RHOW*(RKAPVK*UCN/LOG(ESTAR*YPLUSN))**2
     1                    *RE/FSMACH
C
C Case 2: Y+_N < YV ==> u+ = y+
C
               ELSE
                  IF ( TWSPE ) THEN
                     TW    = TWALL / TINF
                     QW    = TAUW * ( (TW - TN) / (PRT * UPARN + RSMALL) 
     1                              - 0.5E0 * GAMM1 * UPARN ) 
                  ELSE
                     QW    = 0.E0
                     TW    = TN + 0.5E0 * UPARN * UPARN * GAMM1 * PRT
                  ENDIF
                  RHOW = GAMMA * PW / TW
C
C Calculate Tau_w
C
                  TAUWN = (UPARN * RMUN / YNORML)
               ENDIF
C
C Update Tauw
C
               ERROR = TAUWN - TAUW
               TAUW  = (1.E0 - RELAX) * TAUW + RELAX * TAUWN
C
C Converged?
C
               RELERR = ABS(ERROR) / (ABS(TAUW) + 1.E-10)
               IF (RELERR .LT. TOLER) GO TO 305
  300       CONTINUE
            WRITE (IOUT,1010) IBC, JBC, KBC, MAXIT, RELERR, TOLER,YPLUSN
  305       CONTINUE
C
C Pass back the wall shear stress and heat transfer
C Note: this will always develop TAUW >= 0.  The direction of TAUW is
C       assumed to align with the local velocity direction.  SGNM is
C       required here to handle DIM end of coordinate direction.
C
            PROPS (JBC2,KBC2,IBC2,ITQTAU) = TAUW
            PROPS (JBC2,KBC2,IBC2,ITQDOT) = QW
C
C Calculate wall velocity so that du/dy at wall-function cell is from log-law.
C Velocity at the interface of the first cell (wall-function cell) is
C calculated as the average of the values on either side of the interface
C since the solution point in the wall-function cell is equidistant from
C the interface.  This calculates the magnitude of the wall velocity.  Its
C direction is taken as aligned with the velocity parallel to the wall.
C Therefore, UINT, DUDY2, DYINT are all positive quantities.
C
            UINT   = 0.5E0 * ( UPARN + UPARN2 )
            UTAUW  = SQRT (ABS(TAUW) / RHOW * FSMACH / RE)
            DUDY2  = UTAUW / (RKAPVK * YNORML)
            DYINT  = 0.5E0 * (YNORML + PROPS(JN2,KN2,IN2,4))
            UPWALL = UINT - DUDY2 * DYINT
C
C Set wall mean quantities
C
            UPARN  = MAX ( UPARN, RSMALL )
            Q(JBC ,KBC ,IBC ,1) = RHOW
            Q(JBC ,KBC ,IBC ,2) = UPN * UPWALL / UPARN
            Q(JBC ,KBC ,IBC ,3) = VPN * UPWALL / UPARN
            Q(JBC ,KBC ,IBC ,4) = WPN * UPWALL / UPARN
            Q(JBC ,KBC ,IBC ,5) = PW
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
  310    CONTINUE
  320 CONTINUE
C
C 4-a. Boundary conditions for Two-Equation variables
C
C 4-a-i. Boundary conditions for k-epsilon
C
      IF (ITURB .EQ. ITKE) THEN
         DO 420 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 410 I2 = ISTRT(IC2), IEND(IC2)
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
               RHOW   = Q    (JBC ,KBC ,IBC ,1)
C
               TAUW   = PROPS(JBC2,KBC2,IBC2,ITQTAU)
               YNORML = ABS (PROPS(JN  ,KN  ,IN  ,4))
               YNORML = MAX ( YNORML, RSMALL )
C
C Boundary condition for tke     = abs (tau_w) / rho / sqrt(C_mu)
C Boundary condition for epsilon = (sqrt(C_mu) k)^(3/2) / kappa y
C
               CMURT  = SQRT (CMU)
               TKEN   = ABS (TAUW) / RHOW / CMURT * FSMACH / RE
               TKEN   = MAX (TKEN, RSMALL)
C
               TERM   = CMURT * TKEN
               EPSLNN = TERM * SQRT (TERM) / (RKAPVK * YNORML)
               EPSLNN = MAX (EPSLNN, RSMALL)
               EPSLNW = EPSLNN
C
               Q(JN  ,KN  ,IN  ,6) = TKEN
               Q(JN  ,KN  ,IN  ,7) = EPSLNW
C
               Q(JBC ,KBC ,IBC ,6) = 0.E0
               Q(JBC ,KBC ,IBC ,7) = EPSLNW
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  410       CONTINUE
  420    CONTINUE
C
C 4-a-ii. Boundary conditions for k-omega
C
      ELSE IF (ITURB .EQ. ITKW) THEN
      write (iout,'(a)') ' Wall functions not implemented for k-w'
      stop
         DO 440 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 430 I2 = ISTRT(IC2), IEND(IC2)
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
               IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               RHOW   = Q    (JBC,KBC,IBC,1)
               RMUW   = PROPS(JBC,KBC,IBC,1)
C
               YNORML = ABS (PROPS(JN ,KN ,IN ,4))
C
  430       CONTINUE
  440    CONTINUE
C
C 4-b. Boundary conditions for Reynolds stress variables
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 460 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 450 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC   = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC   = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC   = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2  = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2  = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2  = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IN    = IBCN  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JN    = IBCN  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KN    = IBCN  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IN2   = IBCN2 *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JN2   = IBCN2 *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KN2   = IBCN2 *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IM    = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JM    = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KM    = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
C Set up normals for local coordinate system (s=SP, n=RN, t) where
C n is normal to the surface, s is parallel to the surface aligned with
C the flow, and t is perpendicular to s, n (and is not needed)
C
               RNX    = S(JM,KM,IM,1,IDIR) * SGNM
               RNY    = S(JM,KM,IM,2,IDIR) * SGNM
               RNZ    = S(JM,KM,IM,3,IDIR) * SGNM
C
               UDOTN  = Q(JN,KN,IN,2) * RNX + Q(JN,KN,IN,3) * RNY +
     1                  Q(JN,KN,IN,4) * RNZ
               UPN    = Q(JN,KN,IN,2) - UDOTN * RNX
               VPN    = Q(JN,KN,IN,3) - UDOTN * RNY
               WPN    = Q(JN,KN,IN,4) - UDOTN * RNZ
               UPARN  = SQRT( UPN*UPN + VPN*VPN + WPN*WPN )
               IF (UPARN .EQ. 0.E0) THEN
                  SPX    = SHATX
                  SPY    = SHATY
                  SPZ    = SHATZ
               ELSE
                  SPX    = UPN / UPARN
                  SPY    = VPN / UPARN
                  SPZ    = WPN / UPARN
               ENDIF
C
               RHOW   = Q    (JBC ,KBC ,IBC ,1)
C
               TAUW   = PROPS(JBC2,KBC2,IBC2,ITQTAU)
               YNORML = ABS (PROPS(JN ,KN ,IN ,4))
               YNORML = MAX ( YNORML, RSMALL )
C
               RHO    = Q(JN,KN,IN, 1)
               TAUXX  = Q(JN,KN,IN, 6)
               TAUYY  = Q(JN,KN,IN, 7)
               TAUZZ  = Q(JN,KN,IN, 8)
               TAUXY  = Q(JN,KN,IN, 9)
               TAUXZ  = Q(JN,KN,IN,10)
               TAUYZ  = Q(JN,KN,IN,11)
               TKEN   = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
               UTAUW  = SQRT (ABS(TAUW) / RHOW * FSMACH / RE)
               DUPDN  = UTAUW / (RKAPVK * YNORML)
C
C Transform stresses to local streamwise, normal coordinates
C
               TAUNN  =   RNX*RNX*TAUXX + RNY*RNY*TAUYY + RNZ*RNZ*TAUZZ 
     1            + 2.E0*(RNX*RNY*TAUXY + RNX*RNZ*TAUXZ + RNY*RNZ*TAUYZ)
               TAUSN  =   SPX*RNX*TAUXX + SPY*RNY*TAUYY + SPZ*RNZ*TAUZZ 
     1              + (SPX*RNY+SPY*RNX)*TAUXY + (SPX*RNZ+SPZ*RNX)*TAUXZ 
     2              + (SPY*RNZ+SPZ*RNY)*TAUYZ
C
C Set sign of du_parallel/dy_normal as opposite of sign of tausn
C
               DUPDN  = SIGN (DUPDN, (-TAUSN))
C
C Calculate production = - rho * du_i/dx_j * Tau_ij
C
               PROD   = - RHO * TAUSN * DUPDN
C
C Boundary condition for epsilon from Younis' thesis as restriction of
C the epsilon transport equation
C
               CRAT   = CEPS1 / CEPS2 * PROD / RHO
               TAUNN  = ABS (TAUNN)
               EPSLNN = 0.5E0 * (CRAT + SQRT (CRAT*CRAT + 
     1              4.E0*CEPSLN/CEPS2*TKEN*TKEN*TAUNN/(YNORML*YNORML)) )
               EPSLNW = EPSLNN
C
               Q(JN  ,KN  ,IN  ,12) = EPSLNW
C
               Q(JBC ,KBC ,IBC , 6) = 0.E0
               Q(JBC ,KBC ,IBC , 7) = 0.E0
               Q(JBC ,KBC ,IBC , 8) = 0.E0
               Q(JBC ,KBC ,IBC , 9) = 0.E0
               Q(JBC ,KBC ,IBC ,10) = 0.E0
               Q(JBC ,KBC ,IBC ,11) = 0.E0
               Q(JBC ,KBC ,IBC ,12) = EPSLNW
C
               Q(JBC2,KBC2,IBC2, 6) = Q(JBC ,KBC ,IBC , 6)
               Q(JBC2,KBC2,IBC2, 7) = Q(JBC ,KBC ,IBC , 7)
               Q(JBC2,KBC2,IBC2, 8) = Q(JBC ,KBC ,IBC , 8)
               Q(JBC2,KBC2,IBC2, 9) = Q(JBC ,KBC ,IBC , 9)
               Q(JBC2,KBC2,IBC2,10) = Q(JBC ,KBC ,IBC ,10)
               Q(JBC2,KBC2,IBC2,11) = Q(JBC ,KBC ,IBC ,11)
               Q(JBC2,KBC2,IBC2,12) = Q(JBC ,KBC ,IBC ,12)
  450       CONTINUE
  460    CONTINUE
C
C No additional bcs
C
      ENDIF
C
C
C Finished with TURBULENT WALL FUNCTIONS boundary condition
C
      RETURN
 1010 FORMAT (' ','BCWALF: Tau_w not converged at (I:J:K) ', 3I5, 
     1            ' after ', I5, ' iterations.',
     2       /' ','  Rel. Error ', E15.7, 5X, 'Tolerance ', E15.7,
     3            ' Y+ = ',E15.7)
      END
