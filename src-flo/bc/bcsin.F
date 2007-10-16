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
      SUBROUTINE BCSIN  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to do SUBSONIC INFLOW boundary condition.  This routine is based
C on Jeff White's recommendation to use the following: specify total pressure
C and total Temperature, specify flow direction, extrapolate velocities from
C the interior, solve for the speed of sound, get M, calculate static
C conditions using isentropic relations.  Then specify turbulence quantities.
C
C This routine should be used for nozzles, tunnels, etc.
C
C NOTE: This routine cannot currently be used on a boundary that adjoins
C       a WALL or WALL FUNCTION boundary.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for the mean variables
C     3. Set additional boundary conditions for turbulence variables:
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
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
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
         ISTRT(3) = 2
         IEND (3) = 2
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
         IBCM   = 1
         SGNM   = 1.E0
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCM   = IJKDIM(IDIR)
         SGNM   = - 1.E0
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C Conditions at infinity
C
      SINF  = PREF / RHOINF**GAMMA
C
C 2. Boundary conditions for Mean Variables
C
      DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
         DO 110 I2 = ISTRT(IC2), IEND(IC2)
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
C Set the flow direction to use.  Jeff recommends three options:
C     i. Normal to boundary face - default
C    ii. Aligned with local grid - better for highly curved nozzles
C   iii. User input              - to handle other cases
C Current implementation is for i.
C
C The use of the SGNM variable sets up one set of logic
C to handle both the inner and outer boundary cases.
C
            SX    = S(JM,KM,IM,1,IDIR) * SGNM
            SY    = S(JM,KM,IM,2,IDIR) * SGNM
            SZ    = S(JM,KM,IM,3,IDIR) * SGNM
C
            QEXT  = Q(JN,KN,IN,2) * SX + Q(JN,KN,IN,3) * SY 
     1            + Q(JN,KN,IN,4) * SZ
            RHON  = Q(JN,KN,IN,1)
            PN    = Q(JN,KN,IN,5)
C
C One Dimensional Riemann invariant
C - the sign of U dot n handles the sign change in the Riemann invariants
C
            AEXT  = SQRT (GAMMA * PN / RHON)
            REXT  = QEXT - 2.E0 * AEXT / GAMM1
            HTEXT = (AEXT*AEXT/GAMM1 + 0.5E0*QEXT*QEXT)
            A     = ( - REXT + SQRT(REXT*REXT 
     1                          +(GAMMA+1.E0)*(HTEXT-0.5E0*REXT*REXT)) )
     2              * GAMM1 / (GAMMA + 1.E0)
            QDOTN = REXT + 2.E0 * A / GAMM1
            RMACH = QDOTN / A
C
C Isentropic relationships
C
            ARG   = 1.E0 + 0.5E0 * GAMM1 * RMACH * RMACH
            P     = PTOTBC / ARG**(GAMMA/GAMM1)
            T     = TTOTBC / ARG
            RHO   = GAMMA * P / T
C
C Mean Variables
C
            Q(JBC ,KBC ,IBC ,1) = RHO
            Q(JBC ,KBC ,IBC ,2) = QDOTN * SX
            Q(JBC ,KBC ,IBC ,3) = QDOTN * SY
            Q(JBC ,KBC ,IBC ,4) = QDOTN * SZ
            Q(JBC ,KBC ,IBC ,5) = P
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
C
  110    CONTINUE
  120 CONTINUE
C
C 3-a. Boundary conditions for Two-Equation variables
C
C 3-a-i. Boundary conditions for k-epsilon
C
      IF (ITURB .EQ. ITKE) THEN
         DO 220 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 210 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               Q(JBC ,KBC ,IBC ,6) = TKEINF
               Q(JBC ,KBC ,IBC ,7) = EPSINF
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  210       CONTINUE
  220    CONTINUE
C
C 3-a-ii. Boundary conditions for k-omega
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         DO 240 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
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
               Q(JBC ,KBC ,IBC ,6) = TKEINF
               Q(JBC ,KBC ,IBC ,7) = OMEINF
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  230       CONTINUE
  240    CONTINUE
C
C 3-b. Boundary conditions for Reynolds stress variables
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 260 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
            DO 250 I2 = ISTRT(IC2), IEND(IC2)
C
               IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               IBC2 = IBCPT2*IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JBC2 = IBCPT2*IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KBC2 = IBCPT2*IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               Q(JBC ,KBC ,IBC , 6) = TAUINF(1)
               Q(JBC ,KBC ,IBC , 7) = TAUINF(2)
               Q(JBC ,KBC ,IBC , 8) = TAUINF(3)
               Q(JBC ,KBC ,IBC , 9) = TAUINF(4)
               Q(JBC ,KBC ,IBC ,10) = TAUINF(5)
               Q(JBC ,KBC ,IBC ,11) = TAUINF(6)
               Q(JBC ,KBC ,IBC ,12) = EPSINF
C
               Q(JBC2,KBC2,IBC2, 6) = Q(JBC ,KBC ,IBC , 6)
               Q(JBC2,KBC2,IBC2, 7) = Q(JBC ,KBC ,IBC , 7)
               Q(JBC2,KBC2,IBC2, 8) = Q(JBC ,KBC ,IBC , 8)
               Q(JBC2,KBC2,IBC2, 9) = Q(JBC ,KBC ,IBC , 9)
               Q(JBC2,KBC2,IBC2,10) = Q(JBC ,KBC ,IBC ,10)
               Q(JBC2,KBC2,IBC2,11) = Q(JBC ,KBC ,IBC ,11)
               Q(JBC2,KBC2,IBC2,12) = Q(JBC ,KBC ,IBC ,12)
  250       CONTINUE
  260    CONTINUE
C
C No additional bcs
C
      ENDIF
C
C Finished with SUBSONIC INFLOW boundary condition
C
      RETURN
      END
