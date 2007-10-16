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
      SUBROUTINE BCAXI  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to set the axisymmetric boundary condition.
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
C 0. Constants for rotation (2 theta) as the grid is rotated theta and the
C    bc data point is a ghost point.
C
      PI        = ACOS(-1.E0)
      STHETA    = SIN (2.E0 * AXIDTH * PI / 180.E0)
      CTHETA    = COS (2.E0 * AXIDTH * PI / 180.E0)
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
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCPT2 is the boundary condition location outside the block from IBCPT
C      IBCN   is the point to the inside of the block from IBCPT
C      IBCM   is the location of the metrics
C      SGNM   is sign for the rotation (negative for K = 1    boundary, 
C                                       positive for K = KDIM boundary)
C
      IF (INROUT .EQ. 1) THEN
         IBCPT  = 1
         IBCPT2 = 0
         IBCN   = 2
         IBCM   = 1
         SGNM   = - 1.E0
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCM   = IJKDIM(IDIR)
         SGNM   = 1.E0
      ENDIF
C
      STHETA = STHETA * SGNM 
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
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
            U    = Q(JN, KN, IN, 2)
            V    = Q(JN, KN, IN, 3)
            W    = Q(JN, KN, IN, 4)
C
            UAXI = U
            VAXI = V * CTHETA 
            WAXI = V * STHETA
C
            Q(JBC ,KBC ,IBC ,1) = Q(JN  ,KN  ,IN  ,1)
            Q(JBC ,KBC ,IBC ,2) = UAXI
            Q(JBC ,KBC ,IBC ,3) = VAXI
            Q(JBC ,KBC ,IBC ,4) = WAXI
            Q(JBC ,KBC ,IBC ,5) = Q(JN  ,KN  ,IN  ,5)
C
            Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
            Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
            Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
            Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
            Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
  110    CONTINUE
  120 CONTINUE
C
C 3-a. Boundary conditions for Two-Equation variables
C
C 3-a-i. Boundary conditions for k-epsilon, k-omega
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
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
               IN   = IBCN  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JN   = IBCN  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KN   = IBCN  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               Q(JBC ,KBC ,IBC ,6) = Q(JN  ,KN  ,IN  ,6)
               Q(JBC ,KBC ,IBC ,7) = Q(JN  ,KN  ,IN  ,7)
C
               Q(JBC2,KBC2,IBC2,6) = Q(JBC ,KBC ,IBC ,6)
               Q(JBC2,KBC2,IBC2,7) = Q(JBC ,KBC ,IBC ,7)
  210       CONTINUE
  220    CONTINUE
C
C 3-b. Boundary conditions for Reynolds stress variables
C
      ELSE IF (ITURB .EQ. ITRS) THEN
C
C Stresses in cylindrical coordinates are invariant in theta
C => Tau_ij_bc = Tau_ij_n (cylindrical coordinates) 
C    Note also that Tau_r,theta = Tau_theta,z = 0
C
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
               IN   = IBCN  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JN   = IBCN  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KN   = IBCN  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               TAUXXN = Q(JN  ,KN  ,IN  , 6)
               TAUYYN = Q(JN  ,KN  ,IN  , 7)
               TAUZZN = Q(JN  ,KN  ,IN  , 8)
               TAUXYN = Q(JN  ,KN  ,IN  , 9)
               TAUXZN = Q(JN  ,KN  ,IN  ,10)
               TAUYZN = Q(JN  ,KN  ,IN  ,11)
C
C Cylindrical coordinates (r,theta,z correspond to R,T(theta),S(streamwise))
C
               TAURR  = TAUYYN
               TAUTT  = TAUZZN
               TAUSS  = TAUXXN
               TAURS  = TAUXYN
C
C Back to Cartesian
C
               TAUXX  =   TAUSS
               TAUYY  =   CTHETA*CTHETA*TAURR + STHETA*STHETA*TAUTT 
               TAUZZ  =   STHETA*STHETA*TAURR + CTHETA*CTHETA*TAUTT 
               TAUXY  =   CTHETA * TAURS
               TAUXZ  = - STHETA * TAURS
               TAUYZ  =   CTHETA*STHETA*(TAUTT - TAURR)
C
               Q(JBC ,KBC ,IBC , 6) = TAUXX
               Q(JBC ,KBC ,IBC , 7) = TAUYY
               Q(JBC ,KBC ,IBC , 8) = TAUZZ
               Q(JBC ,KBC ,IBC , 9) = TAUXY
               Q(JBC ,KBC ,IBC ,10) = TAUXZ
               Q(JBC ,KBC ,IBC ,11) = TAUYZ
               Q(JBC ,KBC ,IBC ,12) = Q(JN  ,KN  ,IN  ,12)
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
C
C Finished with FIX boundary condition
C
      RETURN
      END
