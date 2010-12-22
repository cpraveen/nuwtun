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
      SUBROUTINE BCTAN  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to do TANGENCY (inviscid wall) boundary condition.
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
         IBCN2  = 3
         IBCM   = 1
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCN2  = IJKDIM(IDIR) - 1
         IBCM   = IJKDIM(IDIR)
      ENDIF
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
            CALL BCTINV(IDIM, JDIM, KDIM, IDIR, IBCPT, IBCPT2, IBCN, 
     1                  IBCN2, IBCM, IC1, IC2, IC3, I2, I3, Q, S)
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
               IN   = IBCN  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JN   = IBCN  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KN   = IBCN  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               TKEN = Q(JN,KN,IN,6)
               EPSN = Q(JN,KN,IN,7)
               TKEW = TKEN
               EPSW = EPSN
C
               Q(JBC ,KBC ,IBC ,6) = TKEW
               Q(JBC ,KBC ,IBC ,7) = EPSW
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
               IN   = IBCN  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JN   = IBCN  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KN   = IBCN  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               TKEN = Q(JN,KN,IN,6)
               OMEN = Q(JN,KN,IN,7)
               TKEW = TKEN
               OMEW = OMEN
C
               Q(JBC ,KBC ,IBC ,6) = TKEW
               Q(JBC ,KBC ,IBC ,7) = OMEW
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
               IN   = IBCN  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JN   = IBCN  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KN   = IBCN  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               TXXN = Q(JN,KN,IN, 6)
               TYYN = Q(JN,KN,IN, 7)
               TZZN = Q(JN,KN,IN, 8)
               TXYN = Q(JN,KN,IN, 9)
               TXZN = Q(JN,KN,IN,10)
               TYZN = Q(JN,KN,IN,11)
               EPSN = Q(JN,KN,IN,12)
               TXXW = TXXN
               TYYW = TYYN
               TZZW = TZZN
               TXYW = TXYN
               TXZW = TXZN
               TYZW = TYZN
               EPSW = EPSN
C
               Q(JBC ,KBC ,IBC , 6) = TXXW
               Q(JBC ,KBC ,IBC , 7) = TYYW
               Q(JBC ,KBC ,IBC , 8) = TZZW
               Q(JBC ,KBC ,IBC , 9) = TXYW
               Q(JBC ,KBC ,IBC ,10) = TXZW
               Q(JBC ,KBC ,IBC ,11) = TYZW
               Q(JBC ,KBC ,IBC ,12) = EPSW
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
C Finished with TANGENCY boundary condition
C
      RETURN
      END

C
C     Apply inviscid wall boundary condition for given face
C
      SUBROUTINE BCTINV(IDIM, JDIM, KDIM, IDIR, IBCPT, IBCPT2, IBCN, 
     1                  IBCN2, IBCM, IC1, IC2, IC3, I2, I3, Q, S)
      INCLUDE '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
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
      RHON  = Q(JN ,KN ,IN ,1)
      PN    = Q(JN ,KN ,IN ,5)
      PN2   = Q(JN2,KN2,IN2,5)
C
c-1st order PW    = PN
      PW    = (9.E0 * PN - PN2) / 8.E0
      RHOW  = RHON * ( PW / PN )**(1.E0/GAMMA)
C
      UBAR  = Q(JN,KN,IN,2) * S(JM,KM,IM,1,IDIR)
     1      + Q(JN,KN,IN,3) * S(JM,KM,IM,2,IDIR)
     2      + Q(JN,KN,IN,4) * S(JM,KM,IM,3,IDIR)
C
      Q(JBC ,KBC ,IBC ,1) = RHOW
      Q(JBC ,KBC ,IBC ,2) = Q(JN,KN,IN,2)-UBAR*S(JM,KM,IM,1,IDIR)
      Q(JBC ,KBC ,IBC ,3) = Q(JN,KN,IN,3)-UBAR*S(JM,KM,IM,2,IDIR)
      Q(JBC ,KBC ,IBC ,4) = Q(JN,KN,IN,4)-UBAR*S(JM,KM,IM,3,IDIR)
      Q(JBC ,KBC ,IBC ,5) = PW
C
      Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)
      Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
      Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
      Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
      Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)
C
      RETURN
      END
C
C     Apply inviscid wall boundary condition for given face
C     using normal momentum equation
C
      SUBROUTINE BCTINV2(IDIM, JDIM, KDIM, IDIR, IBCPT, IBCPT2, IBCN, 
     1                  IBCN2, IBCM, IC1, IC2, IC3, I2, I3, Q, S)
      INCLUDE '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION UXI   (3),
     1          UET   (3),
     2          UZE   (3)
      REAL      RHS, LHS
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
      RHON  = Q(JN ,KN ,IN ,1)
      PN    = Q(JN ,KN ,IN ,5)
      PN2   = Q(JN2,KN2,IN2,5)
C
C Unit normal
C
      SX = S(JM,KM,IM,1,IDIR)
      SY = S(JM,KM,IM,2,IDIR)
      SZ = S(JM,KM,IM,3,IDIR)
C
C Normal velocity
C
      UBAR  = Q(JN,KN,IN,2) * SX
     1      + Q(JN,KN,IN,3) * SY
     2      + Q(JN,KN,IN,4) * SZ
C
      Q(JBC ,KBC ,IBC ,2) = Q(JN,KN,IN,2)-UBAR*SX
      Q(JBC ,KBC ,IBC ,3) = Q(JN,KN,IN,3)-UBAR*SY
      Q(JBC ,KBC ,IBC ,4) = Q(JN,KN,IN,4)-UBAR*SZ
C
      Q(JBC2,KBC2,IBC2,2) = Q(JBC ,KBC ,IBC ,2)
      Q(JBC2,KBC2,IBC2,3) = Q(JBC ,KBC ,IBC ,3)
      Q(JBC2,KBC2,IBC2,4) = Q(JBC ,KBC ,IBC ,4)
C
C Contravariant velocities
C      
      UT = (Q(JBC,KBC,IBC,2) * S(JN,KN,IN,1,1)
     1   +  Q(JBC,KBC,IBC,3) * S(JN,KN,IN,2,1)
     2   +  Q(JBC,KBC,IBC,4) * S(JN,KN,IN,3,1)) * S(JN,KN,IN,4,1)
      VT = (Q(JBC,KBC,IBC,2) * S(JN,KN,IN,1,2)
     1   +  Q(JBC,KBC,IBC,3) * S(JN,KN,IN,2,2)
     2   +  Q(JBC,KBC,IBC,4) * S(JN,KN,IN,3,2)) * S(JN,KN,IN,4,2)
      WT = (Q(JBC,KBC,IBC,2) * S(JN,KN,IN,1,3)
     1   +  Q(JBC,KBC,IBC,3) * S(JN,KN,IN,2,3)
     2   +  Q(JBC,KBC,IBC,4) * S(JN,KN,IN,3,3)) * S(JN,KN,IN,4,3)
C      
C du/d(xi)
C      
      UXI(1) = (Q(JN,KN,IN+1,2) - Q(JN,KN,IN-1,2))/2.0
      UXI(2) = (Q(JN,KN,IN+1,3) - Q(JN,KN,IN-1,3))/2.0
      UXI(3) = (Q(JN,KN,IN+1,4) - Q(JN,KN,IN-1,4))/2.0
C      
C du/d(eta)
C      
      UET(1) = (Q(JN+1,KN,IN,2) - Q(JN-1,KN,IN,2))/2.0
      UET(2) = (Q(JN+1,KN,IN,3) - Q(JN-1,KN,IN,3))/2.0
      UET(3) = (Q(JN+1,KN,IN,4) - Q(JN-1,KN,IN,4))/2.0
C      
C du/d(zeta)
C      
      UZE(1) = (Q(JN,KN+1,IN,2) - Q(JN,KN-1,IN,2))/2.0
      UZE(2) = (Q(JN,KN+1,IN,3) - Q(JN,KN-1,IN,3))/2.0
      UZE(3) = (Q(JN,KN+1,IN,4) - Q(JN,KN-1,IN,4))/2.0
C
C Dot product with normal vector
C
      UXIN = UXI(1) * SX + UXI(2) * SY + UXI(3) * SZ
      UETN = UET(1) * SX + UET(2) * SY + UET(3) * SZ
      UZEN = UZE(1) * SX + UZE(2) * SY + UZE(3) * SZ
      RHS = -RHON * ( UT*UXIN + VT*UETN + WT*UZEN )

      IF(IDIR.EQ.1)THEN
         PET = (Q(JN+1,KN,IN,5) - Q(JN-1,KN,IN,5))/2.0
         PZE = (Q(JN,KN+1,IN,5) - Q(JN,KN-1,IN,5))/2.0
         LHS = (SX*S(JN,KN,IN,1,2) 
     1       +  SY*S(JN,KN,IN,2,2)
     2       +  SZ*S(JN,KN,IN,3,2)) * S(JN,KN,IN,4,2) * PET
     3       + (SX*S(JN,KN,IN,1,3) 
     4       +  SY*S(JN,KN,IN,2,3)
     5       +  SZ*S(JN,KN,IN,3,3)) * S(JN,KN,IN,4,3) * PZE
         PXI = (RHS - LHS)/S(JN,KN,IN,4,1)
         PW  = (9.0*PN - PN2 - 3.0*PXI)/8.0
      ELSE IF(IDIR.EQ.2)THEN
         PXI = (Q(JN,KN,IN+1,5) - Q(JN,KN,IN-1,5))/2.0
         PZE = (Q(JN,KN+1,IN,5) - Q(JN,KN-1,IN,5))/2.0
         LHS = (SX*S(JN,KN,IN,1,3) 
     1       +  SY*S(JN,KN,IN,2,3)
     2       +  SZ*S(JN,KN,IN,3,3)) * S(JN,KN,IN,4,3) * PZE
     3       + (SX*S(JN,KN,IN,1,1) 
     4       +  SY*S(JN,KN,IN,2,1)
     5       +  SZ*S(JN,KN,IN,3,1)) * S(JN,KN,IN,4,1) * PXI
         PET = (RHS - LHS)/S(JN,KN,IN,4,2)
         PW  = (9.0*PN - PN2 - 3.0*PET)/8.0
      ELSE
         PXI = (Q(JN,KN,IN+1,5) - Q(JN,KN,IN-1,5))/2.0
         PET = (Q(JN+1,KN,IN,5) - Q(JN-1,KN,IN,5))/2.0
         LHS = (SX*S(JN,KN,IN,1,1) 
     1       +  SY*S(JN,KN,IN,2,1)
     2       +  SZ*S(JN,KN,IN,3,1)) * S(JN,KN,IN,4,1) * PXI
     3       + (SX*S(JN,KN,IN,1,2) 
     4       +  SY*S(JN,KN,IN,2,2)
     5       +  SZ*S(JN,KN,IN,3,2)) * S(JN,KN,IN,4,2) * PET
         PZE = (RHS - LHS)/S(JN,KN,IN,4,3)
         PW  = (9.0*PN - PN2 - 3.0*PZE)/8.0
      ENDIF

      Q(JBC ,KBC ,IBC ,5) = PW
      Q(JBC2,KBC2,IBC2,5) = Q(JBC ,KBC ,IBC ,5)

      RHOW  = RHON * ( PW / PN )**(1.E0/GAMMA)
      Q(JBC ,KBC ,IBC ,1) = RHOW
      Q(JBC2,KBC2,IBC2,1) = Q(JBC ,KBC ,IBC ,1)

      RETURN
      END
