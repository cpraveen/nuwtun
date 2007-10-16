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
      SUBROUTINE BCWALL (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT)
C
C Routine to do VISCOUS WALL boundary condition.
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
C Set control for higher order evaluation of PW
C     IHIORD = 0 for first  order
C     IHIORD = 1 for higher order
C
      IF (FOURTH) THEN
         IHIORD = 1
      ELSE
         IHIORD = 0
      ENDIF
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
         IBCN3  = 4
         IBCM   = 1
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCN2  = IJKDIM(IDIR) - 1
         IBCN3  = IJKDIM(IDIR) - 2
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
            IN3  = IBCN3 *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JN3  = IBCN3 *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KN3  = IBCN3 *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
c-1st order PW    = Q(JN,KN,IN,5)
            PN    = Q(JN ,KN ,IN ,5)
            PN2   = Q(JN2,KN2,IN2,5)
            PN3   = Q(JN3,KN3,IN3,5)
c-2nd order PW    = (9.E0 * PN - PN2) / 8.E0
c-3rd order PW    = (225.E0 * PN - 50.E0 * PN2 + 9.E0 * PN3) / 184.E0
            PW    = (1 - IHIORD) * PN +
     1              IHIORD * (225.E0*PN - 50.E0*PN2 + 9.E0*PN3) / 184.E0
            IF ( TWSPE ) THEN
               RHOW = GAMMA * PW / ( TWALL / TINF )
            ELSE
               RHOW = Q(JN,KN,IN,1)
            ENDIF
            Q(JBC ,KBC ,IBC ,1) = RHOW
            Q(JBC ,KBC ,IBC ,2) = 0.E0
            Q(JBC ,KBC ,IBC ,3) = 0.E0
            Q(JBC ,KBC ,IBC ,4) = 0.E0
            Q(JBC ,KBC ,IBC ,5) = PW
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
               RHOW   = Q    (JBC,KBC,IBC,1)
               RMUW   = PROPS(JBC,KBC,IBC,1)
C
               TKEN   = Q    (JN ,KN ,IN ,6)
               YNORML = ABS (PROPS(JN ,KN ,IN ,4))
C
C Wall boundary condition for epsilon from Speziale
C epsilon+(wall) = 2 NU K+(1) / Y+(1)**2
C
               EPSLNW = 2.E0 * RMUW / RHOW * TKEN
     1                  / (YNORML * YNORML + RSMASQ) * FSMACH / RE
C
               Q(JBC ,KBC ,IBC ,6) = 0.E0
               Q(JBC ,KBC ,IBC ,7) = EPSLNW
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
               IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
               JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
               KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
               RHOW   = Q    (JBC,KBC,IBC,1)
               RMUW   = PROPS(JBC,KBC,IBC,1)
C
               YNORML = ABS (PROPS(JN ,KN ,IN ,4))
C
C-----Begin wall BC based on Wilcox's rough wall approach
C     Calculate the velocity parallel to the wall at the cell center off 
C     the wall.  U(parallel) = U - UBAR where U is the total velocity vector
C     and UBAR is the velocity normal to the wall.  (Velocity on the wall 
C     should be identically zero)
C     VELN is magnitude of the parallel velocity at point N
C
               UBARN = S(JM,KM,IM,1,IDIR) * Q(JN,KN,IN,2)
     1               + S(JM,KM,IM,2,IDIR) * Q(JN,KN,IN,3)
     2               + S(JM,KM,IM,3,IDIR) * Q(JN,KN,IN,4)
               VELN  = SQRT((Q(JN,KN,IN,2)-UBARN*S(JM,KM,IM,1,IDIR))**2
     1                     +(Q(JN,KN,IN,3)-UBARN*S(JM,KM,IM,2,IDIR))**2
     2                     +(Q(JN,KN,IN,4)-UBARN*S(JM,KM,IM,3,IDIR))**2)
C
C     Tau(wall) = mu * du/dn
C
               DVELDN = VELN / (YNORML + RSMALL)
               TAUW   = RMUW * DVELDN
C
C     Kr+         = utau *  rho(wall) * Kr / mu(wall)
C     Omega(wall) = utau**2 rho(wall) * Sr / mu(wall)
C     If Kr = 0, set Kr+ = 5 -> hydraulically smooth
C
               IF (ROUGHK .GT. 0.E0) THEN
                  RKRPLS = SQRT (RHOW * TAUW) * ROUGHK / RMUW
               ELSE
                  RKRPLS = 5.E0
               ENDIF
C
               IF (RKRPLS .GE. 25.0E0) THEN
                  SR =  100.E0 / RKRPLS
               ELSE
                  SR = 2500.E0 / (RKRPLS * RKRPLS)
               ENDIF
C
               OMEGAW = DVELDN * SR
C-----End   wall BC based on Wilcox's rough wall approach
C-----Begin wall BC based on Menter (NASA TM 103975)
c-             OMEGAW = 60.E0*RMUW / (BKW*RHOW*YNORML*YNORML+RSMASQ)
c-   1                  *FSMACH/RE
C-----End   wall BC based on Menter (NASA TM 103975)
C
               Q(JBC ,KBC ,IBC ,6) = 0.E0
               Q(JBC ,KBC ,IBC ,7) = OMEGAW
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
               RHOW  = Q    (JBC,KBC,IBC,1)
               RMUW  = PROPS(JBC,KBC,IBC,1)
C
               TKEN  = 0.5E0*(Q(JN,KN,IN,6)+Q(JN,KN,IN,7)+Q(JN,KN,IN,8))
               YNORML= ABS (PROPS(JN ,KN ,IN ,4))
C
C Wall boundary condition for epsilon from Speziale
C epsilon+(wall) = 2 NU K+(1) / Y+(1)**2
C
               EPSLNW = 2.E0 * RMUW / RHOW * TKEN
     1                  / (YNORML * YNORML + RSMASQ) * FSMACH / RE
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
  250       CONTINUE
  260    CONTINUE
C
C No additional bcs
C
      ENDIF
C
C
C Finished with WALL boundary condition
C
      RETURN
      END
