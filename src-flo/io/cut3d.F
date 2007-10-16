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
      SUBROUTINE CUT3D (IDIM, JDIM, KDIM, Q, NV, ICUTS)
C
C     This routine averages the value of the interior and boundary 
C     condition flow variable at the cell center for the CUT boundary
C     condition type to the face to coincide with the PLOT3D output.
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NV)
C
      DIMENSION ICUTS (21)
C
C     Local variables
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3), 
     2          IEND  (3),
     3          ISGN  (3),
     4          LEN   (3) 
C
C     Store IDIM, JDIM, KDIM in IJKDIM for looping
C
      IJKDIM(1) = IDIM
      IJKDIM(2) = JDIM
      IJKDIM(3) = KDIM
C
C     Store boundary condition data in mnemonic names
C
      IDIR     = ICUTS( 2)
      INOUT    = ICUTS( 3)
      ISTRT(1) = ICUTS( 4)
      IEND (1) = ICUTS( 5)
      ISTRT(2) = ICUTS( 6)
      IEND (2) = ICUTS( 7)
      ISTRT(3) = ICUTS( 8)
      IEND (3) = ICUTS( 9)
C
C     Set needed variables depending on whether the boundary is
C     the inner boundary (INOUT = 1) or the outer boundary (INOUT > 1)
C          INR    = interior point of block
C          IBCPT  = boundary point of block
C
      IF (INOUT .EQ. 1) THEN
         INR    = 2
         IBCPT  = 1
      ELSE
         INR    = IJKDIM(IDIR)
         IBCPT  = IJKDIM(IDIR) + 1
      ENDIF
C
C     Find the length of the two outer loops and loop over these using
C     offsets and delta function to set the two cut data points to the
C     two interior data points of block 2.
C
      DO 100 L = 1, 3
C
         LEN(L) = ABS ( IEND(L) - ISTRT(L) )
C
C     Increment/Decrement 
C
         IF ( IEND(L) .GT. ISTRT(L) ) THEN
            ISGN(L) =   1
         ELSE
            ISGN(L) = - 1
         ENDIF
  100 CONTINUE
C
C     II1 first  boundary point of block 1
C     II2 second boundary point of block 1
C     IN1 first  interior point of block 2
C     IN2 second interior point of block 2
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
      DO 800 I3 = 0, LEN(IC3)
         DO 700 I2 = 0, LEN(IC2)
            IBC = IBCPT                       * IKD (IC1, 1) +
     1            (ISTRT(IC2) + ISGN(IC2)*I2) * IKD (IC2, 1) +
     2            (ISTRT(IC3) + ISGN(IC3)*I3) * IKD (IC3, 1) 
            JBC = IBCPT                       * IKD (IC1, 2) +
     1            (ISTRT(IC2) + ISGN(IC2)*I2) * IKD (IC2, 2) +
     2            (ISTRT(IC3) + ISGN(IC3)*I3) * IKD (IC3, 2) 
            KBC = IBCPT                       * IKD (IC1, 3) +
     1            (ISTRT(IC2) + ISGN(IC2)*I2) * IKD (IC2, 3) +
     2            (ISTRT(IC3) + ISGN(IC3)*I3) * IKD (IC3, 3) 
C
            IN  = INR                         * IKD (IC1, 1) +
     1            (ISTRT(IC2) + ISGN(IC2)*I2) * IKD (IC2, 1) +
     2            (ISTRT(IC3) + ISGN(IC3)*I3) * IKD (IC3, 1) 
            JN  = INR                         * IKD (IC1, 2) +
     1            (ISTRT(IC2) + ISGN(IC2)*I2) * IKD (IC2, 2) +
     2            (ISTRT(IC3) + ISGN(IC3)*I3) * IKD (IC3, 2) 
            KN  = INR                         * IKD (IC1, 3) +
     1            (ISTRT(IC2) + ISGN(IC2)*I2) * IKD (IC2, 3) +
     2            (ISTRT(IC3) + ISGN(IC3)*I3) * IKD (IC3, 3) 
C
CDIR$ SHORTLOOP
            DO 600 L = 1, NV
               Q(JBC,KBC,IBC,L) = 0.5E0 * (Q(JN,KN,IN,L) +
     1                                     Q(JBC,KBC,IBC,L))
 600        CONTINUE
 700     CONTINUE
 800  CONTINUE
C
C     Finished with Cut boundary condition fix to PLOT3D
C
      RETURN
      END
