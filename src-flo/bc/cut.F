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
      SUBROUTINE CUT    (IDIM1, JDIM1, KDIM1, Q1,
     1                   IDIM2, JDIM2, KDIM2, Q2, ICUTS, NV)
C
C Routine to do cut boundary condition.
C Flow data stored in Q1 are updated from data in the interior of Q2.
C This routine works for Q or PROPS and needs the size of the array.
C
C IDIM1,JDIM1,KDIM1 : Block 1 dimensions 
C Q1                : Block 1 data (Q or PROPS)
C IDIM2,JDIM2,KDIM2 : Block 2 dimensions 
C Q2                : Block 2 data (Q or PROPS)
C ICUTS             : Array controlling beginning/ending of cut
C NV                : Number of variables (NV = NQ for Q, NP for PROPS)
C
      include '../header/common.h'
      DIMENSION Q1   (0:JDIM1+2,0:KDIM1+2,0:IDIM1+2,NV),
     1          Q2   (0:JDIM2+2,0:KDIM2+2,0:IDIM2+2,NV)
C
      DIMENSION ICUTS(21)
C
      DIMENSION IJKDM1(3), IJKDM2(3)
      DIMENSION ISTR1(3),  IEND1(3), 
     1          ISTR2(3),  IEND2(3), 
     2          ISGN1(3),  ISGN2(3),
     3          IORD(3),   LEN(3)
C
C     Store IDIM, JDIM, KDIM in IJKDIM for looping
C
      IJKDM1(1) = IDIM1
      IJKDM1(2) = JDIM1
      IJKDM1(3) = KDIM1
      IJKDM2(1) = IDIM2
      IJKDM2(2) = JDIM2
      IJKDM2(3) = KDIM2
C
C     Store boundary condition data in mnemonic names
C
      IDIR1    = ICUTS( 2)
      INOUT1   = ICUTS( 3)
      ISTR1(1) = ICUTS( 4)
      IEND1(1) = ICUTS( 5)
      ISTR1(2) = ICUTS( 6)
      IEND1(2) = ICUTS( 7)
      ISTR1(3) = ICUTS( 8)
      IEND1(3) = ICUTS( 9)
C
      IBLK2    = ICUTS(10)
      IDIR2    = ICUTS(11)
      INOUT2   = ICUTS(12)
      ISTR2(1) = ICUTS(13)
      IEND2(1) = ICUTS(14)
      ISTR2(2) = ICUTS(15)
      IEND2(2) = ICUTS(16)
      ISTR2(3) = ICUTS(17)
      IEND2(3) = ICUTS(18)
C
      IORD(1)  = ICUTS(19)
      IORD(2)  = ICUTS(20)
      IORD(3)  = ICUTS(21)
C
C     Set needed variables depending on whether the boundary is
C     the inner boundary (INOUT1 = 1) or the outer boundary (INOUT1 > 1)
C          IBCPT  = boundary point of block 1
C          IBCINC = increment to second boundary point of block 1
C                   IBCPT2 = IBCPT + IBCINC
C          INR    = interior point of block 2.
C          INRINC = increment to second interior point of block 2
C                   INR2 = INR + INRINC
C
      IF (INOUT1 .EQ. 1) THEN
         IBCPT  =   1
         IBCINC = - 1
      ELSE
         IBCPT  = IJKDM1(IDIR1) + 1
         IBCINC = 1
      ENDIF
C
      IF (INOUT2 .EQ. 1) THEN
         INR    =   2
         INRINC =   1
      ELSE
         INR    =   IJKDM2(IDIR2)
         INRINC = - 1
      ENDIF
C
C     Find the length of the two outer loops and loop over these using
C     offsets and delta function to set the two cut data points to the
C     two interior data points of block 2.
C
      DO 100 L = 1, 3
C
         LEN(L) = ABS ( IEND1(L) - ISTR1(L) )
C
C     Increment/Decrement 
C
         IF ( IEND1(L) .GT. ISTR1(L) ) THEN
            ISGN1(L) =   1
         ELSE
            ISGN1(L) = - 1
         ENDIF
C
C     Increment/Decrement 
C
         IF ( IEND2(L) .GT. ISTR2(L) ) THEN
            ISGN2(L) =   1
         ELSE
            ISGN2(L) = - 1
         ENDIF
  100 CONTINUE
C
C     II1 first  boundary point of block 1
C     II2 second boundary point of block 1
C     IN1 first  interior point of block 2
C     IN2 second interior point of block 2
C
      IC1 = ICY (IDIR1, 1)
      IC2 = ICY (IDIR1, 2)
      IC3 = ICY (IDIR1, 3)
C
      JC1 = IORD (IC1)
      JC2 = IORD (IC2)
      JC3 = IORD (IC3)
C
      DO 800 I3 = 0, LEN(IC3)
         DO 700 I2 = 0, LEN(IC2)
            II1 = IBCPT                        * IKD (IC1, 1) +
     1            (ISTR1(IC2) + ISGN1(IC2)*I2) * IKD (IC2, 1) +
     2            (ISTR1(IC3) + ISGN1(IC3)*I3) * IKD (IC3, 1) 
            JJ1 = IBCPT                        * IKD (IC1, 2) +
     1            (ISTR1(IC2) + ISGN1(IC2)*I2) * IKD (IC2, 2) +
     2            (ISTR1(IC3) + ISGN1(IC3)*I3) * IKD (IC3, 2) 
            KK1 = IBCPT                        * IKD (IC1, 3) +
     1            (ISTR1(IC2) + ISGN1(IC2)*I2) * IKD (IC2, 3) +
     2            (ISTR1(IC3) + ISGN1(IC3)*I3) * IKD (IC3, 3) 
C
            II2 = II1 + IBCINC * IKD (IC1, 1)
            JJ2 = JJ1 + IBCINC * IKD (IC1, 2)
            KK2 = KK1 + IBCINC * IKD (IC1, 3)
C
            IN1 = INR                          * IKD (JC1, 1) +
     1            (ISTR2(JC2) + ISGN2(JC2)*I2) * IKD (JC2, 1) +
     2            (ISTR2(JC3) + ISGN2(JC3)*I3) * IKD (JC3, 1) 
            JN1 = INR                          * IKD (JC1, 2) +
     1            (ISTR2(JC2) + ISGN2(JC2)*I2) * IKD (JC2, 2) +
     2            (ISTR2(JC3) + ISGN2(JC3)*I3) * IKD (JC3, 2) 
            KN1 = INR                          * IKD (JC1, 3) +
     1            (ISTR2(JC2) + ISGN2(JC2)*I2) * IKD (JC2, 3) +
     2            (ISTR2(JC3) + ISGN2(JC3)*I3) * IKD (JC3, 3) 
C
            IN2 = IN1 + INRINC * IKD (JC1, 1)
            JN2 = JN1 + INRINC * IKD (JC1, 2)
            KN2 = KN1 + INRINC * IKD (JC1, 3)
C
CDIR$ SHORTLOOP
            DO 600 L = 1, NV
               Q1(JJ1,KK1,II1,L) = Q2(JN1,KN1,IN1,L)
               Q1(JJ2,KK2,II2,L) = Q2(JN2,KN2,IN2,L)
 600        CONTINUE
 700     CONTINUE
 800  CONTINUE
C
C     Finished with Cut boundary condition
C
      RETURN
      END
