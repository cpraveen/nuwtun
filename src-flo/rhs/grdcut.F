      SUBROUTINE GRDCUT (IDIM1, JDIM1, KDIM1, R1,
     1                   IDIM2, JDIM2, KDIM2, R2, ICUTS)
C
C Routine to do cut boundary condition on grid points R.
C Grid points stored in R1 are updated from grid points in R2.
C
C IDIM1,JDIM1,KDIM1 : Dimensions of current block 1
C R1                : Position vector of grid vertices (x,y,z) for block 1
C IDIM2,JDIM2,KDIM2 : Dimensions of current block 2
C R2                : Position vector of grid vertices (x,y,z) for block 2
C ICUTS             : Data controlling transfer from block 2 to block 1
C
      include '../header/common.h'
      DIMENSION R1   (0:JDIM1+1,0:KDIM1+1,0:IDIM1+1,3),
     1          R2   (0:JDIM2+1,0:KDIM2+1,0:IDIM2+1,3)
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
C     Translate ICUTS data from cell-centered storage to grid point storage.
C     This requires decrementing the minimum ('starting') index. Do this only 
C     for the indices that the cut loops over, not the idir direction.
C
      DO 10 L = 1, 3
C
         IF (L .NE. IDIR1) THEN
            IF (IEND1(L) .GE. ISTR1(L)) THEN
               ISTR1(L) = ISTR1(L) - 1
            ELSE
               IEND1(L) = IEND1(L) - 1
            ENDIF
         ENDIF
C
         IF (L .NE. IDIR2) THEN
            IF (IEND2(L) .GE. ISTR2(L)) THEN
               ISTR2(L) = ISTR2(L) - 1
            ELSE
               IEND2(L) = IEND2(L) - 1
            ENDIF
         ENDIF
   10 CONTINUE
C
C     Set needed variables depending on whether the boundary is
C     the inner boundary (INOUT1 = 1) or the outer boundary (INOUT1 > 1)
C          IBCPT  = boundary grid point of block 1
C          INR    = corresponding point of block 2
C
      IF (INOUT1 .EQ. 1) THEN
         IBCPT  = 0
      ELSE
         IBCPT  = IJKDM1(IDIR1) + 1
      ENDIF
C
      IF (INOUT2 .EQ. 1) THEN
         INR    = 1
      ELSE
         INR    = IJKDM2(IDIR2)
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
            II = IBCPT                        * IKD (IC1, 1) +
     1           (ISTR1(IC2) + ISGN1(IC2)*I2) * IKD (IC2, 1) +
     2           (ISTR1(IC3) + ISGN1(IC3)*I3) * IKD (IC3, 1) 
            JJ = IBCPT                        * IKD (IC1, 2) +
     1           (ISTR1(IC2) + ISGN1(IC2)*I2) * IKD (IC2, 2) +
     2           (ISTR1(IC3) + ISGN1(IC3)*I3) * IKD (IC3, 2) 
            KK = IBCPT                        * IKD (IC1, 3) +
     1           (ISTR1(IC2) + ISGN1(IC2)*I2) * IKD (IC2, 3) +
     2           (ISTR1(IC3) + ISGN1(IC3)*I3) * IKD (IC3, 3) 
C
            IN = INR                          * IKD (JC1, 1) +
     1           (ISTR2(JC2) + ISGN2(JC2)*I2) * IKD (JC2, 1) +
     2           (ISTR2(JC3) + ISGN2(JC3)*I3) * IKD (JC3, 1) 
            JN = INR                          * IKD (JC1, 2) +
     1           (ISTR2(JC2) + ISGN2(JC2)*I2) * IKD (JC2, 2) +
     2           (ISTR2(JC3) + ISGN2(JC3)*I3) * IKD (JC3, 2) 
            KN = INR                          * IKD (JC1, 3) +
     1           (ISTR2(JC2) + ISGN2(JC2)*I2) * IKD (JC2, 3) +
     2           (ISTR2(JC3) + ISGN2(JC3)*I3) * IKD (JC3, 3) 
C
 1000 format (' ',7i5)
            R1(JJ,KK,II,1) = R2(JN,KN,IN,1)
            R1(JJ,KK,II,2) = R2(JN,KN,IN,2)
            R1(JJ,KK,II,3) = R2(JN,KN,IN,3)
 700     CONTINUE
 800  CONTINUE
C
C     Finished with Grid Cut to set image points for metric evaluation
C
      RETURN
      END
