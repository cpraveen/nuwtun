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
      SUBROUTINE BCPRFL (IDIM, JDIM, KDIM, IPSTRT, IPEND, NPTS,
     1                   Q, QBC, S, IBCDAT, IERROR)
C
C Routine to set boundary condition from input profile
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C NPTS           : Number of points in QBC array for this bc segment
C Q              : Primitive variables at cell centers
C QBC            : Boundary data from the PROFILE data
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
C IERROR         : Error code
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QBC   (NPTS,NQ),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION IBCDAT(10)
C
C     Store boundary condition data in mnemonic names
C     N.B. Do not modify ISTRT or IEND as they are used to calculate indices
C
      IBCTYP = IBCDAT(1)
      IDIR   = IBCDAT(2)
      INROUT = IBCDAT(3)
      ISTRT  = IBCDAT(4)
      IEND   = IBCDAT(5)
      JSTRT  = IBCDAT(6)
      JEND   = IBCDAT(7)
      KSTRT  = IBCDAT(8)
      KEND   = IBCDAT(9)
C
C     Set needed variables depending on whether the boundary is
C     the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C          IBC  is the boundary condition location
C          IBC2 is the boundary condition location outside the block from IBC
C          IN   is the point to the inside of the block from IBC
C
C
C     I Direction
C
      IF (IDIR .EQ. 1) THEN
         IF (INROUT .EQ. 1) THEN
            IBC  = 1
            IBC2 = 0
            IN   = 2
         ELSE
            IBC  = IDIM + 1
            IBC2 = IDIM + 2
            IN   = IDIM
         ENDIF
         DO 130 L = 1, NQ
            DO 120 K = KSTRT, KEND
               DO 110 J = JSTRT, JEND
                  JKBCPT = (J-JSTRT+1) + (K-KSTRT)*(JEND-JSTRT+1) 
                  Q(J,K,IBC,L)  = QBC(JKBCPT,L)
                  Q(J,K,IBC2,L) = Q(J,K,IBC,L)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
C
C     J Direction
C
      ELSE IF (IDIR .EQ. 2) THEN
         IF (INROUT .EQ. 1) THEN
            JBC  = 1
            JBC2 = 0
            JN   = 2
         ELSE
            JBC  = JDIM + 1
            JBC2 = JDIM + 2
            JN   = JDIM
         ENDIF
         DO 230 L = 1, NQ
            DO 220 I = MAX(ISTRT, IPSTRT), MIN(IEND, IPEND)
               DO 210 K = KSTRT, KEND
                  IKBCPT = (K-KSTRT+1) + (I-ISTRT)*(KEND-KSTRT+1) 
                  Q(JBC,K,I,L)  = QBC(IKBCPT,L)
                  Q(JBC2,K,I,L) = Q(JBC,K,I,L)
  210          CONTINUE
  220       CONTINUE
  230    CONTINUE
C
C     K Direction
C
      ELSE IF (IDIR .EQ. 3) THEN
         IF (INROUT .EQ. 1) THEN
            KBC  = 1
            KBC2 = 0
            KN   = 2
         ELSE
            KBC  = KDIM + 1
            KBC2 = KDIM + 2
            KN   = KDIM
         ENDIF
         DO 330 L = 1, NQ
            DO 320 I = MAX(ISTRT, IPSTRT), MIN(IEND, IPEND)
               DO 310 J = JSTRT, JEND
                  JIBCPT = (J-JSTRT+1) + (I-ISTRT)*(JEND-JSTRT+1) 
                  Q(J,KBC,I,L)  = QBC(JIBCPT,L)
                  Q(J,KBC2,I,L) = Q(J,KBC,I,L)
  310          CONTINUE
  320       CONTINUE
  330    CONTINUE
      ELSE
         WRITE (IOUT,1000) IDIR
         IERROR = IYES
      ENDIF
C
C     Finished with PROFILE boundary condition
C
      RETURN
 1000 FORMAT (' ','BC    : ERROR-> Illegal direction = ',I3,
     1            ' in BCPRFL.')
      END
