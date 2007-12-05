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


      SUBROUTINE QBCMOD2 (IDIM, JDIM, KDIM, Q, NBCS, IBCDAT)
C
C Routine to extrapolate Q to the ghost cell center using values from
C the first two cells inside the computational domain to then use in the 
C interpolation (prolongation) operator. This is modified from
C subroutine QBCMOD
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Loop over all boundary conditions
C        If boundary condition is a half cell, then
C        a. Replace cell face data with equivalent ghostcell data
C
C IDIM,JDIM,KDIM : Dimensions of current block
C Q              : Primitive variables at cell centers
C NBCS           : Number of boundary conditions for current block
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
      include '../../src-flo/header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION IBCDAT(10,NBCS)
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
C 2. Loop over boundary conditions
C
      DO 200 IBCS = 1, NBCS
         IBCTYP   = IBCDAT(1,IBCS)
C
         IF (IBCTYP .EQ. ITAN .OR. IBCTYP .EQ. IWALL .OR. 
     1       IBCTYP .EQ. IWALFN) THEN
C
C Store boundary condition data in mnemonic names
C
            IDIR     = IBCDAT(2,IBCS)
            INROUT   = IBCDAT(3,IBCS)
            ISTRT(1) = IBCDAT(4,IBCS)
            IEND (1) = IBCDAT(5,IBCS)
            ISTRT(2) = IBCDAT(6,IBCS)
            IEND (2) = IBCDAT(7,IBCS)
            ISTRT(3) = IBCDAT(8,IBCS)
            IEND (3) = IBCDAT(9,IBCS)
C
C Modify beginning, ending indices to extend boundary condition to edge/corner
C
            DO 10 L = 1, 3
               IF (L .NE. IDIR) THEN
                  IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
                  IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
               END IF
   10       CONTINUE
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCN   is the point to the inside of the block from IBCPT
C
            IF (INROUT .EQ. 1) THEN
               IBCPT  = 1
               IBCN   = 2
               IBCN2  = 3
            ELSE
               IBCPT  = IJKDIM(IDIR) + 1
               IBCN   = IJKDIM(IDIR)
               IBCN2  = IJKDIM(IDIR) - 1
            ENDIF
C
            IC1 = ICY (IDIR, 1)
            IC2 = ICY (IDIR, 2)
            IC3 = ICY (IDIR, 3)
C
C 2-a. Replace cell face data with equivalent ghostcell data
C
            DO 130 L = 1, NQ
               DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
                  DO 110 I2 = ISTRT(IC2), IEND(IC2)
C
                     IBC  = IBCPT*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JBC  = IBCPT*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KBC  = IBCPT*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                     IN   = IBCN *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JN   = IBCN *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KN   = IBCN *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                     IN2  = IBCN2*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JN2  = IBCN2*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KN2  = IBCN2*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                     DELTAQ = Q(JN2,KN2,IN2,L) - Q(JN ,KN ,IN ,L)
                     Q(JBC,KBC,IBC,L) = Q(JN,KN,IBC,L) - DELTAQ
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
         ENDIF
  200 CONTINUE
C
C Finished modifying the boundary conditions for interpolation (prolongation)
C
      RETURN
      END
