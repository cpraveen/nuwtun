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
      SUBROUTINE RCENTR (IDIM, JDIM, KDIM, R, RC, IBCDAT, NBCS)
C
C Routine to average the positions of the vertices to get a cell centered 
C value of the grid for output purposes.
C
C This routine is organized as follows:
C      1. Calculate cell centers for (1:dim+1) in all directions
C         a. Corners
C         b. Edges
C         c. Faces
C         d. Interior
C      2. For wall functions, loop over all boundary conditions and
C         correct the cell center calculated above to account for
C         wall function spacing
C
C IDIM,JDIM,KDIM : Dimensions of current block
C R              : Position vector of grid vertices (x,y,z)
C RC             : Calculated position of grid at cell centers
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
C NBCS           : Number of boundary conditions for current block
C
      include '../header/common.h'
C
      DIMENSION R     (0:JDIM+1,0:KDIM+1,0:IDIM+1,3),
     1          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),
     2          IEND  (3)
C
C R  is the grid input at vertices
C RC is the grid at cell centers to be sent to PLOT3D
C
C 1. Calculate cell-centered grid
C
C 1-a. Corners
C
      DO 10 IDIR = 1, 3
         RC(1,1,1,IDIR)                = R(1,1,1,IDIR)
         RC(JDIM+1,1,1,IDIR)           = R(JDIM,1,1,IDIR)
         RC(1,1,IDIM+1,IDIR)           = R(1,1,IDIM,IDIR)
         RC(JDIM+1,1,IDIM+1,IDIR)      = R(JDIM,1,IDIM,IDIR)
         RC(1,KDIM+1,1,IDIR)           = R(1,KDIM,1,IDIR)
         RC(JDIM+1,KDIM+1,1,IDIR)      = R(JDIM,KDIM,1,IDIR)
         RC(1,KDIM+1,IDIM+1,IDIR)      = R(1,KDIM,IDIM,IDIR)
         RC(JDIM+1,KDIM+1,IDIM+1,IDIR) = R(JDIM,KDIM,IDIM,IDIR)
   10 CONTINUE
C
C 1-b. Edges
C
      DO 60 IDIR = 1, 3
         DO 30 I = 2, IDIM
            RC(1,1,I,IDIR)           = 0.5E0 * ( R(1,1,I-1,IDIR)
     1                                         + R(1,1,I,IDIR) )
            RC(JDIM+1,1,I,IDIR)      = 0.5E0 * ( R(JDIM,1,I-1,IDIR)
     1                                         + R(JDIM,1,I,IDIR) )
            RC(1,KDIM+1,I,IDIR)      = 0.5E0 * ( R(1,KDIM,I-1,IDIR)
     1                                         + R(1,KDIM,I,IDIR) )
            RC(JDIM+1,KDIM+1,I,IDIR) = 0.5E0 * ( R(JDIM,KDIM,I-1,IDIR)
     1                                         + R(JDIM,KDIM,I,IDIR) )
   30    CONTINUE
         DO 40 J = 2, JDIM
            RC(J,1,1,IDIR)           = 0.5E0 * ( R(J-1,1,1,IDIR)
     1                                         + R(J,1,1,IDIR) )
            RC(J,1,IDIM+1,IDIR)      = 0.5E0 * ( R(J-1,1,IDIM,IDIR)
     1                                         + R(J,1,IDIM,IDIR) )
            RC(J,KDIM+1,1,IDIR)      = 0.5E0 * ( R(J-1,KDIM,1,IDIR)
     1                                         + R(J,KDIM,1,IDIR) )
            RC(J,KDIM+1,IDIM+1,IDIR) = 0.5E0 * ( R(J-1,KDIM,IDIM,IDIR)
     1                                         + R(J,KDIM,IDIM,IDIR) )
   40    CONTINUE
         DO 50 K = 2, KDIM
            RC(1,K,1,IDIR)           = 0.5E0 * ( R(1,K-1,1,IDIR)
     1                                         + R(1,K,1,IDIR) )
            RC(1,K,IDIM+1,IDIR)      = 0.5E0 * ( R(1,K-1,IDIM,IDIR)
     1                                         + R(1,K,IDIM,IDIR) )
            RC(JDIM+1,K,1,IDIR)      = 0.5E0 * ( R(JDIM,K-1,1,IDIR)
     1                                         + R(JDIM,K,1,IDIR) )
            RC(JDIM+1,K,IDIM+1,IDIR) = 0.5E0 * ( R(JDIM,K-1,IDIM,IDIR)
     1                                         + R(JDIM,K,IDIM,IDIR) )
   50    CONTINUE
   60 CONTINUE
C
C 1-c. J = 1, J = JDIM+1 Planes
C
      DO 120 IDIR = 1, 3
         DO 110 I = 2, IDIM
            DO 100 K = 2, KDIM
               RC(1,K,I,IDIR)      = 0.25E0*
     1                  ( R(1,K,I,IDIR)   + R(1,K,I-1,IDIR)
     2                  + R(1,K-1,I,IDIR) + R(1,K-1,I-1,IDIR) )
               RC(JDIM+1,K,I,IDIR) = 0.25E0*
     1                  ( R(JDIM,K,I,IDIR)   + R(JDIM,K,I-1,IDIR)
     2                  + R(JDIM,K-1,I,IDIR) + R(JDIM,K-1,I-1,IDIR) )
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE
C
C      K = 1, K = KDIM+1 Planes
C
      DO 150 IDIR = 1, 3
         DO 140 I = 2, IDIM
            DO 130 J = 2, JDIM
               RC(J,1,I,IDIR)      = 0.25E0*
     1                  ( R(J,1,I,IDIR)   + R(J,1,I-1,IDIR)
     2                  + R(J-1,1,I,IDIR) + R(J-1,1,I-1,IDIR) )
               RC(J,KDIM+1,I,IDIR) = 0.25E0*
     1                  ( R(J,KDIM,I,IDIR)   + R(J,KDIM,I-1,IDIR)
     2                  + R(J-1,KDIM,I,IDIR) + R(J-1,KDIM,I-1,IDIR) )
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
C
C      I = 1, I = IDIM+1 Planes
C
      DO 180 IDIR = 1, 3
         DO 170 K = 2, KDIM
            DO 160 J = 2, JDIM
               RC(J,K,1,IDIR)      = 0.25E0*
     1                  ( R(J,K,1,IDIR)   + R(J-1,K,1,IDIR)
     2                  + R(J,K-1,1,IDIR) + R(J-1,K-1,1,IDIR) )
               RC(J,K,IDIM+1,IDIR) = 0.25E0*
     1                  ( R(J,K,IDIM,IDIR)   + R(J-1,K,IDIM,IDIR)
     2                  + R(J,K-1,IDIM,IDIR) + R(J-1,K-1,IDIM,IDIR) )
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
C
C 1-d. The interior of the domain
C
      DO 240 IDIR = 1, 3
         DO 230 I = 2, IDIM
            DO 220 K = 2, KDIM
               DO 210 J = 2, JDIM
                  RC(J,K,I,IDIR) = 0.125E0*
     1                     (  R(J,K,I,IDIR)     + R(J,K,I-1,IDIR)
     2                      + R(J-1,K,I,IDIR)   + R(J-1,K,I-1,IDIR)
     3                      + R(J,K-1,I,IDIR)   + R(J,K-1,I-1,IDIR)
     4                      + R(J-1,K-1,I,IDIR) + R(J-1,K-1,I-1,IDIR) )
  210          CONTINUE
  220       CONTINUE
  230    CONTINUE
  240 CONTINUE
C
C 2. Modify cell centers of wall function points to reflect proper location
C
      IF (IFWALF) THEN
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
         IJKDIM(1) = IDIM
         IJKDIM(2) = JDIM
         IJKDIM(3) = KDIM
C
C Loop over all bcs
C
         DO 400 IBC = 1, NBCS
            IBCTYP   = IBCDAT(1,IBC)
            IF (IBCTYP .EQ. IWALFN) THEN
C
C Store boundary condition data in mnemonic names
C
               IDIR     = IBCDAT(2,IBC)
               INROUT   = IBCDAT(3,IBC)
               ISTRT(1) = IBCDAT(4,IBC)
               IEND (1) = IBCDAT(5,IBC)
               ISTRT(2) = IBCDAT(6,IBC)
               IEND (2) = IBCDAT(7,IBC)
               ISTRT(3) = IBCDAT(8,IBC)
               IEND (3) = IBCDAT(9,IBC)
C
C Modify beginning, ending indices to extend boundary condition to edge/corner
C
               DO 300 L = 1, 3
                  IF (L .NE. IDIR) THEN
                     IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
                     IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L)+1
                  END IF
  300          CONTINUE
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
                  IBCN   = 2
                  IBCN2  = 3
                  SGNM   = 1.E0
               ELSE
                  IBCPT  = IJKDIM(IDIR) + 1
                  IBCN   = IJKDIM(IDIR) 
                  IBCN2  = IJKDIM(IDIR) - 1
                  SGNM   = - 1.E0
               ENDIF
C
               IC1 = ICY (IDIR, 1)
               IC2 = ICY (IDIR, 2)
               IC3 = ICY (IDIR, 3)
C
C Reset solution point location for wall functions
C
               DO 320 I3 = ISTRT(IC3), IEND(IC3)
                  DO 310 I2 = ISTRT(IC2), IEND(IC2)
C
                     IW   = IBCPT*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JW   = IBCPT*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KW   = IBCPT*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                     IN   = IBCN *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JN   = IBCN *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KN   = IBCN *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                     IN2  = IBCN2*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JN2  = IBCN2*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KN2  = IBCN2*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
C Calculate distance from 2 to wall, 3 to 2, 3 to interface
C
                     DS2TOW =SQRT((RC(JN ,KN ,IN ,1)-RC(JW,KW,IW,1))**2
     1                           +(RC(JN ,KN ,IN ,2)-RC(JW,KW,IW,2))**2
     2                           +(RC(JN ,KN ,IN ,3)-RC(JW,KW,IW,3))**2)
                     DS3TO2 =SQRT((RC(JN2,KN2,IN2,1)-RC(JN,KN,IN,1))**2
     1                           +(RC(JN2,KN2,IN2,2)-RC(JN,KN,IN,2))**2
     2                           +(RC(JN2,KN2,IN2,3)-RC(JN,KN,IN,3))**2)
                     DS3TOI = (DS2TOW + DS3TO2) - 2.E0 * DS2TOW
C
C Linearly interpolate cell center of first cell to new location
C
                     FACT   = 2.E0 * DS3TOI / DS3TO2
                     OMFACT = 1.E0 - FACT
                     XC     = RC(JN,KN,IN,1)
                     YC     = RC(JN,KN,IN,2)
                     ZC     = RC(JN,KN,IN,3)
                     RC(JN,KN,IN,1) = FACT   * XC
     1                              + OMFACT * RC(JN2,KN2,IN2,1)
                     RC(JN,KN,IN,2) = FACT   * YC
     1                              + OMFACT * RC(JN2,KN2,IN2,2)
                     RC(JN,KN,IN,3) = FACT   * ZC
     1                              + OMFACT * RC(JN2,KN2,IN2,3)
  310             CONTINUE
  320          CONTINUE
            ENDIF
  400    CONTINUE
      ENDIF
C
C Finished with grid averaging
C
      RETURN
      END
