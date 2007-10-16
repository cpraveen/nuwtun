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
      SUBROUTINE RDPRFL (NPTS, QBC, IBCDAT, IERRCD)
C
C Routine to read in a profile and store it into the QBC array to
C be used in the boundary condition routine BCPROF.
C
C NPTS           : Number of points in the QBC array for current segment
C QBC            : Primitive variables read from a file
C IBCDAT         : Data controlling boundary conditions:
C                  IBCDAT( 1,ibc) = Boundary Condition Type
C                          2      = Direction (1=i, 2=j, 3=k)
C                          3      = i(1 or idim)
C                          4      = istrt
C                          5      = iend
C                          6      = jstrt
C                          7      = jend
C                          8      = kstrt
C                          9      = kend
C                         10      = profile identifier
C IERRCD         : Error Code
C                  IERRCD = INO  if NO storage error or error reading file
C                  IERRCD = IYES if    storage error or error reading file
C
      include '../header/common.h'
C
      DIMENSION QBC   (NPTS,NQ)
C
      DIMENSION IBCDAT(10)
C
C     Store boundary condition data in mnemonic names
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
C     I Direction
C
      IF (IDIR .EQ. 1) THEN
         WRITE (IOUT,1010) 'K', KSTRT, KEND
         WRITE (IOUT,1010) 'J', JSTRT, JEND
         DO 120 K = KSTRT, KEND
            DO 110 J = JSTRT, JEND
               JKBCPT = (J-JSTRT+1) + (K-KSTRT)*(JEND-JSTRT+1)
               READ (IRDBC,*,END=900,ERR=910) (QBC(JKBCPT,L),L=1,NQ)
  110       CONTINUE
  120    CONTINUE
C
C     J Direction
C
      ELSE IF (IDIR .EQ. 2) THEN
         WRITE (IOUT,1010) 'I', ISTRT, IEND
         WRITE (IOUT,1010) 'K', KSTRT, KEND
         DO 220 I = ISTRT, IEND
            DO 210 K = KSTRT, KEND
               IKBCPT = (K-KSTRT+1) + (I-ISTRT)*(KEND-KSTRT+1)
               READ (IRDBC,*,END=900,ERR=910) (QBC(IKBCPT,L),L=1,NQ)
  210       CONTINUE
  220    CONTINUE
C
C     K Direction
C
      ELSE IF (IDIR .EQ. 3) THEN
         WRITE (IOUT,1010) 'I', ISTRT, IEND
         WRITE (IOUT,1010) 'J', JSTRT, JEND
         DO 320 I = ISTRT, IEND
            DO 310 J = JSTRT, JEND
               IJBCPT = (J-JSTRT+1) + (I-ISTRT)*(JEND-JSTRT+1)
               READ (IRDBC,*,END=900,ERR=910) (QBC(IJBCPT,L),L=1,NQ)
  310       CONTINUE
  320    CONTINUE
      ELSE
         IERRCD = IYES
         WRITE (IOUT,1200) IDIR
         RETURN
      ENDIF
C
C     Finished reading profile boundary data
C
      RETURN
C
C 2.   Error handling section
C 2-a. End-Of-File Handler
C
  900 CONTINUE
      IERRCD = IYES
      WRITE (IOUT,1205)
      RETURN
C
C 2-b. Error Handler
C
  910 CONTINUE
      IERRCD = IYES
      WRITE (IOUT,1210)
      RETURN
C
 1010 FORMAT (' ','RDPRFL: Looping over ',A1,' from ',I4,' to ',I4)
 1200 FORMAT (' ','RDPRFL: ERROR-> Invalid direction specified ',I4,
     1       /' ', 8X,'RUN ABORTING!')
 1205 FORMAT (' ','RDPRFL: ERROR-> Unexpected End-Of-File Encountered ',
     1       /' ', 8X,'RUN ABORTING!')
 1210 FORMAT (' ','RDPRFL: ERROR-> Unknown Error Encountered ',
     1       /' ', 8X,'RUN ABORTING!')
C
      END
