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
      SUBROUTINE RDP3D  (GRDBIN, IDIM, JDIM, KDIM, R, IEOF, IERRCD)
C
C Routine to read in one block of the grid in PLOT3D format (either formatted
C or binary).  The loops run over kdim ( jdim ( idim ) ) ).  
C
C This routine is organized as follows:
C     1. Read grid
C        a. Binary    (3D or 2D grid)
C        b. Formatted (3D or 2D grid)
C     2. If two-dimensional grid then copy input plane to plane 2
C     3. Error handling section
C        a. End-Of-File Handler
C        b. Error Handler
C
C GRDBIN         : Logical for binary (.TRUE.) or formatted (.FALSE.) file
C IDIM,JDIM,KDIM : Dimensions of current block
C R              : Position vector of grid vertices (x,y,z)
C IEOF           : End-of-file switch
C                  IEOF   = INO  if End-of-file NOT reached
C                  IEOF   = IYES if End-of-file WAS reached
C IERRCD         : Error Code
C                  IERRCD = INO  if NO error reading file
C                  IERRCD = IYES if    error reading file
C
      include '../header/common.h'
C
      DIMENSION R     (0:JDIM+1,0:KDIM+1,0:IDIM+1,3)
C
      LOGICAL GRDBIN
C
C 1.   Read grid
C 1-a. Binary   PLOT3D grid
C
      IF (GRDBIN) THEN
         IF (THREED) THEN
            READ (IGRID,END=900,ERR=910) 
     1              ((((R(J,K,I,L),I=1,IDIM),J=1,JDIM),K=1,KDIM),L=1,3)
         ELSE
            READ (IGRID,END=900,ERR=910) 
     1               (((R(J,1,I,L),I=1,IDIM),J=1,JDIM),L=1,2)
         ENDIF
C
C 1-b. Formatted PLOT3D grid
C
      ELSE
         IF (THREED) THEN
            READ (IGRID,*,END=900,ERR=910) 
     1              ((((R(J,K,I,L),I=1,IDIM),J=1,JDIM),K=1,KDIM),L=1,3)
         ELSE
            READ (IGRID,*,END=900,ERR=910) 
     1               (((R(J,1,I,L),I=1,IDIM),J=1,JDIM),L=1,2)
         ENDIF
      ENDIF
C
C 2.   Axisymmetric or Two-Dimensional Grid Modifications
C 2-a. Axisymmetric: Rotate K=1, K=2 grid planes by AXIDTH degrees
C
      IF (AXISYM) THEN
         PI     = ACOS(-1.E0)
         STHETA = SIN (AXIDTH * PI / 180.E0)
         DO 220 I = 1, IDIM
            DO 210 J = 1, JDIM
               Z          =   R(J,1,I,2) * STHETA
               R(J,1,I,3) = - Z
               R(J,2,I,3) =   Z
               R(J,2,I,1) =   R(J,1,I,1)
               R(J,2,I,2) =   R(J,1,I,2)
 210        CONTINUE
 220     CONTINUE
C
C 2-b. Two-Dimensional: Copy input plane to plane 2 
C
      ELSE IF (.NOT. THREED) THEN
         DO 240 I = 1, IDIM
            DO 230 J = 1, JDIM
               R(J,1,I,3) = 0.E0
               R(J,2,I,3) = 1.E0
               R(J,2,I,1) = R(J,1,I,1)
               R(J,2,I,2) = R(J,1,I,2)
 230        CONTINUE
 240     CONTINUE
      ENDIF
C
C     Succesful completion reading block grid - return
C
      RETURN
C
C 3.   Error handling section
C 3-a. End-Of-File Handler
C
  900 CONTINUE
      IEOF   = IYES
      WRITE (IOUT,1000) 
      RETURN
C
C 3-b. Error Handler
C
  910 CONTINUE
      IERRCD = IYES
      WRITE (IOUT,1005) 
      RETURN
C
 1000 FORMAT (' ','RDP3D : ERROR-> Unexpected End-Of-File Encountered ',
     1       /' ', 8X,'RUN ABORTING!')
 1005 FORMAT (' ','RDP3D : ERROR-> Unknown Error Encountered ',
     1       /' ', 8X,'RUN ABORTING!')
C
      END
