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
      SUBROUTINE RCCOPY (IDIM, JDIM, KDIM, RCIN, RCOUT)
C
C Routine to copy the input cell centered grid (RCIN) to the output cell
C centered grid (RCOUT).
C
C IDIM,JDIM,KDIM : Dimensions of current block
C RCIN           : Input  position vector of grid vertices at cell centers
C RCOUT          : Output position vector of grid vertices at cell centers
C RC             : Calculated position of grid at cell centers
C
      include '../header/common.h'
C
      DIMENSION RCIN  (JDIM+1,KDIM+1,IDIM+1,3),
     1          RCOUT (JDIM+1,KDIM+1,IDIM+1,3)
C
C Loop and copy
C
      DO 130 K = 1, KDIM+1
         DO 120 I = 1, IDIM+1
            DO 110 J = 1, JDIM+1
               RCOUT(J,K,I,1) = RCIN (J,K,I,1)
               RCOUT(J,K,I,2) = RCIN (J,K,I,2)
               RCOUT(J,K,I,3) = RCIN (J,K,I,3)
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C Finished copying
C
      RETURN
      END
