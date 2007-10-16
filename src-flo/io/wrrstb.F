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
      SUBROUTINE WRRSTB (ITURB, ITER, IBLOCK, ILVL, IDIM, JDIM, KDIM, 
     1                   Q, QN, PROPS, WRPROP, WRUNST)
C
C Routine to write the restart data for a block.
C
C N.B.  QN is converted to PRIMITIVE variables for restart file.
C
C This routine is organized as follows:
C     1. Write Flow Variables and Properties
C        a. Write Flow Variables
C        b. Write Properties
C        c. Write flow variables at additional time levels for unsteady flow
C
C ITURB          : Turbulence model
C ITER           : Number of iterations to write residual history
C IBLOCK         : Block number
C ILVL           : Grid level
C IDIM,JDIM,KDIM : Dimensions of current block
C Q              : Primitive variables stored at cell centers
C QN             : Primitive variables at time levels n, n-1, etc.
C                  N.B. Only NTIME-1 time levels are written/read on restart
C PROPS          : Properties at cell centers
C WRPROP         : Logical controlling whether to write property data to file
C                  WRPROP = .true.  Write property data
C                         = .false. Do not write property data
C WRUNST         : Logical controlling whether to write unsteady data to file
C                  WRUNST = .true.  Write unsteady data
C                         = .false. Do not write unsteady data
C
      include '../header/common.h'
      include '../header/histry.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QN    (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ,NTIME-1),
     2          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP)
C
C Control of whether or not to write properties data and additional time levels
C
      LOGICAL   WRPROP, WRUNST
C
C 1.   Write Flow Variables and Properties
C 1-a. Write Flow Variables
C
      WRITE (IWRRST) Q
C
C 1-b. Write Properties
C
      IF (WRPROP) THEN
         WRITE (IWRRST) PROPS
      ENDIF
C
C 1-c. Write flow variables at additional time levels for unsteady flow
C
      IF (WRUNST) THEN
         WRITE (IWRRST) QN
      ENDIF
C
      RETURN
      END
