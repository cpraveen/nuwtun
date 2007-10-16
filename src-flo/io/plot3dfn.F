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
      SUBROUTINE PLOT3DFN (ITURB)
C
C This routine writes the PLOT3D function name file.
C
C ITURB : Turbulence model
C
C This routine is organized as follows:
C     1. Write function names for turbulence quantities
C        a. K-epsilon
C        b. K-omega
C        c. Reynolds stress
C     2. Write function names for properties
C
      include '../header/common.h'
C
C Local storage
C
      INTEGER      L
C
      CHARACTER*25 FUNNAME (5)
C
C 1.   Write function names for turbulence quantities
C 1-a. K-epsilon
C
      IF (ITURB .EQ. ITKE) THEN
         WRITE (IPLT3FN,1000) 'K                        ',
     1                        'Epsilon                  '
C
C 1-b. K-omega
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         WRITE (IPLT3FN,1000) 'K                        ',
     1                        'Omega                    '
C
C 1-c. Reynolds stress
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         WRITE (IPLT3FN,1000) 'Tau_xx                   ',
     1                        'Tau_yy                   ',
     2                        'Tau_zz                   ',
     3                        'Tau_xy                   ',
     4                        'Tau_xz                   ',
     5                        'Tau_yz                   ',
     6                        'Epsilon                  '
C
      ENDIF
C
C 2.   Write function names for properties
C
      FUNNAME(1) = 'mu'
      FUNNAME(2) = 'mu_t'
      FUNNAME(3) = 'y+'
      FUNNAME(4) = 'y_normal'
      FUNNAME(5) = 'C_mu*'
C
      DO L = 1, NP
         WRITE (IPLT3FN,1000) FUNNAME(L)
      END DO
C
C Finished PLOT3DFN writing function name file
C
      RETURN
 1000 FORMAT (A25)
      END
