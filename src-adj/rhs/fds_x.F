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
      SUBROUTINE FDS_X (IFLUX, ITURB, NPTS, QP, QM, S, SB, F, FB)
C
C     Routine to calculate interface fluxes by calling the 
C     appropriate Flux Difference Splitting scheme.
C
C
C     IFLUX          : Flux evaluation scheme
C     ITURB          : Turbulence model
C     NPTS           : Size of array to calculate fluxes on.
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     S              : Metrics
C     F              : Flux at the interface
C
      include '../../src-flo/header/common.h'
C
      DIMENSION QP    (NPTS,NQ),
     1          QM    (NPTS,NQ),
     2          S     (NPTS,3),
     3          F     (NPTS,NF),
     4          SB    (NPTS,3),
     5          FB    (NPTS,NF)
C
C     Local variables
C
      DIMENSION QP1(NQ), QM1(NQ), STMP(3), FTMP(NF)
      DIMENSION STMPB(3), FTMPB(NF)

      IF (IFLUX .EQ. IFROE) THEN
         DO I=1,NPTS
            QP1(:)   = QP(I,:)
            QM1(:)   = QM(I,:)
            STMP(:)  = S(I,:)
            STMPB(:) = 0.0
            FTMPB(:) = FB(I,:)
            CALL FLUROE_BX(ITURB, QP1, QM1, STMP, STMPB, FTMP, FTMPB)
            SB(I,:)  = SB(I,:) + STMPB(:)
         ENDDO
      ENDIF
C
C     Finished calculating flux at the interface 
C
      RETURN
      END
