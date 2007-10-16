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
      SUBROUTINE FDS (IFLUX, ITURB, NPTS, QP, QM, S, F)
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
      include '../header/common.h'
C
      DIMENSION QP    (NPTS,NQ),
     1          QM    (NPTS,NQ),
     2          S     (NPTS,3),
     3          F     (NPTS,NF)

C     Local variables
      DIMENSION QP1(NQ), QM1(NQ), STMP(3), FTMP(NF)

      IF (IFLUX .EQ. IFROE) THEN
         DO I=1,NPTS
            DO J=1,NQ
               QP1(J) = QP(I,J)
               QM1(J) = QM(I,J)
            ENDDO
            DO J=1,3
               STMP(J) = S(I,J)
            ENDDO
            CALL FLUROE(ITURB, QP1, QM1, STMP, FTMP)
            DO J=1,NF
               F(I,J) = FTMP(J)
            ENDDO
         ENDDO
      ENDIF
C
C     Finished calculating flux at the interface 
C
      RETURN
      END
