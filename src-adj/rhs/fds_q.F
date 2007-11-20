      SUBROUTINE FDS_Q (IFLUX, ITURB, NPTS, QP, QPB, QM, QMB, S, F, FB)
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
     4          QPB   (NPTS,NQ),
     5          QMB   (NPTS,NQ),
     6          FB    (NPTS,NF)

C     Local variables
      DIMENSION QP1(NQ), QM1(NQ), STMP(3), FTMP(NF)
      DIMENSION QP1B(NQ), QM1B(NQ), FTMPB(NF)

      IF (IFLUX .EQ. IFROE) THEN
         DO I=1,NPTS
            DO J=1,NQ
               QP1(J) = QP(I,J)
               QM1(J) = QM(I,J)
               QP1B(J)= 0.0
               QM1B(J)= 0.0
            ENDDO
            DO J=1,3
               STMP(J) = S(I,J)
            ENDDO
            DO J=1,NF
               FTMPB(J) = FB(I,J)
            ENDDO
            CALL FLUROE_BQ(ITURB, QP1, QP1B, QM1, QM1B, STMP, 
     1                     FTMP, FTMPB)
            DO J=1,NQ
               QMB(I, J) = QMB(I, J) + QM1B(J)
               QPB(I, J) = QPB(I, J) + QP1B(J)
            ENDDO
         ENDDO
      ENDIF
C
C     Finished calculating flux at the interface 
C
      RETURN
      END
