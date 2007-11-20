C
      SUBROUTINE INTQA1( IDIM, JDIM, KDIM, QA )
      include '../../src-flo/header/common.h'
      DIMENSION QA (0:JDIM+2, 0:KDIM+2, 0:IDIM+2, NQ)

      DO L=1,NQ
         DO J=2,JDIM
            DO K=2,KDIM
               DO I=2,IDIM
                  QA(J,K,I,L) = 1.0
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
