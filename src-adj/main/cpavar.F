      SUBROUTINE CPAVAR(IDIM, JDIM, KDIM, QA, RESB)
      include '../../src-flo/header/common.h'
      DIMENSION QA   (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          RESB (JDIM+1,KDIM+1,IDIM+1,NF)

      DO L=1,NQ
         DO J=2,JDIM
            DO K=2,KDIM
               DO I=2,IDIM
                  RESB(J,K,I,L) = QA(J,K,I,L)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END


