      SUBROUTINE ADDDPT(IDIM, JDIM, KDIM, QD, RD, QB, RB, QDQB, RDRB)
      include '../../src-flo/header/common.h'
      DIMENSION QD (0:JDIM+2, 0:KDIM+2, 0:IDIM+2, NQ),
     1          RD (0:JDIM+2, 0:KDIM+2, 0:IDIM+2, NQ),
     2          QB (0:JDIM+2, 0:KDIM+2, 0:IDIM+2, NQ),
     3          RB (0:JDIM+2, 0:KDIM+2, 0:IDIM+2, NQ)

      DO L=1,NQ
         DO J=2,JDIM
            DO K=2,KDIM
               DO I=2,IDIM
                  QDQB = QDQB + QD(J,K,I,L)*QB(J,K,I,L)
                  RDRB = RDRB + RD(J,K,I,L)*RB(J,K,I,L)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
