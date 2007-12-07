      SUBROUTINE sumgrad(IDIM, JDIM, KDIM, R, RB, grad)
      DIMENSION R  (0:JDIM+1,0:KDIM+1,0:IDIM+1,3),
     1          RB (0:JDIM+1,0:KDIM+1,0:IDIM+1,3)

C     dJ/d(alpha) = y * dJ/dx - x * dJ/dy
      DO I=1,IDIM
         DO J=1,JDIM
            DO K=1,KDIM
               grad = grad + R(J,K,I,2)*RB(J,K,I,1)
     1                     - R(J,K,I,1)*RB(J,K,I,2)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
