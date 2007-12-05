      SUBROUTINE GDXTND (IDIM, JDIM, KDIM, R)
C
C Routine to extend the grid by copying the information on I=1,IDIM
C J=1,JDIM and K=1,KDIM grid faces into additional storage.
C
C IDIM,JDIM,KDIM : Dimensions of current block
C R              : Position vector of grid vertices (x,y,z)
C
      include '../header/common.h'
C
      DIMENSION R     (0:JDIM+1,0:KDIM+1,0:IDIM+1,3)
C
C 1. Extend K face
C
      DO 130 L = 1, 3
         DO 120 I = 1, IDIM
            DO 110 J = 1, JDIM
               R(J,0     ,I,L) = R(J,1   ,I,L)
               R(J,KDIM+1,I,L) = R(J,KDIM,I,L)
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C 2. Extend J face
C
      DO 230 L = 1, 3
         DO 220 I = 1, IDIM
            DO 210 K = 0, KDIM+1
               R(0     ,K,I,L) = R(1   ,K,I,L)
               R(JDIM+1,K,I,L) = R(JDIM,K,I,L)
  210       CONTINUE
  220    CONTINUE
  230 CONTINUE
C
C 3. Extend I face
C
      DO 330 L = 1, 3
         DO 320 K = 0, KDIM+1
            DO 310 J = 0, JDIM+1
               R(J,K,0     ,L) = R(J,K,1   ,L)
               R(J,K,IDIM+1,L) = R(J,K,IDIM,L)
  310       CONTINUE
  320    CONTINUE
  330 CONTINUE
C
C Finished extending grid.
C
      RETURN
      END
