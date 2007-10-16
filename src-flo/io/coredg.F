C
C     This file is part of NuWTun, see <http://nuwtun.berlios.de> for
C     licensing conditions. This program is distributed WITHOUT ANY WARRANTY;
C     without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
C     PARTICULAR PURPOSE. 
C
      SUBROUTINE COREDG (IDIM, JDIM, KDIM, Q, NV)
C
C     Assign solution to corner point of a block for PLOT3D output. It
C     takes the solution in the cell adjacent to the corner point and
C     puts it for the corner point, basically zeroth order
C     interpolation. This is not continuous since the same corner point
C     in an adjacent block may have a slightly different solution.
C     For points along the block edges, the two adjacent face values are
C     averaged.
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NV)
C
      IF (THREED) THEN
         KBEG = 1
         KEND = KDIM + 1
         KDEL = KDIM
      ELSE
C        For 2-d, we only need solution in K=2 plane
         KBEG = 2
         KEND = 2
         KDEL = 1
      ENDIF
C
C     (IN,JN,KN) is the cell sharing the corner point (I,J,K)
C
      DO I = 1,IDIM+1,IDIM
         DO J = 1,JDIM+1,JDIM
            DO K = KBEG,KEND,KDEL
               IF(I.EQ.1)THEN
                  IN = I + 1
               ELSE
                  IN = I - 1
               ENDIF
               IF(J.EQ.1)THEN
                  JN = J + 1
               ELSE
                  JN = J - 1
               ENDIF
               IF (THREED) THEN
                  IF(K.EQ.1)THEN
                     KN = K + 1
                  ELSE
                     KN = K - 1
                  ENDIF
               ELSE
                  KN = 2
               ENDIF
               DO L=1,NV
                  Q(J,K,I,L) = Q(JN,KN,IN,L)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      IF (.NOT. THREED) RETURN
C
C     Edges: For 3-D, for each point on edge of block, we average the 
C     two adjacent face values
C
C     Edges along K
C
      DO I = 1,IDIM+1,IDIM
         DO J = 1,JDIM+1,JDIM
            DO K = 2,KDIM
               IF(I.EQ.1)THEN
                  IN = I + 1
               ELSE
                  IN = I - 1
               ENDIF
               IF(J.EQ.1)THEN
                  JN = J + 1
               ELSE
                  JN = J - 1
               ENDIF
               DO L=1,NV
                  Q(J,K,I,L) = 0.5E0 * (Q(J,K,IN,L) + Q(JN,K,I,L))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     Edges along J
C
      DO I = 1,IDIM+1,IDIM
         DO K = 1,KDIM+1,KDIM
            DO J = 2,JDIM
               IF(I.EQ.1)THEN
                  IN = I + 1
               ELSE
                  IN = I - 1
               ENDIF
               IF(K.EQ.1)THEN
                  KN = K + 1
               ELSE
                  KN = K - 1
               ENDIF
               DO L=1,NV
                  Q(J,K,I,L) = 0.5E0 * (Q(J,K,IN,L) + Q(J,KN,I,L))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     Edges along I
C
      DO J = 1,JDIM+1,JDIM
         DO K = 1,KDIM+1,KDIM
            DO I = 2,IDIM
               IF(J.EQ.1)THEN
                  JN = J + 1
               ELSE
                  JN = J - 1
               ENDIF
               IF(K.EQ.1)THEN
                  KN = K + 1
               ELSE
                  KN = K - 1
               ENDIF
               DO L=1,NV
                  Q(J,K,I,L) = 0.5E0 * (Q(JN,K,I,L) + Q(J,KN,I,L))
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
