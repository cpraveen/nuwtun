      SUBROUTINE CVTC2P (IDIM, JDIM, KDIM, QC, QP)
C
C     Convert primitive to conserved variables. Note that this
C     subroutine is not called from anywhere else. It is only used to
C     generate CVTP2C_Q.
C
C     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
C     Later, we should include TKE in the energy term
C     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
C
      include '../../src-flo/header/common.h'
C
      DIMENSION QC    (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QP    (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C 
C 1. Convert mean flow from primitive to 'conserved' variables 
C    NOTE: do not include TKE in conversion from pressure to total energy
C          to match the equation of state that is assumed in FAST/PLOT3D
C 
      DO 130 I = 2, IDIM
         DO 120 K = 2, KDIM
            DO 110 J = 2, JDIM
               QP(J,K,I,1) = QC(J,K,I,1)
               QP(J,K,I,2) = QC(J,K,I,2)/QC(J,K,I,1)
               QP(J,K,I,3) = QC(J,K,I,3)/QC(J,K,I,1)
               QP(J,K,I,4) = QC(J,K,I,4)/QC(J,K,I,1)
               QP(J,K,I,5) = GAMM1*(QC(J,K,I,5) - 
     1                               0.5*(QC(J,K,I,2)**2 +
     2                                    QC(J,K,I,3)**2 +
     3                                    QC(J,K,I,4)**2)/QC(J,K,I,1))
  110       CONTINUE 
  120    CONTINUE 
  130 CONTINUE 

      RETURN
      END
