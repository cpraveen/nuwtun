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

      SUBROUTINE VBTRI (NPTS, NSIM, IRANK, A, B, C, R) 
C
C     Routine to perform the Block Tridiagonal Inversion
C     a block tridiagonal system where A, B and C contain
C     the (i-1), (i) and (i+1) of the block tridiagonal system
C     and R is the right hand side which is returned as the solution
C     vector.  A, B and C are overwritten.
C
C     NPTS:      The number of 
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     IRANK:     The rank of the system to invert
C     A:         The supra-diagonal block 7x7 of L
C     B:         The lu decomposition of the diagonal block 7x7 of L
C     C:         The super-diagonal block 7x7 of U
C     R:         The right hand side
C
      include '../header/common.h'
C
      DIMENSION A (NSIM,NPTS,NRANK,NRANK),
     1          B (NSIM,NPTS,NRANK,NRANK),
     2          C (NSIM,NPTS,NRANK,NRANK),
     3          R (NSIM,NPTS,NRANK)
C
      IF (IRANK .EQ. 5) THEN
         CALL VB5TRI (NPTS, NSIM, A, B, C, R)
      ELSE IF (IRANK .EQ. 7) THEN
         CALL VB7TRI (NPTS, NSIM, NRANK, A, B, C, R)
      ELSE IF (IRANK .EQ. 12) THEN
         CALL VB12TR (NPTS, NSIM, NRANK, A, B, C, R)
      ELSE 
         WRITE (IOUT,1000) IRANK, IRANK
         STOP
      ENDIF
C
C     Finished 
C
      RETURN
 1000 FORMAT (' ','Block solver required for a ',I3,' x ',I3,' block.',
     1       /' ','This solver does not exist.  Run Aborting.')
      END






      SUBROUTINE VB5TRI (NPTS, NSIM, A, B, C, R)
C
C     Routine to invert a 5x5 block tridiagonal system
C     Inversion is by Upper Triangularization with a
C     unity diagonal followed by a backsubstitution step.
C     Vectorization is attained by inverting multiple systems
C     at one time with the recursive inverrsion algorithm.
C     The A, B, C arrays are destroyed in the solution
C     procedure.
C
C     NPTS:      The number of 
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     A:         The supra-diagonal block 5x5
C     B:         The       diagonal block 5x5
C     C:         The super-diagonal block 5x5
C     R:         The right hand side
C
C     R contains the solution vector at return
C
C     The storage of the system is such that vectorization occurs over
C     the first index to provide a unit stride for maximum speed.
C
      DIMENSION A (NSIM,NPTS,5,5),
     1          B (NSIM,NPTS,5,5),
     2          C (NSIM,NPTS,5,5),
     3          R (NSIM,NPTS,5)
C
C     Storage for the LU decomposition of a block 5x5
C
      REAL L11, 
     1     L21, L22, 
     2     L31, L32, L33, 
     3     L41, L42, L43, L44, 
     4     L51, L52, L53, L54, L55
C
C     Forward Substitution Loop - Upper Triangularization
C     Loop only to NPTS-1 to optimize vectorization while avoiding unnecessary
C     calculations of the L terms
C
      DO 200 J = 1, NPTS - 1
C
C     Eliminate A
C
C     B = B - A*C(j-1)
C
         IF (J .GT. 1) THEN
            DO 160 L = 1, 5
               DO 140 M = 1, 5
                  DO 130 I = 1, NSIM
                     B(I,J,L,M) = B(I,J,L,M)
     1                          - A(I,J,L,1)*C(I,J-1,1,M) 
     2                          - A(I,J,L,2)*C(I,J-1,2,M)
     3                          - A(I,J,L,3)*C(I,J-1,3,M) 
     4                          - A(I,J,L,4)*C(I,J-1,4,M)
     5                          - A(I,J,L,5)*C(I,J-1,5,M)
  130          CONTINUE
  140       CONTINUE
C
C     R = R - A*R(j-1)
C
               DO 150 I = 1, NSIM
                  R(I,J,L) = R(I,J,L) 
     1                     - A(I,J,L,1)*R(I,J-1,1) 
     2                     - A(I,J,L,2)*R(I,J-1,2) 
     3                     - A(I,J,L,3)*R(I,J-1,3) 
     4                     - A(I,J,L,4)*R(I,J-1,4) 
     5                     - A(I,J,L,5)*R(I,J-1,5)
  150          CONTINUE
  160       CONTINUE
         ENDIF
C
C     Normalize row by multiplication of B(-1) to get unity diagonal
C
C     Form the LU of B(1)
C
         DO 190 I = 1, NSIM
            L11 =  B(I,J,1,1)
            U12 =  B(I,J,1,2) / L11
            U13 =  B(I,J,1,3) / L11
            U14 =  B(I,J,1,4) / L11
            U15 =  B(I,J,1,5) / L11
C
            L21 =  B(I,J,2,1) 
            L22 =  B(I,J,2,2) - L21*U12
            U23 = (B(I,J,2,3) - L21*U13) / L22
            U24 = (B(I,J,2,4) - L21*U14) / L22
            U25 = (B(I,J,2,5) - L21*U15) / L22
C
            L31 =  B(I,J,3,1)
            L32 =  B(I,J,3,2) - L31*U12
            L33 =  B(I,J,3,3) - L31*U13 - L32*U23
            U34 = (B(I,J,3,4) - L31*U14 - L32*U24) / L33
            U35 = (B(I,J,3,5) - L31*U15 - L32*U25) / L33
C
            L41 =  B(I,J,4,1)
            L42 =  B(I,J,4,2) - L41*U12
            L43 =  B(I,J,4,3) - L41*U13 - L42*U23
            L44 =  B(I,J,4,4) - L41*U14 - L42*U24 - L43*U34
            U45 = (B(I,J,4,5) - L41*U15 - L42*U25 - L43*U35) / L44
C
            L51 =  B(I,J,5,1)
            L52 =  B(I,J,5,2) - L51*U12
            L53 =  B(I,J,5,3) - L51*U13 - L52*U23
            L54 =  B(I,J,5,4) - L51*U14 - L52*U24 - L53*U34
            L55 =  B(I,J,5,5) - L51*U15 - L52*U25 - L53*U35 - L54*U45
C
C     Calculate B(-1) * R
C
            R1       =  R(I,J,1) / L11
            R2       = (R(I,J,2) - L21*R1) / L22
            R3       = (R(I,J,3) - L31*R1 - L32*R2) / L33
            R4       = (R(I,J,4) - L41*R1 - L42*R2 - L43*R3) / L44
            R(I,J,5) = (R(I,J,5) - L51*R1 - L52*R2 - L53*R3 - L54*R4)
     1                 / L55
            R(I,J,4) =  R4-U45*R(I,J,5)
            R(I,J,3) =  R3-U35*R(I,J,5)-U34*R(I,J,4)
            R(I,J,2) =  R2-U25*R(I,J,5)-U24*R(I,J,4)-U23*R(I,J,3)
            R(I,J,1) =  R1-U15*R(I,J,5)-U14*R(I,J,4)-U13*R(I,J,3)
     1                    -U12*R(I,J,2)
C
C     Calculate B(-1) * C
C
            C1         =  C(I,J,1,1) / L11
            C2         = (C(I,J,2,1) - L21*C1) / L22
            C3         = (C(I,J,3,1) - L31*C1 - L32*C2) / L33
            C4         = (C(I,J,4,1) - L41*C1 - L42*C2 - L43*C3)/L44
            C(I,J,5,1) = (C(I,J,5,1) - L51*C1 - L52*C2 - L53*C3 -
     1                                                   L54*C4)/L55
            C(I,J,4,1) =  C4-U45*C(I,J,5,1)
            C(I,J,3,1) =  C3-U35*C(I,J,5,1)-U34*C(I,J,4,1)
            C(I,J,2,1) =  C2-U25*C(I,J,5,1)-U24*C(I,J,4,1)
     1                      -U23*C(I,J,3,1)
            C(I,J,1,1) =  C1-U15*C(I,J,5,1)-U14*C(I,J,4,1)
     1                      -U13*C(I,J,3,1)-U12*C(I,J,2,1)
C
            C1         =  C(I,J,1,2) / L11
            C2         = (C(I,J,2,2) - L21*C1) / L22
            C3         = (C(I,J,3,2) - L31*C1 - L32*C2) / L33
            C4         = (C(I,J,4,2) - L41*C1 - L42*C2 - L43*C3)/L44
            C(I,J,5,2) = (C(I,J,5,2) - L51*C1 - L52*C2 - L53*C3 -
     1                                                   L54*C4)/L55
            C(I,J,4,2) =  C4-U45*C(I,J,5,2)
            C(I,J,3,2) =  C3-U35*C(I,J,5,2)-U34*C(I,J,4,2)
            C(I,J,2,2) =  C2-U25*C(I,J,5,2)-U24*C(I,J,4,2)
     1                      -U23*C(I,J,3,2)
            C(I,J,1,2) =  C1-U15*C(I,J,5,2)-U14*C(I,J,4,2)
     1                      -U13*C(I,J,3,2)-U12*C(I,J,2,2)
C
            C1         =  C(I,J,1,3) / L11
            C2         = (C(I,J,2,3) - L21*C1) / L22
            C3         = (C(I,J,3,3) - L31*C1 - L32*C2) / L33
            C4         = (C(I,J,4,3) - L41*C1 - L42*C2 - L43*C3)/L44
            C(I,J,5,3) = (C(I,J,5,3) - L51*C1 - L52*C2 - L53*C3 -
     1                                                   L54*C4)/L55
            C(I,J,4,3) =  C4-U45*C(I,J,5,3)
            C(I,J,3,3) =  C3-U35*C(I,J,5,3)-U34*C(I,J,4,3)
            C(I,J,2,3) =  C2-U25*C(I,J,5,3)-U24*C(I,J,4,3)
     1                      -U23*C(I,J,3,3)
            C(I,J,1,3) =  C1-U15*C(I,J,5,3)-U14*C(I,J,4,3)
     1                      -U13*C(I,J,3,3)-U12*C(I,J,2,3)
C
            C1         =  C(I,J,1,4) / L11
            C2         = (C(I,J,2,4) - L21*C1) / L22
            C3         = (C(I,J,3,4) - L31*C1 - L32*C2) / L33
            C4         = (C(I,J,4,4) - L41*C1 - L42*C2 - L43*C3)/L44
            C(I,J,5,4) = (C(I,J,5,4) - L51*C1 - L52*C2 - L53*C3 -
     1                                                   L54*C4)/L55
            C(I,J,4,4) =  C4-U45*C(I,J,5,4)
            C(I,J,3,4) =  C3-U35*C(I,J,5,4)-U34*C(I,J,4,4)
            C(I,J,2,4) =  C2-U25*C(I,J,5,4)-U24*C(I,J,4,4)
     1                      -U23*C(I,J,3,4)
            C(I,J,1,4) =  C1-U15*C(I,J,5,4)-U14*C(I,J,4,4)
     1                      -U13*C(I,J,3,4)-U12*C(I,J,2,4)
C
            C1         =  C(I,J,1,5) / L11
            C2         = (C(I,J,2,5) - L21*C1) / L22
            C3         = (C(I,J,3,5) - L31*C1 - L32*C2) / L33
            C4         = (C(I,J,4,5) - L41*C1 - L42*C2 - L43*C3)/L44
            C(I,J,5,5) = (C(I,J,5,5) - L51*C1 - L52*C2 - L53*C3 -
     1                                                   L54*C4)/L55
            C(I,J,4,5) =  C4-U45*C(I,J,5,5)
            C(I,J,3,5) =  C3-U35*C(I,J,5,5)-U34*C(I,J,4,5)
            C(I,J,2,5) =  C2-U25*C(I,J,5,5)-U24*C(I,J,4,5)
     1                      -U23*C(I,J,3,5)
            C(I,J,1,5) =  C1-U15*C(I,J,5,5)-U14*C(I,J,4,5)
     1                      -U13*C(I,J,3,5)-U12*C(I,J,2,5)
  190    CONTINUE
  200 CONTINUE
C
C     J = NPTS Elimination
C     Do not do this loop if NPTS = 1
C
C     Eliminate A
C
C     B = B - A*C(j-1)
C
      J = NPTS
      IF (NPTS .GT. 1) THEN
         DO 260 L = 1, 5
            DO 240 M = 1, 5
               DO 230 I = 1, NSIM
                  B(I,J,L,M) = B(I,J,L,M)
     1                       - A(I,J,L,1)*C(I,J-1,1,M) 
     2                       - A(I,J,L,2)*C(I,J-1,2,M)
     3                       - A(I,J,L,3)*C(I,J-1,3,M) 
     4                       - A(I,J,L,4)*C(I,J-1,4,M)
     5                       - A(I,J,L,5)*C(I,J-1,5,M)
  230          CONTINUE
  240       CONTINUE
C
C     R = R - A*R(j-1)
C
            DO 250 I = 1, NSIM
               R(I,J,L) = R(I,J,L) 
     1                  - A(I,J,L,1)*R(I,J-1,1) 
     2                  - A(I,J,L,2)*R(I,J-1,2) 
     3                  - A(I,J,L,3)*R(I,J-1,3) 
     4                  - A(I,J,L,4)*R(I,J-1,4) 
     5                  - A(I,J,L,5)*R(I,J-1,5)
  250       CONTINUE
  260    CONTINUE
      ENDIF
C
C     Normalize row by multiplication of B(-1) to get unity diagonal
C
C     Form the LU of B(1)
C
      DO 300 I = 1, NSIM
         L11 =  B(I,J,1,1)
         U12 =  B(I,J,1,2) / L11
         U13 =  B(I,J,1,3) / L11
         U14 =  B(I,J,1,4) / L11
         U15 =  B(I,J,1,5) / L11
C
         L21 =  B(I,J,2,1) 
         L22 =  B(I,J,2,2) - L21*U12
         U23 = (B(I,J,2,3) - L21*U13) / L22
         U24 = (B(I,J,2,4) - L21*U14) / L22
         U25 = (B(I,J,2,5) - L21*U15) / L22
C
         L31 =  B(I,J,3,1)
         L32 =  B(I,J,3,2) - L31*U12
         L33 =  B(I,J,3,3) - L31*U13 - L32*U23
         U34 = (B(I,J,3,4) - L31*U14 - L32*U24) / L33
         U35 = (B(I,J,3,5) - L31*U15 - L32*U25) / L33
C
         L41 =  B(I,J,4,1)
         L42 =  B(I,J,4,2) - L41*U12
         L43 =  B(I,J,4,3) - L41*U13 - L42*U23
         L44 =  B(I,J,4,4) - L41*U14 - L42*U24 - L43*U34
         U45 = (B(I,J,4,5) - L41*U15 - L42*U25 - L43*U35) / L44
C
         L51 =  B(I,J,5,1)
         L52 =  B(I,J,5,2) - L51*U12
         L53 =  B(I,J,5,3) - L51*U13 - L52*U23
         L54 =  B(I,J,5,4) - L51*U14 - L52*U24 - L53*U34
         L55 =  B(I,J,5,5) - L51*U15 - L52*U25 - L53*U35 - L54*U45
C
C     Calculate B(-1) * R
C
         R1       =  R(I,J,1) / L11
         R2       = (R(I,J,2) - L21*R1) / L22
         R3       = (R(I,J,3) - L31*R1 - L32*R2) / L33
         R4       = (R(I,J,4) - L41*R1 - L42*R2 - L43*R3) / L44
         R(I,J,5) = (R(I,J,5) - L51*R1 - L52*R2 - L53*R3 - L54*R4)
     1              / L55
         R(I,J,4) =  R4-U45*R(I,J,5)
         R(I,J,3) =  R3-U35*R(I,J,5)-U34*R(I,J,4)
         R(I,J,2) =  R2-U25*R(I,J,5)-U24*R(I,J,4)-U23*R(I,J,3)
         R(I,J,1) =  R1-U15*R(I,J,5)-U14*R(I,J,4)-U13*R(I,J,3)
     1                    -U12*R(I,J,2)
  300 CONTINUE
C
C     Backward Substitution Sweep
C
      DO 400 J = NPTS-1, 1, -1
         DO 390 L = 1, 5
            DO 380 I = 1, NSIM
               R(I,J,L) = R(I,J,L)
     1                  - C(I,J,L,1)*R(I,J+1,1) 
     2                  - C(I,J,L,2)*R(I,J+1,2)
     3                  - C(I,J,L,3)*R(I,J+1,3) 
     4                  - C(I,J,L,4)*R(I,J+1,4)
     5                  - C(I,J,L,5)*R(I,J+1,5)
  380       CONTINUE
  390    CONTINUE
  400 CONTINUE
C
C     Finished
C
      RETURN
      END






      SUBROUTINE VB7TRI (NPTS, NSIM, NRANK, A, B, C, R)
C
C     Routine to invert a 7x7 block tridiagonal system
C     Inversion is by Upper Triangularization with a
C     unity diagonal followed by a backsubstitution step.
C     Vectorization is attained by inverting multiple systems
C     at one time with the recursive inverrsion algorithm.
C     The A, B, C arrays are destroyed in the solution
C     procedure.
C
C     NPTS:      The number of 
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     NRANK:     The actual dimension of the storage
C     A:         The supra-diagonal block 7x7
C     B:         The       diagonal block 7x7
C     C:         The super-diagonal block 7x7
C     R:         The right hand side
C
C     R contains the solution vector at return
C
C     The storage of the system is such that vectorization occurs over
C     the first index to provide a unit stride for maximum speed.
C
      DIMENSION A (NSIM,NPTS,NRANK,NRANK),
     1          B (NSIM,NPTS,NRANK,NRANK),
     2          C (NSIM,NPTS,NRANK,NRANK),
     3          R (NSIM,NPTS,NRANK)
C
C     Storage for the LU decomposition of a block 7x7
C
      REAL L11, 
     1     L21, L22, 
     2     L31, L32, L33, 
     3     L41, L42, L43, L44, 
     4     L51, L52, L53, L54, L55,
     5     L61, L62, L63, L64, L65, L66,
     6     L71, L72, L73, L74, L75, L76, L77
C
C     Forward Substitution Loop - Upper Triangularization
C     Loop only to NPTS-1 to optimize vectorization while avoiding unnecessary
C     calculations of the L terms
C
      DO 200 J = 1, NPTS - 1
C
C     Eliminate A
C
C     B = B - A*C(j-1)
C
         IF (J .GT. 1) THEN
            DO 160 L = 1, 7
               DO 140 M = 1, 7
                  DO 130 I = 1, NSIM
                     B(I,J,L,M) = B(I,J,L,M)
     1                          - A(I,J,L,1)*C(I,J-1,1,M) 
     2                          - A(I,J,L,2)*C(I,J-1,2,M)
     3                          - A(I,J,L,3)*C(I,J-1,3,M) 
     4                          - A(I,J,L,4)*C(I,J-1,4,M)
     5                          - A(I,J,L,5)*C(I,J-1,5,M)
     6                          - A(I,J,L,6)*C(I,J-1,6,M)
     7                          - A(I,J,L,7)*C(I,J-1,7,M)
  130             CONTINUE
  140          CONTINUE
C
C     R = R - A*R(j-1)
C
               DO 150 I = 1, NSIM
                  R(I,J,L) = R(I,J,L) 
     1                     - A(I,J,L,1)*R(I,J-1,1) 
     2                     - A(I,J,L,2)*R(I,J-1,2) 
     3                     - A(I,J,L,3)*R(I,J-1,3) 
     4                     - A(I,J,L,4)*R(I,J-1,4) 
     5                     - A(I,J,L,5)*R(I,J-1,5) 
     6                     - A(I,J,L,6)*R(I,J-1,6) 
     7                     - A(I,J,L,7)*R(I,J-1,7)
  150          CONTINUE
  160       CONTINUE
         ENDIF
C
C     Normalize row by multiplication of B(-1) to get unity diagonal
C
C     Form the LU of B(1)
C
         DO 190 I = 1, NSIM
            L11 =  B(I,J,1,1)
            U12 =  B(I,J,1,2) / L11
            U13 =  B(I,J,1,3) / L11
            U14 =  B(I,J,1,4) / L11
            U15 =  B(I,J,1,5) / L11
            U16 =  B(I,J,1,6) / L11
            U17 =  B(I,J,1,7) / L11
C
            L21 =  B(I,J,2,1) 
            L22 =  B(I,J,2,2) - L21*U12
            U23 = (B(I,J,2,3) - L21*U13) / L22
            U24 = (B(I,J,2,4) - L21*U14) / L22
            U25 = (B(I,J,2,5) - L21*U15) / L22
            U26 = (B(I,J,2,6) - L21*U16) / L22
            U27 = (B(I,J,2,7) - L21*U17) / L22
C
            L31 =  B(I,J,3,1)
            L32 =  B(I,J,3,2) - L31*U12
            L33 =  B(I,J,3,3) - L31*U13 - L32*U23
            U34 = (B(I,J,3,4) - L31*U14 - L32*U24) / L33
            U35 = (B(I,J,3,5) - L31*U15 - L32*U25) / L33
            U36 = (B(I,J,3,6) - L31*U16 - L32*U26) / L33
            U37 = (B(I,J,3,7) - L31*U17 - L32*U27) / L33
C
            L41 =  B(I,J,4,1)
            L42 =  B(I,J,4,2) - L41*U12
            L43 =  B(I,J,4,3) - L41*U13 - L42*U23
            L44 =  B(I,J,4,4) - L41*U14 - L42*U24 - L43*U34
            U45 = (B(I,J,4,5) - L41*U15 - L42*U25 - L43*U35) / L44
            U46 = (B(I,J,4,6) - L41*U16 - L42*U26 - L43*U36) / L44
            U47 = (B(I,J,4,7) - L41*U17 - L42*U27 - L43*U37) / L44
C
            L51 =  B(I,J,5,1)
            L52 =  B(I,J,5,2) - L51*U12
            L53 =  B(I,J,5,3) - L51*U13 - L52*U23
            L54 =  B(I,J,5,4) - L51*U14 - L52*U24 - L53*U34
            L55 =  B(I,J,5,5) - L51*U15 - L52*U25 - L53*U35 - L54*U45
            U56 = (B(I,J,5,6) - L51*U16 - L52*U26 - L53*U36 - L54*U46) 
     1            / L55
            U57 = (B(I,J,5,7) - L51*U17 - L52*U27 - L53*U37 - L54*U47) 
     1            / L55
C
            L61 =  B(I,J,6,1)
            L62 =  B(I,J,6,2) - L61*U12
            L63 =  B(I,J,6,3) - L61*U13 - L62*U23
            L64 =  B(I,J,6,4) - L61*U14 - L62*U24 - L63*U34
            L65 =  B(I,J,6,5) - L61*U15 - L62*U25 - L63*U35 - L64*U45
            L66 =  B(I,J,6,6) - L61*U16 - L62*U26 - L63*U36 - L64*U46
     1                        - L65*U56
            U67 = (B(I,J,6,7) - L61*U17 - L62*U27 - L63*U37 - L64*U47
     1                        - L65*U57) / L66
C
            L71 =  B(I,J,7,1)
            L72 =  B(I,J,7,2) - L71*U12
            L73 =  B(I,J,7,3) - L71*U13 - L72*U23
            L74 =  B(I,J,7,4) - L71*U14 - L72*U24 - L73*U34
            L75 =  B(I,J,7,5) - L71*U15 - L72*U25 - L73*U35 - L74*U45
            L76 =  B(I,J,7,6) - L71*U16 - L72*U26 - L73*U36 - L74*U46
     1                        - L75*U56
            L77 =  B(I,J,7,7) - L71*U17 - L72*U27 - L73*U37 - L74*U47
     1                        - L75*U57 - L76*U67
C
C     Calculate B(-1) * R
C
            R1       =  R(I,J,1) / L11
            R2       = (R(I,J,2) - L21*R1) / L22
            R3       = (R(I,J,3) - L31*R1 - L32*R2) / L33
            R4       = (R(I,J,4) - L41*R1 - L42*R2 - L43*R3) / L44
            R5       = (R(I,J,5) - L51*R1 - L52*R2 - L53*R3 - L54*R4)
     1                 / L55
            R6       = (R(I,J,6) - L61*R1 - L62*R2 - L63*R3 - L64*R4 
     1                           - L65*R5) / L66
            R(I,J,7) = (R(I,J,7) - L71*R1 - L72*R2 - L73*R3 - L74*R4
     1                           - L75*R5 - L76*R6) / L77
            R(I,J,6) =  R6-U67*R(I,J,7)
            R(I,J,5) =  R5-U57*R(I,J,7)-U56*R(I,J,6)
            R(I,J,4) =  R4-U47*R(I,J,7)-U46*R(I,J,6)-U45*R(I,J,5)
            R(I,J,3) =  R3-U37*R(I,J,7)-U36*R(I,J,6)-U35*R(I,J,5)
     1                    -U34*R(I,J,4)
            R(I,J,2) =  R2-U27*R(I,J,7)-U26*R(I,J,6)-U25*R(I,J,5)
     1                    -U24*R(I,J,4)-U23*R(I,J,3)
            R(I,J,1) =  R1-U17*R(I,J,7)-U16*R(I,J,6)-U15*R(I,J,5)
     1                    -U14*R(I,J,4)-U13*R(I,J,3)-U12*R(I,J,2)
C
C     Calculate B(-1) * C
C
            C1         =  C(I,J,1,1) / L11
            C2         = (C(I,J,2,1) - L21*C1) / L22
            C3         = (C(I,J,3,1) - L31*C1 - L32*C2) / L33
            C4         = (C(I,J,4,1) - L41*C1 - L42*C2 - L43*C3)/L44
            C5         = (C(I,J,5,1) - L51*C1 - L52*C2 - L53*C3
     1                               - L54*C4)/L55
            C6         = (C(I,J,6,1) - L61*C1 - L62*C2 - L63*C3
     1                               - L64*C4 - L65*C5)/L66
            C(I,J,7,1) = (C(I,J,7,1) - L71*C1 - L72*C2 - L73*C3
     1                               - L74*C4 - L75*C5 - L76*C6)/L77
            C(I,J,6,1) = C6-U67*C(I,J,7,1) 
            C(I,J,5,1) = C5-U57*C(I,J,7,1)-U56*C(I,J,6,1)
            C(I,J,4,1) = C4-U47*C(I,J,7,1)-U46*C(I,J,6,1)-U45*C(I,J,5,1)
            C(I,J,3,1) = C3-U37*C(I,J,7,1)-U36*C(I,J,6,1)-U35*C(I,J,5,1)
     1                     -U34*C(I,J,4,1)
            C(I,J,2,1) = C2-U27*C(I,J,7,1)-U26*C(I,J,6,1)-U25*C(I,J,5,1)
     1                     -U24*C(I,J,4,1)-U23*C(I,J,3,1)
            C(I,J,1,1) = C1-U17*C(I,J,7,1)-U16*C(I,J,6,1)-U15*C(I,J,5,1)
     1                     -U14*C(I,J,4,1)-U13*C(I,J,3,1)-U12*C(I,J,2,1)
C
            C1         =  C(I,J,1,2) / L11
            C2         = (C(I,J,2,2) - L21*C1) / L22
            C3         = (C(I,J,3,2) - L31*C1 - L32*C2) / L33
            C4         = (C(I,J,4,2) - L41*C1 - L42*C2 - L43*C3)/L44
            C5         = (C(I,J,5,2) - L51*C1 - L52*C2 - L53*C3
     1                               - L54*C4)/L55
            C6         = (C(I,J,6,2) - L61*C1 - L62*C2 - L63*C3
     1                               - L64*C4 - L65*C5)/L66
            C(I,J,7,2) = (C(I,J,7,2) - L71*C1 - L72*C2 - L73*C3
     1                               - L74*C4 - L75*C5 - L76*C6)/L77
            C(I,J,6,2) = C6-U67*C(I,J,7,2) 
            C(I,J,5,2) = C5-U57*C(I,J,7,2)-U56*C(I,J,6,2)
            C(I,J,4,2) = C4-U47*C(I,J,7,2)-U46*C(I,J,6,2)-U45*C(I,J,5,2)
            C(I,J,3,2) = C3-U37*C(I,J,7,2)-U36*C(I,J,6,2)-U35*C(I,J,5,2)
     1                     -U34*C(I,J,4,2)
            C(I,J,2,2) = C2-U27*C(I,J,7,2)-U26*C(I,J,6,2)-U25*C(I,J,5,2)
     1                     -U24*C(I,J,4,2)-U23*C(I,J,3,2)
            C(I,J,1,2) = C1-U17*C(I,J,7,2)-U16*C(I,J,6,2)-U15*C(I,J,5,2)
     1                     -U14*C(I,J,4,2)-U13*C(I,J,3,2)-U12*C(I,J,2,2)
C
            C1         =  C(I,J,1,3) / L11
            C2         = (C(I,J,2,3) - L21*C1) / L22
            C3         = (C(I,J,3,3) - L31*C1 - L32*C2) / L33
            C4         = (C(I,J,4,3) - L41*C1 - L42*C2 - L43*C3)/L44
            C5         = (C(I,J,5,3) - L51*C1 - L52*C2 - L53*C3
     1                               - L54*C4)/L55
            C6         = (C(I,J,6,3) - L61*C1 - L62*C2 - L63*C3
     1                               - L64*C4 - L65*C5)/L66
            C(I,J,7,3) = (C(I,J,7,3) - L71*C1 - L72*C2 - L73*C3
     1                               - L74*C4 - L75*C5 - L76*C6)/L77
            C(I,J,6,3) = C6-U67*C(I,J,7,3) 
            C(I,J,5,3) = C5-U57*C(I,J,7,3)-U56*C(I,J,6,3)
            C(I,J,4,3) = C4-U47*C(I,J,7,3)-U46*C(I,J,6,3)-U45*C(I,J,5,3)
            C(I,J,3,3) = C3-U37*C(I,J,7,3)-U36*C(I,J,6,3)-U35*C(I,J,5,3)
     1                     -U34*C(I,J,4,3)
            C(I,J,2,3) = C2-U27*C(I,J,7,3)-U26*C(I,J,6,3)-U25*C(I,J,5,3)
     1                     -U24*C(I,J,4,3)-U23*C(I,J,3,3)
            C(I,J,1,3) = C1-U17*C(I,J,7,3)-U16*C(I,J,6,3)-U15*C(I,J,5,3)
     1                     -U14*C(I,J,4,3)-U13*C(I,J,3,3)-U12*C(I,J,2,3)
C
            C1         =  C(I,J,1,4) / L11
            C2         = (C(I,J,2,4) - L21*C1) / L22
            C3         = (C(I,J,3,4) - L31*C1 - L32*C2) / L33
            C4         = (C(I,J,4,4) - L41*C1 - L42*C2 - L43*C3)/L44
            C5         = (C(I,J,5,4) - L51*C1 - L52*C2 - L53*C3
     1                               - L54*C4)/L55
            C6         = (C(I,J,6,4) - L61*C1 - L62*C2 - L63*C3
     1                               - L64*C4 - L65*C5)/L66
            C(I,J,7,4) = (C(I,J,7,4) - L71*C1 - L72*C2 - L73*C3
     1                               - L74*C4 - L75*C5 - L76*C6)/L77
            C(I,J,6,4) = C6-U67*C(I,J,7,4) 
            C(I,J,5,4) = C5-U57*C(I,J,7,4)-U56*C(I,J,6,4)
            C(I,J,4,4) = C4-U47*C(I,J,7,4)-U46*C(I,J,6,4)-U45*C(I,J,5,4)
            C(I,J,3,4) = C3-U37*C(I,J,7,4)-U36*C(I,J,6,4)-U35*C(I,J,5,4)
     1                     -U34*C(I,J,4,4)
            C(I,J,2,4) = C2-U27*C(I,J,7,4)-U26*C(I,J,6,4)-U25*C(I,J,5,4)
     1                     -U24*C(I,J,4,4)-U23*C(I,J,3,4)
            C(I,J,1,4) = C1-U17*C(I,J,7,4)-U16*C(I,J,6,4)-U15*C(I,J,5,4)
     1                     -U14*C(I,J,4,4)-U13*C(I,J,3,4)-U12*C(I,J,2,4)
C
            C1         =  C(I,J,1,5) / L11
            C2         = (C(I,J,2,5) - L21*C1) / L22
            C3         = (C(I,J,3,5) - L31*C1 - L32*C2) / L33
            C4         = (C(I,J,4,5) - L41*C1 - L42*C2 - L43*C3)/L44
            C5         = (C(I,J,5,5) - L51*C1 - L52*C2 - L53*C3
     1                               - L54*C4)/L55
            C6         = (C(I,J,6,5) - L61*C1 - L62*C2 - L63*C3
     1                               - L64*C4 - L65*C5)/L66
            C(I,J,7,5) = (C(I,J,7,5) - L71*C1 - L72*C2 - L73*C3
     1                               - L74*C4 - L75*C5 - L76*C6)/L77
            C(I,J,6,5) = C6-U67*C(I,J,7,5) 
            C(I,J,5,5) = C5-U57*C(I,J,7,5)-U56*C(I,J,6,5)
            C(I,J,4,5) = C4-U47*C(I,J,7,5)-U46*C(I,J,6,5)-U45*C(I,J,5,5)
            C(I,J,3,5) = C3-U37*C(I,J,7,5)-U36*C(I,J,6,5)-U35*C(I,J,5,5)
     1                     -U34*C(I,J,4,5)
            C(I,J,2,5) = C2-U27*C(I,J,7,5)-U26*C(I,J,6,5)-U25*C(I,J,5,5)
     1                     -U24*C(I,J,4,5)-U23*C(I,J,3,5)
            C(I,J,1,5) = C1-U17*C(I,J,7,5)-U16*C(I,J,6,5)-U15*C(I,J,5,5)
     1                     -U14*C(I,J,4,5)-U13*C(I,J,3,5)-U12*C(I,J,2,5)
C
            C1         =  C(I,J,1,6) / L11
            C2         = (C(I,J,2,6) - L21*C1) / L22
            C3         = (C(I,J,3,6) - L31*C1 - L32*C2) / L33
            C4         = (C(I,J,4,6) - L41*C1 - L42*C2 - L43*C3)/L44
            C5         = (C(I,J,5,6) - L51*C1 - L52*C2 - L53*C3
     1                               - L54*C4)/L55
            C6         = (C(I,J,6,6) - L61*C1 - L62*C2 - L63*C3
     1                               - L64*C4 - L65*C5)/L66
            C(I,J,7,6) = (C(I,J,7,6) - L71*C1 - L72*C2 - L73*C3
     1                               - L74*C4 - L75*C5 - L76*C6)/L77
            C(I,J,6,6) = C6-U67*C(I,J,7,6) 
            C(I,J,5,6) = C5-U57*C(I,J,7,6)-U56*C(I,J,6,6)
            C(I,J,4,6) = C4-U47*C(I,J,7,6)-U46*C(I,J,6,6)-U45*C(I,J,5,6)
            C(I,J,3,6) = C3-U37*C(I,J,7,6)-U36*C(I,J,6,6)-U35*C(I,J,5,6)
     1                     -U34*C(I,J,4,6)
            C(I,J,2,6) = C2-U27*C(I,J,7,6)-U26*C(I,J,6,6)-U25*C(I,J,5,6)
     1                     -U24*C(I,J,4,6)-U23*C(I,J,3,6)
            C(I,J,1,6) = C1-U17*C(I,J,7,6)-U16*C(I,J,6,6)-U15*C(I,J,5,6)
     1                     -U14*C(I,J,4,6)-U13*C(I,J,3,6)-U12*C(I,J,2,6)
C
            C1         =  C(I,J,1,7) / L11
            C2         = (C(I,J,2,7) - L21*C1) / L22
            C3         = (C(I,J,3,7) - L31*C1 - L32*C2) / L33
            C4         = (C(I,J,4,7) - L41*C1 - L42*C2 - L43*C3)/L44
            C5         = (C(I,J,5,7) - L51*C1 - L52*C2 - L53*C3
     1                               - L54*C4)/L55
            C6         = (C(I,J,6,7) - L61*C1 - L62*C2 - L63*C3
     1                               - L64*C4 - L65*C5)/L66
            C(I,J,7,7) = (C(I,J,7,7) - L71*C1 - L72*C2 - L73*C3
     1                               - L74*C4 - L75*C5 - L76*C6)/L77
            C(I,J,6,7) = C6-U67*C(I,J,7,7) 
            C(I,J,5,7) = C5-U57*C(I,J,7,7)-U56*C(I,J,6,7)
            C(I,J,4,7) = C4-U47*C(I,J,7,7)-U46*C(I,J,6,7)-U45*C(I,J,5,7)
            C(I,J,3,7) = C3-U37*C(I,J,7,7)-U36*C(I,J,6,7)-U35*C(I,J,5,7)
     1                     -U34*C(I,J,4,7)
            C(I,J,2,7) = C2-U27*C(I,J,7,7)-U26*C(I,J,6,7)-U25*C(I,J,5,7)
     1                     -U24*C(I,J,4,7)-U23*C(I,J,3,7)
            C(I,J,1,7) = C1-U17*C(I,J,7,7)-U16*C(I,J,6,7)-U15*C(I,J,5,7)
     1                     -U14*C(I,J,4,7)-U13*C(I,J,3,7)-U12*C(I,J,2,7)
  190    CONTINUE
  200 CONTINUE
C
C     J = NPTS Elimination
C     Do not do this loop if NPTS = 1
C
C     Eliminate A
C
C     B = B - A*C(j-1)
C
      J = NPTS
      IF (NPTS .GT. 1) THEN
         DO 260 L = 1, 7
            DO 240 M = 1, 7
               DO 230 I = 1, NSIM
                  B(I,J,L,M) = B(I,J,L,M)
     1                       - A(I,J,L,1)*C(I,J-1,1,M) 
     2                       - A(I,J,L,2)*C(I,J-1,2,M)
     3                       - A(I,J,L,3)*C(I,J-1,3,M) 
     4                       - A(I,J,L,4)*C(I,J-1,4,M)
     5                       - A(I,J,L,5)*C(I,J-1,5,M)
     6                       - A(I,J,L,6)*C(I,J-1,6,M)
     7                       - A(I,J,L,7)*C(I,J-1,7,M)
  230          CONTINUE
  240       CONTINUE
C
C     R = R - A*R(j-1)
C
            DO 250 I = 1, NSIM
               R(I,J,L) = R(I,J,L) 
     1                  - A(I,J,L,1)*R(I,J-1,1) 
     2                  - A(I,J,L,2)*R(I,J-1,2) 
     3                  - A(I,J,L,3)*R(I,J-1,3) 
     4                  - A(I,J,L,4)*R(I,J-1,4) 
     5                  - A(I,J,L,5)*R(I,J-1,5)
     6                  - A(I,J,L,6)*R(I,J-1,6)
     7                  - A(I,J,L,7)*R(I,J-1,7)
  250       CONTINUE
  260    CONTINUE
      ENDIF
C
C     Normalize row by multiplication of B(-1) to get unity diagonal
C
C     Form the LU of B(1)
C
      DO 300 I = 1, NSIM
         L11 =  B(I,J,1,1)
         U12 =  B(I,J,1,2) / L11
         U13 =  B(I,J,1,3) / L11
         U14 =  B(I,J,1,4) / L11
         U15 =  B(I,J,1,5) / L11
         U16 =  B(I,J,1,6) / L11
         U17 =  B(I,J,1,7) / L11
C
         L21 =  B(I,J,2,1) 
         L22 =  B(I,J,2,2) - L21*U12
         U23 = (B(I,J,2,3) - L21*U13) / L22
         U24 = (B(I,J,2,4) - L21*U14) / L22
         U25 = (B(I,J,2,5) - L21*U15) / L22
         U26 = (B(I,J,2,6) - L21*U16) / L22
         U27 = (B(I,J,2,7) - L21*U17) / L22
C
         L31 =  B(I,J,3,1)
         L32 =  B(I,J,3,2) - L31*U12
         L33 =  B(I,J,3,3) - L31*U13 - L32*U23
         U34 = (B(I,J,3,4) - L31*U14 - L32*U24) / L33
         U35 = (B(I,J,3,5) - L31*U15 - L32*U25) / L33
         U36 = (B(I,J,3,6) - L31*U16 - L32*U26) / L33
         U37 = (B(I,J,3,7) - L31*U17 - L32*U27) / L33
C
         L41 =  B(I,J,4,1)
         L42 =  B(I,J,4,2) - L41*U12
         L43 =  B(I,J,4,3) - L41*U13 - L42*U23
         L44 =  B(I,J,4,4) - L41*U14 - L42*U24 - L43*U34
         U45 = (B(I,J,4,5) - L41*U15 - L42*U25 - L43*U35) / L44
         U46 = (B(I,J,4,6) - L41*U16 - L42*U26 - L43*U36) / L44
         U47 = (B(I,J,4,7) - L41*U17 - L42*U27 - L43*U37) / L44
C
         L51 =  B(I,J,5,1)
         L52 =  B(I,J,5,2) - L51*U12
         L53 =  B(I,J,5,3) - L51*U13 - L52*U23
         L54 =  B(I,J,5,4) - L51*U14 - L52*U24 - L53*U34
         L55 =  B(I,J,5,5) - L51*U15 - L52*U25 - L53*U35 - L54*U45
         U56 = (B(I,J,5,6) - L51*U16 - L52*U26 - L53*U36 - L54*U46) 
     1         / L55
         U57 = (B(I,J,5,7) - L51*U17 - L52*U27 - L53*U37 - L54*U47) 
     1         / L55
C
         L61 =  B(I,J,6,1)
         L62 =  B(I,J,6,2) - L61*U12
         L63 =  B(I,J,6,3) - L61*U13 - L62*U23
         L64 =  B(I,J,6,4) - L61*U14 - L62*U24 - L63*U34
         L65 =  B(I,J,6,5) - L61*U15 - L62*U25 - L63*U35 - L64*U45
         L66 =  B(I,J,6,6) - L61*U16 - L62*U26 - L63*U36 - L64*U46
     1                     - L65*U56
         U67 = (B(I,J,6,7) - L61*U17 - L62*U27 - L63*U37 - L64*U47
     1                     - L65*U57) / L66
C
         L71 =  B(I,J,7,1)
         L72 =  B(I,J,7,2) - L71*U12
         L73 =  B(I,J,7,3) - L71*U13 - L72*U23
         L74 =  B(I,J,7,4) - L71*U14 - L72*U24 - L73*U34
         L75 =  B(I,J,7,5) - L71*U15 - L72*U25 - L73*U35 - L74*U45
         L76 =  B(I,J,7,6) - L71*U16 - L72*U26 - L73*U36 - L74*U46
     1                     - L75*U56
         L77 =  B(I,J,7,7) - L71*U17 - L72*U27 - L73*U37 - L74*U47
     1                     - L75*U57 - L76*U67
C
C     Calculate B(-1) * R
C
         R1       =  R(I,J,1) / L11
         R2       = (R(I,J,2) - L21*R1) / L22
         R3       = (R(I,J,3) - L31*R1 - L32*R2) / L33
         R4       = (R(I,J,4) - L41*R1 - L42*R2 - L43*R3) / L44
         R5       = (R(I,J,5) - L51*R1 - L52*R2 - L53*R3 - L54*R4)
     1              / L55
         R6       = (R(I,J,6) - L61*R1 - L62*R2 - L63*R3 - L64*R4 
     1                        - L65*R5) / L66
         R(I,J,7) = (R(I,J,7) - L71*R1 - L72*R2 - L73*R3 - L74*R4
     1                        - L75*R5 - L76*R6) / L77
         R(I,J,6) =  R6-U67*R(I,J,7)
         R(I,J,5) =  R5-U57*R(I,J,7)-U56*R(I,J,6)
         R(I,J,4) =  R4-U47*R(I,J,7)-U46*R(I,J,6)-U45*R(I,J,5)
         R(I,J,3) =  R3-U37*R(I,J,7)-U36*R(I,J,6)-U35*R(I,J,5)
     1                 -U34*R(I,J,4)
         R(I,J,2) =  R2-U27*R(I,J,7)-U26*R(I,J,6)-U25*R(I,J,5)
     1                 -U24*R(I,J,4)-U23*R(I,J,3)
         R(I,J,1) =  R1-U17*R(I,J,7)-U16*R(I,J,6)-U15*R(I,J,5)
     1                 -U14*R(I,J,4)-U13*R(I,J,3)-U12*R(I,J,2)
  300 CONTINUE
C
C     Backward Substitution Sweep
C
      DO 400 J = NPTS-1, 1, -1
         DO 390 L = 1, 7
            DO 380 I = 1, NSIM
               R(I,J,L) = R(I,J,L)
     1                  - C(I,J,L,1)*R(I,J+1,1) 
     2                  - C(I,J,L,2)*R(I,J+1,2)
     3                  - C(I,J,L,3)*R(I,J+1,3) 
     4                  - C(I,J,L,4)*R(I,J+1,4)
     5                  - C(I,J,L,5)*R(I,J+1,5)
     6                  - C(I,J,L,6)*R(I,J+1,6)
     7                  - C(I,J,L,7)*R(I,J+1,7)
  380       CONTINUE
  390    CONTINUE
  400 CONTINUE
C
C     Finished
C
      RETURN
      END






      SUBROUTINE VB12TR (NPTS, NSIM, NRANK, A, B, C, R)
C
C     Routine to invert a 12x12 block tridiagonal system
C     Inversion is by Upper Triangularization with a
C     unity diagonal followed by a backsubstitution step.
C     Vectorization is attained by inverting multiple systems
C     at one time with the recursive inversion algorithm.
C     The A, B, C arrays are destroyed in the solution
C     procedure.
C
C     NPTS:      The number of 
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     NRANK:     The actual dimension of the storage
C     A:         The supra-diagonal block 12x12
C     B:         The       diagonal block 12x12
C     C:         The super-diagonal block 12x12
C     R:         The right hand side
C
C     R contains the solution vector at return
C
C     The storage of the system is such that vectorization occurs over
C     the first index to provide a unit stride for maximum speed.
C
      DIMENSION A (NSIM,NPTS,NRANK,NRANK),
     1          B (NSIM,NPTS,NRANK,NRANK),
     2          C (NSIM,NPTS,NRANK,NRANK),
     3          R (NSIM,NPTS,NRANK)
C
C     Storage for the LU decomposition of a block 12x12
C
      REAL L11, 
     1     L21, L22, 
     2     L31, L32, L33, 
     3     L41, L42, L43, L44, 
     4     L51, L52, L53, L54, L55,
     5     L61, L62, L63, L64, L65, L66,
     6     L71, L72, L73, L74, L75, L76, L77,
     7     L81, L82, L83, L84, L85, L86, L87, L88,
     8     L91, L92, L93, L94, L95, L96, L97, L98, L99,
     9     LA1, LA2, LA3, LA4, LA5, LA6, LA7, LA8, LA9, LAA,
     A     LB1, LB2, LB3, LB4, LB5, LB6, LB7, LB8, LB9, LBA, LBB,
     B     LC1, LC2, LC3, LC4, LC5, LC6, LC7, LC8, LC9, LCA, LCB, LCC
C
C     Forward Substitution Loop - Upper Triangularization
C     Loop only to NPTS-1 to optimize vectorization while avoiding unnecessary
C     calculations of the L terms
C
      DO 200 J = 1, NPTS - 1
C
C     Eliminate A
C
C     B = B - A*C(j-1)
C
         IF (J .GT. 1) THEN
            DO 160 L = 1, 12
               DO 140 M = 1, 12
                  DO 130 I = 1, NSIM
                     B(I,J,L,M) = B(I,J,L,M)
     1                          - A(I,J,L, 1)*C(I,J-1, 1,M) 
     2                          - A(I,J,L, 2)*C(I,J-1, 2,M)
     3                          - A(I,J,L, 3)*C(I,J-1, 3,M) 
     4                          - A(I,J,L, 4)*C(I,J-1, 4,M)
     5                          - A(I,J,L, 5)*C(I,J-1, 5,M)
     6                          - A(I,J,L, 6)*C(I,J-1, 6,M)
     7                          - A(I,J,L, 7)*C(I,J-1, 7,M)
     8                          - A(I,J,L, 8)*C(I,J-1, 8,M)
     9                          - A(I,J,L, 9)*C(I,J-1, 9,M)
     A                          - A(I,J,L,10)*C(I,J-1,10,M)
     B                          - A(I,J,L,11)*C(I,J-1,11,M)
     C                          - A(I,J,L,12)*C(I,J-1,12,M)
  130          CONTINUE
  140       CONTINUE
C
C     R = R - A*R(j-1)
C
               DO 150 I = 1, NSIM
                  R(I,J,L) = R(I,J,L) 
     1                     - A(I,J,L, 1)*R(I,J-1, 1) 
     2                     - A(I,J,L, 2)*R(I,J-1, 2) 
     3                     - A(I,J,L, 3)*R(I,J-1, 3) 
     4                     - A(I,J,L, 4)*R(I,J-1, 4) 
     5                     - A(I,J,L, 5)*R(I,J-1, 5) 
     6                     - A(I,J,L, 6)*R(I,J-1, 6) 
     7                     - A(I,J,L, 7)*R(I,J-1, 7)
     8                     - A(I,J,L, 8)*R(I,J-1, 8)
     9                     - A(I,J,L, 9)*R(I,J-1, 9)
     A                     - A(I,J,L,10)*R(I,J-1,10)
     B                     - A(I,J,L,11)*R(I,J-1,11)
     C                     - A(I,J,L,12)*R(I,J-1,12)
  150          CONTINUE
  160       CONTINUE
         ENDIF
C
C     Normalize row by multiplication of B(-1) to get unity diagonal
C
C     Form the LU of B(1)
C
C     The IVDEP compiler directive is necessary to overcome a difficulty with
C     overflowing dependency checking tables in the Cray compiler 5.0.4.7
CDIR$ IVDEP
         DO 190 I = 1, NSIM
            L11 =  B(I,J,1, 1)
            U12 =  B(I,J,1, 2) / L11
            U13 =  B(I,J,1, 3) / L11
            U14 =  B(I,J,1, 4) / L11
            U15 =  B(I,J,1, 5) / L11
            U16 =  B(I,J,1, 6) / L11
            U17 =  B(I,J,1, 7) / L11
            U18 =  B(I,J,1, 8) / L11
            U19 =  B(I,J,1, 9) / L11
            U1A =  B(I,J,1,10) / L11
            U1B =  B(I,J,1,11) / L11
            U1C =  B(I,J,1,12) / L11
C
            L21 =  B(I,J,2, 1) 
            L22 =  B(I,J,2, 2) - L21*U12
            U23 = (B(I,J,2, 3) - L21*U13) / L22
            U24 = (B(I,J,2, 4) - L21*U14) / L22
            U25 = (B(I,J,2, 5) - L21*U15) / L22
            U26 = (B(I,J,2, 6) - L21*U16) / L22
            U27 = (B(I,J,2, 7) - L21*U17) / L22
            U28 = (B(I,J,2, 8) - L21*U18) / L22
            U29 = (B(I,J,2, 9) - L21*U19) / L22
            U2A = (B(I,J,2,10) - L21*U1A) / L22
            U2B = (B(I,J,2,11) - L21*U1B) / L22
            U2C = (B(I,J,2,12) - L21*U1C) / L22
C
            L31 =  B(I,J,3, 1)
            L32 =  B(I,J,3, 2) - L31*U12
            L33 =  B(I,J,3, 3) - L31*U13 - L32*U23
            U34 = (B(I,J,3, 4) - L31*U14 - L32*U24) / L33
            U35 = (B(I,J,3, 5) - L31*U15 - L32*U25) / L33
            U36 = (B(I,J,3, 6) - L31*U16 - L32*U26) / L33
            U37 = (B(I,J,3, 7) - L31*U17 - L32*U27) / L33
            U38 = (B(I,J,3, 8) - L31*U18 - L32*U28) / L33
            U39 = (B(I,J,3, 9) - L31*U19 - L32*U29) / L33
            U3A = (B(I,J,3,10) - L31*U1A - L32*U2A) / L33
            U3B = (B(I,J,3,11) - L31*U1B - L32*U2B) / L33
            U3C = (B(I,J,3,12) - L31*U1C - L32*U2C) / L33
C
            L41 =  B(I,J,4, 1)
            L42 =  B(I,J,4, 2) - L41*U12
            L43 =  B(I,J,4, 3) - L41*U13 - L42*U23
            L44 =  B(I,J,4, 4) - L41*U14 - L42*U24 - L43*U34
            U45 = (B(I,J,4, 5) - L41*U15 - L42*U25 - L43*U35) / L44
            U46 = (B(I,J,4, 6) - L41*U16 - L42*U26 - L43*U36) / L44
            U47 = (B(I,J,4, 7) - L41*U17 - L42*U27 - L43*U37) / L44
            U48 = (B(I,J,4, 8) - L41*U18 - L42*U28 - L43*U38) / L44
            U49 = (B(I,J,4, 9) - L41*U19 - L42*U29 - L43*U39) / L44
            U4A = (B(I,J,4,10) - L41*U1A - L42*U2A - L43*U3A) / L44
            U4B = (B(I,J,4,11) - L41*U1B - L42*U2B - L43*U3B) / L44
            U4C = (B(I,J,4,12) - L41*U1C - L42*U2C - L43*U3C) / L44
C
            L51 =  B(I,J,5, 1)
            L52 =  B(I,J,5, 2) - L51*U12
            L53 =  B(I,J,5, 3) - L51*U13 - L52*U23
            L54 =  B(I,J,5, 4) - L51*U14 - L52*U24 - L53*U34
            L55 =  B(I,J,5, 5) - L51*U15 - L52*U25 - L53*U35 - L54*U45
            U56 = (B(I,J,5, 6) - L51*U16 - L52*U26 - L53*U36 - L54*U46) 
     1            / L55
            U57 = (B(I,J,5, 7) - L51*U17 - L52*U27 - L53*U37 - L54*U47) 
     1            / L55
            U58 = (B(I,J,5, 8) - L51*U18 - L52*U28 - L53*U38 - L54*U48) 
     1            / L55
            U59 = (B(I,J,5, 9) - L51*U19 - L52*U29 - L53*U39 - L54*U49) 
     1            / L55
            U5A = (B(I,J,5,10) - L51*U1A - L52*U2A - L53*U3A - L54*U4A) 
     1            / L55
            U5B = (B(I,J,5,11) - L51*U1B - L52*U2B - L53*U3B - L54*U4B) 
     1            / L55
            U5C = (B(I,J,5,12) - L51*U1C - L52*U2C - L53*U3C - L54*U4C) 
     1            / L55
C
            L61 =  B(I,J,6, 1)
            L62 =  B(I,J,6, 2) - L61*U12
            L63 =  B(I,J,6, 3) - L61*U13 - L62*U23
            L64 =  B(I,J,6, 4) - L61*U14 - L62*U24 - L63*U34
            L65 =  B(I,J,6, 5) - L61*U15 - L62*U25 - L63*U35 - L64*U45
            L66 =  B(I,J,6, 6) - L61*U16 - L62*U26 - L63*U36 - L64*U46
     1                         - L65*U56
            U67 = (B(I,J,6, 7) - L61*U17 - L62*U27 - L63*U37 - L64*U47
     1                         - L65*U57) / L66
            U68 = (B(I,J,6, 8) - L61*U18 - L62*U28 - L63*U38 - L64*U48
     1                         - L65*U58) / L66
            U69 = (B(I,J,6, 9) - L61*U19 - L62*U29 - L63*U39 - L64*U49
     1                         - L65*U59) / L66
            U6A = (B(I,J,6,10) - L61*U1A - L62*U2A - L63*U3A - L64*U4A
     1                         - L65*U5A) / L66
            U6B = (B(I,J,6,11) - L61*U1B - L62*U2B - L63*U3B - L64*U4B
     1                         - L65*U5B) / L66
            U6C = (B(I,J,6,12) - L61*U1C - L62*U2C - L63*U3C - L64*U4C
     1                         - L65*U5C) / L66
C
            L71 =  B(I,J,7, 1)
            L72 =  B(I,J,7, 2) - L71*U12
            L73 =  B(I,J,7, 3) - L71*U13 - L72*U23
            L74 =  B(I,J,7, 4) - L71*U14 - L72*U24 - L73*U34
            L75 =  B(I,J,7, 5) - L71*U15 - L72*U25 - L73*U35 - L74*U45
            L76 =  B(I,J,7, 6) - L71*U16 - L72*U26 - L73*U36 - L74*U46
     1                         - L75*U56
            L77 =  B(I,J,7, 7) - L71*U17 - L72*U27 - L73*U37 - L74*U47
     1                         - L75*U57 - L76*U67
            U78 = (B(I,J,7, 8) - L71*U18 - L72*U28 - L73*U38 - L74*U48
     1                         - L75*U58 - L76*U68) / L77
            U79 = (B(I,J,7, 9) - L71*U19 - L72*U29 - L73*U39 - L74*U49
     1                         - L75*U59 - L76*U69) / L77
            U7A = (B(I,J,7,10) - L71*U1A - L72*U2A - L73*U3A - L74*U4A
     1                         - L75*U5A - L76*U6A) / L77
            U7B = (B(I,J,7,11) - L71*U1B - L72*U2B - L73*U3B - L74*U4B
     1                         - L75*U5B - L76*U6B) / L77
            U7C = (B(I,J,7,12) - L71*U1C - L72*U2C - L73*U3C - L74*U4C
     1                         - L75*U5C - L76*U6C) / L77
C
            L81 =  B(I,J,8, 1)
            L82 =  B(I,J,8, 2) - L81*U12
            L83 =  B(I,J,8, 3) - L81*U13 - L82*U23
            L84 =  B(I,J,8, 4) - L81*U14 - L82*U24 - L83*U34
            L85 =  B(I,J,8, 5) - L81*U15 - L82*U25 - L83*U35 - L84*U45
            L86 =  B(I,J,8, 6) - L81*U16 - L82*U26 - L83*U36 - L84*U46
     1                         - L85*U56
            L87 =  B(I,J,8, 7) - L81*U17 - L82*U27 - L83*U37 - L84*U47
     1                         - L85*U57 - L86*U67
            L88 =  B(I,J,8, 8) - L81*U18 - L82*U28 - L83*U38 - L84*U48
     1                         - L85*U58 - L86*U68 - L87*U78
            U89 = (B(I,J,8, 9) - L81*U19 - L82*U29 - L83*U39 - L84*U49
     1                         - L85*U59 - L86*U69 - L87*U79) / L88
            U8A = (B(I,J,8,10) - L81*U1A - L82*U2A - L83*U3A - L84*U4A
     1                         - L85*U5A - L86*U6A - L87*U7A) / L88
            U8B = (B(I,J,8,11) - L81*U1B - L82*U2B - L83*U3B - L84*U4B
     1                         - L85*U5B - L86*U6B - L87*U7B) / L88
            U8C = (B(I,J,8,12) - L81*U1C - L82*U2C - L83*U3C - L84*U4C
     1                         - L85*U5C - L86*U6C - L87*U7C) / L88
C
            L91 =  B(I,J,9, 1)
            L92 =  B(I,J,9, 2) - L91*U12
            L93 =  B(I,J,9, 3) - L91*U13 - L92*U23
            L94 =  B(I,J,9, 4) - L91*U14 - L92*U24 - L93*U34
            L95 =  B(I,J,9, 5) - L91*U15 - L92*U25 - L93*U35 - L94*U45
            L96 =  B(I,J,9, 6) - L91*U16 - L92*U26 - L93*U36 - L94*U46
     1                         - L95*U56
            L97 =  B(I,J,9, 7) - L91*U17 - L92*U27 - L93*U37 - L94*U47
     1                         - L95*U57 - L96*U67
            L98 =  B(I,J,9, 8) - L91*U18 - L92*U28 - L93*U38 - L94*U48
     1                         - L95*U58 - L96*U68 - L97*U78
            L99 =  B(I,J,9, 9) - L91*U19 - L92*U29 - L93*U39 - L94*U49
     1                         - L95*U59 - L96*U69 - L97*U79 - L98*U89
            U9A = (B(I,J,9,10) - L91*U1A - L92*U2A - L93*U3A - L94*U4A
     1                         - L95*U5A - L96*U6A - L97*U7A - L98*U8A) 
     2            / L99
            U9B = (B(I,J,9,11) - L91*U1B - L92*U2B - L93*U3B - L94*U4B
     1                         - L95*U5B - L96*U6B - L97*U7B - L98*U8B) 
     2            / L99
            U9C = (B(I,J,9,12) - L91*U1C - L92*U2C - L93*U3C - L94*U4C
     1                         - L95*U5C - L96*U6C - L97*U7C - L98*U8C) 
     2            / L99
C
            LA1 =  B(I,J,10, 1)
            LA2 =  B(I,J,10, 2) - LA1*U12
            LA3 =  B(I,J,10, 3) - LA1*U13 - LA2*U23
            LA4 =  B(I,J,10, 4) - LA1*U14 - LA2*U24 - LA3*U34
            LA5 =  B(I,J,10, 5) - LA1*U15 - LA2*U25 - LA3*U35 - LA4*U45
            LA6 =  B(I,J,10, 6) - LA1*U16 - LA2*U26 - LA3*U36 - LA4*U46
     1                          - LA5*U56
            LA7 =  B(I,J,10, 7) - LA1*U17 - LA2*U27 - LA3*U37 - LA4*U47
     1                          - LA5*U57 - LA6*U67
            LA8 =  B(I,J,10, 8) - LA1*U18 - LA2*U28 - LA3*U38 - LA4*U48
     1                          - LA5*U58 - LA6*U68 - LA7*U78
            LA9 =  B(I,J,10, 9) - LA1*U19 - LA2*U29 - LA3*U39 - LA4*U49
     1                          - LA5*U59 - LA6*U69 - LA7*U79 - LA8*U89
            LAA =  B(I,J,10,10) - LA1*U1A - LA2*U2A - LA3*U3A - LA4*U4A
     1                          - LA5*U5A - LA6*U6A - LA7*U7A - LA8*U8A 
     2                          - LA9*U9A
            UAB = (B(I,J,10,11) - LA1*U1B - LA2*U2B - LA3*U3B - LA4*U4B
     1                          - LA5*U5B - LA6*U6B - LA7*U7B - LA8*U8B
     2                          - LA9*U9B) / LAA
            UAC = (B(I,J,10,12) - LA1*U1C - LA2*U2C - LA3*U3C - LA4*U4C
     1                          - LA5*U5C - LA6*U6C - LA7*U7C - LA8*U8C
     2                          - LA9*U9C) / LAA
C
            LB1 =  B(I,J,11, 1)
            LB2 =  B(I,J,11, 2) - LB1*U12
            LB3 =  B(I,J,11, 3) - LB1*U13 - LB2*U23
            LB4 =  B(I,J,11, 4) - LB1*U14 - LB2*U24 - LB3*U34
            LB5 =  B(I,J,11, 5) - LB1*U15 - LB2*U25 - LB3*U35 - LB4*U45
            LB6 =  B(I,J,11, 6) - LB1*U16 - LB2*U26 - LB3*U36 - LB4*U46
     1                          - LB5*U56
            LB7 =  B(I,J,11, 7) - LB1*U17 - LB2*U27 - LB3*U37 - LB4*U47
     1                          - LB5*U57 - LB6*U67
            LB8 =  B(I,J,11, 8) - LB1*U18 - LB2*U28 - LB3*U38 - LB4*U48
     1                          - LB5*U58 - LB6*U68 - LB7*U78
            LB9 =  B(I,J,11, 9) - LB1*U19 - LB2*U29 - LB3*U39 - LB4*U49
     1                          - LB5*U59 - LB6*U69 - LB7*U79 - LB8*U89
            LBA =  B(I,J,11,10) - LB1*U1A - LB2*U2A - LB3*U3A - LB4*U4A
     1                          - LB5*U5A - LB6*U6A - LB7*U7A - LB8*U8A 
     2                          - LB9*U9A
            LBB =  B(I,J,11,11) - LB1*U1B - LB2*U2B - LB3*U3B - LB4*U4B
     1                          - LB5*U5B - LB6*U6B - LB7*U7B - LB8*U8B
     2                          - LB9*U9B - LBA*UAB
            UBC = (B(I,J,11,12) - LB1*U1C - LB2*U2C - LB3*U3C - LB4*U4C
     1                          - LB5*U5C - LB6*U6C - LB7*U7C - LB8*U8C
     2                          - LB9*U9C - LBA*UAC) / LAA
C
            LC1 =  B(I,J,12, 1)
            LC2 =  B(I,J,12, 2) - LC1*U12
            LC3 =  B(I,J,12, 3) - LC1*U13 - LC2*U23
            LC4 =  B(I,J,12, 4) - LC1*U14 - LC2*U24 - LC3*U34
            LC5 =  B(I,J,12, 5) - LC1*U15 - LC2*U25 - LC3*U35 - LC4*U45
            LC6 =  B(I,J,12, 6) - LC1*U16 - LC2*U26 - LC3*U36 - LC4*U46
     1                          - LC5*U56
            LC7 =  B(I,J,12, 7) - LC1*U17 - LC2*U27 - LC3*U37 - LC4*U47
     1                          - LC5*U57 - LC6*U67
            LC8 =  B(I,J,12, 8) - LC1*U18 - LC2*U28 - LC3*U38 - LC4*U48
     1                          - LC5*U58 - LC6*U68 - LC7*U78
            LC9 =  B(I,J,12, 9) - LC1*U19 - LC2*U29 - LC3*U39 - LC4*U49
     1                          - LC5*U59 - LC6*U69 - LC7*U79 - LC8*U89
            LCA =  B(I,J,12,10) - LC1*U1A - LC2*U2A - LC3*U3A - LC4*U4A
     1                          - LC5*U5A - LC6*U6A - LC7*U7A - LC8*U8A 
     2                          - LC9*U9A
            LCB =  B(I,J,12,11) - LC1*U1B - LC2*U2B - LC3*U3B - LC4*U4B
     1                          - LC5*U5B - LC6*U6B - LC7*U7B - LC8*U8B
     2                          - LC9*U9B - LCA*UAB
            LCC =  B(I,J,12,12) - LC1*U1C - LC2*U2C - LC3*U3C - LC4*U4C
     1                          - LC5*U5C - LC6*U6C - LC7*U7C - LC8*U8C
     2                          - LC9*U9C - LCA*UAC - LCB*UBC
C
C     Calculate B(-1) * R
C
            R1        =  R(I,J, 1) / L11
            R2        = (R(I,J, 2) - L21*R1) / L22
            R3        = (R(I,J, 3) - L31*R1 - L32*R2) / L33
            R4        = (R(I,J, 4) - L41*R1 - L42*R2 - L43*R3) / L44
            R5        = (R(I,J, 5) - L51*R1 - L52*R2 - L53*R3 - L54*R4)
     1                  / L55
            R6        = (R(I,J, 6) - L61*R1 - L62*R2 - L63*R3 - L64*R4 
     1                             - L65*R5) / L66
            R7        = (R(I,J, 7) - L71*R1 - L72*R2 - L73*R3 - L74*R4
     1                             - L75*R5 - L76*R6) / L77
            R8        = (R(I,J, 8) - L81*R1 - L82*R2 - L83*R3 - L84*R4
     1                             - L85*R5 - L86*R6 - L87*R7) / L88
            R9        = (R(I,J, 9) - L91*R1 - L92*R2 - L93*R3 - L94*R4
     1                             - L95*R5 - L96*R6 - L97*R7 - L98*R8) 
     2                  / L99
            RA        = (R(I,J,10) - LA1*R1 - LA2*R2 - LA3*R3 - LA4*R4
     1                             - LA5*R5 - LA6*R6 - LA7*R7 - LA8*R8 
     2                             - LA9*R9) / LAA
            RB        = (R(I,J,11) - LB1*R1 - LB2*R2 - LB3*R3 - LB4*R4
     1                             - LB5*R5 - LB6*R6 - LB7*R7 - LB8*R8 
     2                             - LB9*R9 - LBA*RA) / LBB
            R(I,J,12) = (R(I,J,12) - LC1*R1 - LC2*R2 - LC3*R3 - LC4*R4
     1                             - LC5*R5 - LC6*R6 - LC7*R7 - LC8*R8 
     2                             - LC9*R9 - LCA*RA - LCB*RB) / LCC
            R(I,J,11) = RB-UBC*R(I,J,12)
            R(I,J,10) = RA-UAC*R(I,J,12)-UAB*R(I,J,11)
            R(I,J, 9) = R9-U9C*R(I,J,12)-U9B*R(I,J,11)-U9A*R(I,J,10)
            R(I,J, 8) = R8-U8C*R(I,J,12)-U8B*R(I,J,11)-U8A*R(I,J,10)
     1                    -U89*R(I,J, 9)
            R(I,J, 7) = R7-U7C*R(I,J,12)-U7B*R(I,J,11)-U7A*R(I,J,10)
     1                    -U79*R(I,J, 9)-U78*R(I,J, 8)
            R(I,J, 6) = R6-U6C*R(I,J,12)-U6B*R(I,J,11)-U6A*R(I,J,10)
     1                    -U69*R(I,J, 9)-U68*R(I,J, 8)-U67*R(I,J, 7)
            R(I,J, 5) = R5-U5C*R(I,J,12)-U5B*R(I,J,11)-U5A*R(I,J,10)
     1                    -U59*R(I,J, 9)-U58*R(I,J, 8)-U57*R(I,J, 7)
     2                    -U56*R(I,J, 6)
            R(I,J, 4) = R4-U4C*R(I,J,12)-U4B*R(I,J,11)-U4A*R(I,J,10)
     1                    -U49*R(I,J, 9)-U48*R(I,J, 8)-U47*R(I,J, 7)
     2                    -U46*R(I,J, 6)-U45*R(I,J, 5)
            R(I,J, 3) = R3-U3C*R(I,J,12)-U3B*R(I,J,11)-U3A*R(I,J,10)
     1                    -U39*R(I,J, 9)-U38*R(I,J, 8)-U37*R(I,J, 7)
     2                    -U36*R(I,J, 6)-U35*R(I,J, 5)-U34*R(I,J, 4)
            R(I,J, 2) = R2-U2C*R(I,J,12)-U2B*R(I,J,11)-U2A*R(I,J,10)
     1                    -U29*R(I,J, 9)-U28*R(I,J, 8)-U27*R(I,J, 7)
     2                    -U26*R(I,J, 6)-U25*R(I,J, 5)-U24*R(I,J, 4)
     3                    -U23*R(I,J, 3)
            R(I,J, 1) = R1-U1C*R(I,J,12)-U1B*R(I,J,11)-U1A*R(I,J,10)
     1                    -U19*R(I,J, 9)-U18*R(I,J, 8)-U17*R(I,J, 7)
     2                    -U16*R(I,J, 6)-U15*R(I,J, 5)-U14*R(I,J, 4)
     3                    -U13*R(I,J, 3)-U12*R(I,J, 2)
C
C     Calculate B(-1) * C
C
C     Column 1
C
            C1           =  C(I,J, 1, 1) / L11
            C2           = (C(I,J, 2, 1)-L21*C1) / L22
            C3           = (C(I,J, 3, 1)-L31*C1-L32*C2) / L33
            C4           = (C(I,J, 4, 1)-L41*C1-L42*C2-L43*C3) / L44
            C5           = (C(I,J, 5, 1)-L51*C1-L52*C2-L53*C3-L54*C4)
     1                     / L55
            C6           = (C(I,J, 6, 1)-L61*C1-L62*C2-L63*C3-L64*C4 
     1                                  -L65*C5) / L66
            C7           = (C(I,J, 7, 1)-L71*C1-L72*C2-L73*C3-L74*C4
     1                                  -L75*C5-L76*C6) / L77
            C8           = (C(I,J, 8, 1)-L81*C1-L82*C2-L83*C3-L84*C4
     1                                  -L85*C5-L86*C6-L87*C7) / L88
            C9           = (C(I,J, 9, 1)-L91*C1-L92*C2-L93*C3-L94*C4
     1                                  -L95*C5-L96*C6-L97*C7-L98*C8) 
     2                     / L99
            CA           = (C(I,J,10, 1)-LA1*C1-LA2*C2-LA3*C3-LA4*C4
     1                                  -LA5*C5-LA6*C6-LA7*C7-LA8*C8 
     2                                  -LA9*C9) / LAA
            CB           = (C(I,J,11, 1)-LB1*C1-LB2*C2-LB3*C3-LB4*C4
     1                                  -LB5*C5-LB6*C6-LB7*C7-LB8*C8 
     2                                  -LB9*C9-LBA*CA) / LBB
            C(I,J,12, 1) = (C(I,J,12, 1)-LC1*C1-LC2*C2-LC3*C3-LC4*C4
     1                                  -LC5*C5-LC6*C6-LC7*C7-LC8*C8 
     2                                  -LC9*C9-LCA*CA-LCB*CB) / LCC
            C(I,J,11, 1) = CB-UBC*C(I,J,12, 1)
            C(I,J,10, 1) = CA-UAC*C(I,J,12, 1)-UAB*C(I,J,11, 1)
            C(I,J, 9, 1) = C9-U9C*C(I,J,12, 1)-U9B*C(I,J,11, 1)
     1                       -U9A*C(I,J,10, 1)
            C(I,J, 8, 1) = C8-U8C*C(I,J,12, 1)-U8B*C(I,J,11, 1)
     1                       -U8A*C(I,J,10, 1)-U89*C(I,J, 9, 1)
            C(I,J, 7, 1) = C7-U7C*C(I,J,12, 1)-U7B*C(I,J,11, 1)
     1                       -U7A*C(I,J,10, 1)-U79*C(I,J, 9, 1)
     2                       -U78*C(I,J, 8, 1)
            C(I,J, 6, 1) = C6-U6C*C(I,J,12, 1)-U6B*C(I,J,11, 1)
     1                       -U6A*C(I,J,10, 1)-U69*C(I,J, 9, 1)
     2                       -U68*C(I,J, 8, 1)-U67*C(I,J, 7, 1)
            C(I,J, 5, 1) = C5-U5C*C(I,J,12, 1)-U5B*C(I,J,11, 1)
     1                       -U5A*C(I,J,10, 1)-U59*C(I,J, 9, 1)
     2                       -U58*C(I,J, 8, 1)-U57*C(I,J, 7, 1)
     3                       -U56*C(I,J, 6, 1)
            C(I,J, 4, 1) = C4-U4C*C(I,J,12, 1)-U4B*C(I,J,11, 1)
     1                       -U4A*C(I,J,10, 1)-U49*C(I,J, 9, 1)
     2                       -U48*C(I,J, 8, 1)-U47*C(I,J, 7, 1)
     3                       -U46*C(I,J, 6, 1)-U45*C(I,J, 5, 1)
            C(I,J, 3, 1) = C3-U3C*C(I,J,12, 1)-U3B*C(I,J,11, 1)
     1                       -U3A*C(I,J,10, 1)-U39*C(I,J, 9, 1)
     2                       -U38*C(I,J, 8, 1)-U37*C(I,J, 7, 1)
     3                       -U36*C(I,J, 6, 1)-U35*C(I,J, 5, 1)
     4                       -U34*C(I,J, 4, 1)
            C(I,J, 2, 1) = C2-U2C*C(I,J,12, 1)-U2B*C(I,J,11, 1)
     1                       -U2A*C(I,J,10, 1)-U29*C(I,J, 9, 1)
     2                       -U28*C(I,J, 8, 1)-U27*C(I,J, 7, 1)
     3                       -U26*C(I,J, 6, 1)-U25*C(I,J, 5, 1)
     4                       -U24*C(I,J, 4, 1)-U23*C(I,J, 3, 1)
            C(I,J, 1, 1) = C1-U1C*C(I,J,12, 1)-U1B*C(I,J,11, 1)
     1                       -U1A*C(I,J,10, 1)-U19*C(I,J, 9, 1)
     2                       -U18*C(I,J, 8, 1)-U17*C(I,J, 7, 1)
     3                       -U16*C(I,J, 6, 1)-U15*C(I,J, 5, 1)
     4                       -U14*C(I,J, 4, 1)-U13*C(I,J, 3, 1)
     5                       -U12*C(I,J, 2, 1)
C
C     Column 2
C
            C1           =  C(I,J, 1, 2) / L11
            C2           = (C(I,J, 2, 2)-L21*C1) / L22
            C3           = (C(I,J, 3, 2)-L31*C1-L32*C2) / L33
            C4           = (C(I,J, 4, 2)-L41*C1-L42*C2-L43*C3) / L44
            C5           = (C(I,J, 5, 2)-L51*C1-L52*C2-L53*C3-L54*C4)
     1                     / L55
            C6           = (C(I,J, 6, 2)-L61*C1-L62*C2-L63*C3-L64*C4 
     1                                  -L65*C5) / L66
            C7           = (C(I,J, 7, 2)-L71*C1-L72*C2-L73*C3-L74*C4
     1                                  -L75*C5-L76*C6) / L77
            C8           = (C(I,J, 8, 2)-L81*C1-L82*C2-L83*C3-L84*C4
     1                                  -L85*C5-L86*C6-L87*C7) / L88
            C9           = (C(I,J, 9, 2)-L91*C1-L92*C2-L93*C3-L94*C4
     1                                  -L95*C5-L96*C6-L97*C7-L98*C8) 
     2                     / L99
            CA           = (C(I,J,10, 2)-LA1*C1-LA2*C2-LA3*C3-LA4*C4
     1                                  -LA5*C5-LA6*C6-LA7*C7-LA8*C8 
     2                                  -LA9*C9) / LAA
            CB           = (C(I,J,11, 2)-LB1*C1-LB2*C2-LB3*C3-LB4*C4
     1                                  -LB5*C5-LB6*C6-LB7*C7-LB8*C8 
     2                                  -LB9*C9-LBA*CA) / LBB
            C(I,J,12, 2) = (C(I,J,12, 2)-LC1*C1-LC2*C2-LC3*C3-LC4*C4
     1                                  -LC5*C5-LC6*C6-LC7*C7-LC8*C8 
     2                                  -LC9*C9-LCA*CA-LCB*CB) / LCC
            C(I,J,11, 2) = CB-UBC*C(I,J,12, 2)
            C(I,J,10, 2) = CA-UAC*C(I,J,12, 2)-UAB*C(I,J,11, 2)
            C(I,J, 9, 2) = C9-U9C*C(I,J,12, 2)-U9B*C(I,J,11, 2)
     1                       -U9A*C(I,J,10, 2)
            C(I,J, 8, 2) = C8-U8C*C(I,J,12, 2)-U8B*C(I,J,11, 2)
     1                       -U8A*C(I,J,10, 2)-U89*C(I,J, 9, 2)
            C(I,J, 7, 2) = C7-U7C*C(I,J,12, 2)-U7B*C(I,J,11, 2)
     1                       -U7A*C(I,J,10, 2)-U79*C(I,J, 9, 2)
     2                       -U78*C(I,J, 8, 2)
            C(I,J, 6, 2) = C6-U6C*C(I,J,12, 2)-U6B*C(I,J,11, 2)
     1                       -U6A*C(I,J,10, 2)-U69*C(I,J, 9, 2)
     2                       -U68*C(I,J, 8, 2)-U67*C(I,J, 7, 2)
            C(I,J, 5, 2) = C5-U5C*C(I,J,12, 2)-U5B*C(I,J,11, 2)
     1                       -U5A*C(I,J,10, 2)-U59*C(I,J, 9, 2)
     2                       -U58*C(I,J, 8, 2)-U57*C(I,J, 7, 2)
     3                       -U56*C(I,J, 6, 2)
            C(I,J, 4, 2) = C4-U4C*C(I,J,12, 2)-U4B*C(I,J,11, 2)
     1                       -U4A*C(I,J,10, 2)-U49*C(I,J, 9, 2)
     2                       -U48*C(I,J, 8, 2)-U47*C(I,J, 7, 2)
     3                       -U46*C(I,J, 6, 2)-U45*C(I,J, 5, 2)
            C(I,J, 3, 2) = C3-U3C*C(I,J,12, 2)-U3B*C(I,J,11, 2)
     1                       -U3A*C(I,J,10, 2)-U39*C(I,J, 9, 2)
     2                       -U38*C(I,J, 8, 2)-U37*C(I,J, 7, 2)
     3                       -U36*C(I,J, 6, 2)-U35*C(I,J, 5, 2)
     4                       -U34*C(I,J, 4, 2)
            C(I,J, 2, 2) = C2-U2C*C(I,J,12, 2)-U2B*C(I,J,11, 2)
     1                       -U2A*C(I,J,10, 2)-U29*C(I,J, 9, 2)
     2                       -U28*C(I,J, 8, 2)-U27*C(I,J, 7, 2)
     3                       -U26*C(I,J, 6, 2)-U25*C(I,J, 5, 2)
     4                       -U24*C(I,J, 4, 2)-U23*C(I,J, 3, 2)
            C(I,J, 1, 2) = C1-U1C*C(I,J,12, 2)-U1B*C(I,J,11, 2)
     1                       -U1A*C(I,J,10, 2)-U19*C(I,J, 9, 2)
     2                       -U18*C(I,J, 8, 2)-U17*C(I,J, 7, 2)
     3                       -U16*C(I,J, 6, 2)-U15*C(I,J, 5, 2)
     4                       -U14*C(I,J, 4, 2)-U13*C(I,J, 3, 2)
     5                       -U12*C(I,J, 2, 2)
C
C     Column 3
C
            C1           =  C(I,J, 1, 3) / L11
            C2           = (C(I,J, 2, 3)-L21*C1) / L22
            C3           = (C(I,J, 3, 3)-L31*C1-L32*C2) / L33
            C4           = (C(I,J, 4, 3)-L41*C1-L42*C2-L43*C3) / L44
            C5           = (C(I,J, 5, 3)-L51*C1-L52*C2-L53*C3-L54*C4)
     1                     / L55
            C6           = (C(I,J, 6, 3)-L61*C1-L62*C2-L63*C3-L64*C4 
     1                                  -L65*C5) / L66
            C7           = (C(I,J, 7, 3)-L71*C1-L72*C2-L73*C3-L74*C4
     1                                  -L75*C5-L76*C6) / L77
            C8           = (C(I,J, 8, 3)-L81*C1-L82*C2-L83*C3-L84*C4
     1                                  -L85*C5-L86*C6-L87*C7) / L88
            C9           = (C(I,J, 9, 3)-L91*C1-L92*C2-L93*C3-L94*C4
     1                                  -L95*C5-L96*C6-L97*C7-L98*C8) 
     2                     / L99
            CA           = (C(I,J,10, 3)-LA1*C1-LA2*C2-LA3*C3-LA4*C4
     1                                  -LA5*C5-LA6*C6-LA7*C7-LA8*C8 
     2                                  -LA9*C9) / LAA
            CB           = (C(I,J,11, 3)-LB1*C1-LB2*C2-LB3*C3-LB4*C4
     1                                  -LB5*C5-LB6*C6-LB7*C7-LB8*C8 
     2                                  -LB9*C9-LBA*CA) / LBB
            C(I,J,12, 3) = (C(I,J,12, 3)-LC1*C1-LC2*C2-LC3*C3-LC4*C4
     1                                  -LC5*C5-LC6*C6-LC7*C7-LC8*C8 
     2                                  -LC9*C9-LCA*CA-LCB*CB) / LCC
            C(I,J,11, 3) = CB-UBC*C(I,J,12, 3)
            C(I,J,10, 3) = CA-UAC*C(I,J,12, 3)-UAB*C(I,J,11, 3)
            C(I,J, 9, 3) = C9-U9C*C(I,J,12, 3)-U9B*C(I,J,11, 3)
     1                       -U9A*C(I,J,10, 3)
            C(I,J, 8, 3) = C8-U8C*C(I,J,12, 3)-U8B*C(I,J,11, 3)
     1                       -U8A*C(I,J,10, 3)-U89*C(I,J, 9, 3)
            C(I,J, 7, 3) = C7-U7C*C(I,J,12, 3)-U7B*C(I,J,11, 3)
     1                       -U7A*C(I,J,10, 3)-U79*C(I,J, 9, 3)
     2                       -U78*C(I,J, 8, 3)
            C(I,J, 6, 3) = C6-U6C*C(I,J,12, 3)-U6B*C(I,J,11, 3)
     1                       -U6A*C(I,J,10, 3)-U69*C(I,J, 9, 3)
     2                       -U68*C(I,J, 8, 3)-U67*C(I,J, 7, 3)
            C(I,J, 5, 3) = C5-U5C*C(I,J,12, 3)-U5B*C(I,J,11, 3)
     1                       -U5A*C(I,J,10, 3)-U59*C(I,J, 9, 3)
     2                       -U58*C(I,J, 8, 3)-U57*C(I,J, 7, 3)
     3                       -U56*C(I,J, 6, 3)
            C(I,J, 4, 3) = C4-U4C*C(I,J,12, 3)-U4B*C(I,J,11, 3)
     1                       -U4A*C(I,J,10, 3)-U49*C(I,J, 9, 3)
     2                       -U48*C(I,J, 8, 3)-U47*C(I,J, 7, 3)
     3                       -U46*C(I,J, 6, 3)-U45*C(I,J, 5, 3)
            C(I,J, 3, 3) = C3-U3C*C(I,J,12, 3)-U3B*C(I,J,11, 3)
     1                       -U3A*C(I,J,10, 3)-U39*C(I,J, 9, 3)
     2                       -U38*C(I,J, 8, 3)-U37*C(I,J, 7, 3)
     3                       -U36*C(I,J, 6, 3)-U35*C(I,J, 5, 3)
     4                       -U34*C(I,J, 4, 3)
            C(I,J, 2, 3) = C2-U2C*C(I,J,12, 3)-U2B*C(I,J,11, 3)
     1                       -U2A*C(I,J,10, 3)-U29*C(I,J, 9, 3)
     2                       -U28*C(I,J, 8, 3)-U27*C(I,J, 7, 3)
     3                       -U26*C(I,J, 6, 3)-U25*C(I,J, 5, 3)
     4                       -U24*C(I,J, 4, 3)-U23*C(I,J, 3, 3)
            C(I,J, 1, 3) = C1-U1C*C(I,J,12, 3)-U1B*C(I,J,11, 3)
     1                       -U1A*C(I,J,10, 3)-U19*C(I,J, 9, 3)
     2                       -U18*C(I,J, 8, 3)-U17*C(I,J, 7, 3)
     3                       -U16*C(I,J, 6, 3)-U15*C(I,J, 5, 3)
     4                       -U14*C(I,J, 4, 3)-U13*C(I,J, 3, 3)
     5                       -U12*C(I,J, 2, 3)
C
C     Column 4
C
            C1           =  C(I,J, 1, 4) / L11
            C2           = (C(I,J, 2, 4)-L21*C1) / L22
            C3           = (C(I,J, 3, 4)-L31*C1-L32*C2) / L33
            C4           = (C(I,J, 4, 4)-L41*C1-L42*C2-L43*C3) / L44
            C5           = (C(I,J, 5, 4)-L51*C1-L52*C2-L53*C3-L54*C4)
     1                     / L55
            C6           = (C(I,J, 6, 4)-L61*C1-L62*C2-L63*C3-L64*C4 
     1                                  -L65*C5) / L66
            C7           = (C(I,J, 7, 4)-L71*C1-L72*C2-L73*C3-L74*C4
     1                                  -L75*C5-L76*C6) / L77
            C8           = (C(I,J, 8, 4)-L81*C1-L82*C2-L83*C3-L84*C4
     1                                  -L85*C5-L86*C6-L87*C7) / L88
            C9           = (C(I,J, 9, 4)-L91*C1-L92*C2-L93*C3-L94*C4
     1                                  -L95*C5-L96*C6-L97*C7-L98*C8) 
     2                     / L99
            CA           = (C(I,J,10, 4)-LA1*C1-LA2*C2-LA3*C3-LA4*C4
     1                                  -LA5*C5-LA6*C6-LA7*C7-LA8*C8 
     2                                  -LA9*C9) / LAA
            CB           = (C(I,J,11, 4)-LB1*C1-LB2*C2-LB3*C3-LB4*C4
     1                                  -LB5*C5-LB6*C6-LB7*C7-LB8*C8 
     2                                  -LB9*C9-LBA*CA) / LBB
            C(I,J,12, 4) = (C(I,J,12, 4)-LC1*C1-LC2*C2-LC3*C3-LC4*C4
     1                                  -LC5*C5-LC6*C6-LC7*C7-LC8*C8 
     2                                  -LC9*C9-LCA*CA-LCB*CB) / LCC
            C(I,J,11, 4) = CB-UBC*C(I,J,12, 4)
            C(I,J,10, 4) = CA-UAC*C(I,J,12, 4)-UAB*C(I,J,11, 4)
            C(I,J, 9, 4) = C9-U9C*C(I,J,12, 4)-U9B*C(I,J,11, 4)
     1                       -U9A*C(I,J,10, 4)
            C(I,J, 8, 4) = C8-U8C*C(I,J,12, 4)-U8B*C(I,J,11, 4)
     1                       -U8A*C(I,J,10, 4)-U89*C(I,J, 9, 4)
            C(I,J, 7, 4) = C7-U7C*C(I,J,12, 4)-U7B*C(I,J,11, 4)
     1                       -U7A*C(I,J,10, 4)-U79*C(I,J, 9, 4)
     2                       -U78*C(I,J, 8, 4)
            C(I,J, 6, 4) = C6-U6C*C(I,J,12, 4)-U6B*C(I,J,11, 4)
     1                       -U6A*C(I,J,10, 4)-U69*C(I,J, 9, 4)
     2                       -U68*C(I,J, 8, 4)-U67*C(I,J, 7, 4)
            C(I,J, 5, 4) = C5-U5C*C(I,J,12, 4)-U5B*C(I,J,11, 4)
     1                       -U5A*C(I,J,10, 4)-U59*C(I,J, 9, 4)
     2                       -U58*C(I,J, 8, 4)-U57*C(I,J, 7, 4)
     3                       -U56*C(I,J, 6, 4)
            C(I,J, 4, 4) = C4-U4C*C(I,J,12, 4)-U4B*C(I,J,11, 4)
     1                       -U4A*C(I,J,10, 4)-U49*C(I,J, 9, 4)
     2                       -U48*C(I,J, 8, 4)-U47*C(I,J, 7, 4)
     3                       -U46*C(I,J, 6, 4)-U45*C(I,J, 5, 4)
            C(I,J, 3, 4) = C3-U3C*C(I,J,12, 4)-U3B*C(I,J,11, 4)
     1                       -U3A*C(I,J,10, 4)-U39*C(I,J, 9, 4)
     2                       -U38*C(I,J, 8, 4)-U37*C(I,J, 7, 4)
     3                       -U36*C(I,J, 6, 4)-U35*C(I,J, 5, 4)
     4                       -U34*C(I,J, 4, 4)
            C(I,J, 2, 4) = C2-U2C*C(I,J,12, 4)-U2B*C(I,J,11, 4)
     1                       -U2A*C(I,J,10, 4)-U29*C(I,J, 9, 4)
     2                       -U28*C(I,J, 8, 4)-U27*C(I,J, 7, 4)
     3                       -U26*C(I,J, 6, 4)-U25*C(I,J, 5, 4)
     4                       -U24*C(I,J, 4, 4)-U23*C(I,J, 3, 4)
            C(I,J, 1, 4) = C1-U1C*C(I,J,12, 4)-U1B*C(I,J,11, 4)
     1                       -U1A*C(I,J,10, 4)-U19*C(I,J, 9, 4)
     2                       -U18*C(I,J, 8, 4)-U17*C(I,J, 7, 4)
     3                       -U16*C(I,J, 6, 4)-U15*C(I,J, 5, 4)
     4                       -U14*C(I,J, 4, 4)-U13*C(I,J, 3, 4)
     5                       -U12*C(I,J, 2, 4)
C
C     Column 5
C
            C1           =  C(I,J, 1, 5) / L11
            C2           = (C(I,J, 2, 5)-L21*C1) / L22
            C3           = (C(I,J, 3, 5)-L31*C1-L32*C2) / L33
            C4           = (C(I,J, 4, 5)-L41*C1-L42*C2-L43*C3) / L44
            C5           = (C(I,J, 5, 5)-L51*C1-L52*C2-L53*C3-L54*C4)
     1                     / L55
            C6           = (C(I,J, 6, 5)-L61*C1-L62*C2-L63*C3-L64*C4 
     1                                  -L65*C5) / L66
            C7           = (C(I,J, 7, 5)-L71*C1-L72*C2-L73*C3-L74*C4
     1                                  -L75*C5-L76*C6) / L77
            C8           = (C(I,J, 8, 5)-L81*C1-L82*C2-L83*C3-L84*C4
     1                                  -L85*C5-L86*C6-L87*C7) / L88
            C9           = (C(I,J, 9, 5)-L91*C1-L92*C2-L93*C3-L94*C4
     1                                  -L95*C5-L96*C6-L97*C7-L98*C8) 
     2                     / L99
            CA           = (C(I,J,10, 5)-LA1*C1-LA2*C2-LA3*C3-LA4*C4
     1                                  -LA5*C5-LA6*C6-LA7*C7-LA8*C8 
     2                                  -LA9*C9) / LAA
            CB           = (C(I,J,11, 5)-LB1*C1-LB2*C2-LB3*C3-LB4*C4
     1                                  -LB5*C5-LB6*C6-LB7*C7-LB8*C8 
     2                                  -LB9*C9-LBA*CA) / LBB
            C(I,J,12, 5) = (C(I,J,12, 5)-LC1*C1-LC2*C2-LC3*C3-LC4*C4
     1                                  -LC5*C5-LC6*C6-LC7*C7-LC8*C8 
     2                                  -LC9*C9-LCA*CA-LCB*CB) / LCC
            C(I,J,11, 5) = CB-UBC*C(I,J,12, 5)
            C(I,J,10, 5) = CA-UAC*C(I,J,12, 5)-UAB*C(I,J,11, 5)
            C(I,J, 9, 5) = C9-U9C*C(I,J,12, 5)-U9B*C(I,J,11, 5)
     1                       -U9A*C(I,J,10, 5)
            C(I,J, 8, 5) = C8-U8C*C(I,J,12, 5)-U8B*C(I,J,11, 5)
     1                       -U8A*C(I,J,10, 5)-U89*C(I,J, 9, 5)
            C(I,J, 7, 5) = C7-U7C*C(I,J,12, 5)-U7B*C(I,J,11, 5)
     1                       -U7A*C(I,J,10, 5)-U79*C(I,J, 9, 5)
     2                       -U78*C(I,J, 8, 5)
            C(I,J, 6, 5) = C6-U6C*C(I,J,12, 5)-U6B*C(I,J,11, 5)
     1                       -U6A*C(I,J,10, 5)-U69*C(I,J, 9, 5)
     2                       -U68*C(I,J, 8, 5)-U67*C(I,J, 7, 5)
            C(I,J, 5, 5) = C5-U5C*C(I,J,12, 5)-U5B*C(I,J,11, 5)
     1                       -U5A*C(I,J,10, 5)-U59*C(I,J, 9, 5)
     2                       -U58*C(I,J, 8, 5)-U57*C(I,J, 7, 5)
     3                       -U56*C(I,J, 6, 5)
            C(I,J, 4, 5) = C4-U4C*C(I,J,12, 5)-U4B*C(I,J,11, 5)
     1                       -U4A*C(I,J,10, 5)-U49*C(I,J, 9, 5)
     2                       -U48*C(I,J, 8, 5)-U47*C(I,J, 7, 5)
     3                       -U46*C(I,J, 6, 5)-U45*C(I,J, 5, 5)
            C(I,J, 3, 5) = C3-U3C*C(I,J,12, 5)-U3B*C(I,J,11, 5)
     1                       -U3A*C(I,J,10, 5)-U39*C(I,J, 9, 5)
     2                       -U38*C(I,J, 8, 5)-U37*C(I,J, 7, 5)
     3                       -U36*C(I,J, 6, 5)-U35*C(I,J, 5, 5)
     4                       -U34*C(I,J, 4, 5)
            C(I,J, 2, 5) = C2-U2C*C(I,J,12, 5)-U2B*C(I,J,11, 5)
     1                       -U2A*C(I,J,10, 5)-U29*C(I,J, 9, 5)
     2                       -U28*C(I,J, 8, 5)-U27*C(I,J, 7, 5)
     3                       -U26*C(I,J, 6, 5)-U25*C(I,J, 5, 5)
     4                       -U24*C(I,J, 4, 5)-U23*C(I,J, 3, 5)
            C(I,J, 1, 5) = C1-U1C*C(I,J,12, 5)-U1B*C(I,J,11, 5)
     1                       -U1A*C(I,J,10, 5)-U19*C(I,J, 9, 5)
     2                       -U18*C(I,J, 8, 5)-U17*C(I,J, 7, 5)
     3                       -U16*C(I,J, 6, 5)-U15*C(I,J, 5, 5)
     4                       -U14*C(I,J, 4, 5)-U13*C(I,J, 3, 5)
     5                       -U12*C(I,J, 2, 5)
C
C     Column 6
C
            C1           =  C(I,J, 1, 6) / L11
            C2           = (C(I,J, 2, 6)-L21*C1) / L22
            C3           = (C(I,J, 3, 6)-L31*C1-L32*C2) / L33
            C4           = (C(I,J, 4, 6)-L41*C1-L42*C2-L43*C3) / L44
            C5           = (C(I,J, 5, 6)-L51*C1-L52*C2-L53*C3-L54*C4)
     1                     / L55
            C6           = (C(I,J, 6, 6)-L61*C1-L62*C2-L63*C3-L64*C4 
     1                                  -L65*C5) / L66
            C7           = (C(I,J, 7, 6)-L71*C1-L72*C2-L73*C3-L74*C4
     1                                  -L75*C5-L76*C6) / L77
            C8           = (C(I,J, 8, 6)-L81*C1-L82*C2-L83*C3-L84*C4
     1                                  -L85*C5-L86*C6-L87*C7) / L88
            C9           = (C(I,J, 9, 6)-L91*C1-L92*C2-L93*C3-L94*C4
     1                                  -L95*C5-L96*C6-L97*C7-L98*C8) 
     2                     / L99
            CA           = (C(I,J,10, 6)-LA1*C1-LA2*C2-LA3*C3-LA4*C4
     1                                  -LA5*C5-LA6*C6-LA7*C7-LA8*C8 
     2                                  -LA9*C9) / LAA
            CB           = (C(I,J,11, 6)-LB1*C1-LB2*C2-LB3*C3-LB4*C4
     1                                  -LB5*C5-LB6*C6-LB7*C7-LB8*C8 
     2                                  -LB9*C9-LBA*CA) / LBB
            C(I,J,12, 6) = (C(I,J,12, 6)-LC1*C1-LC2*C2-LC3*C3-LC4*C4
     1                                  -LC5*C5-LC6*C6-LC7*C7-LC8*C8 
     2                                  -LC9*C9-LCA*CA-LCB*CB) / LCC
            C(I,J,11, 6) = CB-UBC*C(I,J,12, 6)
            C(I,J,10, 6) = CA-UAC*C(I,J,12, 6)-UAB*C(I,J,11, 6)
            C(I,J, 9, 6) = C9-U9C*C(I,J,12, 6)-U9B*C(I,J,11, 6)
     1                       -U9A*C(I,J,10, 6)
            C(I,J, 8, 6) = C8-U8C*C(I,J,12, 6)-U8B*C(I,J,11, 6)
     1                       -U8A*C(I,J,10, 6)-U89*C(I,J, 9, 6)
            C(I,J, 7, 6) = C7-U7C*C(I,J,12, 6)-U7B*C(I,J,11, 6)
     1                       -U7A*C(I,J,10, 6)-U79*C(I,J, 9, 6)
     2                       -U78*C(I,J, 8, 6)
            C(I,J, 6, 6) = C6-U6C*C(I,J,12, 6)-U6B*C(I,J,11, 6)
     1                       -U6A*C(I,J,10, 6)-U69*C(I,J, 9, 6)
     2                       -U68*C(I,J, 8, 6)-U67*C(I,J, 7, 6)
            C(I,J, 5, 6) = C5-U5C*C(I,J,12, 6)-U5B*C(I,J,11, 6)
     1                       -U5A*C(I,J,10, 6)-U59*C(I,J, 9, 6)
     2                       -U58*C(I,J, 8, 6)-U57*C(I,J, 7, 6)
     3                       -U56*C(I,J, 6, 6)
            C(I,J, 4, 6) = C4-U4C*C(I,J,12, 6)-U4B*C(I,J,11, 6)
     1                       -U4A*C(I,J,10, 6)-U49*C(I,J, 9, 6)
     2                       -U48*C(I,J, 8, 6)-U47*C(I,J, 7, 6)
     3                       -U46*C(I,J, 6, 6)-U45*C(I,J, 5, 6)
            C(I,J, 3, 6) = C3-U3C*C(I,J,12, 6)-U3B*C(I,J,11, 6)
     1                       -U3A*C(I,J,10, 6)-U39*C(I,J, 9, 6)
     2                       -U38*C(I,J, 8, 6)-U37*C(I,J, 7, 6)
     3                       -U36*C(I,J, 6, 6)-U35*C(I,J, 5, 6)
     4                       -U34*C(I,J, 4, 6)
            C(I,J, 2, 6) = C2-U2C*C(I,J,12, 6)-U2B*C(I,J,11, 6)
     1                       -U2A*C(I,J,10, 6)-U29*C(I,J, 9, 6)
     2                       -U28*C(I,J, 8, 6)-U27*C(I,J, 7, 6)
     3                       -U26*C(I,J, 6, 6)-U25*C(I,J, 5, 6)
     4                       -U24*C(I,J, 4, 6)-U23*C(I,J, 3, 6)
            C(I,J, 1, 6) = C1-U1C*C(I,J,12, 6)-U1B*C(I,J,11, 6)
     1                       -U1A*C(I,J,10, 6)-U19*C(I,J, 9, 6)
     2                       -U18*C(I,J, 8, 6)-U17*C(I,J, 7, 6)
     3                       -U16*C(I,J, 6, 6)-U15*C(I,J, 5, 6)
     4                       -U14*C(I,J, 4, 6)-U13*C(I,J, 3, 6)
     5                       -U12*C(I,J, 2, 6)
C
C     Column 7
C
            C1           =  C(I,J, 1, 7) / L11
            C2           = (C(I,J, 2, 7)-L21*C1) / L22
            C3           = (C(I,J, 3, 7)-L31*C1-L32*C2) / L33
            C4           = (C(I,J, 4, 7)-L41*C1-L42*C2-L43*C3) / L44
            C5           = (C(I,J, 5, 7)-L51*C1-L52*C2-L53*C3-L54*C4)
     1                     / L55
            C6           = (C(I,J, 6, 7)-L61*C1-L62*C2-L63*C3-L64*C4 
     1                                  -L65*C5) / L66
            C7           = (C(I,J, 7, 7)-L71*C1-L72*C2-L73*C3-L74*C4
     1                                  -L75*C5-L76*C6) / L77
            C8           = (C(I,J, 8, 7)-L81*C1-L82*C2-L83*C3-L84*C4
     1                                  -L85*C5-L86*C6-L87*C7) / L88
            C9           = (C(I,J, 9, 7)-L91*C1-L92*C2-L93*C3-L94*C4
     1                                  -L95*C5-L96*C6-L97*C7-L98*C8) 
     2                     / L99
            CA           = (C(I,J,10, 7)-LA1*C1-LA2*C2-LA3*C3-LA4*C4
     1                                  -LA5*C5-LA6*C6-LA7*C7-LA8*C8 
     2                                  -LA9*C9) / LAA
            CB           = (C(I,J,11, 7)-LB1*C1-LB2*C2-LB3*C3-LB4*C4
     1                                  -LB5*C5-LB6*C6-LB7*C7-LB8*C8 
     2                                  -LB9*C9-LBA*CA) / LBB
            C(I,J,12, 7) = (C(I,J,12, 7)-LC1*C1-LC2*C2-LC3*C3-LC4*C4
     1                                  -LC5*C5-LC6*C6-LC7*C7-LC8*C8 
     2                                  -LC9*C9-LCA*CA-LCB*CB) / LCC
            C(I,J,11, 7) = CB-UBC*C(I,J,12, 7)
            C(I,J,10, 7) = CA-UAC*C(I,J,12, 7)-UAB*C(I,J,11, 7)
            C(I,J, 9, 7) = C9-U9C*C(I,J,12, 7)-U9B*C(I,J,11, 7)
     1                       -U9A*C(I,J,10, 7)
            C(I,J, 8, 7) = C8-U8C*C(I,J,12, 7)-U8B*C(I,J,11, 7)
     1                       -U8A*C(I,J,10, 7)-U89*C(I,J, 9, 7)
            C(I,J, 7, 7) = C7-U7C*C(I,J,12, 7)-U7B*C(I,J,11, 7)
     1                       -U7A*C(I,J,10, 7)-U79*C(I,J, 9, 7)
     2                       -U78*C(I,J, 8, 7)
            C(I,J, 6, 7) = C6-U6C*C(I,J,12, 7)-U6B*C(I,J,11, 7)
     1                       -U6A*C(I,J,10, 7)-U69*C(I,J, 9, 7)
     2                       -U68*C(I,J, 8, 7)-U67*C(I,J, 7, 7)
            C(I,J, 5, 7) = C5-U5C*C(I,J,12, 7)-U5B*C(I,J,11, 7)
     1                       -U5A*C(I,J,10, 7)-U59*C(I,J, 9, 7)
     2                       -U58*C(I,J, 8, 7)-U57*C(I,J, 7, 7)
     3                       -U56*C(I,J, 6, 7)
            C(I,J, 4, 7) = C4-U4C*C(I,J,12, 7)-U4B*C(I,J,11, 7)
     1                       -U4A*C(I,J,10, 7)-U49*C(I,J, 9, 7)
     2                       -U48*C(I,J, 8, 7)-U47*C(I,J, 7, 7)
     3                       -U46*C(I,J, 6, 7)-U45*C(I,J, 5, 7)
            C(I,J, 3, 7) = C3-U3C*C(I,J,12, 7)-U3B*C(I,J,11, 7)
     1                       -U3A*C(I,J,10, 7)-U39*C(I,J, 9, 7)
     2                       -U38*C(I,J, 8, 7)-U37*C(I,J, 7, 7)
     3                       -U36*C(I,J, 6, 7)-U35*C(I,J, 5, 7)
     4                       -U34*C(I,J, 4, 7)
            C(I,J, 2, 7) = C2-U2C*C(I,J,12, 7)-U2B*C(I,J,11, 7)
     1                       -U2A*C(I,J,10, 7)-U29*C(I,J, 9, 7)
     2                       -U28*C(I,J, 8, 7)-U27*C(I,J, 7, 7)
     3                       -U26*C(I,J, 6, 7)-U25*C(I,J, 5, 7)
     4                       -U24*C(I,J, 4, 7)-U23*C(I,J, 3, 7)
            C(I,J, 1, 7) = C1-U1C*C(I,J,12, 7)-U1B*C(I,J,11, 7)
     1                       -U1A*C(I,J,10, 7)-U19*C(I,J, 9, 7)
     2                       -U18*C(I,J, 8, 7)-U17*C(I,J, 7, 7)
     3                       -U16*C(I,J, 6, 7)-U15*C(I,J, 5, 7)
     4                       -U14*C(I,J, 4, 7)-U13*C(I,J, 3, 7)
     5                       -U12*C(I,J, 2, 7)
C
C     Column 8
C
            C1           =  C(I,J, 1, 8) / L11
            C2           = (C(I,J, 2, 8)-L21*C1) / L22
            C3           = (C(I,J, 3, 8)-L31*C1-L32*C2) / L33
            C4           = (C(I,J, 4, 8)-L41*C1-L42*C2-L43*C3) / L44
            C5           = (C(I,J, 5, 8)-L51*C1-L52*C2-L53*C3-L54*C4)
     1                     / L55
            C6           = (C(I,J, 6, 8)-L61*C1-L62*C2-L63*C3-L64*C4 
     1                                  -L65*C5) / L66
            C7           = (C(I,J, 7, 8)-L71*C1-L72*C2-L73*C3-L74*C4
     1                                  -L75*C5-L76*C6) / L77
            C8           = (C(I,J, 8, 8)-L81*C1-L82*C2-L83*C3-L84*C4
     1                                  -L85*C5-L86*C6-L87*C7) / L88
            C9           = (C(I,J, 9, 8)-L91*C1-L92*C2-L93*C3-L94*C4
     1                                  -L95*C5-L96*C6-L97*C7-L98*C8) 
     2                     / L99
            CA           = (C(I,J,10, 8)-LA1*C1-LA2*C2-LA3*C3-LA4*C4
     1                                  -LA5*C5-LA6*C6-LA7*C7-LA8*C8 
     2                                  -LA9*C9) / LAA
            CB           = (C(I,J,11, 8)-LB1*C1-LB2*C2-LB3*C3-LB4*C4
     1                                  -LB5*C5-LB6*C6-LB7*C7-LB8*C8 
     2                                  -LB9*C9-LBA*CA) / LBB
            C(I,J,12, 8) = (C(I,J,12, 8)-LC1*C1-LC2*C2-LC3*C3-LC4*C4
     1                                  -LC5*C5-LC6*C6-LC7*C7-LC8*C8 
     2                                  -LC9*C9-LCA*CA-LCB*CB) / LCC
            C(I,J,11, 8) = CB-UBC*C(I,J,12, 8)
            C(I,J,10, 8) = CA-UAC*C(I,J,12, 8)-UAB*C(I,J,11, 8)
            C(I,J, 9, 8) = C9-U9C*C(I,J,12, 8)-U9B*C(I,J,11, 8)
     1                       -U9A*C(I,J,10, 8)
            C(I,J, 8, 8) = C8-U8C*C(I,J,12, 8)-U8B*C(I,J,11, 8)
     1                       -U8A*C(I,J,10, 8)-U89*C(I,J, 9, 8)
            C(I,J, 7, 8) = C7-U7C*C(I,J,12, 8)-U7B*C(I,J,11, 8)
     1                       -U7A*C(I,J,10, 8)-U79*C(I,J, 9, 8)
     2                       -U78*C(I,J, 8, 8)
            C(I,J, 6, 8) = C6-U6C*C(I,J,12, 8)-U6B*C(I,J,11, 8)
     1                       -U6A*C(I,J,10, 8)-U69*C(I,J, 9, 8)
     2                       -U68*C(I,J, 8, 8)-U67*C(I,J, 7, 8)
            C(I,J, 5, 8) = C5-U5C*C(I,J,12, 8)-U5B*C(I,J,11, 8)
     1                       -U5A*C(I,J,10, 8)-U59*C(I,J, 9, 8)
     2                       -U58*C(I,J, 8, 8)-U57*C(I,J, 7, 8)
     3                       -U56*C(I,J, 6, 8)
            C(I,J, 4, 8) = C4-U4C*C(I,J,12, 8)-U4B*C(I,J,11, 8)
     1                       -U4A*C(I,J,10, 8)-U49*C(I,J, 9, 8)
     2                       -U48*C(I,J, 8, 8)-U47*C(I,J, 7, 8)
     3                       -U46*C(I,J, 6, 8)-U45*C(I,J, 5, 8)
            C(I,J, 3, 8) = C3-U3C*C(I,J,12, 8)-U3B*C(I,J,11, 8)
     1                       -U3A*C(I,J,10, 8)-U39*C(I,J, 9, 8)
     2                       -U38*C(I,J, 8, 8)-U37*C(I,J, 7, 8)
     3                       -U36*C(I,J, 6, 8)-U35*C(I,J, 5, 8)
     4                       -U34*C(I,J, 4, 8)
            C(I,J, 2, 8) = C2-U2C*C(I,J,12, 8)-U2B*C(I,J,11, 8)
     1                       -U2A*C(I,J,10, 8)-U29*C(I,J, 9, 8)
     2                       -U28*C(I,J, 8, 8)-U27*C(I,J, 7, 8)
     3                       -U26*C(I,J, 6, 8)-U25*C(I,J, 5, 8)
     4                       -U24*C(I,J, 4, 8)-U23*C(I,J, 3, 8)
            C(I,J, 1, 8) = C1-U1C*C(I,J,12, 8)-U1B*C(I,J,11, 8)
     1                       -U1A*C(I,J,10, 8)-U19*C(I,J, 9, 8)
     2                       -U18*C(I,J, 8, 8)-U17*C(I,J, 7, 8)
     3                       -U16*C(I,J, 6, 8)-U15*C(I,J, 5, 8)
     4                       -U14*C(I,J, 4, 8)-U13*C(I,J, 3, 8)
     5                       -U12*C(I,J, 2, 8)
C
C     Column 9
C
            C1           =  C(I,J, 1, 9) / L11
            C2           = (C(I,J, 2, 9)-L21*C1) / L22
            C3           = (C(I,J, 3, 9)-L31*C1-L32*C2) / L33
            C4           = (C(I,J, 4, 9)-L41*C1-L42*C2-L43*C3) / L44
            C5           = (C(I,J, 5, 9)-L51*C1-L52*C2-L53*C3-L54*C4)
     1                     / L55
            C6           = (C(I,J, 6, 9)-L61*C1-L62*C2-L63*C3-L64*C4 
     1                                  -L65*C5) / L66
            C7           = (C(I,J, 7, 9)-L71*C1-L72*C2-L73*C3-L74*C4
     1                                  -L75*C5-L76*C6) / L77
            C8           = (C(I,J, 8, 9)-L81*C1-L82*C2-L83*C3-L84*C4
     1                                  -L85*C5-L86*C6-L87*C7) / L88
            C9           = (C(I,J, 9, 9)-L91*C1-L92*C2-L93*C3-L94*C4
     1                                  -L95*C5-L96*C6-L97*C7-L98*C8) 
     2                     / L99
            CA           = (C(I,J,10, 9)-LA1*C1-LA2*C2-LA3*C3-LA4*C4
     1                                  -LA5*C5-LA6*C6-LA7*C7-LA8*C8 
     2                                  -LA9*C9) / LAA
            CB           = (C(I,J,11, 9)-LB1*C1-LB2*C2-LB3*C3-LB4*C4
     1                                  -LB5*C5-LB6*C6-LB7*C7-LB8*C8 
     2                                  -LB9*C9-LBA*CA) / LBB
            C(I,J,12, 9) = (C(I,J,12, 9)-LC1*C1-LC2*C2-LC3*C3-LC4*C4
     1                                  -LC5*C5-LC6*C6-LC7*C7-LC8*C8 
     2                                  -LC9*C9-LCA*CA-LCB*CB) / LCC
            C(I,J,11, 9) = CB-UBC*C(I,J,12, 9)
            C(I,J,10, 9) = CA-UAC*C(I,J,12, 9)-UAB*C(I,J,11, 9)
            C(I,J, 9, 9) = C9-U9C*C(I,J,12, 9)-U9B*C(I,J,11, 9)
     1                       -U9A*C(I,J,10, 9)
            C(I,J, 8, 9) = C8-U8C*C(I,J,12, 9)-U8B*C(I,J,11, 9)
     1                       -U8A*C(I,J,10, 9)-U89*C(I,J, 9, 9)
            C(I,J, 7, 9) = C7-U7C*C(I,J,12, 9)-U7B*C(I,J,11, 9)
     1                       -U7A*C(I,J,10, 9)-U79*C(I,J, 9, 9)
     2                       -U78*C(I,J, 8, 9)
            C(I,J, 6, 9) = C6-U6C*C(I,J,12, 9)-U6B*C(I,J,11, 9)
     1                       -U6A*C(I,J,10, 9)-U69*C(I,J, 9, 9)
     2                       -U68*C(I,J, 8, 9)-U67*C(I,J, 7, 9)
            C(I,J, 5, 9) = C5-U5C*C(I,J,12, 9)-U5B*C(I,J,11, 9)
     1                       -U5A*C(I,J,10, 9)-U59*C(I,J, 9, 9)
     2                       -U58*C(I,J, 8, 9)-U57*C(I,J, 7, 9)
     3                       -U56*C(I,J, 6, 9)
            C(I,J, 4, 9) = C4-U4C*C(I,J,12, 9)-U4B*C(I,J,11, 9)
     1                       -U4A*C(I,J,10, 9)-U49*C(I,J, 9, 9)
     2                       -U48*C(I,J, 8, 9)-U47*C(I,J, 7, 9)
     3                       -U46*C(I,J, 6, 9)-U45*C(I,J, 5, 9)
            C(I,J, 3, 9) = C3-U3C*C(I,J,12, 9)-U3B*C(I,J,11, 9)
     1                       -U3A*C(I,J,10, 9)-U39*C(I,J, 9, 9)
     2                       -U38*C(I,J, 8, 9)-U37*C(I,J, 7, 9)
     3                       -U36*C(I,J, 6, 9)-U35*C(I,J, 5, 9)
     4                       -U34*C(I,J, 4, 9)
            C(I,J, 2, 9) = C2-U2C*C(I,J,12, 9)-U2B*C(I,J,11, 9)
     1                       -U2A*C(I,J,10, 9)-U29*C(I,J, 9, 9)
     2                       -U28*C(I,J, 8, 9)-U27*C(I,J, 7, 9)
     3                       -U26*C(I,J, 6, 9)-U25*C(I,J, 5, 9)
     4                       -U24*C(I,J, 4, 9)-U23*C(I,J, 3, 9)
            C(I,J, 1, 9) = C1-U1C*C(I,J,12, 9)-U1B*C(I,J,11, 9)
     1                       -U1A*C(I,J,10, 9)-U19*C(I,J, 9, 9)
     2                       -U18*C(I,J, 8, 9)-U17*C(I,J, 7, 9)
     3                       -U16*C(I,J, 6, 9)-U15*C(I,J, 5, 9)
     4                       -U14*C(I,J, 4, 9)-U13*C(I,J, 3, 9)
     5                       -U12*C(I,J, 2, 9)
C
C     Column 10
C
            C1           =  C(I,J, 1,10) / L11
            C2           = (C(I,J, 2,10)-L21*C1) / L22
            C3           = (C(I,J, 3,10)-L31*C1-L32*C2) / L33
            C4           = (C(I,J, 4,10)-L41*C1-L42*C2-L43*C3) / L44
            C5           = (C(I,J, 5,10)-L51*C1-L52*C2-L53*C3-L54*C4)
     1                     / L55
            C6           = (C(I,J, 6,10)-L61*C1-L62*C2-L63*C3-L64*C4 
     1                                  -L65*C5) / L66
            C7           = (C(I,J, 7,10)-L71*C1-L72*C2-L73*C3-L74*C4
     1                                  -L75*C5-L76*C6) / L77
            C8           = (C(I,J, 8,10)-L81*C1-L82*C2-L83*C3-L84*C4
     1                                  -L85*C5-L86*C6-L87*C7) / L88
            C9           = (C(I,J, 9,10)-L91*C1-L92*C2-L93*C3-L94*C4
     1                                  -L95*C5-L96*C6-L97*C7-L98*C8) 
     2                     / L99
            CA           = (C(I,J,10,10)-LA1*C1-LA2*C2-LA3*C3-LA4*C4
     1                                  -LA5*C5-LA6*C6-LA7*C7-LA8*C8 
     2                                  -LA9*C9) / LAA
            CB           = (C(I,J,11,10)-LB1*C1-LB2*C2-LB3*C3-LB4*C4
     1                                  -LB5*C5-LB6*C6-LB7*C7-LB8*C8 
     2                                  -LB9*C9-LBA*CA) / LBB
            C(I,J,12,10) = (C(I,J,12,10)-LC1*C1-LC2*C2-LC3*C3-LC4*C4
     1                                  -LC5*C5-LC6*C6-LC7*C7-LC8*C8 
     2                                  -LC9*C9-LCA*CA-LCB*CB) / LCC
            C(I,J,11,10) = CB-UBC*C(I,J,12,10)
            C(I,J,10,10) = CA-UAC*C(I,J,12,10)-UAB*C(I,J,11,10)
            C(I,J, 9,10) = C9-U9C*C(I,J,12,10)-U9B*C(I,J,11,10)
     1                       -U9A*C(I,J,10,10)
            C(I,J, 8,10) = C8-U8C*C(I,J,12,10)-U8B*C(I,J,11,10)
     1                       -U8A*C(I,J,10,10)-U89*C(I,J, 9,10)
            C(I,J, 7,10) = C7-U7C*C(I,J,12,10)-U7B*C(I,J,11,10)
     1                       -U7A*C(I,J,10,10)-U79*C(I,J, 9,10)
     2                       -U78*C(I,J, 8,10)
            C(I,J, 6,10) = C6-U6C*C(I,J,12,10)-U6B*C(I,J,11,10)
     1                       -U6A*C(I,J,10,10)-U69*C(I,J, 9,10)
     2                       -U68*C(I,J, 8,10)-U67*C(I,J, 7,10)
            C(I,J, 5,10) = C5-U5C*C(I,J,12,10)-U5B*C(I,J,11,10)
     1                       -U5A*C(I,J,10,10)-U59*C(I,J, 9,10)
     2                       -U58*C(I,J, 8,10)-U57*C(I,J, 7,10)
     3                       -U56*C(I,J, 6,10)
            C(I,J, 4,10) = C4-U4C*C(I,J,12,10)-U4B*C(I,J,11,10)
     1                       -U4A*C(I,J,10,10)-U49*C(I,J, 9,10)
     2                       -U48*C(I,J, 8,10)-U47*C(I,J, 7,10)
     3                       -U46*C(I,J, 6,10)-U45*C(I,J, 5,10)
            C(I,J, 3,10) = C3-U3C*C(I,J,12,10)-U3B*C(I,J,11,10)
     1                       -U3A*C(I,J,10,10)-U39*C(I,J, 9,10)
     2                       -U38*C(I,J, 8,10)-U37*C(I,J, 7,10)
     3                       -U36*C(I,J, 6,10)-U35*C(I,J, 5,10)
     4                       -U34*C(I,J, 4,10)
            C(I,J, 2,10) = C2-U2C*C(I,J,12,10)-U2B*C(I,J,11,10)
     1                       -U2A*C(I,J,10,10)-U29*C(I,J, 9,10)
     2                       -U28*C(I,J, 8,10)-U27*C(I,J, 7,10)
     3                       -U26*C(I,J, 6,10)-U25*C(I,J, 5,10)
     4                       -U24*C(I,J, 4,10)-U23*C(I,J, 3,10)
            C(I,J, 1,10) = C1-U1C*C(I,J,12,10)-U1B*C(I,J,11,10)
     1                       -U1A*C(I,J,10,10)-U19*C(I,J, 9,10)
     2                       -U18*C(I,J, 8,10)-U17*C(I,J, 7,10)
     3                       -U16*C(I,J, 6,10)-U15*C(I,J, 5,10)
     4                       -U14*C(I,J, 4,10)-U13*C(I,J, 3,10)
     5                       -U12*C(I,J, 2,10)
C
C     Column 11
C
            C1           =  C(I,J, 1,11) / L11
            C2           = (C(I,J, 2,11)-L21*C1) / L22
            C3           = (C(I,J, 3,11)-L31*C1-L32*C2) / L33
            C4           = (C(I,J, 4,11)-L41*C1-L42*C2-L43*C3) / L44
            C5           = (C(I,J, 5,11)-L51*C1-L52*C2-L53*C3-L54*C4)
     1                     / L55
            C6           = (C(I,J, 6,11)-L61*C1-L62*C2-L63*C3-L64*C4 
     1                                  -L65*C5) / L66
            C7           = (C(I,J, 7,11)-L71*C1-L72*C2-L73*C3-L74*C4
     1                                  -L75*C5-L76*C6) / L77
            C8           = (C(I,J, 8,11)-L81*C1-L82*C2-L83*C3-L84*C4
     1                                  -L85*C5-L86*C6-L87*C7) / L88
            C9           = (C(I,J, 9,11)-L91*C1-L92*C2-L93*C3-L94*C4
     1                                  -L95*C5-L96*C6-L97*C7-L98*C8) 
     2                     / L99
            CA           = (C(I,J,10,11)-LA1*C1-LA2*C2-LA3*C3-LA4*C4
     1                                  -LA5*C5-LA6*C6-LA7*C7-LA8*C8 
     2                                  -LA9*C9) / LAA
            CB           = (C(I,J,11,11)-LB1*C1-LB2*C2-LB3*C3-LB4*C4
     1                                  -LB5*C5-LB6*C6-LB7*C7-LB8*C8 
     2                                  -LB9*C9-LBA*CA) / LBB
            C(I,J,12,11) = (C(I,J,12,11)-LC1*C1-LC2*C2-LC3*C3-LC4*C4
     1                                  -LC5*C5-LC6*C6-LC7*C7-LC8*C8 
     2                                  -LC9*C9-LCA*CA-LCB*CB) / LCC
            C(I,J,11,11) = CB-UBC*C(I,J,12,11)
            C(I,J,10,11) = CA-UAC*C(I,J,12,11)-UAB*C(I,J,11,11)
            C(I,J, 9,11) = C9-U9C*C(I,J,12,11)-U9B*C(I,J,11,11)
     1                       -U9A*C(I,J,10,11)
            C(I,J, 8,11) = C8-U8C*C(I,J,12,11)-U8B*C(I,J,11,11)
     1                       -U8A*C(I,J,10,11)-U89*C(I,J, 9,11)
            C(I,J, 7,11) = C7-U7C*C(I,J,12,11)-U7B*C(I,J,11,11)
     1                       -U7A*C(I,J,10,11)-U79*C(I,J, 9,11)
     2                       -U78*C(I,J, 8,11)
            C(I,J, 6,11) = C6-U6C*C(I,J,12,11)-U6B*C(I,J,11,11)
     1                       -U6A*C(I,J,10,11)-U69*C(I,J, 9,11)
     2                       -U68*C(I,J, 8,11)-U67*C(I,J, 7,11)
            C(I,J, 5,11) = C5-U5C*C(I,J,12,11)-U5B*C(I,J,11,11)
     1                       -U5A*C(I,J,10,11)-U59*C(I,J, 9,11)
     2                       -U58*C(I,J, 8,11)-U57*C(I,J, 7,11)
     3                       -U56*C(I,J, 6,11)
            C(I,J, 4,11) = C4-U4C*C(I,J,12,11)-U4B*C(I,J,11,11)
     1                       -U4A*C(I,J,10,11)-U49*C(I,J, 9,11)
     2                       -U48*C(I,J, 8,11)-U47*C(I,J, 7,11)
     3                       -U46*C(I,J, 6,11)-U45*C(I,J, 5,11)
            C(I,J, 3,11) = C3-U3C*C(I,J,12,11)-U3B*C(I,J,11,11)
     1                       -U3A*C(I,J,10,11)-U39*C(I,J, 9,11)
     2                       -U38*C(I,J, 8,11)-U37*C(I,J, 7,11)
     3                       -U36*C(I,J, 6,11)-U35*C(I,J, 5,11)
     4                       -U34*C(I,J, 4,11)
            C(I,J, 2,11) = C2-U2C*C(I,J,12,11)-U2B*C(I,J,11,11)
     1                       -U2A*C(I,J,10,11)-U29*C(I,J, 9,11)
     2                       -U28*C(I,J, 8,11)-U27*C(I,J, 7,11)
     3                       -U26*C(I,J, 6,11)-U25*C(I,J, 5,11)
     4                       -U24*C(I,J, 4,11)-U23*C(I,J, 3,11)
            C(I,J, 1,11) = C1-U1C*C(I,J,12,11)-U1B*C(I,J,11,11)
     1                       -U1A*C(I,J,10,11)-U19*C(I,J, 9,11)
     2                       -U18*C(I,J, 8,11)-U17*C(I,J, 7,11)
     3                       -U16*C(I,J, 6,11)-U15*C(I,J, 5,11)
     4                       -U14*C(I,J, 4,11)-U13*C(I,J, 3,11)
     5                       -U12*C(I,J, 2,11)
C
C     Column 12
C
            C1           =  C(I,J, 1,12) / L11
            C2           = (C(I,J, 2,12)-L21*C1) / L22
            C3           = (C(I,J, 3,12)-L31*C1-L32*C2) / L33
            C4           = (C(I,J, 4,12)-L41*C1-L42*C2-L43*C3) / L44
            C5           = (C(I,J, 5,12)-L51*C1-L52*C2-L53*C3-L54*C4)
     1                     / L55
            C6           = (C(I,J, 6,12)-L61*C1-L62*C2-L63*C3-L64*C4 
     1                                  -L65*C5) / L66
            C7           = (C(I,J, 7,12)-L71*C1-L72*C2-L73*C3-L74*C4
     1                                  -L75*C5-L76*C6) / L77
            C8           = (C(I,J, 8,12)-L81*C1-L82*C2-L83*C3-L84*C4
     1                                  -L85*C5-L86*C6-L87*C7) / L88
            C9           = (C(I,J, 9,12)-L91*C1-L92*C2-L93*C3-L94*C4
     1                                  -L95*C5-L96*C6-L97*C7-L98*C8) 
     2                     / L99
            CA           = (C(I,J,10,12)-LA1*C1-LA2*C2-LA3*C3-LA4*C4
     1                                  -LA5*C5-LA6*C6-LA7*C7-LA8*C8 
     2                                  -LA9*C9) / LAA
            CB           = (C(I,J,11,12)-LB1*C1-LB2*C2-LB3*C3-LB4*C4
     1                                  -LB5*C5-LB6*C6-LB7*C7-LB8*C8 
     2                                  -LB9*C9-LBA*CA) / LBB
            C(I,J,12,12) = (C(I,J,12,12)-LC1*C1-LC2*C2-LC3*C3-LC4*C4
     1                                  -LC5*C5-LC6*C6-LC7*C7-LC8*C8 
     2                                  -LC9*C9-LCA*CA-LCB*CB) / LCC
            C(I,J,11,12) = CB-UBC*C(I,J,12,12)
            C(I,J,10,12) = CA-UAC*C(I,J,12,12)-UAB*C(I,J,11,12)
            C(I,J, 9,12) = C9-U9C*C(I,J,12,12)-U9B*C(I,J,11,12)
     1                       -U9A*C(I,J,10,12)
            C(I,J, 8,12) = C8-U8C*C(I,J,12,12)-U8B*C(I,J,11,12)
     1                       -U8A*C(I,J,10,12)-U89*C(I,J, 9,12)
            C(I,J, 7,12) = C7-U7C*C(I,J,12,12)-U7B*C(I,J,11,12)
     1                       -U7A*C(I,J,10,12)-U79*C(I,J, 9,12)
     2                       -U78*C(I,J, 8,12)
            C(I,J, 6,12) = C6-U6C*C(I,J,12,12)-U6B*C(I,J,11,12)
     1                       -U6A*C(I,J,10,12)-U69*C(I,J, 9,12)
     2                       -U68*C(I,J, 8,12)-U67*C(I,J, 7,12)
            C(I,J, 5,12) = C5-U5C*C(I,J,12,12)-U5B*C(I,J,11,12)
     1                       -U5A*C(I,J,10,12)-U59*C(I,J, 9,12)
     2                       -U58*C(I,J, 8,12)-U57*C(I,J, 7,12)
     3                       -U56*C(I,J, 6,12)
            C(I,J, 4,12) = C4-U4C*C(I,J,12,12)-U4B*C(I,J,11,12)
     1                       -U4A*C(I,J,10,12)-U49*C(I,J, 9,12)
     2                       -U48*C(I,J, 8,12)-U47*C(I,J, 7,12)
     3                       -U46*C(I,J, 6,12)-U45*C(I,J, 5,12)
            C(I,J, 3,12) = C3-U3C*C(I,J,12,12)-U3B*C(I,J,11,12)
     1                       -U3A*C(I,J,10,12)-U39*C(I,J, 9,12)
     2                       -U38*C(I,J, 8,12)-U37*C(I,J, 7,12)
     3                       -U36*C(I,J, 6,12)-U35*C(I,J, 5,12)
     4                       -U34*C(I,J, 4,12)
            C(I,J, 2,12) = C2-U2C*C(I,J,12,12)-U2B*C(I,J,11,12)
     1                       -U2A*C(I,J,10,12)-U29*C(I,J, 9,12)
     2                       -U28*C(I,J, 8,12)-U27*C(I,J, 7,12)
     3                       -U26*C(I,J, 6,12)-U25*C(I,J, 5,12)
     4                       -U24*C(I,J, 4,12)-U23*C(I,J, 3,12)
            C(I,J, 1,12) = C1-U1C*C(I,J,12,12)-U1B*C(I,J,11,12)
     1                       -U1A*C(I,J,10,12)-U19*C(I,J, 9,12)
     2                       -U18*C(I,J, 8,12)-U17*C(I,J, 7,12)
     3                       -U16*C(I,J, 6,12)-U15*C(I,J, 5,12)
     4                       -U14*C(I,J, 4,12)-U13*C(I,J, 3,12)
     5                       -U12*C(I,J, 2,12)
C
  190    CONTINUE
  200 CONTINUE
C
C     J = NPTS Elimination
C     Do not do this loop if NPTS = 1
C
C     Eliminate A
C
C     B = B - A*C(j-1)
C
      J = NPTS
      IF (NPTS .GT. 1) THEN
         DO 260 L = 1, 12
            DO 240 M = 1, 12
               DO 230 I = 1, NSIM
                  B(I,J,L,M) = B(I,J,L,M)
     1                       - A(I,J,L, 1)*C(I,J-1, 1,M) 
     2                       - A(I,J,L, 2)*C(I,J-1, 2,M)
     3                       - A(I,J,L, 3)*C(I,J-1, 3,M) 
     4                       - A(I,J,L, 4)*C(I,J-1, 4,M)
     5                       - A(I,J,L, 5)*C(I,J-1, 5,M)
     6                       - A(I,J,L, 6)*C(I,J-1, 6,M)
     7                       - A(I,J,L, 7)*C(I,J-1, 7,M)
     8                       - A(I,J,L, 8)*C(I,J-1, 8,M)
     9                       - A(I,J,L, 9)*C(I,J-1, 9,M)
     A                       - A(I,J,L,10)*C(I,J-1,10,M)
     B                       - A(I,J,L,11)*C(I,J-1,11,M)
     C                       - A(I,J,L,12)*C(I,J-1,12,M)
  230          CONTINUE
  240       CONTINUE
C
C     R = R - A*R(j-1)
C
            DO 250 I = 1, NSIM
               R(I,J,L) = R(I,J,L) 
     1                  - A(I,J,L, 1)*R(I,J-1, 1) 
     2                  - A(I,J,L, 2)*R(I,J-1, 2) 
     3                  - A(I,J,L, 3)*R(I,J-1, 3) 
     4                  - A(I,J,L, 4)*R(I,J-1, 4) 
     5                  - A(I,J,L, 5)*R(I,J-1, 5)
     6                  - A(I,J,L, 6)*R(I,J-1, 6)
     7                  - A(I,J,L, 7)*R(I,J-1, 7)
     8                  - A(I,J,L, 8)*R(I,J-1, 8)
     9                  - A(I,J,L, 9)*R(I,J-1, 9)
     A                  - A(I,J,L,10)*R(I,J-1,10)
     B                  - A(I,J,L,11)*R(I,J-1,11)
     C                  - A(I,J,L,12)*R(I,J-1,12)
  250       CONTINUE
  260    CONTINUE
      ENDIF
C
C     Normalize row by multiplication of B(-1) to get unity diagonal
C
C     Form the LU of B(1)
C
      DO 300 I = 1, NSIM
         L11 =  B(I,J,1, 1)
         U12 =  B(I,J,1, 2) / L11
         U13 =  B(I,J,1, 3) / L11
         U14 =  B(I,J,1, 4) / L11
         U15 =  B(I,J,1, 5) / L11
         U16 =  B(I,J,1, 6) / L11
         U17 =  B(I,J,1, 7) / L11
         U18 =  B(I,J,1, 8) / L11
         U19 =  B(I,J,1, 9) / L11
         U1A =  B(I,J,1,10) / L11
         U1B =  B(I,J,1,11) / L11
         U1C =  B(I,J,1,12) / L11
C
         L21 =  B(I,J,2, 1) 
         L22 =  B(I,J,2, 2) - L21*U12
         U23 = (B(I,J,2, 3) - L21*U13) / L22
         U24 = (B(I,J,2, 4) - L21*U14) / L22
         U25 = (B(I,J,2, 5) - L21*U15) / L22
         U26 = (B(I,J,2, 6) - L21*U16) / L22
         U27 = (B(I,J,2, 7) - L21*U17) / L22
         U28 = (B(I,J,2, 8) - L21*U18) / L22
         U29 = (B(I,J,2, 9) - L21*U19) / L22
         U2A = (B(I,J,2,10) - L21*U1A) / L22
         U2B = (B(I,J,2,11) - L21*U1B) / L22
         U2C = (B(I,J,2,12) - L21*U1C) / L22
C
         L31 =  B(I,J,3, 1)
         L32 =  B(I,J,3, 2) - L31*U12
         L33 =  B(I,J,3, 3) - L31*U13 - L32*U23
         U34 = (B(I,J,3, 4) - L31*U14 - L32*U24) / L33
         U35 = (B(I,J,3, 5) - L31*U15 - L32*U25) / L33
         U36 = (B(I,J,3, 6) - L31*U16 - L32*U26) / L33
         U37 = (B(I,J,3, 7) - L31*U17 - L32*U27) / L33
         U38 = (B(I,J,3, 8) - L31*U18 - L32*U28) / L33
         U39 = (B(I,J,3, 9) - L31*U19 - L32*U29) / L33
         U3A = (B(I,J,3,10) - L31*U1A - L32*U2A) / L33
         U3B = (B(I,J,3,11) - L31*U1B - L32*U2B) / L33
         U3C = (B(I,J,3,12) - L31*U1C - L32*U2C) / L33
C
         L41 =  B(I,J,4, 1)
         L42 =  B(I,J,4, 2) - L41*U12
         L43 =  B(I,J,4, 3) - L41*U13 - L42*U23
         L44 =  B(I,J,4, 4) - L41*U14 - L42*U24 - L43*U34
         U45 = (B(I,J,4, 5) - L41*U15 - L42*U25 - L43*U35) / L44
         U46 = (B(I,J,4, 6) - L41*U16 - L42*U26 - L43*U36) / L44
         U47 = (B(I,J,4, 7) - L41*U17 - L42*U27 - L43*U37) / L44
         U48 = (B(I,J,4, 8) - L41*U18 - L42*U28 - L43*U38) / L44
         U49 = (B(I,J,4, 9) - L41*U19 - L42*U29 - L43*U39) / L44
         U4A = (B(I,J,4,10) - L41*U1A - L42*U2A - L43*U3A) / L44
         U4B = (B(I,J,4,11) - L41*U1B - L42*U2B - L43*U3B) / L44
         U4C = (B(I,J,4,12) - L41*U1C - L42*U2C - L43*U3C) / L44
C
         L51 =  B(I,J,5, 1)
         L52 =  B(I,J,5, 2) - L51*U12
         L53 =  B(I,J,5, 3) - L51*U13 - L52*U23
         L54 =  B(I,J,5, 4) - L51*U14 - L52*U24 - L53*U34
         L55 =  B(I,J,5, 5) - L51*U15 - L52*U25 - L53*U35 - L54*U45
         U56 = (B(I,J,5, 6) - L51*U16 - L52*U26 - L53*U36 - L54*U46) 
     1         / L55
         U57 = (B(I,J,5, 7) - L51*U17 - L52*U27 - L53*U37 - L54*U47) 
     1         / L55
         U58 = (B(I,J,5, 8) - L51*U18 - L52*U28 - L53*U38 - L54*U48) 
     1         / L55
         U59 = (B(I,J,5, 9) - L51*U19 - L52*U29 - L53*U39 - L54*U49) 
     1         / L55
         U5A = (B(I,J,5,10) - L51*U1A - L52*U2A - L53*U3A - L54*U4A) 
     1         / L55
         U5B = (B(I,J,5,11) - L51*U1B - L52*U2B - L53*U3B - L54*U4B) 
     1         / L55
         U5C = (B(I,J,5,12) - L51*U1C - L52*U2C - L53*U3C - L54*U4C) 
     1         / L55
C
         L61 =  B(I,J,6, 1)
         L62 =  B(I,J,6, 2) - L61*U12
         L63 =  B(I,J,6, 3) - L61*U13 - L62*U23
         L64 =  B(I,J,6, 4) - L61*U14 - L62*U24 - L63*U34
         L65 =  B(I,J,6, 5) - L61*U15 - L62*U25 - L63*U35 - L64*U45
         L66 =  B(I,J,6, 6) - L61*U16 - L62*U26 - L63*U36 - L64*U46
     1                      - L65*U56
         U67 = (B(I,J,6, 7) - L61*U17 - L62*U27 - L63*U37 - L64*U47
     1                      - L65*U57) / L66
         U68 = (B(I,J,6, 8) - L61*U18 - L62*U28 - L63*U38 - L64*U48
     1                      - L65*U58) / L66
         U69 = (B(I,J,6, 9) - L61*U19 - L62*U29 - L63*U39 - L64*U49
     1                      - L65*U59) / L66
         U6A = (B(I,J,6,10) - L61*U1A - L62*U2A - L63*U3A - L64*U4A
     1                      - L65*U5A) / L66
         U6B = (B(I,J,6,11) - L61*U1B - L62*U2B - L63*U3B - L64*U4B
     1                      - L65*U5B) / L66
         U6C = (B(I,J,6,12) - L61*U1C - L62*U2C - L63*U3C - L64*U4C
     1                      - L65*U5C) / L66
C
         L71 =  B(I,J,7, 1)
         L72 =  B(I,J,7, 2) - L71*U12
         L73 =  B(I,J,7, 3) - L71*U13 - L72*U23
         L74 =  B(I,J,7, 4) - L71*U14 - L72*U24 - L73*U34
         L75 =  B(I,J,7, 5) - L71*U15 - L72*U25 - L73*U35 - L74*U45
         L76 =  B(I,J,7, 6) - L71*U16 - L72*U26 - L73*U36 - L74*U46
     1                      - L75*U56
         L77 =  B(I,J,7, 7) - L71*U17 - L72*U27 - L73*U37 - L74*U47
     1                      - L75*U57 - L76*U67
         U78 = (B(I,J,7, 8) - L71*U18 - L72*U28 - L73*U38 - L74*U48
     1                      - L75*U58 - L76*U68) / L77
         U79 = (B(I,J,7, 9) - L71*U19 - L72*U29 - L73*U39 - L74*U49
     1                      - L75*U59 - L76*U69) / L77
         U7A = (B(I,J,7,10) - L71*U1A - L72*U2A - L73*U3A - L74*U4A
     1                      - L75*U5A - L76*U6A) / L77
         U7B = (B(I,J,7,11) - L71*U1B - L72*U2B - L73*U3B - L74*U4B
     1                      - L75*U5B - L76*U6B) / L77
         U7C = (B(I,J,7,12) - L71*U1C - L72*U2C - L73*U3C - L74*U4C
     1                      - L75*U5C - L76*U6C) / L77
C
         L81 =  B(I,J,8, 1)
         L82 =  B(I,J,8, 2) - L81*U12
         L83 =  B(I,J,8, 3) - L81*U13 - L82*U23
         L84 =  B(I,J,8, 4) - L81*U14 - L82*U24 - L83*U34
         L85 =  B(I,J,8, 5) - L81*U15 - L82*U25 - L83*U35 - L84*U45
         L86 =  B(I,J,8, 6) - L81*U16 - L82*U26 - L83*U36 - L84*U46
     1                      - L85*U56
         L87 =  B(I,J,8, 7) - L81*U17 - L82*U27 - L83*U37 - L84*U47
     1                      - L85*U57 - L86*U67
         L88 =  B(I,J,8, 8) - L81*U18 - L82*U28 - L83*U38 - L84*U48
     1                      - L85*U58 - L86*U68 - L87*U78
         U89 = (B(I,J,8, 9) - L81*U19 - L82*U29 - L83*U39 - L84*U49
     1                      - L85*U59 - L86*U69 - L87*U79) / L88
         U8A = (B(I,J,8,10) - L81*U1A - L82*U2A - L83*U3A - L84*U4A
     1                      - L85*U5A - L86*U6A - L87*U7A) / L88
         U8B = (B(I,J,8,11) - L81*U1B - L82*U2B - L83*U3B - L84*U4B
     1                      - L85*U5B - L86*U6B - L87*U7B) / L88
         U8C = (B(I,J,8,12) - L81*U1C - L82*U2C - L83*U3C - L84*U4C
     1                      - L85*U5C - L86*U6C - L87*U7C) / L88
C
         L91 =  B(I,J,9, 1)
         L92 =  B(I,J,9, 2) - L91*U12
         L93 =  B(I,J,9, 3) - L91*U13 - L92*U23
         L94 =  B(I,J,9, 4) - L91*U14 - L92*U24 - L93*U34
         L95 =  B(I,J,9, 5) - L91*U15 - L92*U25 - L93*U35 - L94*U45
         L96 =  B(I,J,9, 6) - L91*U16 - L92*U26 - L93*U36 - L94*U46
     1                      - L95*U56
         L97 =  B(I,J,9, 7) - L91*U17 - L92*U27 - L93*U37 - L94*U47
     1                      - L95*U57 - L96*U67
         L98 =  B(I,J,9, 8) - L91*U18 - L92*U28 - L93*U38 - L94*U48
     1                      - L95*U58 - L96*U68 - L97*U78
         L99 =  B(I,J,9, 9) - L91*U19 - L92*U29 - L93*U39 - L94*U49
     1                      - L95*U59 - L96*U69 - L97*U79 - L98*U89
         U9A = (B(I,J,9,10) - L91*U1A - L92*U2A - L93*U3A - L94*U4A
     1                      - L95*U5A - L96*U6A - L97*U7A - L98*U8A) 
     2         / L99
         U9B = (B(I,J,9,11) - L91*U1B - L92*U2B - L93*U3B - L94*U4B
     1                      - L95*U5B - L96*U6B - L97*U7B - L98*U8B) 
     2         / L99
         U9C = (B(I,J,9,12) - L91*U1C - L92*U2C - L93*U3C - L94*U4C
     1                      - L95*U5C - L96*U6C - L97*U7C - L98*U8C) 
     2         / L99
C
         LA1 =  B(I,J,10, 1)
         LA2 =  B(I,J,10, 2) - LA1*U12
         LA3 =  B(I,J,10, 3) - LA1*U13 - LA2*U23
         LA4 =  B(I,J,10, 4) - LA1*U14 - LA2*U24 - LA3*U34
         LA5 =  B(I,J,10, 5) - LA1*U15 - LA2*U25 - LA3*U35 - LA4*U45
         LA6 =  B(I,J,10, 6) - LA1*U16 - LA2*U26 - LA3*U36 - LA4*U46
     1                       - LA5*U56
         LA7 =  B(I,J,10, 7) - LA1*U17 - LA2*U27 - LA3*U37 - LA4*U47
     1                       - LA5*U57 - LA6*U67
         LA8 =  B(I,J,10, 8) - LA1*U18 - LA2*U28 - LA3*U38 - LA4*U48
     1                       - LA5*U58 - LA6*U68 - LA7*U78
         LA9 =  B(I,J,10, 9) - LA1*U19 - LA2*U29 - LA3*U39 - LA4*U49
     1                       - LA5*U59 - LA6*U69 - LA7*U79 - LA8*U89
         LAA =  B(I,J,10,10) - LA1*U1A - LA2*U2A - LA3*U3A - LA4*U4A
     1                       - LA5*U5A - LA6*U6A - LA7*U7A - LA8*U8A 
     2                       - LA9*U9A
         UAB = (B(I,J,10,11) - LA1*U1B - LA2*U2B - LA3*U3B - LA4*U4B
     1                       - LA5*U5B - LA6*U6B - LA7*U7B - LA8*U8B
     2                       - LA9*U9B) / LAA
         UAC = (B(I,J,10,12) - LA1*U1C - LA2*U2C - LA3*U3C - LA4*U4C
     1                       - LA5*U5C - LA6*U6C - LA7*U7C - LA8*U8C
     2                       - LA9*U9C) / LAA
C
         LB1 =  B(I,J,11, 1)
         LB2 =  B(I,J,11, 2) - LB1*U12
         LB3 =  B(I,J,11, 3) - LB1*U13 - LB2*U23
         LB4 =  B(I,J,11, 4) - LB1*U14 - LB2*U24 - LB3*U34
         LB5 =  B(I,J,11, 5) - LB1*U15 - LB2*U25 - LB3*U35 - LB4*U45
         LB6 =  B(I,J,11, 6) - LB1*U16 - LB2*U26 - LB3*U36 - LB4*U46
     1                       - LB5*U56
         LB7 =  B(I,J,11, 7) - LB1*U17 - LB2*U27 - LB3*U37 - LB4*U47
     1                       - LB5*U57 - LB6*U67
         LB8 =  B(I,J,11, 8) - LB1*U18 - LB2*U28 - LB3*U38 - LB4*U48
     1                       - LB5*U58 - LB6*U68 - LB7*U78
         LB9 =  B(I,J,11, 9) - LB1*U19 - LB2*U29 - LB3*U39 - LB4*U49
     1                       - LB5*U59 - LB6*U69 - LB7*U79 - LB8*U89
         LBA =  B(I,J,11,10) - LB1*U1A - LB2*U2A - LB3*U3A - LB4*U4A
     1                       - LB5*U5A - LB6*U6A - LB7*U7A - LB8*U8A 
     2                       - LB9*U9A
         LBB =  B(I,J,11,11) - LB1*U1B - LB2*U2B - LB3*U3B - LB4*U4B
     1                       - LB5*U5B - LB6*U6B - LB7*U7B - LB8*U8B
     2                       - LB9*U9B - LBA*UAB
         UBC = (B(I,J,11,12) - LB1*U1C - LB2*U2C - LB3*U3C - LB4*U4C
     1                       - LB5*U5C - LB6*U6C - LB7*U7C - LB8*U8C
     2                       - LB9*U9C - LBA*UAC) / LAA
C
         LC1 =  B(I,J,12, 1)
         LC2 =  B(I,J,12, 2) - LC1*U12
         LC3 =  B(I,J,12, 3) - LC1*U13 - LC2*U23
         LC4 =  B(I,J,12, 4) - LC1*U14 - LC2*U24 - LC3*U34
         LC5 =  B(I,J,12, 5) - LC1*U15 - LC2*U25 - LC3*U35 - LC4*U45
         LC6 =  B(I,J,12, 6) - LC1*U16 - LC2*U26 - LC3*U36 - LC4*U46
     1                       - LC5*U56
         LC7 =  B(I,J,12, 7) - LC1*U17 - LC2*U27 - LC3*U37 - LC4*U47
     1                       - LC5*U57 - LC6*U67
         LC8 =  B(I,J,12, 8) - LC1*U18 - LC2*U28 - LC3*U38 - LC4*U48
     1                       - LC5*U58 - LC6*U68 - LC7*U78
         LC9 =  B(I,J,12, 9) - LC1*U19 - LC2*U29 - LC3*U39 - LC4*U49
     1                       - LC5*U59 - LC6*U69 - LC7*U79 - LC8*U89
         LCA =  B(I,J,12,10) - LC1*U1A - LC2*U2A - LC3*U3A - LC4*U4A
     1                       - LC5*U5A - LC6*U6A - LC7*U7A - LC8*U8A 
     2                       - LC9*U9A
         LCB =  B(I,J,12,11) - LC1*U1B - LC2*U2B - LC3*U3B - LC4*U4B
     1                       - LC5*U5B - LC6*U6B - LC7*U7B - LC8*U8B
     2                       - LC9*U9B - LCA*UAB
         LCC =  B(I,J,12,12) - LC1*U1C - LC2*U2C - LC3*U3C - LC4*U4C
     1                       - LC5*U5C - LC6*U6C - LC7*U7C - LC8*U8C
     2                       - LC9*U9C - LCA*UAC - LCB*UBC
C
C     Calculate B(-1) * R
C
         R1        =  R(I,J, 1) / L11
         R2        = (R(I,J, 2) - L21*R1) / L22
         R3        = (R(I,J, 3) - L31*R1 - L32*R2) / L33
         R4        = (R(I,J, 4) - L41*R1 - L42*R2 - L43*R3) / L44
         R5        = (R(I,J, 5) - L51*R1 - L52*R2 - L53*R3 - L54*R4)
     1               / L55
         R6        = (R(I,J, 6) - L61*R1 - L62*R2 - L63*R3 - L64*R4 
     1                          - L65*R5) / L66
         R7        = (R(I,J, 7) - L71*R1 - L72*R2 - L73*R3 - L74*R4
     1                          - L75*R5 - L76*R6) / L77
         R8        = (R(I,J, 8) - L81*R1 - L82*R2 - L83*R3 - L84*R4
     1                          - L85*R5 - L86*R6 - L87*R7) / L88
         R9        = (R(I,J, 9) - L91*R1 - L92*R2 - L93*R3 - L94*R4
     1                          - L95*R5 - L96*R6 - L97*R7 - L98*R8) 
     2               / L99
         RA        = (R(I,J,10) - LA1*R1 - LA2*R2 - LA3*R3 - LA4*R4
     1                          - LA5*R5 - LA6*R6 - LA7*R7 - LA8*R8 
     2                          - LA9*R9) / LAA
         RB        = (R(I,J,11) - LB1*R1 - LB2*R2 - LB3*R3 - LB4*R4
     1                          - LB5*R5 - LB6*R6 - LB7*R7 - LB8*R8 
     2                          - LB9*R9 - LBA*RA) / LBB
         R(I,J,12) = (R(I,J,12) - LC1*R1 - LC2*R2 - LC3*R3 - LC4*R4
     1                          - LC5*R5 - LC6*R6 - LC7*R7 - LC8*R8 
     2                          - LC9*R9 - LCA*RA - LCB*RB) / LCC
         R(I,J,11) = RB-UBC*R(I,J,12)
         R(I,J,10) = RA-UAC*R(I,J,12)-UAB*R(I,J,11)
         R(I,J, 9) = R9-U9C*R(I,J,12)-U9B*R(I,J,11)-U9A*R(I,J,10)
         R(I,J, 8) = R8-U8C*R(I,J,12)-U8B*R(I,J,11)-U8A*R(I,J,10)
     1                 -U89*R(I,J, 9)
         R(I,J, 7) = R7-U7C*R(I,J,12)-U7B*R(I,J,11)-U7A*R(I,J,10)
     1                 -U79*R(I,J, 9)-U78*R(I,J, 8)
         R(I,J, 6) = R6-U6C*R(I,J,12)-U6B*R(I,J,11)-U6A*R(I,J,10)
     1                 -U69*R(I,J, 9)-U68*R(I,J, 8)-U67*R(I,J, 7)
         R(I,J, 5) = R5-U5C*R(I,J,12)-U5B*R(I,J,11)-U5A*R(I,J,10)
     1                 -U59*R(I,J, 9)-U58*R(I,J, 8)-U57*R(I,J, 7)
     2                 -U56*R(I,J, 6)
         R(I,J, 4) = R4-U4C*R(I,J,12)-U4B*R(I,J,11)-U4A*R(I,J,10)
     1                 -U49*R(I,J, 9)-U48*R(I,J, 8)-U47*R(I,J, 7)
     2                 -U46*R(I,J, 6)-U45*R(I,J, 5)
         R(I,J, 3) = R3-U3C*R(I,J,12)-U3B*R(I,J,11)-U3A*R(I,J,10)
     1                 -U39*R(I,J, 9)-U38*R(I,J, 8)-U37*R(I,J, 7)
     2                 -U36*R(I,J, 6)-U35*R(I,J, 5)-U34*R(I,J, 4)
         R(I,J, 2) = R2-U2C*R(I,J,12)-U2B*R(I,J,11)-U2A*R(I,J,10)
     1                 -U29*R(I,J, 9)-U28*R(I,J, 8)-U27*R(I,J, 7)
     2                 -U26*R(I,J, 6)-U25*R(I,J, 5)-U24*R(I,J, 4)
     3                 -U23*R(I,J, 3)
         R(I,J, 1) = R1-U1C*R(I,J,12)-U1B*R(I,J,11)-U1A*R(I,J,10)
     1                 -U19*R(I,J, 9)-U18*R(I,J, 8)-U17*R(I,J, 7)
     2                 -U16*R(I,J, 6)-U15*R(I,J, 5)-U14*R(I,J, 4)
     3                 -U13*R(I,J, 3)-U12*R(I,J, 2)
  300 CONTINUE
C
C     Backward Substitution Sweep
C
      DO 400 J = NPTS-1, 1, -1
         DO 390 L = 1, 12
            DO 380 I = 1, NSIM
               R(I,J,L) = R(I,J,L)
     1                  - C(I,J,L, 1)*R(I,J+1, 1) 
     2                  - C(I,J,L, 2)*R(I,J+1, 2)
     3                  - C(I,J,L, 3)*R(I,J+1, 3) 
     4                  - C(I,J,L, 4)*R(I,J+1, 4)
     5                  - C(I,J,L, 5)*R(I,J+1, 5)
     6                  - C(I,J,L, 6)*R(I,J+1, 6)
     7                  - C(I,J,L, 7)*R(I,J+1, 7)
     8                  - C(I,J,L, 8)*R(I,J+1, 8)
     9                  - C(I,J,L, 9)*R(I,J+1, 9)
     A                  - C(I,J,L,10)*R(I,J+1,10)
     B                  - C(I,J,L,11)*R(I,J+1,11)
     C                  - C(I,J,L,12)*R(I,J+1,12)
  380       CONTINUE
  390    CONTINUE
  400 CONTINUE
C
C     Finished
C
      RETURN
      END






      SUBROUTINE VSTRI  (NPTS, NSIM, A, B, C, R) 
C
C     Routine to perform the solve a scalar tridiagonal system 
C
C     NPTS:      The number of 
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     A:         The supra-diagonal element 
C     B:         The diagonal element
C     C:         The super-diagonal element 
C     R:         The right hand side
C
C     A, B and C are destroyed by this routine
C     R contains the solution vector at return
C
C     The storage of the system is such that vectorization occurs over
C     the first index to provide a unit stride for maximum speed.
C
      DIMENSION A (NSIM,NPTS),
     1          B (NSIM,NPTS),
     2          C (NSIM,NPTS),
     3          R (NSIM,NPTS)
C
C     Forward Substitution Loop 
C
C     J = 1
C     Scale first row by B(1)
C
      J = 1
c---CMIC$ DO ALL AUTOSCOPE VECTOR
c---C$DOACROSS LOCAL(I)
      DO 10 I = 1, NSIM
         C(I,J) = C(I,J) / B(I,J)
         R(I,J) = R(I,J) / B(I,J)
   10 CONTINUE
C
C     Loop over the interior of the inversion - 2 to NPTS
C     Eliminate A
C
      DO 200 J = 2, NPTS
C
C     R = (R - L(i,i-1)*R(j-1)) / L(i,i)
C
c---CMIC$ DO ALL AUTOSCOPE VECTOR
c---C$DOACROSS LOCAL(I)
         DO 150 I = 1, NSIM
            B(I,J) = (B(I,J) - A(I,J)*C(I,J-1))
            C(I,J) =  C(I,J) / B(I,J)
            R(I,J) = (R(I,J) - A(I,J)*R(I,J-1)) / B(I,J)
  150    CONTINUE
  200 CONTINUE
C
C     Backward Substitution Sweep 
C
      DO 400 J = NPTS-1, 1, -1
c---CMIC$ DO ALL AUTOSCOPE VECTOR
c---C$DOACROSS LOCAL(I)
         DO 350 I = 1, NSIM
            R(I,J) = R(I,J) - C(I,J)*R(I,J+1) 
  350    CONTINUE
  400 CONTINUE
C
C     Finished
C
      RETURN
      END





      SUBROUTINE VBDIA (NPTS, NSIM, IRANK, B, R) 
C
C     Routine to perform the Block Diagonal Inversion of
C     a block diagonal system where B is the block diagonal and
C     R is the right hand side which is returned as the solution
C     vector.  B is NOT overwritten.
C
C     NPTS:      The number of 
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     IRANK:     The rank of the system to invert
C     B:         The lu decomposition of the diagonal block 7x7 of L
C     R:         The right hand side
C
      include '../header/common.h'
C
      DIMENSION B (NSIM,NPTS,NRANK,NRANK),
     1          R (NSIM,NPTS,NRANK)
C
      IF (IRANK .EQ. 7) THEN
         CALL VB7DIA (NPTS, NSIM, NRANK, B, R)
      ELSE IF (IRANK .EQ. 12) THEN
         CALL VB12DI (NPTS, NSIM, NRANK, B, R)
      ELSE 
         WRITE (IOUT,1000) IRANK, IRANK
         STOP
      ENDIF
C
C     Finished 
C
      RETURN
 1000 FORMAT (' ','Block solver required for a ',I3,' x ',I3,' block.',
     1       /' ','This solver does not exist.  Run Aborting.')
      END






      SUBROUTINE VB7DIA (NPTS, NSIM, NRANK, B, R)
C
C     Routine to invert a 7x7 block diagonal system
C     Inversion is by Upper Triangularization with a
C     unity diagonal followed by a backsubstitution step.
C     Vectorization is attained by inverting multiple systems
C     at one time with the recursive inversion algorithm.
C     The B array is NOT destroyed in the solution procedure.
C
C     NPTS:      The number of 
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     NRANK:     The actual dimension of the storage
C     B:         The diagonal block 7x7
C     R:         The right hand side
C
C     R contains the solution vector at return
C
C     The storage of the system is such that vectorization occurs over
C     the first index to provide a unit stride for maximum speed.
C
      DIMENSION B (NSIM,NPTS,NRANK,NRANK),
     1          R (NSIM,NPTS,NRANK)
C
C     Storage for the LU decomposition of a block 7x7
C
      REAL L11, 
     1     L21, L22, 
     2     L31, L32, L33, 
     3     L41, L42, L43, L44, 
     4     L51, L52, L53, L54, L55,
     5     L61, L62, L63, L64, L65, L66,
     6     L71, L72, L73, L74, L75, L76, L77
C
C     Loop over all points
C
      DO 200 J = 1, NPTS
C
C     Form the LU of B(1)
C
         DO 190 I = 1, NSIM
            L11 =  B(I,J,1,1)
            U12 =  B(I,J,1,2) / L11
            U13 =  B(I,J,1,3) / L11
            U14 =  B(I,J,1,4) / L11
            U15 =  B(I,J,1,5) / L11
            U16 =  B(I,J,1,6) / L11
            U17 =  B(I,J,1,7) / L11
C
            L21 =  B(I,J,2,1) 
            L22 =  B(I,J,2,2) - L21*U12
            U23 = (B(I,J,2,3) - L21*U13) / L22
            U24 = (B(I,J,2,4) - L21*U14) / L22
            U25 = (B(I,J,2,5) - L21*U15) / L22
            U26 = (B(I,J,2,6) - L21*U16) / L22
            U27 = (B(I,J,2,7) - L21*U17) / L22
C
            L31 =  B(I,J,3,1)
            L32 =  B(I,J,3,2) - L31*U12
            L33 =  B(I,J,3,3) - L31*U13 - L32*U23
            U34 = (B(I,J,3,4) - L31*U14 - L32*U24) / L33
            U35 = (B(I,J,3,5) - L31*U15 - L32*U25) / L33
            U36 = (B(I,J,3,6) - L31*U16 - L32*U26) / L33
            U37 = (B(I,J,3,7) - L31*U17 - L32*U27) / L33
C
            L41 =  B(I,J,4,1)
            L42 =  B(I,J,4,2) - L41*U12
            L43 =  B(I,J,4,3) - L41*U13 - L42*U23
            L44 =  B(I,J,4,4) - L41*U14 - L42*U24 - L43*U34
            U45 = (B(I,J,4,5) - L41*U15 - L42*U25 - L43*U35) / L44
            U46 = (B(I,J,4,6) - L41*U16 - L42*U26 - L43*U36) / L44
            U47 = (B(I,J,4,7) - L41*U17 - L42*U27 - L43*U37) / L44
C
            L51 =  B(I,J,5,1)
            L52 =  B(I,J,5,2) - L51*U12
            L53 =  B(I,J,5,3) - L51*U13 - L52*U23
            L54 =  B(I,J,5,4) - L51*U14 - L52*U24 - L53*U34
            L55 =  B(I,J,5,5) - L51*U15 - L52*U25 - L53*U35 - L54*U45
            U56 = (B(I,J,5,6) - L51*U16 - L52*U26 - L53*U36 - L54*U46) 
     1            / L55
            U57 = (B(I,J,5,7) - L51*U17 - L52*U27 - L53*U37 - L54*U47) 
     1            / L55
C
            L61 =  B(I,J,6,1)
            L62 =  B(I,J,6,2) - L61*U12
            L63 =  B(I,J,6,3) - L61*U13 - L62*U23
            L64 =  B(I,J,6,4) - L61*U14 - L62*U24 - L63*U34
            L65 =  B(I,J,6,5) - L61*U15 - L62*U25 - L63*U35 - L64*U45
            L66 =  B(I,J,6,6) - L61*U16 - L62*U26 - L63*U36 - L64*U46
     1                        - L65*U56
            U67 = (B(I,J,6,7) - L61*U17 - L62*U27 - L63*U37 - L64*U47
     1                        - L65*U57) / L66
C
            L71 =  B(I,J,7,1)
            L72 =  B(I,J,7,2) - L71*U12
            L73 =  B(I,J,7,3) - L71*U13 - L72*U23
            L74 =  B(I,J,7,4) - L71*U14 - L72*U24 - L73*U34
            L75 =  B(I,J,7,5) - L71*U15 - L72*U25 - L73*U35 - L74*U45
            L76 =  B(I,J,7,6) - L71*U16 - L72*U26 - L73*U36 - L74*U46
     1                        - L75*U56
            L77 =  B(I,J,7,7) - L71*U17 - L72*U27 - L73*U37 - L74*U47
     1                        - L75*U57 - L76*U67
C
C     Calculate B(-1) * R
C
            R1       =  R(I,J,1) / L11
            R2       = (R(I,J,2) - L21*R1) / L22
            R3       = (R(I,J,3) - L31*R1 - L32*R2) / L33
            R4       = (R(I,J,4) - L41*R1 - L42*R2 - L43*R3) / L44
            R5       = (R(I,J,5) - L51*R1 - L52*R2 - L53*R3 - L54*R4)
     1                 / L55
            R6       = (R(I,J,6) - L61*R1 - L62*R2 - L63*R3 - L64*R4 
     1                           - L65*R5) / L66
            R(I,J,7) = (R(I,J,7) - L71*R1 - L72*R2 - L73*R3 - L74*R4
     1                           - L75*R5 - L76*R6) / L77
            R(I,J,6) =  R6-U67*R(I,J,7)
            R(I,J,5) =  R5-U57*R(I,J,7)-U56*R(I,J,6)
            R(I,J,4) =  R4-U47*R(I,J,7)-U46*R(I,J,6)-U45*R(I,J,5)
            R(I,J,3) =  R3-U37*R(I,J,7)-U36*R(I,J,6)-U35*R(I,J,5)
     1                    -U34*R(I,J,4)
            R(I,J,2) =  R2-U27*R(I,J,7)-U26*R(I,J,6)-U25*R(I,J,5)
     1                    -U24*R(I,J,4)-U23*R(I,J,3)
            R(I,J,1) =  R1-U17*R(I,J,7)-U16*R(I,J,6)-U15*R(I,J,5)
     1                    -U14*R(I,J,4)-U13*R(I,J,3)-U12*R(I,J,2)
  190    CONTINUE
  200 CONTINUE
C
C     Finished
C
      RETURN
      END






      SUBROUTINE VB12DI (NPTS, NSIM, NRANK, B, R)
C
C     Routine to invert a 12x12 block diagonal system
C     Inversion is by Upper Triangularization with a
C     unity diagonal followed by a backsubstitution step.
C     Vectorization is attained by inverting multiple systems
C     at one time with the recursive inversion algorithm.
C     The B array is NOT destroyed in the solution procedure.
C
C     NPTS:      The number of 
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     NRANK:     The actual dimension of the storage
C     B:         The diagonal block 12x12
C     R:         The right hand side
C
C     R contains the solution vector at return
C
C     The storage of the system is such that vectorization occurs over
C     the first index to provide a unit stride for maximum speed.
C
      DIMENSION B (NSIM,NPTS,NRANK,NRANK),
     1          R (NSIM,NPTS,NRANK)
C
C     Storage for the LU decomposition of a block 12x12
C     Note: A = 10, B = 11, C = 12
C
      REAL L11, 
     1     L21, L22, 
     2     L31, L32, L33, 
     3     L41, L42, L43, L44, 
     4     L51, L52, L53, L54, L55,
     5     L61, L62, L63, L64, L65, L66,
     6     L71, L72, L73, L74, L75, L76, L77,
     7     L81, L82, L83, L84, L85, L86, L87, L88,
     8     L91, L92, L93, L94, L95, L96, L97, L98, L99,
     9     LA1, LA2, LA3, LA4, LA5, LA6, LA7, LA8, LA9, LAA,
     A     LB1, LB2, LB3, LB4, LB5, LB6, LB7, LB8, LB9, LBA, LBB,
     B     LC1, LC2, LC3, LC4, LC5, LC6, LC7, LC8, LC9, LCA, LCB, LCC
C
C     Loop over all points
C
      DO 200 J = 1, NPTS
C
C     Form the LU of B(1)
C
C Row 1
C
         DO 190 I = 1, NSIM
            L11 =  B(I,J,1, 1)
            U12 =  B(I,J,1, 2) / L11
            U13 =  B(I,J,1, 3) / L11
            U14 =  B(I,J,1, 4) / L11
            U15 =  B(I,J,1, 5) / L11
            U16 =  B(I,J,1, 6) / L11
            U17 =  B(I,J,1, 7) / L11
            U18 =  B(I,J,1, 8) / L11
            U19 =  B(I,J,1, 9) / L11
            U1A =  B(I,J,1,10) / L11
            U1B =  B(I,J,1,11) / L11
            U1C =  B(I,J,1,12) / L11
C
C Row 2
C
            L21 =  B(I,J,2, 1) 
            L22 =  B(I,J,2, 2) - L21*U12
            U23 = (B(I,J,2, 3) - L21*U13) / L22
            U24 = (B(I,J,2, 4) - L21*U14) / L22
            U25 = (B(I,J,2, 5) - L21*U15) / L22
            U26 = (B(I,J,2, 6) - L21*U16) / L22
            U27 = (B(I,J,2, 7) - L21*U17) / L22
            U28 = (B(I,J,2, 8) - L21*U18) / L22
            U29 = (B(I,J,2, 9) - L21*U19) / L22
            U2A = (B(I,J,2,10) - L21*U1A) / L22
            U2B = (B(I,J,2,11) - L21*U1B) / L22
            U2C = (B(I,J,2,12) - L21*U1C) / L22
C
C Row 3
C
            L31 =  B(I,J,3, 1)
            L32 =  B(I,J,3, 2) - L31*U12
            L33 =  B(I,J,3, 3) - L31*U13 - L32*U23
            U34 = (B(I,J,3, 4) - L31*U14 - L32*U24) / L33
            U35 = (B(I,J,3, 5) - L31*U15 - L32*U25) / L33
            U36 = (B(I,J,3, 6) - L31*U16 - L32*U26) / L33
            U37 = (B(I,J,3, 7) - L31*U17 - L32*U27) / L33
            U38 = (B(I,J,3, 8) - L31*U18 - L32*U28) / L33
            U39 = (B(I,J,3, 9) - L31*U19 - L32*U29) / L33
            U3A = (B(I,J,3,10) - L31*U1A - L32*U2A) / L33
            U3B = (B(I,J,3,11) - L31*U1B - L32*U2B) / L33
            U3C = (B(I,J,3,12) - L31*U1C - L32*U2C) / L33
C
C Row 4
C
            L41 =  B(I,J,4, 1)
            L42 =  B(I,J,4, 2) - L41*U12
            L43 =  B(I,J,4, 3) - L41*U13 - L42*U23
            L44 =  B(I,J,4, 4) - L41*U14 - L42*U24 - L43*U34
            U45 = (B(I,J,4, 5) - L41*U15 - L42*U25 - L43*U35) / L44
            U46 = (B(I,J,4, 6) - L41*U16 - L42*U26 - L43*U36) / L44
            U47 = (B(I,J,4, 7) - L41*U17 - L42*U27 - L43*U37) / L44
            U48 = (B(I,J,4, 8) - L41*U18 - L42*U28 - L43*U38) / L44
            U49 = (B(I,J,4, 9) - L41*U19 - L42*U29 - L43*U39) / L44
            U4A = (B(I,J,4,10) - L41*U1A - L42*U2A - L43*U3A) / L44
            U4B = (B(I,J,4,11) - L41*U1B - L42*U2B - L43*U3B) / L44
            U4C = (B(I,J,4,12) - L41*U1C - L42*U2C - L43*U3C) / L44
C
C Row 5
C
            L51 =  B(I,J,5, 1)
            L52 =  B(I,J,5, 2) - L51*U12
            L53 =  B(I,J,5, 3) - L51*U13 - L52*U23
            L54 =  B(I,J,5, 4) - L51*U14 - L52*U24 - L53*U34
            L55 =  B(I,J,5, 5) - L51*U15 - L52*U25 - L53*U35 - L54*U45
            U56 = (B(I,J,5, 6) - L51*U16 - L52*U26 - L53*U36 - L54*U46) 
     1            / L55
            U57 = (B(I,J,5, 7) - L51*U17 - L52*U27 - L53*U37 - L54*U47) 
     1            / L55
            U58 = (B(I,J,5, 8) - L51*U18 - L52*U28 - L53*U38 - L54*U48) 
     1            / L55
            U59 = (B(I,J,5, 9) - L51*U19 - L52*U29 - L53*U39 - L54*U49) 
     1            / L55
            U5A = (B(I,J,5,10) - L51*U1A - L52*U2A - L53*U3A - L54*U4A) 
     1            / L55
            U5B = (B(I,J,5,11) - L51*U1B - L52*U2B - L53*U3B - L54*U4B) 
     1            / L55
            U5C = (B(I,J,5,12) - L51*U1C - L52*U2C - L53*U3C - L54*U4C) 
     1            / L55
C
C Row 6
C
            L61 =  B(I,J,6, 1)
            L62 =  B(I,J,6, 2) - L61*U12
            L63 =  B(I,J,6, 3) - L61*U13 - L62*U23
            L64 =  B(I,J,6, 4) - L61*U14 - L62*U24 - L63*U34
            L65 =  B(I,J,6, 5) - L61*U15 - L62*U25 - L63*U35 - L64*U45
            L66 =  B(I,J,6, 6) - L61*U16 - L62*U26 - L63*U36 - L64*U46
     1                         - L65*U56
            U67 = (B(I,J,6, 7) - L61*U17 - L62*U27 - L63*U37 - L64*U47
     1                         - L65*U57) / L66
            U68 = (B(I,J,6, 8) - L61*U18 - L62*U28 - L63*U38 - L64*U48
     1                         - L65*U58) / L66
            U69 = (B(I,J,6, 9) - L61*U19 - L62*U29 - L63*U39 - L64*U49
     1                         - L65*U59) / L66
            U6A = (B(I,J,6,10) - L61*U1A - L62*U2A - L63*U3A - L64*U4A
     1                         - L65*U5A) / L66
            U6B = (B(I,J,6,11) - L61*U1B - L62*U2B - L63*U3B - L64*U4B
     1                         - L65*U5B) / L66
            U6C = (B(I,J,6,12) - L61*U1C - L62*U2C - L63*U3C - L64*U4C
     1                         - L65*U5C) / L66
C
C Row 7
C
            L71 =  B(I,J,7, 1)
            L72 =  B(I,J,7, 2) - L71*U12
            L73 =  B(I,J,7, 3) - L71*U13 - L72*U23
            L74 =  B(I,J,7, 4) - L71*U14 - L72*U24 - L73*U34
            L75 =  B(I,J,7, 5) - L71*U15 - L72*U25 - L73*U35 - L74*U45
            L76 =  B(I,J,7, 6) - L71*U16 - L72*U26 - L73*U36 - L74*U46
     1                         - L75*U56
            L77 =  B(I,J,7, 7) - L71*U17 - L72*U27 - L73*U37 - L74*U47
     1                         - L75*U57 - L76*U67
            U78 = (B(I,J,7, 8) - L71*U18 - L72*U28 - L73*U38 - L74*U48
     1                         - L75*U58 - L76*U68) / L77
            U79 = (B(I,J,7, 9) - L71*U19 - L72*U29 - L73*U39 - L74*U49
     1                         - L75*U59 - L76*U69) / L77
            U7A = (B(I,J,7,10) - L71*U1A - L72*U2A - L73*U3A - L74*U4A
     1                         - L75*U5A - L76*U6A) / L77
            U7B = (B(I,J,7,11) - L71*U1B - L72*U2B - L73*U3B - L74*U4B
     1                         - L75*U5B - L76*U6B) / L77
            U7C = (B(I,J,7,12) - L71*U1C - L72*U2C - L73*U3C - L74*U4C
     1                         - L75*U5C - L76*U6C) / L77
C
C Row 8
C
            L81 =  B(I,J,8, 1)
            L82 =  B(I,J,8, 2) - L81*U12
            L83 =  B(I,J,8, 3) - L81*U13 - L82*U23
            L84 =  B(I,J,8, 4) - L81*U14 - L82*U24 - L83*U34
            L85 =  B(I,J,8, 5) - L81*U15 - L82*U25 - L83*U35 - L84*U45
            L86 =  B(I,J,8, 6) - L81*U16 - L82*U26 - L83*U36 - L84*U46
     1                         - L85*U56
            L87 =  B(I,J,8, 7) - L81*U17 - L82*U27 - L83*U37 - L84*U47
     1                         - L85*U57 - L86*U67
            L88 =  B(I,J,8, 8) - L81*U18 - L82*U28 - L83*U38 - L84*U48
     1                         - L85*U58 - L86*U68 - L87*U78
            U89 = (B(I,J,8, 9) - L81*U19 - L82*U29 - L83*U39 - L84*U49
     1                         - L85*U59 - L86*U69 - L87*U79) / L88
            U8A = (B(I,J,8,10) - L81*U1A - L82*U2A - L83*U3A - L84*U4A
     1                         - L85*U5A - L86*U6A - L87*U7A) / L88
            U8B = (B(I,J,8,11) - L81*U1B - L82*U2B - L83*U3B - L84*U4B
     1                         - L85*U5B - L86*U6B - L87*U7B) / L88
            U8C = (B(I,J,8,12) - L81*U1C - L82*U2C - L83*U3C - L84*U4C
     1                         - L85*U5C - L86*U6C - L87*U7C) / L88
C
C Row 9
C
            L91 =  B(I,J,9, 1)
            L92 =  B(I,J,9, 2) - L91*U12
            L93 =  B(I,J,9, 3) - L91*U13 - L92*U23
            L94 =  B(I,J,9, 4) - L91*U14 - L92*U24 - L93*U34
            L95 =  B(I,J,9, 5) - L91*U15 - L92*U25 - L93*U35 - L94*U45
            L96 =  B(I,J,9, 6) - L91*U16 - L92*U26 - L93*U36 - L94*U46
     1                         - L95*U56
            L97 =  B(I,J,9, 7) - L91*U17 - L92*U27 - L93*U37 - L94*U47
     1                         - L95*U57 - L96*U67
            L98 =  B(I,J,9, 8) - L91*U18 - L92*U28 - L93*U38 - L94*U48
     1                         - L95*U58 - L96*U68 - L97*U78
            L99 =  B(I,J,9, 9) - L91*U19 - L92*U29 - L93*U39 - L94*U49
     1                         - L95*U59 - L96*U69 - L97*U79 - L98*U89
            U9A = (B(I,J,9,10) - L91*U1A - L92*U2A - L93*U3A - L94*U4A
     1                         - L95*U5A - L96*U6A - L97*U7A - L98*U8A)
     2            / L99
            U9B = (B(I,J,9,11) - L91*U1B - L92*U2B - L93*U3B - L94*U4B
     1                         - L95*U5B - L96*U6B - L97*U7B - L98*U8B)
     2            / L99
            U9C = (B(I,J,9,12) - L91*U1C - L92*U2C - L93*U3C - L94*U4C
     1                         - L95*U5C - L96*U6C - L97*U7C - L98*U8C)
     2            / L99
C
C Row 10
C
            LA1 =  B(I,J,10, 1)
            LA2 =  B(I,J,10, 2) - LA1*U12
            LA3 =  B(I,J,10, 3) - LA1*U13 - LA2*U23
            LA4 =  B(I,J,10, 4) - LA1*U14 - LA2*U24 - LA3*U34
            LA5 =  B(I,J,10, 5) - LA1*U15 - LA2*U25 - LA3*U35 - LA4*U45
            LA6 =  B(I,J,10, 6) - LA1*U16 - LA2*U26 - LA3*U36 - LA4*U46
     1                          - LA5*U56
            LA7 =  B(I,J,10, 7) - LA1*U17 - LA2*U27 - LA3*U37 - LA4*U47
     1                          - LA5*U57 - LA6*U67
            LA8 =  B(I,J,10, 8) - LA1*U18 - LA2*U28 - LA3*U38 - LA4*U48
     1                          - LA5*U58 - LA6*U68 - LA7*U78
            LA9 =  B(I,J,10, 9) - LA1*U19 - LA2*U29 - LA3*U39 - LA4*U49
     1                          - LA5*U59 - LA6*U69 - LA7*U79 - LA8*U89
            LAA =  B(I,J,10,10) - LA1*U1A - LA2*U2A - LA3*U3A - LA4*U4A
     1                          - LA5*U5A - LA6*U6A - LA7*U7A - LA8*U8A
     2                          - LA9*U9A
            UAB = (B(I,J,10,11) - LA1*U1B - LA2*U2B - LA3*U3B - LA4*U4B
     1                          - LA5*U5B - LA6*U6B - LA7*U7B - LA8*U8B 
     2                          - LA9*U9B) / LAA
            UAC = (B(I,J,10,12) - LA1*U1C - LA2*U2C - LA3*U3C - LA4*U4C
     1                          - LA5*U5C - LA6*U6C - LA7*U7C - LA8*U8C 
     2                          - LA9*U9C) / LAA
C
C Row 11
C
            LB1 =  B(I,J,11, 1)
            LB2 =  B(I,J,11, 2) - LB1*U12
            LB3 =  B(I,J,11, 3) - LB1*U13 - LB2*U23
            LB4 =  B(I,J,11, 4) - LB1*U14 - LB2*U24 - LB3*U34
            LB5 =  B(I,J,11, 5) - LB1*U15 - LB2*U25 - LB3*U35 - LB4*U45
            LB6 =  B(I,J,11, 6) - LB1*U16 - LB2*U26 - LB3*U36 - LB4*U46
     1                          - LB5*U56
            LB7 =  B(I,J,11, 7) - LB1*U17 - LB2*U27 - LB3*U37 - LB4*U47
     1                          - LB5*U57 - LB6*U67
            LB8 =  B(I,J,11, 8) - LB1*U18 - LB2*U28 - LB3*U38 - LB4*U48
     1                          - LB5*U58 - LB6*U68 - LB7*U78
            LB9 =  B(I,J,11, 9) - LB1*U19 - LB2*U29 - LB3*U39 - LB4*U49
     1                          - LB5*U59 - LB6*U69 - LB7*U79 - LB8*U89
            LBA =  B(I,J,11,10) - LB1*U1A - LB2*U2A - LB3*U3A - LB4*U4A
     1                          - LB5*U5A - LB6*U6A - LB7*U7A - LB8*U8A
     2                          - LB9*U9A
            LBB =  B(I,J,11,11) - LB1*U1B - LB2*U2B - LB3*U3B - LB4*U4B
     1                          - LB5*U5B - LB6*U6B - LB7*U7B - LB8*U8B 
     2                          - LB9*U9B - LBA*UAB
            UBC = (B(I,J,11,12) - LB1*U1C - LB2*U2C - LB3*U3C - LB4*U4C
     1                          - LB5*U5C - LB6*U6C - LB7*U7C - LB8*U8C 
     2                          - LB9*U9C - LBA*UAC) / LBB
C
C Row 12
C
            LC1 =  B(I,J,12, 1)
            LC2 =  B(I,J,12, 2) - LC1*U12
            LC3 =  B(I,J,12, 3) - LC1*U13 - LC2*U23
            LC4 =  B(I,J,12, 4) - LC1*U14 - LC2*U24 - LC3*U34
            LC5 =  B(I,J,12, 5) - LC1*U15 - LC2*U25 - LC3*U35 - LC4*U45
            LC6 =  B(I,J,12, 6) - LC1*U16 - LC2*U26 - LC3*U36 - LC4*U46
     1                          - LC5*U56
            LC7 =  B(I,J,12, 7) - LC1*U17 - LC2*U27 - LC3*U37 - LC4*U47
     1                          - LC5*U57 - LC6*U67
            LC8 =  B(I,J,12, 8) - LC1*U18 - LC2*U28 - LC3*U38 - LC4*U48
     1                          - LC5*U58 - LC6*U68 - LC7*U78
            LC9 =  B(I,J,12, 9) - LC1*U19 - LC2*U29 - LC3*U39 - LC4*U49
     1                          - LC5*U59 - LC6*U69 - LC7*U79 - LC8*U89
            LCA =  B(I,J,12,10) - LC1*U1A - LC2*U2A - LC3*U3A - LC4*U4A
     1                          - LC5*U5A - LC6*U6A - LC7*U7A - LC8*U8A
     2                          - LC9*U9A
            LCB =  B(I,J,12,11) - LC1*U1B - LC2*U2B - LC3*U3B - LC4*U4B
     1                          - LC5*U5B - LC6*U6B - LC7*U7B - LC8*U8B 
     2                          - LC9*U9B - LCA*UAB
            LCC =  B(I,J,12,12) - LC1*U1C - LC2*U2C - LC3*U3C - LC4*U4C
     1                          - LC5*U5C - LC6*U6C - LC7*U7C - LC8*U8C 
     2                          - LC9*U9C - LCA*UAC - LCB*UBC
C
C     Calculate B(-1) * R
C
            R1        =  R(I,J, 1) / L11
            R2        = (R(I,J, 2) - L21*R1) / L22
            R3        = (R(I,J, 3) - L31*R1 - L32*R2) / L33
            R4        = (R(I,J, 4) - L41*R1 - L42*R2 - L43*R3) / L44
            R5        = (R(I,J, 5) - L51*R1 - L52*R2 - L53*R3 - L54*R4)
     1                  / L55
            R6        = (R(I,J, 6) - L61*R1 - L62*R2 - L63*R3 - L64*R4 
     1                             - L65*R5) / L66
            R7        = (R(I,J, 7) - L71*R1 - L72*R2 - L73*R3 - L74*R4
     1                             - L75*R5 - L76*R6) / L77
            R8        = (R(I,J, 8) - L81*R1 - L82*R2 - L83*R3 - L84*R4
     1                             - L85*R5 - L86*R6 - L87*R7) / L88
            R9        = (R(I,J, 9) - L91*R1 - L92*R2 - L93*R3 - L94*R4
     1                             - L95*R5 - L96*R6 - L97*R7 - L98*R8) 
     2                  / L99
            RA        = (R(I,J,10) - LA1*R1 - LA2*R2 - LA3*R3 - LA4*R4
     1                             - LA5*R5 - LA6*R6 - LA7*R7 - LA8*R8
     2                             - LA9*R9) / LAA
            RB        = (R(I,J,11) - LB1*R1 - LB2*R2 - LB3*R3 - LB4*R4
     1                             - LB5*R5 - LB6*R6 - LB7*R7 - LB8*R8
     2                             - LB9*R9 - LBA*RA) / LBB
            RC        = (R(I,J,12) - LC1*R1 - LC2*R2 - LC3*R3 - LC4*R4
     1                             - LC5*R5 - LC6*R6 - LC7*R7 - LC8*R8
     2                             - LC9*R9 - LCA*RA - LCB*RB) / LCC
            R(I,J,12) = RC
            R(I,J,11) = RB-UBC*R(I,J,12)
            R(I,J,10) = RA-UAC*R(I,J,12)-UAB*R(I,J,11)
            R(I,J, 9) = R9-U9C*R(I,J,12)-U9B*R(I,J,11)-U9A*R(I,J,10)
            R(I,J, 8) = R8-U8C*R(I,J,12)-U8B*R(I,J,11)-U8A*R(I,J,10)
     1                    -U89*R(I,J, 9)
            R(I,J, 7) = R7-U7C*R(I,J,12)-U7B*R(I,J,11)-U7A*R(I,J,10)
     1                    -U79*R(I,J, 9)-U78*R(I,J, 8)
            R(I,J, 6) = R6-U6C*R(I,J,12)-U6B*R(I,J,11)-U6A*R(I,J,10)
     1                    -U69*R(I,J, 9)-U68*R(I,J, 8)-U67*R(I,J, 7)
            R(I,J, 5) = R5-U5C*R(I,J,12)-U5B*R(I,J,11)-U5A*R(I,J,10)
     1                    -U59*R(I,J, 9)-U58*R(I,J, 8)-U57*R(I,J, 7)
     2                    -U56*R(I,J, 6)
            R(I,J, 4) = R4-U4C*R(I,J,12)-U4B*R(I,J,11)-U4A*R(I,J,10)
     1                    -U49*R(I,J, 9)-U48*R(I,J, 8)-U47*R(I,J, 7)
     2                    -U46*R(I,J, 6)-U45*R(I,J, 5)
            R(I,J, 3) = R3-U3C*R(I,J,12)-U3B*R(I,J,11)-U3A*R(I,J,10)
     1                    -U39*R(I,J, 9)-U38*R(I,J, 8)-U37*R(I,J, 7)
     2                    -U36*R(I,J, 6)-U35*R(I,J, 5)-U34*R(I,J, 4)
            R(I,J, 2) = R2-U2C*R(I,J,12)-U2B*R(I,J,11)-U2A*R(I,J,10)
     1                    -U29*R(I,J, 9)-U28*R(I,J, 8)-U27*R(I,J, 7)
     2                    -U26*R(I,J, 6)-U25*R(I,J, 5)-U24*R(I,J, 4)
     3                    -U23*R(I,J, 3)
            R(I,J, 1) = R1-U1C*R(I,J,12)-U1B*R(I,J,11)-U1A*R(I,J,10)
     1                    -U19*R(I,J, 9)-U18*R(I,J, 8)-U17*R(I,J, 7)
     2                    -U16*R(I,J, 6)-U15*R(I,J, 5)-U14*R(I,J, 4)
     3                    -U13*R(I,J, 3)-U12*R(I,J, 2)
  190    CONTINUE
  200 CONTINUE
C
C     Finished
C
      RETURN
      END






      SUBROUTINE VFBTR (NPTS, NSIM, IRANK, A, B, C, R) 
C
C     Routine to perform the Forward/Backward substitution sweeps on 
C     a block tridiagonal system where A, B and C contain
C     the LU decomposition of the block tridiagonal system.
C
C     VLUTR5 calculates the LU decomposition of the block tridiagonal system
C     VFBTR5 performs the forward/backward sweeps to complete the inversion
C
C     NPTS:      The number of 
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     IRANK:     The rank of the system to invert
C     A:         The supra-diagonal block 7x7 of L
C     B:         The lu decomposition of the diagonal block 7x7 of L
C     C:         The super-diagonal block 7x7 of U
C     R:         The right hand side
C
      include '../header/common.h'
C
      DIMENSION A (NSIM,NPTS,NRANK,NRANK),
     1          B (NSIM,NPTS,NRANK,NRANK),
     2          C (NSIM,NPTS,NRANK,NRANK),
     3          R (NSIM,NPTS,NRANK)
C
      IF (IRANK .EQ. 5) THEN
         CALL VFBTR5 (NPTS, NSIM, A, B, C, R)
      ELSE IF (IRANK .EQ. 7) THEN
         CALL VFBTR7 (NPTS, NSIM, NRANK, A, B, C, R)
      ELSE IF (IRANK .EQ. 12) THEN
         CALL VFBT12 (NPTS, NSIM, NRANK, A, B, C, R)
      ELSE 
         WRITE (IOUT,1000) IRANK, IRANK
         STOP
      ENDIF
C
C     Finished 
C
      RETURN
 1000 FORMAT (' ','Block solver required for a ',I3,' x ',I3,' block.',
     1       /' ','This solver does not exist.  Run Aborting.')
      END






      SUBROUTINE VFBTR5 (NPTS, NSIM, A, B, C, R) 
C
C     Routine to perform the Forward/Backward substitution sweeps on 
C     a 5x5 block tridiagonal system where A, B and C contain
C     the LU decomposition of the block tridiagonal system.
C
C     VLUTR5 calculates the LU decomposition of the block tridiagonal system
C     VFBTR5 performs the forward/backward sweeps to complete the inversion
C
C     NPTS:      The number of 
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     A:         The supra-diagonal block 5x5 of L
C     B:         The lu decomposition of the diagonal block 5x5 of L
C     C:         The super-diagonal block 5x5 of U
C     R:         The right hand side
C
C     A, B and C are unmodified by this routine
C     R contains the solution vector at return
C
C     The storage of the system is such that vectorization occurs over
C     the first index to provide a unit stride for maximum speed.
C
      DIMENSION A (NSIM,NPTS,5,5),
     1          B (NSIM,NPTS,5,5),
     2          C (NSIM,NPTS,5,5),
     3          R (NSIM,NPTS,5)
C
C     Forward Substitution Loop 
C
C     J = 1
C
      J = 1
      DO 10 I = 1, NSIM
C
C     Calculate L[inverse](i,i) * R
C
         R1       =  R(I,J,1) / B(I,J,1,1)
         R2       = (R(I,J,2)-B(I,J,2,1)*R1) / B(I,J,2,2)
         R3       = (R(I,J,3)-B(I,J,3,1)*R1-B(I,J,3,2)*R2) / B(I,J,3,3)
         R4       = (R(I,J,4)-B(I,J,4,1)*R1-B(I,J,4,2)*R2 
     1                       -B(I,J,4,3)*R3)               / B(I,J,4,4)
         R(I,J,5) = (R(I,J,5)-B(I,J,5,1)*R1-B(I,J,5,2)*R2 
     1                       -B(I,J,5,3)*R3-B(I,J,5,4)*R4) / B(I,J,5,5)
         R(I,J,4) =  R4-B(I,J,4,5)*R(I,J,5)
         R(I,J,3) =  R3-B(I,J,3,5)*R(I,J,5)-B(I,J,3,4)*R(I,J,4)
         R(I,J,2) =  R2-B(I,J,2,5)*R(I,J,5)-B(I,J,2,4)*R(I,J,4)
     1                 -B(I,J,2,3)*R(I,J,3)
         R(I,J,1) =  R1-B(I,J,1,5)*R(I,J,5)-B(I,J,1,4)*R(I,J,4)
     1                 -B(I,J,1,3)*R(I,J,3)-B(I,J,1,2)*R(I,J,2)
   10 CONTINUE
C
C     Loop over the interior of the inversion - 2 to NPTS
C
      DO 200 J = 2, NPTS
C
C     R = R - L(i,i-1)*R(j-1)
C
         DO 160 l = 1, 5
            DO 150 I = 1, NSIM
               R(I,J,L) = R(I,J,L) 
     1                  - A(I,J,L,1)*R(I,J-1,1) - A(I,J,L,2)*R(I,J-1,2) 
     2                  - A(I,J,L,3)*R(I,J-1,3) - A(I,J,L,4)*R(I,J-1,4) 
     3                  - A(I,J,L,5)*R(I,J-1,5)
  150       CONTINUE
  160    CONTINUE
C
C     Normalize row by multiplication of L[inverse](i,i) to get unity diagonal
C
         DO 190 I = 1, NSIM
C
C     Calculate L[inverse](i,i) * R
C
            R1       =  R(I,J,1) / B(I,J,1,1)
            R2       = (R(I,J,2)-B(I,J,2,1)*R1) / B(I,J,2,2)
            R3       = (R(I,J,3)-B(I,J,3,1)*R1-B(I,J,3,2)*R2)/B(I,J,3,3)
            R4       = (R(I,J,4)-B(I,J,4,1)*R1-B(I,J,4,2)*R2 
     1                          -B(I,J,4,3)*R3)              /B(I,J,4,4)
            R(I,J,5) = (R(I,J,5)-B(I,J,5,1)*R1-B(I,J,5,2)*R2 
     1                          -B(I,J,5,3)*R3-B(I,J,5,4)*R4)/B(I,J,5,5)
            R(I,J,4) =  R4-B(I,J,4,5)*R(I,J,5)
            R(I,J,3) =  R3-B(I,J,3,5)*R(I,J,5)-B(I,J,3,4)*R(I,J,4)
            R(I,J,2) =  R2-B(I,J,2,5)*R(I,J,5)-B(I,J,2,4)*R(I,J,4)
     1                    -B(I,J,2,3)*R(I,J,3)
            R(I,J,1) =  R1-B(I,J,1,5)*R(I,J,5)-B(I,J,1,4)*R(I,J,4)
     1                    -B(I,J,1,3)*R(I,J,3)-B(I,J,1,2)*R(I,J,2)
  190    CONTINUE
  200 CONTINUE
C
C     Backward Substitution Sweep - Invert U
C
      DO 400 J = NPTS-1, 1, -1
         DO 390 L = 1, 5
            DO 380 I = 1, NSIM
               R(I,J,L) = R(I,J,L)
     1                  - C(I,J,L,1)*R(I,J+1,1) 
     2                  - C(I,J,L,2)*R(I,J+1,2)
     3                  - C(I,J,L,3)*R(I,J+1,3) 
     4                  - C(I,J,L,4)*R(I,J+1,4)
     5                  - C(I,J,L,5)*R(I,J+1,5)
  380       CONTINUE
  390    CONTINUE
  400 CONTINUE
C
C     Finished
C
      RETURN
      END






      SUBROUTINE VFBTR7 (NPTS, NSIM, NRANK, A, B, C, R) 
C
C     Routine to perform the Forward/Backward substitution sweeps on 
C     a 7x7 block tridiagonal system where A, B and C contain
C     the LU decomposition of the block tridiagonal system.
C
C     VLUTR5 calculates the LU decomposition of the block tridiagonal system
C     VFBTR5 performs the forward/backward sweeps to complete the inversion
C
C     NPTS:      The number of 
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     NRANK:     The actual dimension of the storage
C     A:         The supra-diagonal block 7x7 of L
C     B:         The lu decomposition of the diagonal block 7x7 of L
C     C:         The super-diagonal block 7x7 of U
C     R:         The right hand side
C
C     A, B and C are unmodified by this routine
C     R contains the solution vector at return
C
C     The storage of the system is such that vectorization occurs over
C     the first index to provide a unit stride for maximum speed.
C
      DIMENSION A (NSIM,NPTS,NRANK,NRANK),
     1          B (NSIM,NPTS,NRANK,NRANK),
     2          C (NSIM,NPTS,NRANK,NRANK),
     3          R (NSIM,NPTS,NRANK)
C
C     Forward Substitution Loop 
C
C     J = 1
C
      J = 1
      DO 10 I = 1, NSIM
C
C     Calculate L[inverse](i,i) * R
C
         R1       =  R(I,J,1) / B(I,J,1,1)
         R2       = (R(I,J,2)-B(I,J,2,1)*R1) / B(I,J,2,2)
         R3       = (R(I,J,3)-B(I,J,3,1)*R1-B(I,J,3,2)*R2)/B(I,J,3,3)
         R4       = (R(I,J,4)-B(I,J,4,1)*R1-B(I,J,4,2)*R2 
     1                       -B(I,J,4,3)*R3)              /B(I,J,4,4)
         R5       = (R(I,J,5)-B(I,J,5,1)*R1-B(I,J,5,2)*R2 
     1                       -B(I,J,5,3)*R3-B(I,J,5,4)*R4)/B(I,J,5,5)
         R6       = (R(I,J,6)-B(I,J,6,1)*R1-B(I,J,6,2)*R2 
     1                       -B(I,J,6,3)*R3-B(I,J,6,4)*R4
     2                       -B(I,J,6,5)*R5              )/B(I,J,6,6)
         R(I,J,7) = (R(I,J,7)-B(I,J,7,1)*R1-B(I,J,7,2)*R2 
     1                       -B(I,J,7,3)*R3-B(I,J,7,4)*R4
     2                       -B(I,J,7,5)*R5-B(I,J,7,6)*R6)/B(I,J,7,7)
         R(I,J,6) =  R6-B(I,J,6,7)*R(I,J,7)
         R(I,J,5) =  R5-B(I,J,5,7)*R(I,J,7)-B(I,J,5,6)*R(I,J,6)
         R(I,J,4) =  R4-B(I,J,4,7)*R(I,J,7)-B(I,J,4,6)*R(I,J,6)
     1                 -B(I,J,4,5)*R(I,J,5)
         R(I,J,3) =  R3-B(I,J,3,7)*R(I,J,7)-B(I,J,3,6)*R(I,J,6)
     1                 -B(I,J,3,5)*R(I,J,5)-B(I,J,3,4)*R(I,J,4)
         R(I,J,2) =  R2-B(I,J,2,7)*R(I,J,7)-B(I,J,2,6)*R(I,J,6)
     1                 -B(I,J,2,5)*R(I,J,5)-B(I,J,2,4)*R(I,J,4)
     2                 -B(I,J,2,3)*R(I,J,3)
         R(I,J,1) =  R1-B(I,J,1,7)*R(I,J,7)-B(I,J,1,6)*R(I,J,6)
     1                 -B(I,J,1,5)*R(I,J,5)-B(I,J,1,4)*R(I,J,4)
     2                 -B(I,J,1,3)*R(I,J,3)-B(I,J,1,2)*R(I,J,2)
   10 CONTINUE
C
C     Loop over the interior of the inversion - 2 to NPTS
C
      DO 200 J = 2, NPTS
C
C     R = R - L(i,i-1)*R(j-1)
C
         DO 160 l = 1, 7
            DO 150 I = 1, NSIM
               R(I,J,L) = R(I,J,L) 
     1                  - A(I,J,L,1)*R(I,J-1,1) - A(I,J,L,2)*R(I,J-1,2) 
     2                  - A(I,J,L,3)*R(I,J-1,3) - A(I,J,L,4)*R(I,J-1,4) 
     3                  - A(I,J,L,5)*R(I,J-1,5) - A(I,J,L,6)*R(I,J-1,6) 
     4                  - A(I,J,L,7)*R(I,J-1,7)
  150       CONTINUE
  160    CONTINUE
C
C     Normalize row by multiplication of L[inverse](i,i) to get unity diagonal
C
         DO 190 I = 1, NSIM
C
C     Calculate L[inverse](i,i) * R
C
            R1       =  R(I,J,1) / B(I,J,1,1)
            R2       = (R(I,J,2)-B(I,J,2,1)*R1) / B(I,J,2,2)
            R3       = (R(I,J,3)-B(I,J,3,1)*R1-B(I,J,3,2)*R2)/B(I,J,3,3)
            R4       = (R(I,J,4)-B(I,J,4,1)*R1-B(I,J,4,2)*R2 
     1                          -B(I,J,4,3)*R3)              /B(I,J,4,4)
            R5       = (R(I,J,5)-B(I,J,5,1)*R1-B(I,J,5,2)*R2 
     1                          -B(I,J,5,3)*R3-B(I,J,5,4)*R4)/B(I,J,5,5)
            R6       = (R(I,J,6)-B(I,J,6,1)*R1-B(I,J,6,2)*R2 
     1                          -B(I,J,6,3)*R3-B(I,J,6,4)*R4
     2                          -B(I,J,6,5)*R5              )/B(I,J,6,6)
            R(I,J,7) = (R(I,J,7)-B(I,J,7,1)*R1-B(I,J,7,2)*R2 
     1                          -B(I,J,7,3)*R3-B(I,J,7,4)*R4
     2                          -B(I,J,7,5)*R5-B(I,J,7,6)*R6)/B(I,J,7,7)
            R(I,J,6) =  R6-B(I,J,6,7)*R(I,J,7)
            R(I,J,5) =  R5-B(I,J,5,7)*R(I,J,7)-B(I,J,5,6)*R(I,J,6)
            R(I,J,4) =  R4-B(I,J,4,7)*R(I,J,7)-B(I,J,4,6)*R(I,J,6)
     1                    -B(I,J,4,5)*R(I,J,5)
            R(I,J,3) =  R3-B(I,J,3,7)*R(I,J,7)-B(I,J,3,6)*R(I,J,6)
     1                    -B(I,J,3,5)*R(I,J,5)-B(I,J,3,4)*R(I,J,4)
            R(I,J,2) =  R2-B(I,J,2,7)*R(I,J,7)-B(I,J,2,6)*R(I,J,6)
     1                    -B(I,J,2,5)*R(I,J,5)-B(I,J,2,4)*R(I,J,4)
     2                    -B(I,J,2,3)*R(I,J,3)
            R(I,J,1) =  R1-B(I,J,1,7)*R(I,J,7)-B(I,J,1,6)*R(I,J,6)
     1                    -B(I,J,1,5)*R(I,J,5)-B(I,J,1,4)*R(I,J,4)
     2                    -B(I,J,1,3)*R(I,J,3)-B(I,J,1,2)*R(I,J,2)
  190    CONTINUE
  200 CONTINUE
C
C     Backward Substitution Sweep - Invert U
C
      DO 400 J = NPTS-1, 1, -1
         DO 390 L = 1, 7
            DO 380 I = 1, NSIM
               R(I,J,L) = R(I,J,L)
     1                  - C(I,J,L,1)*R(I,J+1,1) 
     2                  - C(I,J,L,2)*R(I,J+1,2)
     3                  - C(I,J,L,3)*R(I,J+1,3) 
     4                  - C(I,J,L,4)*R(I,J+1,4)
     5                  - C(I,J,L,5)*R(I,J+1,5)
     6                  - C(I,J,L,6)*R(I,J+1,6)
     7                  - C(I,J,L,7)*R(I,J+1,7)
  380       CONTINUE
  390    CONTINUE
  400 CONTINUE
C
C     Finished
C
      RETURN
      END






      SUBROUTINE VFBT12 (NPTS, NSIM, NRANK, A, B, C, R) 
C
C     Routine to perform the Forward/Backward substitution sweeps on 
C     a 12x12 block tridiagonal system where A, B and C contain
C     the LU decomposition of the block tridiagonal system.
C
C     VLUTR5 calculates the LU decomposition of the block tridiagonal system
C     VFBTR5 performs the forward/backward sweeps to complete the inversion
C
C     NPTS:      The number of 
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     NRANK:     The actual dimension of the storage
C     A:         The supra-diagonal block 12x12 of L
C     B:         The lu decomposition of the diagonal block 12x12 of L
C     C:         The super-diagonal block 12x12 of U
C     R:         The right hand side
C
C     A, B and C are unmodified by this routine
C     R contains the solution vector at return
C
C     The storage of the system is such that vectorization occurs over
C     the first index to provide a unit stride for maximum speed.
C
      DIMENSION A (NSIM,NPTS,NRANK,NRANK),
     1          B (NSIM,NPTS,NRANK,NRANK),
     2          C (NSIM,NPTS,NRANK,NRANK),
     3          R (NSIM,NPTS,NRANK)
C
C     Forward Substitution Loop 
C
C     Loop over the number of points - 1 to NPTS
C
      DO 200 J = 1, NPTS
C
C     R = R - L(i,i-1)*R(j-1)
C
         IF (J .NE. 1) THEN
            DO 160 l = 1, 12
               DO 150 I = 1, NSIM
                  R(I,J,L) = R(I,J,L) 
     1              - A(I,J,L, 1)*R(I,J-1, 1) - A(I,J,L, 2)*R(I,J-1, 2) 
     2              - A(I,J,L, 3)*R(I,J-1, 3) - A(I,J,L, 4)*R(I,J-1, 4) 
     3              - A(I,J,L, 5)*R(I,J-1, 5) - A(I,J,L, 6)*R(I,J-1, 6) 
     4              - A(I,J,L, 7)*R(I,J-1, 7) - A(I,J,L, 8)*R(I,J-1, 8)
     5              - A(I,J,L, 9)*R(I,J-1, 9) - A(I,J,L,10)*R(I,J-1,10)
     6              - A(I,J,L,11)*R(I,J-1,11) - A(I,J,L,12)*R(I,J-1,12)
  150          CONTINUE
  160       CONTINUE
         ENDIF
C
C     Normalize row by multiplication of L[inverse](i,i) to get unity diagonal
C
         DO 190 I = 1, NSIM
C
C     Calculate L[inverse](i,i) * R
C
            R1       =  R(I,J,1) / B(I,J,1,1)
            R2       = (R(I,J,2)-B(I,J,2,1)*R1) / B(I,J,2,2)
            R3       = (R(I,J,3)-B(I,J,3,1)*R1-B(I,J,3,2)*R2)/B(I,J,3,3)
            R4       = (R(I,J,4)-B(I,J,4,1)*R1-B(I,J,4,2)*R2 
     1                          -B(I,J,4,3)*R3)              /B(I,J,4,4)
            R5       = (R(I,J,5)-B(I,J,5,1)*R1-B(I,J,5,2)*R2 
     1                          -B(I,J,5,3)*R3-B(I,J,5,4)*R4)/B(I,J,5,5)
            R6       = (R(I,J,6)-B(I,J,6,1)*R1-B(I,J,6,2)*R2 
     1                          -B(I,J,6,3)*R3-B(I,J,6,4)*R4
     2                          -B(I,J,6,5)*R5              )/B(I,J,6,6)
            R7       = (R(I,J,7)-B(I,J,7,1)*R1-B(I,J,7,2)*R2 
     1                          -B(I,J,7,3)*R3-B(I,J,7,4)*R4
     2                          -B(I,J,7,5)*R5-B(I,J,7,6)*R6)/B(I,J,7,7)
            R8       = (R(I,J,8)-B(I,J,8,1)*R1-B(I,J,8,2)*R2 
     1                          -B(I,J,8,3)*R3-B(I,J,8,4)*R4
     2                          -B(I,J,8,5)*R5-B(I,J,8,6)*R6
     3                          -B(I,J,8,7)*R7              )/B(I,J,8,8)
            R9       = (R(I,J,9)-B(I,J,9,1)*R1-B(I,J,9,2)*R2 
     1                          -B(I,J,9,3)*R3-B(I,J,9,4)*R4
     2                          -B(I,J,9,5)*R5-B(I,J,9,6)*R6
     3                          -B(I,J,9,7)*R7-B(I,J,9,8)*R8)/B(I,J,9,9)
            R10      = (R(I,J,10)-B(I,J,10,1)*R1-B(I,J,10,2)*R2 
     1                           -B(I,J,10,3)*R3-B(I,J,10,4)*R4
     2                           -B(I,J,10,5)*R5-B(I,J,10,6)*R6
     3                           -B(I,J,10,7)*R7-B(I,J,10,8)*R8
     4                           -B(I,J,10,9)*R9               )
     5                                                     /B(I,J,10,10)
            R11      = (R(I,J,11)-B(I,J,11,1)*R1-B(I,J,11,2)*R2 
     1                           -B(I,J,11,3)*R3-B(I,J,11,4)*R4
     2                           -B(I,J,11,5)*R5-B(I,J,11,6)*R6
     3                           -B(I,J,11,7)*R7-B(I,J,11,8)*R8
     4                           -B(I,J,11,9)*R9-B(I,J,11,10)*R10)
     5                                                     /B(I,J,11,11)
            R(I,J,12)= (R(I,J,12)-B(I,J,12,1)*R1-B(I,J,12,2)*R2 
     1                           -B(I,J,12,3)*R3-B(I,J,12,4)*R4
     2                           -B(I,J,12,5)*R5-B(I,J,12,6)*R6
     3                           -B(I,J,12,7)*R7-B(I,J,12,8)*R8
     4                           -B(I,J,12,9)*R9-B(I,J,12,10)*R10
     5                           -B(I,J,12,11)*R11)        /B(I,J,12,12)
            R(I,J,11)= R11-B(I,J,11,12)*R(I,J,12)
            R(I,J,10)= R10-B(I,J,10,12)*R(I,J,12)-B(I,J,10,11)*R(I,J,11) 
            R(I,J,9) = R9 -B(I,J, 9,12)*R(I,J,12)-B(I,J, 9,11)*R(I,J,11)
     1                    -B(I,J, 9,10)*R(I,J,10)
            R(I,J,8) = R8 -B(I,J, 8,12)*R(I,J,12)-B(I,J, 8,11)*R(I,J,11)
     1                    -B(I,J, 8,10)*R(I,J,10)-B(I,J, 8, 9)*R(I,J, 9)
            R(I,J,7) = R7 -B(I,J, 7,12)*R(I,J,12)-B(I,J, 7,11)*R(I,J,11)
     1                    -B(I,J, 7,10)*R(I,J,10)-B(I,J, 7, 9)*R(I,J, 9)
     2                    -B(I,J, 7, 8)*R(I,J, 8)
            R(I,J,6) = R6 -B(I,J, 6,12)*R(I,J,12)-B(I,J, 6,11)*R(I,J,11)
     1                    -B(I,J, 6,10)*R(I,J,10)-B(I,J, 6, 9)*R(I,J, 9)
     2                    -B(I,J, 6, 8)*R(I,J, 8)-B(I,J, 6, 7)*R(I,J, 7)
            R(I,J,5) = R5 -B(I,J, 5,12)*R(I,J,12)-B(I,J, 5,11)*R(I,J,11)
     1                    -B(I,J, 5,10)*R(I,J,10)-B(I,J, 5, 9)*R(I,J, 9)
     2                    -B(I,J, 5, 8)*R(I,J, 8)-B(I,J, 5, 7)*R(I,J, 7)
     3                    -B(I,J, 5, 6)*R(I,J, 6)
            R(I,J,4) = R4 -B(I,J, 4,12)*R(I,J,12)-B(I,J, 4,11)*R(I,J,11)
     1                    -B(I,J, 4,10)*R(I,J,10)-B(I,J, 4, 9)*R(I,J, 9)
     2                    -B(I,J, 4, 8)*R(I,J, 8)-B(I,J, 4, 7)*R(I,J, 7)
     3                    -B(I,J, 4, 6)*R(I,J, 6)-B(I,J, 4, 5)*R(I,J, 5)
            R(I,J,3) = R3 -B(I,J, 3,12)*R(I,J,12)-B(I,J, 3,11)*R(I,J,11)
     1                    -B(I,J, 3,10)*R(I,J,10)-B(I,J, 3, 9)*R(I,J, 9)
     2                    -B(I,J, 3, 8)*R(I,J, 8)-B(I,J, 3, 7)*R(I,J, 7)
     3                    -B(I,J, 3, 6)*R(I,J, 6)-B(I,J, 3, 5)*R(I,J, 5)
     4                    -B(I,J, 3, 4)*R(I,J, 4)
            R(I,J,2) = R2 -B(I,J, 2,12)*R(I,J,12)-B(I,J, 2,11)*R(I,J,11)
     1                    -B(I,J, 2,10)*R(I,J,10)-B(I,J, 2, 9)*R(I,J, 9)
     2                    -B(I,J, 2, 8)*R(I,J, 8)-B(I,J, 2, 7)*R(I,J, 7)
     3                    -B(I,J, 2, 6)*R(I,J, 6)-B(I,J, 2, 5)*R(I,J, 5)
     4                    -B(I,J, 2, 4)*R(I,J, 4)-B(I,J, 2, 3)*R(I,J, 3)
            R(I,J,1) = R1 -B(I,J, 1,12)*R(I,J,12)-B(I,J, 1,11)*R(I,J,11)
     1                    -B(I,J, 1,10)*R(I,J,10)-B(I,J, 1, 9)*R(I,J, 9)
     2                    -B(I,J, 1, 8)*R(I,J, 8)-B(I,J, 1, 7)*R(I,J, 7)
     3                    -B(I,J, 1, 6)*R(I,J, 6)-B(I,J, 1, 5)*R(I,J, 5)
     4                    -B(I,J, 1, 4)*R(I,J, 4)-B(I,J, 1, 3)*R(I,J, 3)
     5                    -B(I,J, 1, 2)*R(I,J, 2)
  190    CONTINUE
  200 CONTINUE
C
C     Backward Substitution Sweep - Invert U
C
      DO 400 J = NPTS-1, 1, -1
         DO 390 L = 1, 12
            DO 380 I = 1, NSIM
               R(I,J,L) = R(I,J,L)
     1                  - C(I,J,L, 1)*R(I,J+1, 1) 
     2                  - C(I,J,L, 2)*R(I,J+1, 2)
     3                  - C(I,J,L, 3)*R(I,J+1, 3) 
     4                  - C(I,J,L, 4)*R(I,J+1, 4)
     5                  - C(I,J,L, 5)*R(I,J+1, 5)
     6                  - C(I,J,L, 6)*R(I,J+1, 6)
     7                  - C(I,J,L, 7)*R(I,J+1, 7)
     8                  - C(I,J,L, 8)*R(I,J+1, 8)
     9                  - C(I,J,L, 9)*R(I,J+1, 9)
     A                  - C(I,J,L,10)*R(I,J+1,10)
     B                  - C(I,J,L,11)*R(I,J+1,11)
     C                  - C(I,J,L,12)*R(I,J+1,12)
  380       CONTINUE
  390    CONTINUE
  400 CONTINUE
C
C     Finished
C
      RETURN
      END






      SUBROUTINE VFBTRS (NPTS, NSIM, A, B, C, R) 
C
C     Routine to perform the Forward/Backward substitution sweeps on 
C     a scalar tridiagonal system where A, B and C contain
C     the LU decomposition of the block tridiagonal system.
C
C     VLUTRS calculates the LU decomposition of the scalar tridiagonal system
C     VFBTRS performs the forward/backward sweeps to complete the inversion
C
C     NPTS:      The number of 
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     A:         The supra-diagonal element of L
C     B:         The lu decomposition of the diagonal element of L
C     C:         The super-diagonal element of U
C     R:         The right hand side
C
C     A, B and C are unmodified by this routine
C     R contains the solution vector at return
C
C     The storage of the system is such that vectorization occurs over
C     the first index to provide a unit stride for maximum speed.
C
      DIMENSION A (NSIM,NPTS),
     1          B (NSIM,NPTS),
     2          C (NSIM,NPTS),
     3          R (NSIM,NPTS)
C
C     Forward Substitution Loop 
C
C     J = 1
C     Calculate L[inverse](i,i) * R
C
      J = 1
      DO 10 I = 1, NSIM
         R(I,J) = R(I,J) / B(I,J)
   10 CONTINUE
C
C     Loop over the interior of the inversion - 2 to NPTS
C
      DO 200 J = 2, NPTS
C
C     R = (R - L(i,i-1)*R(j-1)) / L(i,i)
C
         DO 150 I = 1, NSIM
            R(I,J) = (R(I,J) - A(I,J)*R(I,J-1)) / B(I,J)
  150    CONTINUE
  200 CONTINUE
C
C     Backward Substitution Sweep - Invert U
C
      DO 400 J = NPTS-1, 1, -1
         DO 350 I = 1, NSIM
            R(I,J) = R(I,J) - C(I,J)*R(I,J+1) 
  350    CONTINUE
  400 CONTINUE
C
C     Finished
C
      RETURN
      END






      SUBROUTINE VLUTR (NPTS, NSIM, IRANK, A, B, C)
C
C     Routine to do a vectorized LU decomposition of a 
C     block tridiagonal system of equations.
C     The recursive nature of the LU decomposition leads
C     to a form which is not vectorizable.  Therefore, the
C     approach of this routine is to not attempt vectorization
C     of the LU, but to accomplish the vectorization through
C     simulataneously doing NSIM LU decompositions.
C     This form is with the identity on the diagonal of U.
C
C     NPTS:      The number of elements in the tridiagonal system
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     IRANK:     The rank of the system to invert
C     A:         The supra-diagonal block 7x7 term
C     B:         The       diagonal block 7x7 term
C     C:         The super-diagonal block 7x7 term
C
C     On Exit:
C     A:         L(j,j-1) The supra-diagonal term of L 
C     B:         The non-pivoted lu decomposition of L(i,i)
C     C:         U(j,j+1) The super-diagonal term of U 
C
      include '../header/common.h'
C
      DIMENSION A (NSIM,NPTS,NRANK,NRANK), 
     1          B (NSIM,NPTS,NRANK,NRANK), 
     2          C (NSIM,NPTS,NRANK,NRANK)
C
      IF (IRANK .EQ. 5) THEN
         CALL VLUTR5 (NPTS, NSIM, A, B, C)
      ELSE IF (IRANK .EQ. 7) THEN
         CALL VLUTR7 (NPTS, NSIM, NRANK, A, B, C)
      ELSE IF (IRANK .EQ. 12) THEN
         CALL VLUT12 (NPTS, NSIM, NRANK, A, B, C)
      ELSE 
         WRITE (IOUT,1000) IRANK, IRANK
         STOP
      ENDIF
C
C     Finished forming LU
C
      RETURN
 1000 FORMAT (' ','Block solver required for a ',I3,' x ',I3,' block.',
     1       /' ','This solver does not exist.  Run Aborting.')
      END






      SUBROUTINE VLUTR5 (NPTS, NSIM, A, B, C)
C
C     Routine to do a vectorized LU decomposition of a 
C     5x5 block tridiagonal system of equations.
C     The recursive nature of the LU decomposition leads
C     to a form which is not vectorizable.  Therefore, the
C     approach of this routine is to not attempt vectorization
C     of the LU, but to accomplish the vectorization through
C     simulataneously doing NSIM LU decompositions.
C     This form is with the identity on the diagonal of U.
C
C     NPTS:      The number of elements in the tridiagonal system
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     A:         The supra-diagonal block 5x5 term
C     B:         The       diagonal block 5x5 term
C     C:         The super-diagonal block 5x5 term
C
C     On Exit:
C     A:         L(j,j-1) The supra-diagonal term of L 
C     B:         The non-pivoted lu decomposition of L(i,i)
C                (this is needed both to form the LU and to invert the system)
C
C                One 5x5 block of B is of the form:
C
C                         l(1,1) u(1,2) u(1,3) u(1,4) u(1,5)
C                         l(2,1) l(2,2) u(2,3) u(2,4) u(2,5)
C                B[i,i] = l(3,1) l(3,2) l(3,3) u(3,4) u(3,5)
C                         l(4,1) l(4,2) l(4,3) l(4,4) u(4,5)
C                         l(5,1) l(5,2) l(5,3) l(5,4) l(5,5)
C
C     C:         U(j,j+1) The super-diagonal term of U 
C
C
C     L(i,i-1) = A(i) therefore no operations need to be performed on A.
C     L(i,i)   = B(i) - L(i,i-1)*U(i-1,i)
C     U(i,i+1) = L[inverse](i,i)*C(i)
C
      DIMENSION A (NSIM,NPTS,5,5), 
     1          B (NSIM,NPTS,5,5), 
     2          C (NSIM,NPTS,5,5)
C
C     Form L(1,1) and U(1,2)
C
      J = 1
C
C     Form the lu decomposition of L(1,1)
C
      DO 100 I = 1, NSIM
C
C     Row 1
C
         B(I,J,1,2) = B(I,J,1,2) / B(I,J,1,1)
         B(I,J,1,3) = B(I,J,1,3) / B(I,J,1,1)
         B(I,J,1,4) = B(I,J,1,4) / B(I,J,1,1)
         B(I,J,1,5) = B(I,J,1,5) / B(I,J,1,1)
C
C     Row 2
C
         B(I,J,2,2) =  B(I,J,2,2) - B(I,J,2,1) * B(I,J,1,2)
         B(I,J,2,3) = (B(I,J,2,3) - B(I,J,2,1)*B(I,J,1,3)) / B(I,J,2,2)
         B(I,J,2,4) = (B(I,J,2,4) - B(I,J,2,1)*B(I,J,1,4)) / B(I,J,2,2)
         B(I,J,2,5) = (B(I,J,2,5) - B(I,J,2,1)*B(I,J,1,5)) / B(I,J,2,2)
C
C     Row 3
C
         B(I,J,3,2) =  B(I,J,3,2) - B(I,J,3,1) * B(I,J,1,2)
         B(I,J,3,3) =  B(I,J,3,3) - B(I,J,3,1) * B(I,J,1,3) 
     1                            - B(I,J,3,2) * B(I,J,2,3)
         B(I,J,3,4) = (B(I,J,3,4) - B(I,J,3,1)*B(I,J,1,4) 
     1                            - B(I,J,3,2)*B(I,J,2,4)) / B(I,J,3,3)
         B(I,J,3,5) = (B(I,J,3,5) - B(I,J,3,1)*B(I,J,1,5) 
     1                            - B(I,J,3,2)*B(I,J,2,5)) / B(I,J,3,3)
C
C     Row 4
C
         B(I,J,4,2) =  B(I,J,4,2) - B(I,J,4,1) * B(I,J,1,2)
         B(I,J,4,3) =  B(I,J,4,3) - B(I,J,4,1) * B(I,J,1,3) 
     1                            - B(I,J,4,2) * B(I,J,2,3)
         B(I,J,4,4) =  B(I,J,4,4) - B(I,J,4,1) * B(I,J,1,4) 
     1                            - B(I,J,4,2) * B(I,J,2,4)
     2                            - B(I,J,4,3) * B(I,J,3,4)
         B(I,J,4,5) = (B(I,J,4,5) - B(I,J,4,1)*B(I,J,1,5) 
     1                            - B(I,J,4,2)*B(I,J,2,5)
     2                            - B(I,J,4,3)*B(I,J,3,5)) / B(I,J,4,4)
C
C     Row 5
C
         B(I,J,5,2) =  B(I,J,5,2) - B(I,J,5,1) * B(I,J,1,2)
         B(I,J,5,3) =  B(I,J,5,3) - B(I,J,5,1) * B(I,J,1,3) 
     1                            - B(I,J,5,2) * B(I,J,2,3)
         B(I,J,5,4) =  B(I,J,5,4) - B(I,J,5,1) * B(I,J,1,4) 
     1                            - B(I,J,5,2) * B(I,J,2,4)
     2                            - B(I,J,5,3) * B(I,J,3,4)
         B(I,J,5,5) =  B(I,J,5,5) - B(I,J,5,1) * B(I,J,1,5) 
     1                            - B(I,J,5,2) * B(I,J,2,5)
     2                            - B(I,J,5,3) * B(I,J,3,5)
     3                            - B(I,J,5,4) * B(I,J,4,5)
 100  CONTINUE
C
C     Form U(1,2) = L[inverse](1,1) * C(1)
C
      DO 200 L = 1, 5
         DO 180 I = 1, NSIM
C
C     Forward sweep of lu decomposition to form L[inverse]
C
            C(I,J,1,L) = C(I,J,1,L) / B(I,J,1,1)
            C(I,J,2,L) = ( C(I,J,2,L) - B(I,J,2,1)*C(I,J,1,L) ) 
     1                   / B(I,J,2,2)
            C(I,J,3,L) = ( C(I,J,3,L) - B(I,J,3,1)*C(I,J,1,L) 
     1                                - B(I,J,3,2)*C(I,J,2,L) ) 
     2                   / B(I,J,3,3)
            C(I,J,4,L) = ( C(I,J,4,L) - B(I,J,4,1)*C(I,J,1,L) 
     1                                - B(I,J,4,2)*C(I,J,2,L)  
     2                                - B(I,J,4,3)*C(I,J,3,L) ) 
     3                   / B(I,J,4,4)
            C(I,J,5,L) = ( C(I,J,5,L) - B(I,J,5,1)*C(I,J,1,L) 
     1                                - B(I,J,5,2)*C(I,J,2,L) 
     2                                - B(I,J,5,3)*C(I,J,3,L)  
     3                                - B(I,J,5,4)*C(I,J,4,L) )
     4                   / B(I,J,5,5)
C
C     Backward Sweep of lu
C
            C(I,J,4,L) = C(I,J,4,L) - B(I,J,4,5)*C(I,J,5,L)
            C(I,J,3,L) = C(I,J,3,L) - B(I,J,3,4)*C(I,J,4,L)
     1                              - B(I,J,3,5)*C(I,J,5,L)
            C(I,J,2,L) = C(I,J,2,L) - B(I,J,2,3)*C(I,J,3,L)
     1                              - B(I,J,2,4)*C(I,J,4,L)
     2                              - B(I,J,2,5)*C(I,J,5,L)
            C(I,J,1,L) = C(I,J,1,L) - B(I,J,1,2)*C(I,J,2,L)
     1                              - B(I,J,1,3)*C(I,J,3,L)
     2                              - B(I,J,1,4)*C(I,J,4,L)
     3                              - B(I,J,1,5)*C(I,J,5,L)
 180     CONTINUE
 200  CONTINUE
C
C     Loop over J = 2 to NPTS forming the LU terms
C
      DO 600 J = 2, NPTS
C
C     Form L(i,i) = B(i) - L(i,i-1)*U(i-1,i)
C
         DO 300 L = 1, 5
            DO 280 M = 1, 5
               DO 260 I = 1, NSIM
                  B(I,J,L,M) = B(I,J,L,M) - A(I,J,L,1)*C(I,J-1,1,M)
     1                                    - A(I,J,L,2)*C(I,J-1,2,M)
     2                                    - A(I,J,L,3)*C(I,J-1,3,M)
     3                                    - A(I,J,L,4)*C(I,J-1,4,M)
     4                                    - A(I,J,L,5)*C(I,J-1,5,M)
 260           CONTINUE
 280        CONTINUE
 300     CONTINUE
C
C     Form lu of L(i,i)
C
         DO 400 I = 1, NSIM
C
C     Row 1
C
            B(I,J,1,2) = B(I,J,1,2) / B(I,J,1,1)
            B(I,J,1,3) = B(I,J,1,3) / B(I,J,1,1)
            B(I,J,1,4) = B(I,J,1,4) / B(I,J,1,1)
            B(I,J,1,5) = B(I,J,1,5) / B(I,J,1,1)
C
C     Row 2
C
            B(I,J,2,2) =  B(I,J,2,2) - B(I,J,2,1) * B(I,J,1,2)
            B(I,J,2,3) = (B(I,J,2,3) - B(I,J,2,1)*B(I,J,1,3))/B(I,J,2,2)
            B(I,J,2,4) = (B(I,J,2,4) - B(I,J,2,1)*B(I,J,1,4))/B(I,J,2,2)
            B(I,J,2,5) = (B(I,J,2,5) - B(I,J,2,1)*B(I,J,1,5))/B(I,J,2,2)
C
C     Row 3
C
            B(I,J,3,2) =  B(I,J,3,2) - B(I,J,3,1) * B(I,J,1,2)
            B(I,J,3,3) =  B(I,J,3,3) - B(I,J,3,1) * B(I,J,1,3) 
     1                               - B(I,J,3,2) * B(I,J,2,3)
            B(I,J,3,4) = (B(I,J,3,4) - B(I,J,3,1)*B(I,J,1,4) 
     1                               - B(I,J,3,2)*B(I,J,2,4))/B(I,J,3,3)
            B(I,J,3,5) = (B(I,J,3,5) - B(I,J,3,1)*B(I,J,1,5) 
     1                               - B(I,J,3,2)*B(I,J,2,5))/B(I,J,3,3)
C
C     Row 4
C
            B(I,J,4,2) =  B(I,J,4,2) - B(I,J,4,1) * B(I,J,1,2)
            B(I,J,4,3) =  B(I,J,4,3) - B(I,J,4,1) * B(I,J,1,3) 
     1                               - B(I,J,4,2) * B(I,J,2,3)
            B(I,J,4,4) =  B(I,J,4,4) - B(I,J,4,1) * B(I,J,1,4) 
     1                               - B(I,J,4,2) * B(I,J,2,4)
     2                               - B(I,J,4,3) * B(I,J,3,4)
            B(I,J,4,5) = (B(I,J,4,5) - B(I,J,4,1)*B(I,J,1,5) 
     1                               - B(I,J,4,2)*B(I,J,2,5)
     2                               - B(I,J,4,3)*B(I,J,3,5))/B(I,J,4,4)
C
C     Row 5
C
            B(I,J,5,2) =  B(I,J,5,2) - B(I,J,5,1) * B(I,J,1,2)
            B(I,J,5,3) =  B(I,J,5,3) - B(I,J,5,1) * B(I,J,1,3) 
     1                               - B(I,J,5,2) * B(I,J,2,3)
            B(I,J,5,4) =  B(I,J,5,4) - B(I,J,5,1) * B(I,J,1,4) 
     1                               - B(I,J,5,2) * B(I,J,2,4)
     2                               - B(I,J,5,3) * B(I,J,3,4)
            B(I,J,5,5) =  B(I,J,5,5) - B(I,J,5,1) * B(I,J,1,5) 
     1                               - B(I,J,5,2) * B(I,J,2,5)
     2                               - B(I,J,5,3) * B(I,J,3,5)
     3                               - B(I,J,5,4) * B(I,J,4,5)
 400     CONTINUE
C
C     Form U(i,i+1) = L[inverse](i,i) * C(i) for all except last term
C
         IF (J .NE. NPTS) THEN
            DO 500 L = 1, 5
               DO 480 I = 1, NSIM
C
C     Forward sweep of lu decomposition to form L[inverse]
C
                  C(I,J,1,L) = C(I,J,1,L) / B(I,J,1,1)
                  C(I,J,2,L) = ( C(I,J,2,L) - B(I,J,2,1)*C(I,J,1,L) ) 
     1                         / B(I,J,2,2)
                  C(I,J,3,L) = ( C(I,J,3,L) - B(I,J,3,1)*C(I,J,1,L) 
     1                                      - B(I,J,3,2)*C(I,J,2,L) ) 
     2                         / B(I,J,3,3)
                  C(I,J,4,L) = ( C(I,J,4,L) - B(I,J,4,1)*C(I,J,1,L) 
     1                                      - B(I,J,4,2)*C(I,J,2,L)  
     2                                      - B(I,J,4,3)*C(I,J,3,L) ) 
     3                         / B(I,J,4,4)
                  C(I,J,5,L) = ( C(I,J,5,L) - B(I,J,5,1)*C(I,J,1,L) 
     1                                      - B(I,J,5,2)*C(I,J,2,L) 
     2                                      - B(I,J,5,3)*C(I,J,3,L)  
     3                                      - B(I,J,5,4)*C(I,J,4,L) )
     4                         / B(I,J,5,5)
C
C     Backward Sweep of lu
C
                  C(I,J,4,L) = C(I,J,4,L) - B(I,J,4,5)*C(I,J,5,L)
                  C(I,J,3,L) = C(I,J,3,L) - B(I,J,3,4)*C(I,J,4,L)
     1                                    - B(I,J,3,5)*C(I,J,5,L)
                  C(I,J,2,L) = C(I,J,2,L) - B(I,J,2,3)*C(I,J,3,L)
     1                                    - B(I,J,2,4)*C(I,J,4,L)
     2                                    - B(I,J,2,5)*C(I,J,5,L)
                  C(I,J,1,L) = C(I,J,1,L) - B(I,J,1,2)*C(I,J,2,L)
     1                                    - B(I,J,1,3)*C(I,J,3,L)
     2                                    - B(I,J,1,4)*C(I,J,4,L)
     3                                    - B(I,J,1,5)*C(I,J,5,L)
 480           CONTINUE
 500        CONTINUE
         ENDIF
 600  CONTINUE
C
C     Finished forming LU
C
      RETURN
      END






      SUBROUTINE VLUTR7 (NPTS, NSIM, NRANK, A, B, C)
C
C     Routine to do a vectorized LU decomposition of a 
C     7x7 block tridiagonal system of equations.
C     The recursive nature of the LU decomposition leads
C     to a form which is not vectorizable.  Therefore, the
C     approach of this routine is to not attempt vectorization
C     of the LU, but to accomplish the vectorization through
C     simulataneously doing NSIM LU decompositions.
C     This form is with the identity on the diagonal of U.
C
C     NPTS:      The number of elements in the tridiagonal system
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     NRANK:     The actual dimension of the storage of the tri-diagonal system
C     A:         The supra-diagonal block 7x7 term
C     B:         The       diagonal block 7x7 term
C     C:         The super-diagonal block 7x7 term
C
C     On Exit:
C     A:         L(j,j-1) The supra-diagonal term of L 
C     B:         The non-pivoted lu decomposition of L(i,i)
C                (this is needed both to form the LU and to invert the system)
C
C                One 7x7 block of B is of the form:
C
C                         l(1,1) u(1,2) u(1,3) u(1,4) u(1,5) u(1,6) u(1,7)
C                         l(2,1) l(2,2) u(2,3) u(2,4) u(2,5) u(2,6) u(2,7)
C                         l(3,1) l(3,2) l(3,3) u(3,4) u(3,5) u(3,6) u(3,7)
C                B[i,i] = l(4,1) l(4,2) l(4,3) l(4,4) u(4,5) u(4,6) u(4,7)
C                         l(5,1) l(5,2) l(5,3) l(5,4) l(5,5) u(5,6) u(5,7)
C                         l(6,1) l(6,2) l(6,3) l(6,4) l(6,5) l(6,6) u(6,7)
C                         l(7,1) l(7,2) l(7,3) l(7,4) l(7,5) l(7,6) l(7,7)
C
C     C:         U(j,j+1) The super-diagonal term of U 
C
C
C     L(i,i-1) = A(i) therefore no operations need to be performed on A.
C     L(i,i)   = B(i) - L(i,i-1)*U(i-1,i)
C     U(i,i+1) = L[inverse](i,i)*C(i)
C
      DIMENSION A (NSIM,NPTS,NRANK,NRANK), 
     1          B (NSIM,NPTS,NRANK,NRANK), 
     2          C (NSIM,NPTS,NRANK,NRANK)
C
C     Form L(1,1) and U(1,2)
C
      J = 1
C
C     Form the lu decomposition of L(1,1)
C
      DO 100 I = 1, NSIM
C
C     Row 1
C
         B(I,J,1,2) = B(I,J,1,2) / B(I,J,1,1)
         B(I,J,1,3) = B(I,J,1,3) / B(I,J,1,1)
         B(I,J,1,4) = B(I,J,1,4) / B(I,J,1,1)
         B(I,J,1,5) = B(I,J,1,5) / B(I,J,1,1)
         B(I,J,1,6) = B(I,J,1,6) / B(I,J,1,1)
         B(I,J,1,7) = B(I,J,1,7) / B(I,J,1,1)
C
C     Row 2
C
         B(I,J,2,2) =  B(I,J,2,2) - B(I,J,2,1) * B(I,J,1,2)
         B(I,J,2,3) = (B(I,J,2,3) - B(I,J,2,1)*B(I,J,1,3)) / B(I,J,2,2)
         B(I,J,2,4) = (B(I,J,2,4) - B(I,J,2,1)*B(I,J,1,4)) / B(I,J,2,2)
         B(I,J,2,5) = (B(I,J,2,5) - B(I,J,2,1)*B(I,J,1,5)) / B(I,J,2,2)
         B(I,J,2,6) = (B(I,J,2,6) - B(I,J,2,1)*B(I,J,1,6)) / B(I,J,2,2)
         B(I,J,2,7) = (B(I,J,2,7) - B(I,J,2,1)*B(I,J,1,7)) / B(I,J,2,2)
C
C     Row 3
C
         B(I,J,3,2) =  B(I,J,3,2) - B(I,J,3,1) * B(I,J,1,2)
         B(I,J,3,3) =  B(I,J,3,3) - B(I,J,3,1) * B(I,J,1,3) 
     1                            - B(I,J,3,2) * B(I,J,2,3)
         B(I,J,3,4) = (B(I,J,3,4) - B(I,J,3,1)*B(I,J,1,4) 
     1                            - B(I,J,3,2)*B(I,J,2,4)) / B(I,J,3,3)
         B(I,J,3,5) = (B(I,J,3,5) - B(I,J,3,1)*B(I,J,1,5) 
     1                            - B(I,J,3,2)*B(I,J,2,5)) / B(I,J,3,3)
         B(I,J,3,6) = (B(I,J,3,6) - B(I,J,3,1)*B(I,J,1,6) 
     1                            - B(I,J,3,2)*B(I,J,2,6)) / B(I,J,3,3)
         B(I,J,3,7) = (B(I,J,3,7) - B(I,J,3,1)*B(I,J,1,7) 
     1                            - B(I,J,3,2)*B(I,J,2,7)) / B(I,J,3,3)
C
C     Row 4
C
         B(I,J,4,2) =  B(I,J,4,2) - B(I,J,4,1) * B(I,J,1,2)
         B(I,J,4,3) =  B(I,J,4,3) - B(I,J,4,1) * B(I,J,1,3) 
     1                            - B(I,J,4,2) * B(I,J,2,3)
         B(I,J,4,4) =  B(I,J,4,4) - B(I,J,4,1) * B(I,J,1,4) 
     1                            - B(I,J,4,2) * B(I,J,2,4)
     2                            - B(I,J,4,3) * B(I,J,3,4)
         B(I,J,4,5) = (B(I,J,4,5) - B(I,J,4,1)*B(I,J,1,5) 
     1                            - B(I,J,4,2)*B(I,J,2,5)
     2                            - B(I,J,4,3)*B(I,J,3,5)) / B(I,J,4,4)
         B(I,J,4,6) = (B(I,J,4,6) - B(I,J,4,1)*B(I,J,1,6) 
     1                            - B(I,J,4,2)*B(I,J,2,6)
     2                            - B(I,J,4,3)*B(I,J,3,6)) / B(I,J,4,4)
         B(I,J,4,7) = (B(I,J,4,7) - B(I,J,4,1)*B(I,J,1,7) 
     1                            - B(I,J,4,2)*B(I,J,2,7)
     2                            - B(I,J,4,3)*B(I,J,3,7)) / B(I,J,4,4)
C
C     Row 5
C
         B(I,J,5,2) =  B(I,J,5,2) - B(I,J,5,1) * B(I,J,1,2)
         B(I,J,5,3) =  B(I,J,5,3) - B(I,J,5,1) * B(I,J,1,3) 
     1                            - B(I,J,5,2) * B(I,J,2,3)
         B(I,J,5,4) =  B(I,J,5,4) - B(I,J,5,1) * B(I,J,1,4) 
     1                            - B(I,J,5,2) * B(I,J,2,4)
     2                            - B(I,J,5,3) * B(I,J,3,4)
         B(I,J,5,5) =  B(I,J,5,5) - B(I,J,5,1) * B(I,J,1,5) 
     1                            - B(I,J,5,2) * B(I,J,2,5)
     2                            - B(I,J,5,3) * B(I,J,3,5)
     3                            - B(I,J,5,4) * B(I,J,4,5)
         B(I,J,5,6) = (B(I,J,5,6) - B(I,J,5,1)*B(I,J,1,6) 
     1                            - B(I,J,5,2)*B(I,J,2,6)
     2                            - B(I,J,5,3)*B(I,J,3,6)
     3                            - B(I,J,5,4)*B(I,J,4,6)) / B(I,J,5,5)
         B(I,J,5,7) = (B(I,J,5,7) - B(I,J,5,1)*B(I,J,1,7) 
     1                            - B(I,J,5,2)*B(I,J,2,7)
     2                            - B(I,J,5,3)*B(I,J,3,7)
     3                            - B(I,J,5,4)*B(I,J,4,7)) / B(I,J,5,5)
C
C     Row 6
C
         B(I,J,6,2) =  B(I,J,6,2) - B(I,J,6,1) * B(I,J,1,2)
         B(I,J,6,3) =  B(I,J,6,3) - B(I,J,6,1) * B(I,J,1,3) 
     1                            - B(I,J,6,2) * B(I,J,2,3)
         B(I,J,6,4) =  B(I,J,6,4) - B(I,J,6,1) * B(I,J,1,4) 
     1                            - B(I,J,6,2) * B(I,J,2,4)
     2                            - B(I,J,6,3) * B(I,J,3,4)
         B(I,J,6,5) =  B(I,J,6,5) - B(I,J,6,1) * B(I,J,1,5) 
     1                            - B(I,J,6,2) * B(I,J,2,5)
     2                            - B(I,J,6,3) * B(I,J,3,5)
     3                            - B(I,J,6,4) * B(I,J,4,5)
         B(I,J,6,6) =  B(I,J,6,6) - B(I,J,6,1) * B(I,J,1,6) 
     1                            - B(I,J,6,2) * B(I,J,2,6)
     2                            - B(I,J,6,3) * B(I,J,3,6)
     3                            - B(I,J,6,4) * B(I,J,4,6)
     4                            - B(I,J,6,5) * B(I,J,5,6)
         B(I,J,6,7) = (B(I,J,6,7) - B(I,J,6,1)*B(I,J,1,7) 
     1                            - B(I,J,6,2)*B(I,J,2,7)
     2                            - B(I,J,6,3)*B(I,J,3,7)
     3                            - B(I,J,6,4)*B(I,J,4,7)
     4                            - B(I,J,6,5)*B(I,J,5,7)) / B(I,J,6,6)
C
C     Row 7
C
         B(I,J,7,2) =  B(I,J,7,2) - B(I,J,7,1) * B(I,J,1,2)
         B(I,J,7,3) =  B(I,J,7,3) - B(I,J,7,1) * B(I,J,1,3) 
     1                            - B(I,J,7,2) * B(I,J,2,3)
         B(I,J,7,4) =  B(I,J,7,4) - B(I,J,7,1) * B(I,J,1,4) 
     1                            - B(I,J,7,2) * B(I,J,2,4)
     2                            - B(I,J,7,3) * B(I,J,3,4)
         B(I,J,7,5) =  B(I,J,7,5) - B(I,J,7,1) * B(I,J,1,5) 
     1                            - B(I,J,7,2) * B(I,J,2,5)
     2                            - B(I,J,7,3) * B(I,J,3,5)
     3                            - B(I,J,7,4) * B(I,J,4,5)
         B(I,J,7,6) =  B(I,J,7,6) - B(I,J,7,1) * B(I,J,1,6) 
     1                            - B(I,J,7,2) * B(I,J,2,6)
     2                            - B(I,J,7,3) * B(I,J,3,6)
     3                            - B(I,J,7,4) * B(I,J,4,6)
     4                            - B(I,J,7,5) * B(I,J,5,6)
         B(I,J,7,7) =  B(I,J,7,7) - B(I,J,7,1) * B(I,J,1,7) 
     1                            - B(I,J,7,2) * B(I,J,2,7)
     2                            - B(I,J,7,3) * B(I,J,3,7)
     3                            - B(I,J,7,4) * B(I,J,4,7)
     4                            - B(I,J,7,5) * B(I,J,5,7)
     5                            - B(I,J,7,6) * B(I,J,6,7)
 100  CONTINUE
C
C     Form U(1,2) = L[inverse](1,1) * C(1)
C
      DO 200 L = 1, 7
         DO 180 I = 1, NSIM
C
C     Forward sweep of lu decomposition to form L[inverse]
C
            C(I,J,1,L) = C(I,J,1,L) / B(I,J,1,1)
            C(I,J,2,L) = ( C(I,J,2,L) - B(I,J,2,1)*C(I,J,1,L) ) 
     1                   / B(I,J,2,2)
            C(I,J,3,L) = ( C(I,J,3,L) - B(I,J,3,1)*C(I,J,1,L) 
     1                                - B(I,J,3,2)*C(I,J,2,L) ) 
     2                   / B(I,J,3,3)
            C(I,J,4,L) = ( C(I,J,4,L) - B(I,J,4,1)*C(I,J,1,L) 
     1                                - B(I,J,4,2)*C(I,J,2,L)  
     2                                - B(I,J,4,3)*C(I,J,3,L) ) 
     3                   / B(I,J,4,4)
            C(I,J,5,L) = ( C(I,J,5,L) - B(I,J,5,1)*C(I,J,1,L) 
     1                                - B(I,J,5,2)*C(I,J,2,L) 
     2                                - B(I,J,5,3)*C(I,J,3,L)  
     3                                - B(I,J,5,4)*C(I,J,4,L) )
     4                   / B(I,J,5,5)
            C(I,J,6,L) = ( C(I,J,6,L) - B(I,J,6,1)*C(I,J,1,L) 
     1                                - B(I,J,6,2)*C(I,J,2,L) 
     2                                - B(I,J,6,3)*C(I,J,3,L)  
     3                                - B(I,J,6,4)*C(I,J,4,L) 
     4                                - B(I,J,6,5)*C(I,J,5,L) )
     5                   / B(I,J,6,6)
            C(I,J,7,L) = ( C(I,J,7,L) - B(I,J,7,1)*C(I,J,1,L) 
     1                                - B(I,J,7,2)*C(I,J,2,L) 
     2                                - B(I,J,7,3)*C(I,J,3,L)  
     3                                - B(I,J,7,4)*C(I,J,4,L) 
     4                                - B(I,J,7,5)*C(I,J,5,L)
     5                                - B(I,J,7,6)*C(I,J,6,L) )
     6                   / B(I,J,7,7)
C
C     Backward Sweep of lu
C
            C(I,J,6,L) = C(I,J,6,L) - B(I,J,6,7)*C(I,J,7,L)
            C(I,J,5,L) = C(I,J,5,L) - B(I,J,5,6)*C(I,J,6,L)
     1                              - B(I,J,5,7)*C(I,J,7,L)
            C(I,J,4,L) = C(I,J,4,L) - B(I,J,4,5)*C(I,J,5,L)
     1                              - B(I,J,4,6)*C(I,J,6,L)
     2                              - B(I,J,4,7)*C(I,J,7,L)
            C(I,J,3,L) = C(I,J,3,L) - B(I,J,3,4)*C(I,J,4,L)
     1                              - B(I,J,3,5)*C(I,J,5,L)
     2                              - B(I,J,3,6)*C(I,J,6,L)
     3                              - B(I,J,3,7)*C(I,J,7,L)
            C(I,J,2,L) = C(I,J,2,L) - B(I,J,2,3)*C(I,J,3,L)
     1                              - B(I,J,2,4)*C(I,J,4,L)
     2                              - B(I,J,2,5)*C(I,J,5,L)
     3                              - B(I,J,2,6)*C(I,J,6,L)
     4                              - B(I,J,2,7)*C(I,J,7,L)
            C(I,J,1,L) = C(I,J,1,L) - B(I,J,1,2)*C(I,J,2,L)
     1                              - B(I,J,1,3)*C(I,J,3,L)
     2                              - B(I,J,1,4)*C(I,J,4,L)
     3                              - B(I,J,1,5)*C(I,J,5,L)
     4                              - B(I,J,1,6)*C(I,J,6,L)
     5                              - B(I,J,1,7)*C(I,J,7,L)
 180     CONTINUE
 200  CONTINUE
C
C     Loop over J = 2 to NPTS forming the LU terms
C
      DO 600 J = 2, NPTS
C
C     Form L(i,i) = B(i) - L(i,i-1)*U(i-1,i)
C
         DO 300 L = 1, 7
            DO 280 M = 1, 7
               DO 260 I = 1, NSIM
                  B(I,J,L,M) = B(I,J,L,M) - A(I,J,L,1)*C(I,J-1,1,M)
     1                                    - A(I,J,L,2)*C(I,J-1,2,M)
     2                                    - A(I,J,L,3)*C(I,J-1,3,M)
     3                                    - A(I,J,L,4)*C(I,J-1,4,M)
     4                                    - A(I,J,L,5)*C(I,J-1,5,M)
     5                                    - A(I,J,L,6)*C(I,J-1,6,M)
     6                                    - A(I,J,L,7)*C(I,J-1,7,M)
 260           CONTINUE
 280        CONTINUE
 300     CONTINUE
C
C     Form lu of L(i,i)
C
         DO 400 I = 1, NSIM
C
C     Row 1
C
            B(I,J,1,2) = B(I,J,1,2) / B(I,J,1,1)
            B(I,J,1,3) = B(I,J,1,3) / B(I,J,1,1)
            B(I,J,1,4) = B(I,J,1,4) / B(I,J,1,1)
            B(I,J,1,5) = B(I,J,1,5) / B(I,J,1,1)
            B(I,J,1,6) = B(I,J,1,6) / B(I,J,1,1)
            B(I,J,1,7) = B(I,J,1,7) / B(I,J,1,1)
C
C     Row 2
C
            B(I,J,2,2) =  B(I,J,2,2) - B(I,J,2,1) * B(I,J,1,2)
            B(I,J,2,3) = (B(I,J,2,3) - B(I,J,2,1)*B(I,J,1,3))/B(I,J,2,2)
            B(I,J,2,4) = (B(I,J,2,4) - B(I,J,2,1)*B(I,J,1,4))/B(I,J,2,2)
            B(I,J,2,5) = (B(I,J,2,5) - B(I,J,2,1)*B(I,J,1,5))/B(I,J,2,2)
            B(I,J,2,6) = (B(I,J,2,6) - B(I,J,2,1)*B(I,J,1,6))/B(I,J,2,2)
            B(I,J,2,7) = (B(I,J,2,7) - B(I,J,2,1)*B(I,J,1,7))/B(I,J,2,2)
C
C     Row 3
C
            B(I,J,3,2) =  B(I,J,3,2) - B(I,J,3,1) * B(I,J,1,2)
            B(I,J,3,3) =  B(I,J,3,3) - B(I,J,3,1) * B(I,J,1,3) 
     1                               - B(I,J,3,2) * B(I,J,2,3)
            B(I,J,3,4) = (B(I,J,3,4) - B(I,J,3,1)*B(I,J,1,4) 
     1                               - B(I,J,3,2)*B(I,J,2,4))/B(I,J,3,3)
            B(I,J,3,5) = (B(I,J,3,5) - B(I,J,3,1)*B(I,J,1,5) 
     1                               - B(I,J,3,2)*B(I,J,2,5))/B(I,J,3,3)
            B(I,J,3,6) = (B(I,J,3,6) - B(I,J,3,1)*B(I,J,1,6) 
     1                               - B(I,J,3,2)*B(I,J,2,6))/B(I,J,3,3)
            B(I,J,3,7) = (B(I,J,3,7) - B(I,J,3,1)*B(I,J,1,7) 
     1                               - B(I,J,3,2)*B(I,J,2,7))/B(I,J,3,3)
C
C     Row 4
C
            B(I,J,4,2) =  B(I,J,4,2) - B(I,J,4,1) * B(I,J,1,2)
            B(I,J,4,3) =  B(I,J,4,3) - B(I,J,4,1) * B(I,J,1,3) 
     1                               - B(I,J,4,2) * B(I,J,2,3)
            B(I,J,4,4) =  B(I,J,4,4) - B(I,J,4,1) * B(I,J,1,4) 
     1                               - B(I,J,4,2) * B(I,J,2,4)
     2                               - B(I,J,4,3) * B(I,J,3,4)
            B(I,J,4,5) = (B(I,J,4,5) - B(I,J,4,1)*B(I,J,1,5) 
     1                               - B(I,J,4,2)*B(I,J,2,5)
     2                               - B(I,J,4,3)*B(I,J,3,5))/B(I,J,4,4)
            B(I,J,4,6) = (B(I,J,4,6) - B(I,J,4,1)*B(I,J,1,6) 
     1                               - B(I,J,4,2)*B(I,J,2,6)
     2                               - B(I,J,4,3)*B(I,J,3,6))/B(I,J,4,4)
            B(I,J,4,7) = (B(I,J,4,7) - B(I,J,4,1)*B(I,J,1,7) 
     1                               - B(I,J,4,2)*B(I,J,2,7)
     2                               - B(I,J,4,3)*B(I,J,3,7))/B(I,J,4,4)
C
C     Row 5
C
            B(I,J,5,2) =  B(I,J,5,2) - B(I,J,5,1) * B(I,J,1,2)
            B(I,J,5,3) =  B(I,J,5,3) - B(I,J,5,1) * B(I,J,1,3) 
     1                               - B(I,J,5,2) * B(I,J,2,3)
            B(I,J,5,4) =  B(I,J,5,4) - B(I,J,5,1) * B(I,J,1,4) 
     1                               - B(I,J,5,2) * B(I,J,2,4)
     2                               - B(I,J,5,3) * B(I,J,3,4)
            B(I,J,5,5) =  B(I,J,5,5) - B(I,J,5,1) * B(I,J,1,5) 
     1                               - B(I,J,5,2) * B(I,J,2,5)
     2                               - B(I,J,5,3) * B(I,J,3,5)
     3                               - B(I,J,5,4) * B(I,J,4,5)
            B(I,J,5,6) = (B(I,J,5,6) - B(I,J,5,1)*B(I,J,1,6) 
     1                               - B(I,J,5,2)*B(I,J,2,6)
     2                               - B(I,J,5,3)*B(I,J,3,6)
     3                               - B(I,J,5,4)*B(I,J,4,6))/B(I,J,5,5)
            B(I,J,5,7) = (B(I,J,5,7) - B(I,J,5,1)*B(I,J,1,7) 
     1                               - B(I,J,5,2)*B(I,J,2,7)
     2                               - B(I,J,5,3)*B(I,J,3,7)
     3                               - B(I,J,5,4)*B(I,J,4,7))/B(I,J,5,5)
C
C     Row 6
C
            B(I,J,6,2) =  B(I,J,6,2) - B(I,J,6,1) * B(I,J,1,2)
            B(I,J,6,3) =  B(I,J,6,3) - B(I,J,6,1) * B(I,J,1,3) 
     1                               - B(I,J,6,2) * B(I,J,2,3)
            B(I,J,6,4) =  B(I,J,6,4) - B(I,J,6,1) * B(I,J,1,4) 
     1                               - B(I,J,6,2) * B(I,J,2,4)
     2                               - B(I,J,6,3) * B(I,J,3,4)
            B(I,J,6,5) =  B(I,J,6,5) - B(I,J,6,1) * B(I,J,1,5) 
     1                               - B(I,J,6,2) * B(I,J,2,5)
     2                               - B(I,J,6,3) * B(I,J,3,5)
     3                               - B(I,J,6,4) * B(I,J,4,5)
            B(I,J,6,6) =  B(I,J,6,6) - B(I,J,6,1) * B(I,J,1,6) 
     1                               - B(I,J,6,2) * B(I,J,2,6)
     2                               - B(I,J,6,3) * B(I,J,3,6)
     3                               - B(I,J,6,4) * B(I,J,4,6)
     4                               - B(I,J,6,5) * B(I,J,5,6)
            B(I,J,6,7) = (B(I,J,6,7) - B(I,J,6,1)*B(I,J,1,7) 
     1                               - B(I,J,6,2)*B(I,J,2,7)
     2                               - B(I,J,6,3)*B(I,J,3,7)
     3                               - B(I,J,6,4)*B(I,J,4,7)
     4                               - B(I,J,6,5)*B(I,J,5,7))/B(I,J,6,6)
C
C     Row 7
C
            B(I,J,7,2) =  B(I,J,7,2) - B(I,J,7,1) * B(I,J,1,2)
            B(I,J,7,3) =  B(I,J,7,3) - B(I,J,7,1) * B(I,J,1,3) 
     1                               - B(I,J,7,2) * B(I,J,2,3)
            B(I,J,7,4) =  B(I,J,7,4) - B(I,J,7,1) * B(I,J,1,4) 
     1                               - B(I,J,7,2) * B(I,J,2,4)
     2                               - B(I,J,7,3) * B(I,J,3,4)
            B(I,J,7,5) =  B(I,J,7,5) - B(I,J,7,1) * B(I,J,1,5) 
     1                               - B(I,J,7,2) * B(I,J,2,5)
     2                               - B(I,J,7,3) * B(I,J,3,5)
     3                               - B(I,J,7,4) * B(I,J,4,5)
            B(I,J,7,6) =  B(I,J,7,6) - B(I,J,7,1) * B(I,J,1,6) 
     1                               - B(I,J,7,2) * B(I,J,2,6)
     2                               - B(I,J,7,3) * B(I,J,3,6)
     3                               - B(I,J,7,4) * B(I,J,4,6)
     4                               - B(I,J,7,5) * B(I,J,5,6)
            B(I,J,7,7) =  B(I,J,7,7) - B(I,J,7,1) * B(I,J,1,7) 
     1                               - B(I,J,7,2) * B(I,J,2,7)
     2                               - B(I,J,7,3) * B(I,J,3,7)
     3                               - B(I,J,7,4) * B(I,J,4,7)
     4                               - B(I,J,7,5) * B(I,J,5,7)
     5                               - B(I,J,7,6) * B(I,J,6,7)
 400     CONTINUE
C
C     Form U(i,i+1) = L[inverse](i,i) * C(i) for all except last term
C
         IF (J .NE. NPTS) THEN
            DO 500 L = 1, 7
               DO 480 I = 1, NSIM
C
C     Forward sweep of lu decomposition to form L[inverse]
C
                  C(I,J,1,L) = C(I,J,1,L) / B(I,J,1,1)
                  C(I,J,2,L) = ( C(I,J,2,L) - B(I,J,2,1)*C(I,J,1,L) ) 
     1                         / B(I,J,2,2)
                  C(I,J,3,L) = ( C(I,J,3,L) - B(I,J,3,1)*C(I,J,1,L) 
     1                                      - B(I,J,3,2)*C(I,J,2,L) ) 
     2                         / B(I,J,3,3)
                  C(I,J,4,L) = ( C(I,J,4,L) - B(I,J,4,1)*C(I,J,1,L) 
     1                                      - B(I,J,4,2)*C(I,J,2,L)  
     2                                      - B(I,J,4,3)*C(I,J,3,L) ) 
     3                         / B(I,J,4,4)
                  C(I,J,5,L) = ( C(I,J,5,L) - B(I,J,5,1)*C(I,J,1,L) 
     1                                      - B(I,J,5,2)*C(I,J,2,L) 
     2                                      - B(I,J,5,3)*C(I,J,3,L)  
     3                                      - B(I,J,5,4)*C(I,J,4,L) )
     4                         / B(I,J,5,5)
                  C(I,J,6,L) = ( C(I,J,6,L) - B(I,J,6,1)*C(I,J,1,L) 
     1                                      - B(I,J,6,2)*C(I,J,2,L) 
     2                                      - B(I,J,6,3)*C(I,J,3,L)  
     3                                      - B(I,J,6,4)*C(I,J,4,L)
     4                                      - B(I,J,6,5)*C(I,J,5,L) )
     5                         / B(I,J,6,6)
                  C(I,J,7,L) = ( C(I,J,7,L) - B(I,J,7,1)*C(I,J,1,L) 
     1                                      - B(I,J,7,2)*C(I,J,2,L) 
     2                                      - B(I,J,7,3)*C(I,J,3,L)  
     3                                      - B(I,J,7,4)*C(I,J,4,L)
     4                                      - B(I,J,7,5)*C(I,J,5,L)
     5                                      - B(I,J,7,6)*C(I,J,6,L) )
     6                         / B(I,J,7,7)
C
C     Backward Sweep of lu
C
                  C(I,J,6,L) = C(I,J,6,L) - B(I,J,6,7)*C(I,J,7,L)
                  C(I,J,5,L) = C(I,J,5,L) - B(I,J,5,6)*C(I,J,6,L)
     1                                    - B(I,J,5,7)*C(I,J,7,L)
                  C(I,J,4,L) = C(I,J,4,L) - B(I,J,4,5)*C(I,J,5,L)
     1                                    - B(I,J,4,6)*C(I,J,6,L)
     2                                    - B(I,J,4,7)*C(I,J,7,L)
                  C(I,J,3,L) = C(I,J,3,L) - B(I,J,3,4)*C(I,J,4,L)
     1                                    - B(I,J,3,5)*C(I,J,5,L)
     2                                    - B(I,J,3,6)*C(I,J,6,L)
     3                                    - B(I,J,3,7)*C(I,J,7,L)
                  C(I,J,2,L) = C(I,J,2,L) - B(I,J,2,3)*C(I,J,3,L)
     1                                    - B(I,J,2,4)*C(I,J,4,L)
     2                                    - B(I,J,2,5)*C(I,J,5,L)
     3                                    - B(I,J,2,6)*C(I,J,6,L)
     4                                    - B(I,J,2,7)*C(I,J,7,L)
                  C(I,J,1,L) = C(I,J,1,L) - B(I,J,1,2)*C(I,J,2,L)
     1                                    - B(I,J,1,3)*C(I,J,3,L)
     2                                    - B(I,J,1,4)*C(I,J,4,L)
     3                                    - B(I,J,1,5)*C(I,J,5,L)
     4                                    - B(I,J,1,6)*C(I,J,6,L)
     5                                    - B(I,J,1,7)*C(I,J,7,L)
 480           CONTINUE
 500        CONTINUE
         ENDIF
 600  CONTINUE
C
C     Finished forming LU
C
      RETURN
      END






      SUBROUTINE VLUT12 (NPTS, NSIM, NRANK, A, B, C)
C
C     Routine to do a vectorized LU decomposition of a 
C     12x12 block tridiagonal system of equations.
C     The recursive nature of the LU decomposition leads
C     to a form which is not vectorizable.  Therefore, the
C     approach of this routine is to not attempt vectorization
C     of the LU, but to accomplish the vectorization through
C     simulataneously doing NSIM LU decompositions.
C     This form is with the identity on the diagonal of U.
C
C     NPTS:      The number of elements in the tridiagonal system
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     NRANK:     The actual dimension of the storage of the tri-diagonal system
C     A:         The supra-diagonal block 12x12 term
C     B:         The       diagonal block 12x12 term
C     C:         The super-diagonal block 12x12 term
C
C     On Exit:
C     A:         L(j,j-1) The supra-diagonal term of L 
C     B:         The non-pivoted lu decomposition of L(i,i)
C                (this is needed both to form the LU and to invert the system)
C
C                One 12x12 block of B is of the form:
C
C                 l(1,1)  u(1,2)  u(1,3)  u(1,4)  u(1,5)  ... u(1,11)  u(1,12)
C                 l(2,1)  l(2,2)  u(2,3)  u(2,4)  u(2,5)  ... u(2,11)  u(2,12)
C                 l(3,1)  l(3,2)  l(3,3)  u(3,4)  u(3,5)  ... u(3,11)  u(3,12)
C        B[i,i] = l(4,1)  l(4,2)  l(4,3)  l(4,4)  u(4,5)  ... u(4,11)  u(4,12)
C                 l(5,1)  l(5,2)  l(5,3)  l(5,4)  l(5,5)  ... u(5,11)  u(5,12)
C                    .       .       .       .       .           .        .
C                    .       .       .       .       .           .        .
C                    .       .       .       .       .           .        .
C                 l(11,1) l(11,2) l(11,3) l(11,4) l(11,5) ... l(11,11) u(11,12)
C                 l(12,1) l(12,2) l(12,3) l(12,4) l(12,5) ... l(12,11) l(12,12)
C
C     C:         U(j,j+1) The super-diagonal term of U 
C
C
C     L(i,i-1) = A(i) therefore no operations need to be performed on A.
C     L(i,i)   = B(i) - L(i,i-1)*U(i-1,i)
C     U(i,i+1) = L[inverse](i,i)*C(i)
C
      DIMENSION A (NSIM,NPTS,NRANK,NRANK), 
     1          B (NSIM,NPTS,NRANK,NRANK), 
     2          C (NSIM,NPTS,NRANK,NRANK)
C
C     Loop over J = 1 to NPTS forming the LU terms
C
      DO 600 J = 1, NPTS
C
C     Form L(i,i) = B(i) - L(i,i-1)*U(i-1,i) for all except the J = 1
C
         IF (J .NE. 1) THEN
            DO 130 L = 1, 12
               DO 120 M = 1, 12
                  DO 110 I = 1, NSIM
                     B(I,J,L,M) = B(I,J,L,M) - A(I,J,L, 1)*C(I,J-1, 1,M)
     1                                       - A(I,J,L, 2)*C(I,J-1, 2,M)
     2                                       - A(I,J,L, 3)*C(I,J-1, 3,M)
     3                                       - A(I,J,L, 4)*C(I,J-1, 4,M)
     4                                       - A(I,J,L, 5)*C(I,J-1, 5,M)
     5                                       - A(I,J,L, 6)*C(I,J-1, 6,M)
     6                                       - A(I,J,L, 7)*C(I,J-1, 7,M)
     7                                       - A(I,J,L, 8)*C(I,J-1, 8,M)
     8                                       - A(I,J,L, 9)*C(I,J-1, 9,M)
     9                                       - A(I,J,L,10)*C(I,J-1,10,M)
     A                                       - A(I,J,L,11)*C(I,J-1,11,M)
     B                                       - A(I,J,L,12)*C(I,J-1,12,M)
 110              CONTINUE
 120           CONTINUE
 130        CONTINUE
         ENDIF
C
C     Form the lu decomposition of L(i,i)
C
         DO 200 I = 1, NSIM
C
C     Row 1
C
            B(I,J,1, 2) = B(I,J,1, 2) / B(I,J,1,1)
            B(I,J,1, 3) = B(I,J,1, 3) / B(I,J,1,1)
            B(I,J,1, 4) = B(I,J,1, 4) / B(I,J,1,1)
            B(I,J,1, 5) = B(I,J,1, 5) / B(I,J,1,1)
            B(I,J,1, 6) = B(I,J,1, 6) / B(I,J,1,1)
            B(I,J,1, 7) = B(I,J,1, 7) / B(I,J,1,1)
            B(I,J,1, 8) = B(I,J,1, 8) / B(I,J,1,1)
            B(I,J,1, 9) = B(I,J,1, 9) / B(I,J,1,1)
            B(I,J,1,10) = B(I,J,1,10) / B(I,J,1,1)
            B(I,J,1,11) = B(I,J,1,11) / B(I,J,1,1)
            B(I,J,1,12) = B(I,J,1,12) / B(I,J,1,1)
C
C     Row 2
C
            B(I,J,2, 2) =  B(I,J,2, 2) - B(I,J,2,1)*B(I,J,1, 2)
            B(I,J,2, 3) = (B(I,J,2, 3) - B(I,J,2,1)*B(I,J,1, 3))
     1                    /B(I,J,2,2)
            B(I,J,2, 4) = (B(I,J,2, 4) - B(I,J,2,1)*B(I,J,1, 4))
     1                    /B(I,J,2,2)
            B(I,J,2, 5) = (B(I,J,2, 5) - B(I,J,2,1)*B(I,J,1, 5))
     1                    /B(I,J,2,2)
            B(I,J,2, 6) = (B(I,J,2, 6) - B(I,J,2,1)*B(I,J,1, 6))
     1                    /B(I,J,2,2)
            B(I,J,2, 7) = (B(I,J,2, 7) - B(I,J,2,1)*B(I,J,1, 7))
     1                    /B(I,J,2,2)
            B(I,J,2, 8) = (B(I,J,2, 8) - B(I,J,2,1)*B(I,J,1, 8))
     1                    /B(I,J,2,2)
            B(I,J,2, 9) = (B(I,J,2, 9) - B(I,J,2,1)*B(I,J,1, 9))
     1                    /B(I,J,2,2)
            B(I,J,2,10) = (B(I,J,2,10) - B(I,J,2,1)*B(I,J,1,10))
     1                    /B(I,J,2,2)
            B(I,J,2,11) = (B(I,J,2,11) - B(I,J,2,1)*B(I,J,1,11))
     1                    /B(I,J,2,2)
            B(I,J,2,12) = (B(I,J,2,12) - B(I,J,2,1)*B(I,J,1,12))
     1                    /B(I,J,2,2)
C
C     Row 3
C
            B(I,J,3, 2) =  B(I,J,3, 2) - B(I,J,3,1)*B(I,J,1, 2)
            B(I,J,3, 3) =  B(I,J,3, 3) - B(I,J,3,1)*B(I,J,1, 3) 
     1                                 - B(I,J,3,2)*B(I,J,2, 3)
            B(I,J,3, 4) = (B(I,J,3, 4) - B(I,J,3,1)*B(I,J,1, 4) 
     1                                 - B(I,J,3,2)*B(I,J,2, 4))
     2                    /B(I,J,3,3)
            B(I,J,3, 5) = (B(I,J,3, 5) - B(I,J,3,1)*B(I,J,1, 5) 
     1                                 - B(I,J,3,2)*B(I,J,2, 5))
     2                    /B(I,J,3,3)
            B(I,J,3, 6) = (B(I,J,3, 6) - B(I,J,3,1)*B(I,J,1, 6) 
     1                                 - B(I,J,3,2)*B(I,J,2, 6))
     2                    /B(I,J,3,3)
            B(I,J,3, 7) = (B(I,J,3, 7) - B(I,J,3,1)*B(I,J,1, 7) 
     1                                 - B(I,J,3,2)*B(I,J,2, 7))
     2                    /B(I,J,3,3)
            B(I,J,3, 8) = (B(I,J,3, 8) - B(I,J,3,1)*B(I,J,1, 8) 
     1                                 - B(I,J,3,2)*B(I,J,2, 8))
     2                    /B(I,J,3,3)
            B(I,J,3, 9) = (B(I,J,3, 9) - B(I,J,3,1)*B(I,J,1, 9) 
     1                                 - B(I,J,3,2)*B(I,J,2, 9))
     2                    /B(I,J,3,3)
            B(I,J,3,10) = (B(I,J,3,10) - B(I,J,3,1)*B(I,J,1,10) 
     1                                 - B(I,J,3,2)*B(I,J,2,10))
     2                    /B(I,J,3,3)
            B(I,J,3,11) = (B(I,J,3,11) - B(I,J,3,1)*B(I,J,1,11) 
     1                                 - B(I,J,3,2)*B(I,J,2,11))
     2                    /B(I,J,3,3)
            B(I,J,3,12) = (B(I,J,3,12) - B(I,J,3,1)*B(I,J,1,12) 
     1                                 - B(I,J,3,2)*B(I,J,2,12))
     2                    /B(I,J,3,3)
C
C     Row 4
C
            B(I,J,4, 2) =  B(I,J,4, 2) - B(I,J,4,1)*B(I,J,1, 2)
            B(I,J,4, 3) =  B(I,J,4, 3) - B(I,J,4,1)*B(I,J,1, 3) 
     1                                 - B(I,J,4,2)*B(I,J,2, 3)
            B(I,J,4, 4) =  B(I,J,4, 4) - B(I,J,4,1)*B(I,J,1, 4) 
     1                                 - B(I,J,4,2)*B(I,J,2, 4)
     2                                 - B(I,J,4,3)*B(I,J,3, 4)
            B(I,J,4, 5) = (B(I,J,4, 5) - B(I,J,4,1)*B(I,J,1, 5) 
     1                                 - B(I,J,4,2)*B(I,J,2, 5)
     2                                 - B(I,J,4,3)*B(I,J,3, 5))
     3                    /B(I,J,4,4)
            B(I,J,4, 6) = (B(I,J,4, 6) - B(I,J,4,1)*B(I,J,1, 6) 
     1                                 - B(I,J,4,2)*B(I,J,2, 6)
     2                                 - B(I,J,4,3)*B(I,J,3, 6))
     3                    /B(I,J,4,4)
            B(I,J,4, 7) = (B(I,J,4, 7) - B(I,J,4,1)*B(I,J,1, 7) 
     1                                 - B(I,J,4,2)*B(I,J,2, 7)
     2                                 - B(I,J,4,3)*B(I,J,3, 7))
     3                    /B(I,J,4,4)
            B(I,J,4, 8) = (B(I,J,4, 8) - B(I,J,4,1)*B(I,J,1, 8) 
     1                                 - B(I,J,4,2)*B(I,J,2, 8)
     2                                 - B(I,J,4,3)*B(I,J,3, 8))
     3                    /B(I,J,4,4)
            B(I,J,4, 9) = (B(I,J,4, 9) - B(I,J,4,1)*B(I,J,1, 9) 
     1                                 - B(I,J,4,2)*B(I,J,2, 9)
     2                                 - B(I,J,4,3)*B(I,J,3, 9))
     3                    /B(I,J,4,4)
            B(I,J,4,10) = (B(I,J,4,10) - B(I,J,4,1)*B(I,J,1,10) 
     1                                 - B(I,J,4,2)*B(I,J,2,10)
     2                                 - B(I,J,4,3)*B(I,J,3,10))
     3                    /B(I,J,4,4)
            B(I,J,4,11) = (B(I,J,4,11) - B(I,J,4,1)*B(I,J,1,11) 
     1                                 - B(I,J,4,2)*B(I,J,2,11)
     2                                 - B(I,J,4,3)*B(I,J,3,11))
     3                    /B(I,J,4,4)
            B(I,J,4,12) = (B(I,J,4,12) - B(I,J,4,1)*B(I,J,1,12) 
     1                                 - B(I,J,4,2)*B(I,J,2,12)
     2                                 - B(I,J,4,3)*B(I,J,3,12))
     3                    /B(I,J,4,4)
C
C     Row 5
C
            B(I,J,5, 2) =  B(I,J,5, 2) - B(I,J,5,1)*B(I,J,1, 2)
            B(I,J,5, 3) =  B(I,J,5, 3) - B(I,J,5,1)*B(I,J,1, 3) 
     1                                 - B(I,J,5,2)*B(I,J,2, 3)
            B(I,J,5, 4) =  B(I,J,5, 4) - B(I,J,5,1)*B(I,J,1, 4) 
     1                                 - B(I,J,5,2)*B(I,J,2, 4)
     2                                 - B(I,J,5,3)*B(I,J,3, 4)
            B(I,J,5, 5) =  B(I,J,5, 5) - B(I,J,5,1)*B(I,J,1, 5) 
     1                                 - B(I,J,5,2)*B(I,J,2, 5)
     2                                 - B(I,J,5,3)*B(I,J,3, 5)
     3                                 - B(I,J,5,4)*B(I,J,4, 5)
            B(I,J,5, 6) = (B(I,J,5, 6) - B(I,J,5,1)*B(I,J,1, 6) 
     1                                 - B(I,J,5,2)*B(I,J,2, 6)
     2                                 - B(I,J,5,3)*B(I,J,3, 6)
     3                                 - B(I,J,5,4)*B(I,J,4, 6))
     4                    /B(I,J,5,5)
            B(I,J,5, 7) = (B(I,J,5, 7) - B(I,J,5,1)*B(I,J,1, 7) 
     1                                 - B(I,J,5,2)*B(I,J,2, 7)
     2                                 - B(I,J,5,3)*B(I,J,3, 7)
     3                                 - B(I,J,5,4)*B(I,J,4, 7))
     4                    /B(I,J,5,5)
            B(I,J,5, 8) = (B(I,J,5, 8) - B(I,J,5,1)*B(I,J,1, 8) 
     1                                 - B(I,J,5,2)*B(I,J,2, 8)
     2                                 - B(I,J,5,3)*B(I,J,3, 8)
     3                                 - B(I,J,5,4)*B(I,J,4, 8))
     4                    /B(I,J,5,5)
            B(I,J,5, 9) = (B(I,J,5, 9) - B(I,J,5,1)*B(I,J,1, 9) 
     1                                 - B(I,J,5,2)*B(I,J,2, 9)
     2                                 - B(I,J,5,3)*B(I,J,3, 9)
     3                                 - B(I,J,5,4)*B(I,J,4, 9))
     4                    /B(I,J,5,5)
            B(I,J,5,10) = (B(I,J,5,10) - B(I,J,5,1)*B(I,J,1,10) 
     1                                 - B(I,J,5,2)*B(I,J,2,10)
     2                                 - B(I,J,5,3)*B(I,J,3,10)
     3                                 - B(I,J,5,4)*B(I,J,4,10))
     4                    /B(I,J,5,5)
            B(I,J,5,11) = (B(I,J,5,11) - B(I,J,5,1)*B(I,J,1,11) 
     1                                 - B(I,J,5,2)*B(I,J,2,11)
     2                                 - B(I,J,5,3)*B(I,J,3,11)
     3                                 - B(I,J,5,4)*B(I,J,4,11))
     4                    /B(I,J,5,5)
            B(I,J,5,12) = (B(I,J,5,12) - B(I,J,5,1)*B(I,J,1,12) 
     1                                 - B(I,J,5,2)*B(I,J,2,12)
     2                                 - B(I,J,5,3)*B(I,J,3,12)
     3                                 - B(I,J,5,4)*B(I,J,4,12))
     4                    /B(I,J,5,5)
C
C     Row 6
C
            B(I,J,6, 2) =  B(I,J,6, 2) - B(I,J,6,1)*B(I,J,1, 2)
            B(I,J,6, 3) =  B(I,J,6, 3) - B(I,J,6,1)*B(I,J,1, 3) 
     1                                 - B(I,J,6,2)*B(I,J,2, 3)
            B(I,J,6, 4) =  B(I,J,6, 4) - B(I,J,6,1)*B(I,J,1, 4) 
     1                                 - B(I,J,6,2)*B(I,J,2, 4)
     2                                 - B(I,J,6,3)*B(I,J,3, 4)
            B(I,J,6, 5) =  B(I,J,6, 5) - B(I,J,6,1)*B(I,J,1, 5) 
     1                                 - B(I,J,6,2)*B(I,J,2, 5)
     2                                 - B(I,J,6,3)*B(I,J,3, 5)
     3                                 - B(I,J,6,4)*B(I,J,4, 5)
            B(I,J,6, 6) =  B(I,J,6, 6) - B(I,J,6,1)*B(I,J,1, 6) 
     1                                 - B(I,J,6,2)*B(I,J,2, 6)
     2                                 - B(I,J,6,3)*B(I,J,3, 6)
     3                                 - B(I,J,6,4)*B(I,J,4, 6)
     4                                 - B(I,J,6,5)*B(I,J,5, 6)
            B(I,J,6, 7) = (B(I,J,6, 7) - B(I,J,6,1)*B(I,J,1, 7) 
     1                                 - B(I,J,6,2)*B(I,J,2, 7)
     2                                 - B(I,J,6,3)*B(I,J,3, 7)
     3                                 - B(I,J,6,4)*B(I,J,4, 7)
     4                                 - B(I,J,6,5)*B(I,J,5, 7))
     5                    /B(I,J,6,6)
            B(I,J,6, 8) = (B(I,J,6, 8) - B(I,J,6,1)*B(I,J,1, 8) 
     1                                 - B(I,J,6,2)*B(I,J,2, 8)
     2                                 - B(I,J,6,3)*B(I,J,3, 8)
     3                                 - B(I,J,6,4)*B(I,J,4, 8)
     4                                 - B(I,J,6,5)*B(I,J,5, 8))
     5                    /B(I,J,6,6)
            B(I,J,6, 9) = (B(I,J,6, 9) - B(I,J,6,1)*B(I,J,1, 9) 
     1                                 - B(I,J,6,2)*B(I,J,2, 9)
     2                                 - B(I,J,6,3)*B(I,J,3, 9)
     3                                 - B(I,J,6,4)*B(I,J,4, 9)
     4                                 - B(I,J,6,5)*B(I,J,5, 9))
     5                    /B(I,J,6,6)
            B(I,J,6,10) = (B(I,J,6,10) - B(I,J,6,1)*B(I,J,1,10) 
     1                                 - B(I,J,6,2)*B(I,J,2,10)
     2                                 - B(I,J,6,3)*B(I,J,3,10)
     3                                 - B(I,J,6,4)*B(I,J,4,10)
     4                                 - B(I,J,6,5)*B(I,J,5,10))
     5                    /B(I,J,6,6)
            B(I,J,6,11) = (B(I,J,6,11) - B(I,J,6,1)*B(I,J,1,11) 
     1                                 - B(I,J,6,2)*B(I,J,2,11)
     2                                 - B(I,J,6,3)*B(I,J,3,11)
     3                                 - B(I,J,6,4)*B(I,J,4,11)
     4                                 - B(I,J,6,5)*B(I,J,5,11))
     5                    /B(I,J,6,6)
            B(I,J,6,12) = (B(I,J,6,12) - B(I,J,6,1)*B(I,J,1,12) 
     1                                 - B(I,J,6,2)*B(I,J,2,12)
     2                                 - B(I,J,6,3)*B(I,J,3,12)
     3                                 - B(I,J,6,4)*B(I,J,4,12)
     4                                 - B(I,J,6,5)*B(I,J,5,12))
     5                    /B(I,J,6,6)
C
C     Row 7
C
            B(I,J,7, 2) =  B(I,J,7, 2) - B(I,J,7,1)*B(I,J,1, 2)
            B(I,J,7, 3) =  B(I,J,7, 3) - B(I,J,7,1)*B(I,J,1, 3) 
     1                                 - B(I,J,7,2)*B(I,J,2, 3)
            B(I,J,7, 4) =  B(I,J,7, 4) - B(I,J,7,1)*B(I,J,1, 4) 
     1                                 - B(I,J,7,2)*B(I,J,2, 4)
     2                                 - B(I,J,7,3)*B(I,J,3, 4)
            B(I,J,7, 5) =  B(I,J,7, 5) - B(I,J,7,1)*B(I,J,1, 5) 
     1                                 - B(I,J,7,2)*B(I,J,2, 5)
     2                                 - B(I,J,7,3)*B(I,J,3, 5)
     3                                 - B(I,J,7,4)*B(I,J,4, 5)
            B(I,J,7, 6) =  B(I,J,7, 6) - B(I,J,7,1)*B(I,J,1, 6) 
     1                                 - B(I,J,7,2)*B(I,J,2, 6)
     2                                 - B(I,J,7,3)*B(I,J,3, 6)
     3                                 - B(I,J,7,4)*B(I,J,4, 6)
     4                                 - B(I,J,7,5)*B(I,J,5, 6)
            B(I,J,7, 7) =  B(I,J,7, 7) - B(I,J,7,1)*B(I,J,1, 7) 
     1                                 - B(I,J,7,2)*B(I,J,2, 7)
     2                                 - B(I,J,7,3)*B(I,J,3, 7)
     3                                 - B(I,J,7,4)*B(I,J,4, 7)
     4                                 - B(I,J,7,5)*B(I,J,5, 7)
     5                                 - B(I,J,7,6)*B(I,J,6, 7)
            B(I,J,7, 8) = (B(I,J,7, 8) - B(I,J,7,1)*B(I,J,1, 8) 
     1                                 - B(I,J,7,2)*B(I,J,2, 8)
     2                                 - B(I,J,7,3)*B(I,J,3, 8)
     3                                 - B(I,J,7,4)*B(I,J,4, 8)
     4                                 - B(I,J,7,5)*B(I,J,5, 8)
     5                                 - B(I,J,7,6)*B(I,J,6, 8))
     6                    /B(I,J,7,7)
            B(I,J,7, 9) = (B(I,J,7, 9) - B(I,J,7,1)*B(I,J,1, 9) 
     1                                 - B(I,J,7,2)*B(I,J,2, 9)
     2                                 - B(I,J,7,3)*B(I,J,3, 9)
     3                                 - B(I,J,7,4)*B(I,J,4, 9)
     4                                 - B(I,J,7,5)*B(I,J,5, 9)
     5                                 - B(I,J,7,6)*B(I,J,6, 9))
     6                    /B(I,J,7,7)
            B(I,J,7,10) = (B(I,J,7,10) - B(I,J,7,1)*B(I,J,1,10) 
     1                                 - B(I,J,7,2)*B(I,J,2,10)
     2                                 - B(I,J,7,3)*B(I,J,3,10)
     3                                 - B(I,J,7,4)*B(I,J,4,10)
     4                                 - B(I,J,7,5)*B(I,J,5,10)
     5                                 - B(I,J,7,6)*B(I,J,6,10))
     6                    /B(I,J,7,7)
            B(I,J,7,11) = (B(I,J,7,11) - B(I,J,7,1)*B(I,J,1,11) 
     1                                 - B(I,J,7,2)*B(I,J,2,11)
     2                                 - B(I,J,7,3)*B(I,J,3,11)
     3                                 - B(I,J,7,4)*B(I,J,4,11)
     4                                 - B(I,J,7,5)*B(I,J,5,11)
     5                                 - B(I,J,7,6)*B(I,J,6,11))
     6                    /B(I,J,7,7)
            B(I,J,7,12) = (B(I,J,7,12) - B(I,J,7,1)*B(I,J,1,12) 
     1                                 - B(I,J,7,2)*B(I,J,2,12)
     2                                 - B(I,J,7,3)*B(I,J,3,12)
     3                                 - B(I,J,7,4)*B(I,J,4,12)
     4                                 - B(I,J,7,5)*B(I,J,5,12)
     5                                 - B(I,J,7,6)*B(I,J,6,12))
     6                    /B(I,J,7,7)
C
C     Row 8
C
            B(I,J,8, 2) =  B(I,J,8, 2) - B(I,J,8,1)*B(I,J,1, 2)
            B(I,J,8, 3) =  B(I,J,8, 3) - B(I,J,8,1)*B(I,J,1, 3) 
     1                                 - B(I,J,8,2)*B(I,J,2, 3)
            B(I,J,8, 4) =  B(I,J,8, 4) - B(I,J,8,1)*B(I,J,1, 4) 
     1                                 - B(I,J,8,2)*B(I,J,2, 4)
     2                                 - B(I,J,8,3)*B(I,J,3, 4)
            B(I,J,8, 5) =  B(I,J,8, 5) - B(I,J,8,1)*B(I,J,1, 5) 
     1                                 - B(I,J,8,2)*B(I,J,2, 5)
     2                                 - B(I,J,8,3)*B(I,J,3, 5)
     3                                 - B(I,J,8,4)*B(I,J,4, 5)
            B(I,J,8, 6) =  B(I,J,8, 6) - B(I,J,8,1)*B(I,J,1, 6) 
     1                                 - B(I,J,8,2)*B(I,J,2, 6)
     2                                 - B(I,J,8,3)*B(I,J,3, 6)
     3                                 - B(I,J,8,4)*B(I,J,4, 6)
     4                                 - B(I,J,8,5)*B(I,J,5, 6)
            B(I,J,8, 7) =  B(I,J,8, 7) - B(I,J,8,1)*B(I,J,1, 7) 
     1                                 - B(I,J,8,2)*B(I,J,2, 7)
     2                                 - B(I,J,8,3)*B(I,J,3, 7)
     3                                 - B(I,J,8,4)*B(I,J,4, 7)
     4                                 - B(I,J,8,5)*B(I,J,5, 7)
     5                                 - B(I,J,8,6)*B(I,J,6, 7)
            B(I,J,8, 8) =  B(I,J,8, 8) - B(I,J,8,1)*B(I,J,1, 8) 
     1                                 - B(I,J,8,2)*B(I,J,2, 8)
     2                                 - B(I,J,8,3)*B(I,J,3, 8)
     3                                 - B(I,J,8,4)*B(I,J,4, 8)
     4                                 - B(I,J,8,5)*B(I,J,5, 8)
     5                                 - B(I,J,8,6)*B(I,J,6, 8)
     6                                 - B(I,J,8,7)*B(I,J,7, 8)
            B(I,J,8, 9) = (B(I,J,8, 9) - B(I,J,8,1)*B(I,J,1, 9) 
     1                                 - B(I,J,8,2)*B(I,J,2, 9)
     2                                 - B(I,J,8,3)*B(I,J,3, 9)
     3                                 - B(I,J,8,4)*B(I,J,4, 9)
     4                                 - B(I,J,8,5)*B(I,J,5, 9)
     5                                 - B(I,J,8,6)*B(I,J,6, 9)
     6                                 - B(I,J,8,7)*B(I,J,7, 9))
     7                    /B(I,J,8,8)
            B(I,J,8,10) = (B(I,J,8,10) - B(I,J,8,1)*B(I,J,1,10) 
     1                                 - B(I,J,8,2)*B(I,J,2,10)
     2                                 - B(I,J,8,3)*B(I,J,3,10)
     3                                 - B(I,J,8,4)*B(I,J,4,10)
     4                                 - B(I,J,8,5)*B(I,J,5,10)
     5                                 - B(I,J,8,6)*B(I,J,6,10)
     6                                 - B(I,J,8,7)*B(I,J,7,10))
     7                    /B(I,J,8,8)
            B(I,J,8,11) = (B(I,J,8,11) - B(I,J,8,1)*B(I,J,1,11) 
     1                                 - B(I,J,8,2)*B(I,J,2,11)
     2                                 - B(I,J,8,3)*B(I,J,3,11)
     3                                 - B(I,J,8,4)*B(I,J,4,11)
     4                                 - B(I,J,8,5)*B(I,J,5,11)
     5                                 - B(I,J,8,6)*B(I,J,6,11)
     6                                 - B(I,J,8,7)*B(I,J,7,11))
     7                    /B(I,J,8,8)
            B(I,J,8,12) = (B(I,J,8,12) - B(I,J,8,1)*B(I,J,1,12) 
     1                                 - B(I,J,8,2)*B(I,J,2,12)
     2                                 - B(I,J,8,3)*B(I,J,3,12)
     3                                 - B(I,J,8,4)*B(I,J,4,12)
     4                                 - B(I,J,8,5)*B(I,J,5,12)
     5                                 - B(I,J,8,6)*B(I,J,6,12)
     6                                 - B(I,J,8,7)*B(I,J,7,12))
     7                    /B(I,J,8,8)
C
C     Row 9
C
            B(I,J,9, 2) =  B(I,J,9, 2) - B(I,J,9,1)*B(I,J,1, 2)
            B(I,J,9, 3) =  B(I,J,9, 3) - B(I,J,9,1)*B(I,J,1, 3) 
     1                                 - B(I,J,9,2)*B(I,J,2, 3)
            B(I,J,9, 4) =  B(I,J,9, 4) - B(I,J,9,1)*B(I,J,1, 4) 
     1                                 - B(I,J,9,2)*B(I,J,2, 4)
     2                                 - B(I,J,9,3)*B(I,J,3, 4)
            B(I,J,9, 5) =  B(I,J,9, 5) - B(I,J,9,1)*B(I,J,1, 5) 
     1                                 - B(I,J,9,2)*B(I,J,2, 5)
     2                                 - B(I,J,9,3)*B(I,J,3, 5)
     3                                 - B(I,J,9,4)*B(I,J,4, 5)
            B(I,J,9, 6) =  B(I,J,9, 6) - B(I,J,9,1)*B(I,J,1, 6) 
     1                                 - B(I,J,9,2)*B(I,J,2, 6)
     2                                 - B(I,J,9,3)*B(I,J,3, 6)
     3                                 - B(I,J,9,4)*B(I,J,4, 6)
     4                                 - B(I,J,9,5)*B(I,J,5, 6)
            B(I,J,9, 7) =  B(I,J,9, 7) - B(I,J,9,1)*B(I,J,1, 7) 
     1                                 - B(I,J,9,2)*B(I,J,2, 7)
     2                                 - B(I,J,9,3)*B(I,J,3, 7)
     3                                 - B(I,J,9,4)*B(I,J,4, 7)
     4                                 - B(I,J,9,5)*B(I,J,5, 7)
     5                                 - B(I,J,9,6)*B(I,J,6, 7)
            B(I,J,9, 8) =  B(I,J,9, 8) - B(I,J,9,1)*B(I,J,1, 8) 
     1                                 - B(I,J,9,2)*B(I,J,2, 8)
     2                                 - B(I,J,9,3)*B(I,J,3, 8)
     3                                 - B(I,J,9,4)*B(I,J,4, 8)
     4                                 - B(I,J,9,5)*B(I,J,5, 8)
     5                                 - B(I,J,9,6)*B(I,J,6, 8)
     6                                 - B(I,J,9,7)*B(I,J,7, 8)
            B(I,J,9, 9) =  B(I,J,9, 9) - B(I,J,9,1)*B(I,J,1, 9) 
     1                                 - B(I,J,9,2)*B(I,J,2, 9)
     2                                 - B(I,J,9,3)*B(I,J,3, 9)
     3                                 - B(I,J,9,4)*B(I,J,4, 9)
     4                                 - B(I,J,9,5)*B(I,J,5, 9)
     5                                 - B(I,J,9,6)*B(I,J,6, 9)
     6                                 - B(I,J,9,7)*B(I,J,7, 9)
     7                                 - B(I,J,9,8)*B(I,J,8, 9)
            B(I,J,9,10) = (B(I,J,9,10) - B(I,J,9,1)*B(I,J,1,10) 
     1                                 - B(I,J,9,2)*B(I,J,2,10)
     2                                 - B(I,J,9,3)*B(I,J,3,10)
     3                                 - B(I,J,9,4)*B(I,J,4,10)
     4                                 - B(I,J,9,5)*B(I,J,5,10)
     5                                 - B(I,J,9,6)*B(I,J,6,10)
     6                                 - B(I,J,9,7)*B(I,J,7,10)
     7                                 - B(I,J,9,8)*B(I,J,8,10))
     8                    /B(I,J,9,9)
            B(I,J,9,11) = (B(I,J,9,11) - B(I,J,9,1)*B(I,J,1,11) 
     1                                 - B(I,J,9,2)*B(I,J,2,11)
     2                                 - B(I,J,9,3)*B(I,J,3,11)
     3                                 - B(I,J,9,4)*B(I,J,4,11)
     4                                 - B(I,J,9,5)*B(I,J,5,11)
     5                                 - B(I,J,9,6)*B(I,J,6,11)
     6                                 - B(I,J,9,7)*B(I,J,7,11)
     7                                 - B(I,J,9,8)*B(I,J,8,11))
     8                    /B(I,J,9,9)
            B(I,J,9,12) = (B(I,J,9,12) - B(I,J,9,1)*B(I,J,1,12) 
     1                                 - B(I,J,9,2)*B(I,J,2,12)
     2                                 - B(I,J,9,3)*B(I,J,3,12)
     3                                 - B(I,J,9,4)*B(I,J,4,12)
     4                                 - B(I,J,9,5)*B(I,J,5,12)
     5                                 - B(I,J,9,6)*B(I,J,6,12)
     6                                 - B(I,J,9,7)*B(I,J,7,12)
     7                                 - B(I,J,9,8)*B(I,J,8,12))
     8                    /B(I,J,9,9)
C
C     Row 10
C
            B(I,J,10, 2)= B(I,J,10, 2) - B(I,J,10,1)*B(I,J,1, 2)
            B(I,J,10, 3)= B(I,J,10, 3) - B(I,J,10,1)*B(I,J,1, 3) 
     1                                 - B(I,J,10,2)*B(I,J,2, 3)
            B(I,J,10, 4)= B(I,J,10, 4) - B(I,J,10,1)*B(I,J,1, 4) 
     1                                 - B(I,J,10,2)*B(I,J,2, 4)
     2                                 - B(I,J,10,3)*B(I,J,3, 4)
            B(I,J,10, 5)= B(I,J,10, 5) - B(I,J,10,1)*B(I,J,1, 5) 
     1                                 - B(I,J,10,2)*B(I,J,2, 5)
     2                                 - B(I,J,10,3)*B(I,J,3, 5)
     3                                 - B(I,J,10,4)*B(I,J,4, 5)
            B(I,J,10, 6)= B(I,J,10, 6) - B(I,J,10,1)*B(I,J,1, 6) 
     1                                 - B(I,J,10,2)*B(I,J,2, 6)
     2                                 - B(I,J,10,3)*B(I,J,3, 6)
     3                                 - B(I,J,10,4)*B(I,J,4, 6)
     4                                 - B(I,J,10,5)*B(I,J,5, 6)
            B(I,J,10, 7)= B(I,J,10, 7) - B(I,J,10,1)*B(I,J,1, 7) 
     1                                 - B(I,J,10,2)*B(I,J,2, 7)
     2                                 - B(I,J,10,3)*B(I,J,3, 7)
     3                                 - B(I,J,10,4)*B(I,J,4, 7)
     4                                 - B(I,J,10,5)*B(I,J,5, 7)
     5                                 - B(I,J,10,6)*B(I,J,6, 7)
            B(I,J,10, 8)= B(I,J,10, 8) - B(I,J,10,1)*B(I,J,1, 8) 
     1                                 - B(I,J,10,2)*B(I,J,2, 8)
     2                                 - B(I,J,10,3)*B(I,J,3, 8)
     3                                 - B(I,J,10,4)*B(I,J,4, 8)
     4                                 - B(I,J,10,5)*B(I,J,5, 8)
     5                                 - B(I,J,10,6)*B(I,J,6, 8)
     6                                 - B(I,J,10,7)*B(I,J,7, 8)
            B(I,J,10, 9)= B(I,J,10, 9) - B(I,J,10,1)*B(I,J,1, 9) 
     1                                 - B(I,J,10,2)*B(I,J,2, 9)
     2                                 - B(I,J,10,3)*B(I,J,3, 9)
     3                                 - B(I,J,10,4)*B(I,J,4, 9)
     4                                 - B(I,J,10,5)*B(I,J,5, 9)
     5                                 - B(I,J,10,6)*B(I,J,6, 9)
     6                                 - B(I,J,10,7)*B(I,J,7, 9)
     7                                 - B(I,J,10,8)*B(I,J,8, 9)
            B(I,J,10,10)= B(I,J,10,10) - B(I,J,10,1)*B(I,J,1,10) 
     1                                 - B(I,J,10,2)*B(I,J,2,10)
     2                                 - B(I,J,10,3)*B(I,J,3,10)
     3                                 - B(I,J,10,4)*B(I,J,4,10)
     4                                 - B(I,J,10,5)*B(I,J,5,10)
     5                                 - B(I,J,10,6)*B(I,J,6,10)
     6                                 - B(I,J,10,7)*B(I,J,7,10)
     7                                 - B(I,J,10,8)*B(I,J,8,10)
     8                                 - B(I,J,10,9)*B(I,J,9,10)
            B(I,J,10,11)=(B(I,J,10,11) - B(I,J,10,1)*B(I,J,1,11) 
     1                                 - B(I,J,10,2)*B(I,J,2,11)
     2                                 - B(I,J,10,3)*B(I,J,3,11)
     3                                 - B(I,J,10,4)*B(I,J,4,11)
     4                                 - B(I,J,10,5)*B(I,J,5,11)
     5                                 - B(I,J,10,6)*B(I,J,6,11)
     6                                 - B(I,J,10,7)*B(I,J,7,11)
     7                                 - B(I,J,10,8)*B(I,J,8,11)
     8                                 - B(I,J,10,9)*B(I,J,9,11))
     9                   /B(I,J,10,10)
            B(I,J,10,12)=(B(I,J,10,12) - B(I,J,10,1)*B(I,J,1,12) 
     1                                 - B(I,J,10,2)*B(I,J,2,12)
     2                                 - B(I,J,10,3)*B(I,J,3,12)
     3                                 - B(I,J,10,4)*B(I,J,4,12)
     4                                 - B(I,J,10,5)*B(I,J,5,12)
     5                                 - B(I,J,10,6)*B(I,J,6,12)
     6                                 - B(I,J,10,7)*B(I,J,7,12)
     7                                 - B(I,J,10,8)*B(I,J,8,12)
     8                                 - B(I,J,10,9)*B(I,J,9,12))
     9                   /B(I,J,10,10)
C
C     Row 11
C
            B(I,J,11, 2)= B(I,J,11, 2) - B(I,J,11,1)*B(I,J,1, 2)
            B(I,J,11, 3)= B(I,J,11, 3) - B(I,J,11,1)*B(I,J,1, 3) 
     1                                 - B(I,J,11,2)*B(I,J,2, 3)
            B(I,J,11, 4)= B(I,J,11, 4) - B(I,J,11,1)*B(I,J,1, 4) 
     1                                 - B(I,J,11,2)*B(I,J,2, 4)
     2                                 - B(I,J,11,3)*B(I,J,3, 4)
            B(I,J,11, 5)= B(I,J,11, 5) - B(I,J,11,1)*B(I,J,1, 5) 
     1                                 - B(I,J,11,2)*B(I,J,2, 5)
     2                                 - B(I,J,11,3)*B(I,J,3, 5)
     3                                 - B(I,J,11,4)*B(I,J,4, 5)
            B(I,J,11, 6)= B(I,J,11, 6) - B(I,J,11,1)*B(I,J,1, 6) 
     1                                 - B(I,J,11,2)*B(I,J,2, 6)
     2                                 - B(I,J,11,3)*B(I,J,3, 6)
     3                                 - B(I,J,11,4)*B(I,J,4, 6)
     4                                 - B(I,J,11,5)*B(I,J,5, 6)
            B(I,J,11, 7)= B(I,J,11, 7) - B(I,J,11,1)*B(I,J,1, 7) 
     1                                 - B(I,J,11,2)*B(I,J,2, 7)
     2                                 - B(I,J,11,3)*B(I,J,3, 7)
     3                                 - B(I,J,11,4)*B(I,J,4, 7)
     4                                 - B(I,J,11,5)*B(I,J,5, 7)
     5                                 - B(I,J,11,6)*B(I,J,6, 7)
            B(I,J,11, 8)= B(I,J,11, 8) - B(I,J,11,1)*B(I,J,1, 8) 
     1                                 - B(I,J,11,2)*B(I,J,2, 8)
     2                                 - B(I,J,11,3)*B(I,J,3, 8)
     3                                 - B(I,J,11,4)*B(I,J,4, 8)
     4                                 - B(I,J,11,5)*B(I,J,5, 8)
     5                                 - B(I,J,11,6)*B(I,J,6, 8)
     6                                 - B(I,J,11,7)*B(I,J,7, 8)
            B(I,J,11, 9)= B(I,J,11, 9) - B(I,J,11,1)*B(I,J,1, 9) 
     1                                 - B(I,J,11,2)*B(I,J,2, 9)
     2                                 - B(I,J,11,3)*B(I,J,3, 9)
     3                                 - B(I,J,11,4)*B(I,J,4, 9)
     4                                 - B(I,J,11,5)*B(I,J,5, 9)
     5                                 - B(I,J,11,6)*B(I,J,6, 9)
     6                                 - B(I,J,11,7)*B(I,J,7, 9)
     7                                 - B(I,J,11,8)*B(I,J,8, 9)
            B(I,J,11,10)= B(I,J,11,10) - B(I,J,11,1)*B(I,J,1,10) 
     1                                 - B(I,J,11,2)*B(I,J,2,10)
     2                                 - B(I,J,11,3)*B(I,J,3,10)
     3                                 - B(I,J,11,4)*B(I,J,4,10)
     4                                 - B(I,J,11,5)*B(I,J,5,10)
     5                                 - B(I,J,11,6)*B(I,J,6,10)
     6                                 - B(I,J,11,7)*B(I,J,7,10)
     7                                 - B(I,J,11,8)*B(I,J,8,10)
     8                                 - B(I,J,11,9)*B(I,J,9,10)
            B(I,J,11,11)= B(I,J,11,11) - B(I,J,11,1)*B(I,J,1,11) 
     1                                 - B(I,J,11,2)*B(I,J,2,11)
     2                                 - B(I,J,11,3)*B(I,J,3,11)
     3                                 - B(I,J,11,4)*B(I,J,4,11)
     4                                 - B(I,J,11,5)*B(I,J,5,11)
     5                                 - B(I,J,11,6)*B(I,J,6,11)
     6                                 - B(I,J,11,7)*B(I,J,7,11)
     7                                 - B(I,J,11,8)*B(I,J,8,11)
     8                                 - B(I,J,11,9)*B(I,J,9,11)
     9                                 - B(I,J,11,10)*B(I,J,10,11)
            B(I,J,11,12)=(B(I,J,11,12) - B(I,J,11,1)*B(I,J,1,12) 
     1                                 - B(I,J,11,2)*B(I,J,2,12)
     2                                 - B(I,J,11,3)*B(I,J,3,12)
     3                                 - B(I,J,11,4)*B(I,J,4,12)
     4                                 - B(I,J,11,5)*B(I,J,5,12)
     5                                 - B(I,J,11,6)*B(I,J,6,12)
     6                                 - B(I,J,11,7)*B(I,J,7,12)
     7                                 - B(I,J,11,8)*B(I,J,8,12)
     8                                 - B(I,J,11,9)*B(I,J,9,12)
     9                                 - B(I,J,11,10)*B(I,J,10,12))
     A                   /B(I,J,11,11)
C
C     Row 12
C
            B(I,J,12, 2)= B(I,J,12, 2) - B(I,J,12,1)*B(I,J,1, 2)
            B(I,J,12, 3)= B(I,J,12, 3) - B(I,J,12,1)*B(I,J,1, 3) 
     1                                 - B(I,J,12,2)*B(I,J,2, 3)
            B(I,J,12, 4)= B(I,J,12, 4) - B(I,J,12,1)*B(I,J,1, 4) 
     1                                 - B(I,J,12,2)*B(I,J,2, 4)
     2                                 - B(I,J,12,3)*B(I,J,3, 4)
            B(I,J,12, 5)= B(I,J,12, 5) - B(I,J,12,1)*B(I,J,1, 5) 
     1                                 - B(I,J,12,2)*B(I,J,2, 5)
     2                                 - B(I,J,12,3)*B(I,J,3, 5)
     3                                 - B(I,J,12,4)*B(I,J,4, 5)
            B(I,J,12, 6)= B(I,J,12, 6) - B(I,J,12,1)*B(I,J,1, 6) 
     1                                 - B(I,J,12,2)*B(I,J,2, 6)
     2                                 - B(I,J,12,3)*B(I,J,3, 6)
     3                                 - B(I,J,12,4)*B(I,J,4, 6)
     4                                 - B(I,J,12,5)*B(I,J,5, 6)
            B(I,J,12, 7)= B(I,J,12, 7) - B(I,J,12,1)*B(I,J,1, 7) 
     1                                 - B(I,J,12,2)*B(I,J,2, 7)
     2                                 - B(I,J,12,3)*B(I,J,3, 7)
     3                                 - B(I,J,12,4)*B(I,J,4, 7)
     4                                 - B(I,J,12,5)*B(I,J,5, 7)
     5                                 - B(I,J,12,6)*B(I,J,6, 7)
            B(I,J,12, 8)= B(I,J,12, 8) - B(I,J,12,1)*B(I,J,1, 8) 
     1                                 - B(I,J,12,2)*B(I,J,2, 8)
     2                                 - B(I,J,12,3)*B(I,J,3, 8)
     3                                 - B(I,J,12,4)*B(I,J,4, 8)
     4                                 - B(I,J,12,5)*B(I,J,5, 8)
     5                                 - B(I,J,12,6)*B(I,J,6, 8)
     6                                 - B(I,J,12,7)*B(I,J,7, 8)
            B(I,J,12, 9)= B(I,J,12, 9) - B(I,J,12,1)*B(I,J,1, 9) 
     1                                 - B(I,J,12,2)*B(I,J,2, 9)
     2                                 - B(I,J,12,3)*B(I,J,3, 9)
     3                                 - B(I,J,12,4)*B(I,J,4, 9)
     4                                 - B(I,J,12,5)*B(I,J,5, 9)
     5                                 - B(I,J,12,6)*B(I,J,6, 9)
     6                                 - B(I,J,12,7)*B(I,J,7, 9)
     7                                 - B(I,J,12,8)*B(I,J,8, 9)
            B(I,J,12,10)= B(I,J,12,10) - B(I,J,12,1)*B(I,J,1,10) 
     1                                 - B(I,J,12,2)*B(I,J,2,10)
     2                                 - B(I,J,12,3)*B(I,J,3,10)
     3                                 - B(I,J,12,4)*B(I,J,4,10)
     4                                 - B(I,J,12,5)*B(I,J,5,10)
     5                                 - B(I,J,12,6)*B(I,J,6,10)
     6                                 - B(I,J,12,7)*B(I,J,7,10)
     7                                 - B(I,J,12,8)*B(I,J,8,10)
     8                                 - B(I,J,12,9)*B(I,J,9,10)
            B(I,J,12,11)= B(I,J,12,11) - B(I,J,12,1)*B(I,J,1,11) 
     1                                 - B(I,J,12,2)*B(I,J,2,11)
     2                                 - B(I,J,12,3)*B(I,J,3,11)
     3                                 - B(I,J,12,4)*B(I,J,4,11)
     4                                 - B(I,J,12,5)*B(I,J,5,11)
     5                                 - B(I,J,12,6)*B(I,J,6,11)
     6                                 - B(I,J,12,7)*B(I,J,7,11)
     7                                 - B(I,J,12,8)*B(I,J,8,11)
     8                                 - B(I,J,12,9)*B(I,J,9,11)
     9                                 - B(I,J,12,10)*B(I,J,10,11)
            B(I,J,12,12)= B(I,J,12,12) - B(I,J,12,1)*B(I,J,1,12) 
     1                                 - B(I,J,12,2)*B(I,J,2,12)
     2                                 - B(I,J,12,3)*B(I,J,3,12)
     3                                 - B(I,J,12,4)*B(I,J,4,12)
     4                                 - B(I,J,12,5)*B(I,J,5,12)
     5                                 - B(I,J,12,6)*B(I,J,6,12)
     6                                 - B(I,J,12,7)*B(I,J,7,12)
     7                                 - B(I,J,12,8)*B(I,J,8,12)
     8                                 - B(I,J,12,9)*B(I,J,9,12)
     9                                 - B(I,J,12,10)*B(I,J,10,12)
     A                                 - B(I,J,12,11)*B(I,J,11,12)
 200     CONTINUE
C
C     Form U(i,i+1) = L[inverse](i,i) * C(i) for all except last term
C
         IF (J .NE. NPTS) THEN
            DO 500 L = 1, 12
               DO 480 I = 1, NSIM
C
C     Forward sweep of lu decomposition to form L[inverse]
C
                  C(I,J,1,L) = C(I,J,1,L) / B(I,J,1,1)
                  C(I,J,2,L) = ( C(I,J,2,L) - B(I,J,2,1)*C(I,J,1,L) ) 
     1                         / B(I,J,2,2)
                  C(I,J,3,L) = ( C(I,J,3,L) - B(I,J,3,1)*C(I,J,1,L) 
     1                                      - B(I,J,3,2)*C(I,J,2,L) ) 
     2                         / B(I,J,3,3)
                  C(I,J,4,L) = ( C(I,J,4,L) - B(I,J,4,1)*C(I,J,1,L) 
     1                                      - B(I,J,4,2)*C(I,J,2,L)  
     2                                      - B(I,J,4,3)*C(I,J,3,L) ) 
     3                         / B(I,J,4,4)
                  C(I,J,5,L) = ( C(I,J,5,L) - B(I,J,5,1)*C(I,J,1,L) 
     1                                      - B(I,J,5,2)*C(I,J,2,L) 
     2                                      - B(I,J,5,3)*C(I,J,3,L)  
     3                                      - B(I,J,5,4)*C(I,J,4,L) )
     4                         / B(I,J,5,5)
                  C(I,J,6,L) = ( C(I,J,6,L) - B(I,J,6,1)*C(I,J,1,L) 
     1                                      - B(I,J,6,2)*C(I,J,2,L) 
     2                                      - B(I,J,6,3)*C(I,J,3,L)  
     3                                      - B(I,J,6,4)*C(I,J,4,L)
     4                                      - B(I,J,6,5)*C(I,J,5,L) )
     5                         / B(I,J,6,6)
                  C(I,J,7,L) = ( C(I,J,7,L) - B(I,J,7,1)*C(I,J,1,L) 
     1                                      - B(I,J,7,2)*C(I,J,2,L) 
     2                                      - B(I,J,7,3)*C(I,J,3,L)  
     3                                      - B(I,J,7,4)*C(I,J,4,L)
     4                                      - B(I,J,7,5)*C(I,J,5,L)
     5                                      - B(I,J,7,6)*C(I,J,6,L) )
     6                         / B(I,J,7,7)
                  C(I,J,8,L) = ( C(I,J,8,L) - B(I,J,8,1)*C(I,J,1,L) 
     1                                      - B(I,J,8,2)*C(I,J,2,L) 
     2                                      - B(I,J,8,3)*C(I,J,3,L)  
     3                                      - B(I,J,8,4)*C(I,J,4,L)
     4                                      - B(I,J,8,5)*C(I,J,5,L)
     5                                      - B(I,J,8,6)*C(I,J,6,L)
     6                                      - B(I,J,8,7)*C(I,J,7,L) )
     7                         / B(I,J,8,8)
                  C(I,J,9,L) = ( C(I,J,9,L) - B(I,J,9,1)*C(I,J,1,L) 
     1                                      - B(I,J,9,2)*C(I,J,2,L) 
     2                                      - B(I,J,9,3)*C(I,J,3,L)  
     3                                      - B(I,J,9,4)*C(I,J,4,L)
     4                                      - B(I,J,9,5)*C(I,J,5,L)
     5                                      - B(I,J,9,6)*C(I,J,6,L)
     6                                      - B(I,J,9,7)*C(I,J,7,L)
     7                                      - B(I,J,9,8)*C(I,J,8,L) )
     8                         / B(I,J,9,9)
                  C(I,J,10,L)=( C(I,J,10,L) - B(I,J,10,1)*C(I,J,1,L) 
     1                                      - B(I,J,10,2)*C(I,J,2,L) 
     2                                      - B(I,J,10,3)*C(I,J,3,L)  
     3                                      - B(I,J,10,4)*C(I,J,4,L)
     4                                      - B(I,J,10,5)*C(I,J,5,L)
     5                                      - B(I,J,10,6)*C(I,J,6,L)
     6                                      - B(I,J,10,7)*C(I,J,7,L)
     7                                      - B(I,J,10,8)*C(I,J,8,L)
     8                                      - B(I,J,10,9)*C(I,J,9,L) )
     9                         / B(I,J,10,10)
                  C(I,J,11,L)=( C(I,J,11,L) - B(I,J,11,1)*C(I,J,1,L) 
     1                                      - B(I,J,11,2)*C(I,J,2,L) 
     2                                      - B(I,J,11,3)*C(I,J,3,L)  
     3                                      - B(I,J,11,4)*C(I,J,4,L)
     4                                      - B(I,J,11,5)*C(I,J,5,L)
     5                                      - B(I,J,11,6)*C(I,J,6,L)
     6                                      - B(I,J,11,7)*C(I,J,7,L)
     7                                      - B(I,J,11,8)*C(I,J,8,L)
     8                                      - B(I,J,11,9)*C(I,J,9,L) 
     9                                      - B(I,J,11,10)*C(I,J,10,L) )
     A                         / B(I,J,11,11)
                  C(I,J,12,L)=( C(I,J,12,L) - B(I,J,12,1)*C(I,J,1,L) 
     1                                      - B(I,J,12,2)*C(I,J,2,L) 
     2                                      - B(I,J,12,3)*C(I,J,3,L)  
     3                                      - B(I,J,12,4)*C(I,J,4,L)
     4                                      - B(I,J,12,5)*C(I,J,5,L)
     5                                      - B(I,J,12,6)*C(I,J,6,L)
     6                                      - B(I,J,12,7)*C(I,J,7,L)
     7                                      - B(I,J,12,8)*C(I,J,8,L)
     8                                      - B(I,J,12,9)*C(I,J,9,L) 
     9                                      - B(I,J,12,10)*C(I,J,10,L)
     A                                      - B(I,J,12,11)*C(I,J,11,L) )
     B                         / B(I,J,12,12)
C
C     Backward Sweep of lu
C
                  C(I,J,11,L) = C(I,J,11,L) - B(I,J,11,12)*C(I,J,12,L)
                  C(I,J,10,L) = C(I,J,10,L) - B(I,J,10,11)*C(I,J,11,L)
     1                                      - B(I,J,10,12)*C(I,J,12,L)
                  C(I,J, 9,L) = C(I,J, 9,L) - B(I,J, 9,10)*C(I,J,10,L)
     1                                      - B(I,J, 9,11)*C(I,J,11,L)
     2                                      - B(I,J, 9,12)*C(I,J,12,L)
                  C(I,J, 8,L) = C(I,J, 8,L) - B(I,J, 8, 9)*C(I,J, 9,L)
     1                                      - B(I,J, 8,10)*C(I,J,10,L)
     2                                      - B(I,J, 8,11)*C(I,J,11,L)
     3                                      - B(I,J, 8,12)*C(I,J,12,L)
                  C(I,J, 7,L) = C(I,J, 7,L) - B(I,J, 7, 8)*C(I,J, 8,L)
     1                                      - B(I,J, 7, 9)*C(I,J, 9,L)
     2                                      - B(I,J, 7,10)*C(I,J,10,L)
     3                                      - B(I,J, 7,11)*C(I,J,11,L)
     4                                      - B(I,J, 7,12)*C(I,J,12,L)
                  C(I,J, 6,L) = C(I,J, 6,L) - B(I,J, 6, 7)*C(I,J, 7,L)
     1                                      - B(I,J, 6, 8)*C(I,J, 8,L)
     2                                      - B(I,J, 6, 9)*C(I,J, 9,L)
     3                                      - B(I,J, 6,10)*C(I,J,10,L)
     4                                      - B(I,J, 6,11)*C(I,J,11,L)
     5                                      - B(I,J, 6,12)*C(I,J,12,L)
                  C(I,J, 5,L) = C(I,J, 5,L) - B(I,J, 5, 6)*C(I,J, 6,L)
     1                                      - B(I,J, 5, 7)*C(I,J, 7,L)
     2                                      - B(I,J, 5, 8)*C(I,J, 8,L)
     3                                      - B(I,J, 5, 9)*C(I,J, 9,L)
     4                                      - B(I,J, 5,10)*C(I,J,10,L)
     5                                      - B(I,J, 5,11)*C(I,J,11,L)
     6                                      - B(I,J, 5,12)*C(I,J,12,L)
                  C(I,J, 4,L) = C(I,J, 4,L) - B(I,J, 4, 5)*C(I,J, 5,L)
     1                                      - B(I,J, 4, 6)*C(I,J, 6,L)
     2                                      - B(I,J, 4, 7)*C(I,J, 7,L)
     3                                      - B(I,J, 4, 8)*C(I,J, 8,L)
     4                                      - B(I,J, 4, 9)*C(I,J, 9,L)
     5                                      - B(I,J, 4,10)*C(I,J,10,L)
     6                                      - B(I,J, 4,11)*C(I,J,11,L)
     7                                      - B(I,J, 4,12)*C(I,J,12,L)
                  C(I,J, 3,L) = C(I,J, 3,L) - B(I,J, 3, 4)*C(I,J, 4,L)
     1                                      - B(I,J, 3, 5)*C(I,J, 5,L)
     2                                      - B(I,J, 3, 6)*C(I,J, 6,L)
     3                                      - B(I,J, 3, 7)*C(I,J, 7,L)
     4                                      - B(I,J, 3, 8)*C(I,J, 8,L)
     5                                      - B(I,J, 3, 9)*C(I,J, 9,L)
     6                                      - B(I,J, 3,10)*C(I,J,10,L)
     7                                      - B(I,J, 3,11)*C(I,J,11,L)
     8                                      - B(I,J, 3,12)*C(I,J,12,L)
                  C(I,J, 2,L) = C(I,J, 2,L) - B(I,J, 2, 3)*C(I,J, 3,L)
     1                                      - B(I,J, 2, 4)*C(I,J, 4,L)
     2                                      - B(I,J, 2, 5)*C(I,J, 5,L)
     3                                      - B(I,J, 2, 6)*C(I,J, 6,L)
     4                                      - B(I,J, 2, 7)*C(I,J, 7,L)
     5                                      - B(I,J, 2, 8)*C(I,J, 8,L)
     6                                      - B(I,J, 2, 9)*C(I,J, 9,L)
     7                                      - B(I,J, 2,10)*C(I,J,10,L)
     8                                      - B(I,J, 2,11)*C(I,J,11,L)
     9                                      - B(I,J, 2,12)*C(I,J,12,L)
                  C(I,J, 1,L) = C(I,J, 1,L) - B(I,J, 1, 2)*C(I,J, 2,L)
     1                                      - B(I,J, 1, 3)*C(I,J, 3,L)
     2                                      - B(I,J, 1, 4)*C(I,J, 4,L)
     3                                      - B(I,J, 1, 5)*C(I,J, 5,L)
     4                                      - B(I,J, 1, 6)*C(I,J, 6,L)
     5                                      - B(I,J, 1, 7)*C(I,J, 7,L)
     6                                      - B(I,J, 1, 8)*C(I,J, 8,L)
     7                                      - B(I,J, 1, 9)*C(I,J, 9,L)
     8                                      - B(I,J, 1,10)*C(I,J,10,L)
     9                                      - B(I,J, 1,11)*C(I,J,11,L)
     A                                      - B(I,J, 1,12)*C(I,J,12,L)
 480           CONTINUE
 500        CONTINUE
         ENDIF
 600  CONTINUE
C
C     Finished forming LU
C
      RETURN
      END






      SUBROUTINE VLUTRS (NPTS, NSIM, A, B, C)
C
C     Routine to do a vectorized LU decomposition of a 
C     scalar tridiagonal system of equations.
C     The recursive nature of the LU decomposition leads
C     to a form which is not vectorizable.  Therefore, the
C     approach of this routine is to not attempt vectorization
C     of the LU, but to accomplish the vectorization through
C     simultaneously doing NSIM LU decompositions.
C     This form is with the identity on the diagonal of U.
C
C     NPTS:      The number of elements in the tridiagonal system
C     NSIM:      The number of tridiagonal systems to simultaneously invert
C     A:         The supra-diagonal scalar term
C     B:         The       diagonal scalar term
C     C:         The super-diagonal scalar term
C
C     On Exit:
C     A:         L(j,j-1) The supra-diagonal term of L 
C     B:         L(j,j)[-1]
C     C:         U(j,j+1) The super-diagonal term of U 
C
C
C     L(i,i-1) = A(i) therefore no operations need to be performed on A.
C     L(i,i)   = B(i) - L(i,i-1)*U(i-1,i)
C     U(i,i+1) = L[inverse](i,i)*C(i)
C
      DIMENSION A (NSIM,NPTS), 
     1          B (NSIM,NPTS), 
     2          C (NSIM,NPTS)
C
C     Form L(1,1) and U(1,2)
C
      J = 1
C
C     L(1,1) = B(1) => Leave B(1) unchanged
C     Form U(1,2) = L[inverse](1,1) * C(1)
C
      DO 200 I = 1, NSIM
         C(I,J) = C(I,J) / B(I,J)
  200 CONTINUE
C
C     Loop over J = 2 to NPTS-1 forming the LU terms
C
      DO 400 J = 2, NPTS-1
C
C     Form L(i,i) = B(i) - L(i,i-1)*U(i-1,i)
C     Form U(i,i+1) = L[inverse](i,i) * C(i) 
C
         DO 300 I = 1, NSIM
            B(I,J) = B(I,J) - A(I,J)*C(I,J-1)
            C(I,J) = C(I,J) / B(I,J)
  300    CONTINUE
  400 CONTINUE
C
C     Form L(N,N) = B(N) - L(N,N-1)*U(N-1,N)
C
      J = NPTS
      DO 500 I = 1, NSIM
         B(I,J) = B(I,J) - A(I,J)*C(I,J-1)
  500 CONTINUE
C
C     Finished forming LU
C
      RETURN
      END
