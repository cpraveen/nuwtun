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
      SUBROUTINE AF3DI  (IDIM, JDIM, KDIM, K,
     1                   Q, PROPS, NPPTS, S, VOL, DTJ, RES,
     2                   QI, PROPSI, STEMP, VOLI, EL, EM, TKE,
     3                   IBCDAT, NBCS, IFLUX, VISCOS, ITURB,
     4                   A, B, C, RT, RTEMP, RTEMP2,
     5                   EIGENP, EIGENM, EIGEN)
C
C     Routine to do the diagonalized I implicit sweep of three factor
C     Approximate Factorisation
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     K              : K plane to perform I implicit sweep on
C     Q              : Primitive variables at cell centers
C     PROPS          : Properties at cell centers
C                      PROPS(1) = molecular viscosity
C                      PROPS(2) = turbulent eddy viscosity 
C                      PROPS(3) = Y+
C     NPPTS          : Number of points at which PROPS are stored.  This allows
C                      the storage of PROPS for viscous flows and the avoidance
C                      of the additional storage for inviscid flows.
C     S              : Metrics
C     VOL            : Cell volumes
C     DTJ            : Time step divided by cell volume
C     RES            : Residual - must be over the entire field for
C                      this time integration scheme
C     QI             : Q at cell center for diagonalized scheme
C     PROPSI         : Properties at cell centers for eigenvalue routine
C     STEMP          : Storage for metrics passed to eigenvalue routine
C                      including the area term
C     VOLI           : Storage for volum,e passed to eigenvalue routine
C     EL, EM         : Unit tangency vectors
C     TKE            : Turbulent kinetic energy
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IFLUX          : Flux evaluation scheme
C     VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C     ITURB          : Turbulence modeling scheme
C     A,B,C          : Tridiagonal systems to invert
C                      A*DQ(I-1) + B*DQ(I) + C*DQ(I+1) = RTEMP
C                      Implicit in the I direction.  Vectorization
C                      by inverting multiple systems at one time to
C                      avoid the recursion of the inversion.  
C     RTEMP          : Storage for the residual array in the form 
C                      appropriate for the inversion routine
C     EIGENP,M       : Storage for eigenvalues(+,-)
C     EIGEN          : Storage for viscous eigenvalue
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          DTJ   (JDIM,KDIM,IDIM)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IDIM+1,NF),
     1          QI    (2:JDIM,1:IDIM+1,NQ), 
     2          PROPSI(2:JDIM,1:IDIM+1,NP),
     3          STEMP (2:JDIM,1:IDIM+1,4),
     4          VOLI  (2:JDIM,1:IDIM+1),
     5          EL    (2:JDIM,1:IDIM+1,3),
     6          EM    (2:JDIM,1:IDIM+1,3),
     7          TKE   (2:JDIM,1:IDIM+1)
C
C     Storage for the implicit system to invert.
C     Store the system only in the internal portion of the field and use 
C     explicit boundary conditions of the form Delta(Q) boundaries = 0.
C
      DIMENSION A     (2:JDIM,2:IDIM),
     1          B     (2:JDIM,2:IDIM),
     2          C     (2:JDIM,2:IDIM),
     3          RT    (2:JDIM,2:IDIM),
     4          RTEMP (2:JDIM,1:IDIM+1,NF),
     5          RTEMP2(2:JDIM,1:IDIM+1,NF)
C
C     Storage for the eigenvalues used to calculate A,B and C
C
      DIMENSION EIGENP(2:JDIM,1:IDIM+1,3),
     1          EIGENM(2:JDIM,1:IDIM+1,3),
     2          EIGEN (2:JDIM,1:IDIM+1,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      LOGICAL   VISCOS
C
C
C     Store Q in QI at cell centers
C
      DO 130 L = 1, NQ
         DO 120 I = 1, IDIM+1
            DO 110 J = 2, JDIM
               QI(J,I,L) = Q(J,K,I,L)
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C     Calculate the metrics for the at the cell centers for the eigenvalues 
C     When prevent division by zero, don't set S=<1,1,1> since TANVEC 
C     calculates the tangency vectors by crossing S with <1,1,1>.
C
      DO 180 I = 1, IDIM+1
         DO 170 J = 2, JDIM
            SX           = 0.5E0*( S(J,K,I  ,1,1)*S(J,K,I  ,4,1)
     1                           + S(J,K,I-1,1,1)*S(J,K,I-1,4,1) )
            SY           = 0.5E0*( S(J,K,I  ,2,1)*S(J,K,I  ,4,1)
     1                           + S(J,K,I-1,2,1)*S(J,K,I-1,4,1) )
            SZ           = 0.5E0*( S(J,K,I  ,3,1)*S(J,K,I  ,4,1)
     1                           + S(J,K,I-1,3,1)*S(J,K,I-1,4,1) )
            STEMP(J,I,4) = SQRT (SX*SX + SY*SY + SZ*SZ)
            IF (STEMP(J,I,4) .GT. 0.E0) THEN
               STEMP(J,I,1) = SX / STEMP(J,I,4)
               STEMP(J,I,2) = SY / STEMP(J,I,4)
               STEMP(J,I,3) = SZ / STEMP(J,I,4)
            ELSE
               STEMP(J,I,1) =   SHATX
               STEMP(J,I,2) = - SHATY
               STEMP(J,I,3) =   SHATZ
            ENDIF
  170    CONTINUE
  180 CONTINUE
C
C     Calculate unit tangency vectors
C
      NPTS = (IDIM + 1) * (JDIM - 1)
      CALL TANVEC (NPTS, STEMP, EL, EM)
C
C     Calculate RTEMP = T[-1] * RES
C
      DO 230 L = 1, NF
         DO 220 I = 1, IDIM+1
            DO 210 J = 2, JDIM
               RTEMP2(J,I,L) = RES(J,K,I,L)
  210       CONTINUE
  220    CONTINUE
  230 CONTINUE
C
      CALL GETTKE (ITURB, NPTS, QI, TKE)
      CALL TINVR  (ITURB, NPTS, QI, STEMP, EL, EM, TKE, RTEMP2, RTEMP)
C
C     Calculate eigenvalues 
C
      CALL EIGENI (ITURB, NPTS, QI, STEMP, EIGENP, EIGENM)
C
C     Viscous eigenvalues
C
      IF (VISCOS) THEN
         DO 320 L = 1, NP
            DO 310 I = 1, IDIM+1
               DO 300 J = 2, JDIM
                  PROPSI(J,I,L) = PROPS(J,K,I,L)
  300          CONTINUE
  310       CONTINUE
  320    CONTINUE
C
         DO 340 I = 1, IDIM+1
            DO 330 J = 2, JDIM
               VOLI(J,I) = VOL(J,K,I)
  330       CONTINUE
  340    CONTINUE
C
         CALL EIGENV (ITURB, NPTS, QI, PROPSI, STEMP, VOLI, EIGEN)
      ENDIF
C
C     Setup and invert each system
C
      DO 490 L = 1, NRANK
         IF (L .EQ. 4) THEN
            M = 2
         ELSE IF (L .EQ. 5) THEN
            M = 3
         ELSE
            M = 1
         ENDIF
C
C     Viscous eigenvalue
C
         IF (L .LE. 5) THEN
            MV = 1
         ELSE IF (L .EQ. NQ) THEN
            MV = 3
         ELSE
            MV = 2
         ENDIF
C
C     Fourth-order central scheme
C
         IF (FOURTH) THEN
            DO 410 I = 2, IDIM
               DO 400 J = 2, JDIM
                  DTTERM  = 1.0E0 /
     1                      ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                      + TDPDEL * VOL(J,K,I) / TIMEDT )
                  A (J,I) = - DTTERM*(EIGENP(J,I-1,M) + EIGENM(J,I-1,M))
                  B (J,I) =   1.0E0  
                  C (J,I) =   DTTERM*(EIGENP(J,I+1,M) + EIGENM(J,I+1,M))
                  RT(J,I) =   RTEMP(J,I,L)
  400          CONTINUE
  410       CONTINUE
C
C     Second-order MUSCL scheme
C
         ELSE
            DO 430 I = 2, IDIM
               DO 420 J = 2, JDIM
                  DTTERM  = 1.0E0 /
     1                      ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                      + TDPDEL * VOL(J,K,I) / TIMEDT )
                  A (J,I) = - DTTERM * EIGENP(J,I-1,M)
                  B (J,I) =   1.0E0  + DTTERM * 
     1                                 (EIGENP(J,I,M) - EIGENM(J,I,M))
                  C (J,I) =   DTTERM * EIGENM(J,I+1,M)
                  RT(J,I) =   RTEMP(J,I,L)
  420          CONTINUE
  430       CONTINUE
         ENDIF
C
C     Viscous contribution
C
         IF (VISCOS) THEN
            DO 450 I = 2, IDIM
               DO 440 J = 2, JDIM
                  DTTERM  = 1.0E0 /
     1                      ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                      + TDPDEL * VOL(J,K,I) / TIMEDT )
                  A (J,I) = A(J,I)-     DTTERM * EIGEN(J,I,MV)
                  B (J,I) = B(J,I)+2.E0*DTTERM * EIGEN(J,I,MV)
                  C (J,I) = C(J,I)-     DTTERM * EIGEN(J,I,MV)
  440          CONTINUE
  450       CONTINUE
         ENDIF
C
C     Invert the tridiagonal systems:
C     A*DQ(I-1) + B*DQ(I) + C*DQ(I+1) = RT
C
         NPTS = IDIM - 1
         NSIM = JDIM - 1
         CALL VSTRI (NPTS, NSIM, A, B, C, RT)
C
C     Store inversion back into RTEMP
C
         DO 470 I = 2, IDIM
            DO 460 J = 2, JDIM
               RTEMP(J,I,L) = RT(J,I)
  460       CONTINUE
  470    CONTINUE
C
C     Next equation
C
  490 CONTINUE
C
C     Calculate RTEMP = T * RTEMP
C
      NPTS = (IDIM + 1) * (JDIM - 1)
      CALL TR     (ITURB, NPTS, QI, STEMP, EL, EM, TKE, RTEMP, RTEMP2)
C
C     Update the residual array
C
      DO 830 L = 1, NF
         DO 820 I = 2, IDIM
            DO 810 J = 2, JDIM
               RES(J,K,I,L) = RTEMP2(J,I,L)
  810       CONTINUE
  820    CONTINUE
  830 CONTINUE
C
C     Reset the residual to zero for the case of a two-dimensional calculation
C     as roundoff error in the transform-solve-transform process can lead to 
C     non-zero residual.
C          For Reynolds stress calculations reset W, Tau_xz, Tau_yz
C              inviscid, laminar, two-equation reset W
C
      IF (.NOT. THREED) THEN
         IF (ITURB .EQ. ITRS) THEN
            DO 920 I = 2, IDIM
               DO 910 J = 2, JDIM
                  RES(J,K,I, 4) = 0.E0
                  RES(J,K,I,10) = 0.E0
                  RES(J,K,I,11) = 0.E0
  910          CONTINUE
  920       CONTINUE
         ELSE
            DO 950 I = 2, IDIM
               DO 940 J = 2, JDIM
                  RES(J,K,I, 4) = 0.E0
  940          CONTINUE
  950       CONTINUE
         ENDIF
      ENDIF
C
C     Finished with the I implicit sweep of diagonalized three factor AF
C
      RETURN
      END





      SUBROUTINE AF3DJ  (IDIM, JDIM, KDIM, K,
     1                   Q, PROPS, NPPTS, S, VOL, DTJ, RES,
     2                   QI, PROPSJ, STEMP, VOLJ, EL, EM, TKE, 
     3                   IBCDAT, NBCS, IFLUX, VISCOS, ITURB,
     4                   A, B, C, RT, RTEMP, RTEMP2, 
     5                   EIGENP, EIGENM, EIGEN)
C
C     Routine to do the diagonalized J implicit sweep of three factor
C     Approximate Factorisation
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     K              : K plane to perform I implicit sweep on
C     Q              : Primitive variables at cell centers
C     PROPS          : Properties at cell centers
C                      PROPS(1) = molecular viscosity
C                      PROPS(2) = turbulent eddy viscosity 
C                      PROPS(3) = Y+
C     NPPTS          : Number of points at which PROPS are stored.  This allows
C                      the storage of PROPS for viscous flows and the avoidance
C                      of the additional storage for inviscid flows.
C     S              : Metrics
C     VOL            : Cell volumes
C     DTJ            : Time step divided by cell volume
C     RES            : Residual - must be over the entire field for
C                      this time integration scheme
C     QI             : Q at cell center for diagonalized scheme
C     PROPSJ         : Properties at cell centers for eigenvalue routine
C     STEMP          : Storage for metrics passed to eigenvalue routine
C                      including the area term
C     VOLJ           : Storage for volum,e passed to eigenvalue routine
C     EL, EM         : Unit tangency vectors
C     TKE            : Turbulent kinetic energy
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IFLUX          : Flux evaluation scheme
C     VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C     ITURB          : Turbulence modeling scheme
C     A,B,C          : Tridiagonal systems to invert
C                      A*DQ(I-1) + B*DQ(I) + C*DQ(I+1) = RTEMP
C                      Implicit in the I direction.  Vectorization
C                      by inverting multiple systems at one time to
C                      avoid the recursion of the inversion.  
C     RTEMP          : Storage for the residual array in the form 
C                      appropriate for the inversion routine
C     EIGENP,M       : Storage for eigenvalues(+,-)
C     EIGEN          : Storage for viscous eigenvalue
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          DTJ   (JDIM,KDIM,IDIM)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IDIM+1,NF),
     1          QI    (2:IDIM,1:JDIM+1,NQ), 
     2          PROPSJ(2:IDIM,1:JDIM+1,NP), 
     3          STEMP (2:IDIM,1:JDIM+1,4),
     4          VOLJ  (2:IDIM,1:JDIM+1),
     5          EL    (2:IDIM,1:JDIM+1,3),
     6          EM    (2:IDIM,1:JDIM+1,3),
     7          TKE   (2:IDIM,1:JDIM+1)
C
C     Storage for the implicit system to invert.
C     Store the system only in the internal portion of the field and use 
C     explicit boundary conditions of the form Delta(Q) boundaries = 0.
C
      DIMENSION A     (2:IDIM,2:JDIM),
     1          B     (2:IDIM,2:JDIM),
     2          C     (2:IDIM,2:JDIM),
     3          RT    (2:IDIM,2:JDIM),
     4          RTEMP (2:IDIM,1:JDIM+1,NF),
     5          RTEMP2(2:IDIM,1:JDIM+1,NF)
C
C     Storage for the eigenvalues used to calculate A,B and C
C
      DIMENSION EIGENP(2:IDIM,1:JDIM+1,3),
     1          EIGENM(2:IDIM,1:JDIM+1,3),
     2          EIGEN (2:IDIM,1:JDIM+1,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      LOGICAL   VISCOS
C
C
C     Store Q in QI at cell centers
C
      DO 130 L = 1, NQ
         DO 120 J = 1, JDIM+1
            DO 110 I = 2, IDIM
               QI(I,J,L) = Q(J,K,I,L)
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C     Calculate the metrics for the at the cell centers for the eigenvalues 
C     When prevent division by zero, don't set S=<1,1,1> since TANVEC 
C     calculates the tangency vectors by crossing S with <1,1,1>.
C
      DO 150 J = 1, JDIM+1
         DO 140 I = 2, IDIM
            SX           = 0.5E0*( S(J  ,K,I,1,2)*S(J  ,K,I,4,2)
     1                           + S(J-1,K,I,1,2)*S(J-1,K,I,4,2) )
            SY           = 0.5E0*( S(J  ,K,I,2,2)*S(J  ,K,I,4,2)
     1                           + S(J-1,K,I,2,2)*S(J-1,K,I,4,2) )
            SZ           = 0.5E0*( S(J  ,K,I,3,2)*S(J  ,K,I,4,2)
     1                           + S(J-1,K,I,3,2)*S(J-1,K,I,4,2) )
            STEMP(I,J,4) = SQRT (SX*SX + SY*SY + SZ*SZ)
            IF (STEMP(I,J,4) .GT. 0.E0) THEN
               STEMP(I,J,1) = SX / STEMP(I,J,4)
               STEMP(I,J,2) = SY / STEMP(I,J,4)
               STEMP(I,J,3) = SZ / STEMP(I,J,4)
            ELSE
               STEMP(I,J,1) =   SHATX
               STEMP(I,J,2) = - SHATY
               STEMP(I,J,3) =   SHATZ
            ENDIF
  140    CONTINUE
  150 CONTINUE
C
C     Calculate unit tangency vectors
C
      NPTS = (JDIM + 1) * (IDIM - 1)
      CALL TANVEC (NPTS, STEMP, EL, EM)
C
C     Calculate RTEMP = T[-1] * RES
C
      DO 230 L = 1, NF
         DO 220 J = 1, JDIM+1
            DO 210 I = 2, IDIM
               RTEMP2(I,J,L) = RES(J,K,I,L)
  210       CONTINUE
  220    CONTINUE
  230 CONTINUE
C
      CALL GETTKE (ITURB, NPTS, QI, TKE)
      CALL TINVR  (ITURB, NPTS, QI, STEMP, EL, EM, TKE, RTEMP2, RTEMP)
C
C     Calculate eigenvalues 
C
      CALL EIGENI (ITURB, NPTS, QI, STEMP, EIGENP, EIGENM)
C
C     Viscous eigenvalues
C
      IF (VISCOS) THEN
         DO 320 L = 1, NP
            DO 310 J = 1, JDIM+1
               DO 300 I = 2, IDIM
                  PROPSJ(I,J,L) = PROPS(J,K,I,L)
  300          CONTINUE
  310       CONTINUE
  320    CONTINUE
C
         DO 340 J = 1, JDIM+1
            DO 330 I = 2, IDIM
               VOLJ(I,J) = VOL(J,K,I)
  330       CONTINUE
  340    CONTINUE
C
         CALL EIGENV (ITURB, NPTS, QI, PROPSJ, STEMP, VOLJ, EIGEN)
      ENDIF
C
C     Setup and invert each system
C
      DO 490 L = 1, NRANK
         IF (L .EQ. 4) THEN
            M = 2
         ELSE IF (L .EQ. 5) THEN
            M = 3
         ELSE
            M = 1
         ENDIF
C
C     Viscous eigenvalue
C
         IF (L .LE. 5) THEN
            MV = 1
         ELSE IF (L .EQ. NQ) THEN
            MV = 3
         ELSE
            MV = 2
         ENDIF
C
C     Fourth-order central scheme
C
         IF (FOURTH) THEN
            DO 410 J = 2, JDIM
               DO 400 I = 2, IDIM
                  DTTERM  = 1.0E0 /
     1                      ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                      + TDPDEL * VOL(J,K,I) / TIMEDT )
                  A (I,J) = - DTTERM*(EIGENP(I,J-1,M) + EIGENM(I,J-1,M))
                  B (I,J) =   1.0E0  
                  C (I,J) =   DTTERM*(EIGENP(I,J+1,M) + EIGENM(I,J+1,M))
                  RT(I,J) =   RTEMP(I,J,L)
  400          CONTINUE
  410       CONTINUE
C
C     Second-order MUSCL scheme
C
         ELSE
            DO 430 J = 2, JDIM
               DO 420 I = 2, IDIM
                  DTTERM  = 1.0E0 /
     1                      ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                      + TDPDEL * VOL(J,K,I) / TIMEDT )
                  A (I,J) = - DTTERM * EIGENP(I,J-1,M)
                  B (I,J) =   1.0E0  + DTTERM *
     1                                 (EIGENP(I,J,M) - EIGENM(I,J,M))
                  C (I,J) =   DTTERM * EIGENM(I,J+1,M)
                  RT(I,J) =   RTEMP(I,J,L)
  420          CONTINUE
  430       CONTINUE
         ENDIF
C
C     Viscous contribution
C
         IF (VISCOS) THEN
            DO 450 J = 2, JDIM
               DO 440 I = 2, IDIM
                  DTTERM  = 1.0E0 /
     1                      ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                      + TDPDEL * VOL(J,K,I) / TIMEDT )
                  A (I,J) = A(I,J)-     DTTERM * EIGEN(I,J,MV)
                  B (I,J) = B(I,J)+2.E0*DTTERM * EIGEN(I,J,MV)
                  C (I,J) = C(I,J)-     DTTERM * EIGEN(I,J,MV)
  440          CONTINUE
  450       CONTINUE
         ENDIF
C
C     Invert the tridiagonal systems:
C     A*DQ(I-1) + B*DQ(I) + C*DQ(I+1) = RT
C
         NPTS = JDIM - 1
         NSIM = IDIM - 1
         CALL VSTRI (NPTS, NSIM, A, B, C, RT)
C
C     Store inversion back into RTEMP
C
         DO 470 J = 2, JDIM
            DO 460 I = 2, IDIM
               RTEMP(I,J,L) = RT(I,J)
  460       CONTINUE
  470    CONTINUE
C
C     Next equation
C
  490 CONTINUE
C
C     Calculate RTEMP = T * RTEMP
C
      NPTS = (JDIM + 1) * (IDIM - 1)
      CALL TR     (ITURB, NPTS, QI, STEMP, EL, EM, TKE, RTEMP, RTEMP2)
C
C     Update the residual array
C
      DO 830 L = 1, NF
         DO 820 J = 2, JDIM
            DO 810 I = 2, IDIM
               RES(J,K,I,L) = RTEMP2(I,J,L)
  810       CONTINUE
  820    CONTINUE
  830 CONTINUE
C
C     Reset the residual to zero for the case of a two-dimensional calculation
C     as roundoff error in the transform-solve-transform process can lead to 
C     non-zero residual.
C          For Reynolds stress calculations reset W, Tau_xz, Tau_yz
C              inviscid, laminar, two-equation reset W
C
      IF (.NOT. THREED) THEN
         IF (ITURB .EQ. ITRS) THEN
            DO 920 I = 2, IDIM
               DO 910 J = 2, JDIM
                  RES(J,K,I, 4) = 0.E0
                  RES(J,K,I,10) = 0.E0
                  RES(J,K,I,11) = 0.E0
  910          CONTINUE
  920       CONTINUE
         ELSE
            DO 950 I = 2, IDIM
               DO 940 J = 2, JDIM
                  RES(J,K,I, 4) = 0.E0
  940          CONTINUE
  950       CONTINUE
         ENDIF
      ENDIF
C
C     Finished with the J implicit sweep of diagonalized three factor AF
C
      RETURN
      END





      SUBROUTINE AF3DK  (IDIM, JDIM, KDIM, I,
     1                   Q, PROPS, NPPTS, S, VOL, DTJ, RES,
     2                   QI, PROPSK, STEMP, VOLK, EL, EM, TKE, 
     3                   IBCDAT, NBCS, IFLUX, VISCOS, ITURB,
     4                   A, B, C, RT, RTEMP, RTEMP2, 
     5                   EIGENP, EIGENM, EIGEN)
C
C     Routine to do the diagonalized K implicit sweep of three factor
C     Approximate Factorisation
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     K              : K plane to perform I implicit sweep on
C     Q              : Primitive variables at cell centers
C     PROPS          : Properties at cell centers
C                      PROPS(1) = molecular viscosity
C                      PROPS(2) = turbulent eddy viscosity 
C                      PROPS(3) = Y+
C     NPPTS          : Number of points at which PROPS are stored.  This allows
C                      the storage of PROPS for viscous flows and the avoidance
C                      of the additional storage for inviscid flows.
C     S              : Metrics
C     VOL            : Cell volumes
C     DTJ            : Time step divided by cell volume
C     RES            : Residual - must be over the entire field for
C                      this time integration scheme
C     QI             : Q at cell center for diagonalized scheme
C     PROPSK         : Properties at cell centers for eigenvalue routine
C     STEMP          : Storage for metrics passed to eigenvalue routine
C                      including the area term
C     VOLK           : Storage for volume passed to eigenvalue routine
C     EL, EM         : Unit tangency vectors
C     TKE            : Turbulent kinetic energy
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IFLUX          : Flux evaluation scheme
C     VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C     ITURB          : Turbulence modeling scheme
C     A,B,C          : Tridiagonal systems to invert
C                      A*DQ(I-1) + B*DQ(I) + C*DQ(I+1) = RTEMP
C                      Implicit in the I direction.  Vectorization
C                      by inverting multiple systems at one time to
C                      avoid the recursion of the inversion.  
C     RTEMP          : Storage for the residual array in the form 
C                      appropriate for the inversion routine
C     EIGENP,M       : Storage for eigenvalues(+,-)
C     EIGEN          : Storage for viscous eigenvalue
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          DTJ   (JDIM,KDIM,IDIM)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IDIM+1,NF),
     1          QI    (2:JDIM,1:KDIM+1,NQ), 
     2          PROPSK(2:JDIM,1:KDIM+1,NP), 
     3          STEMP (2:JDIM,1:KDIM+1,4),
     4          VOLK  (2:JDIM,1:KDIM+1),
     5          EL    (2:JDIM,1:KDIM+1,3),
     6          EM    (2:JDIM,1:KDIM+1,3),
     7          TKE   (2:JDIM,1:KDIM+1)
C
C     Storage for the implicit system to invert.
C     Store the system only in the internal portion of the field and use 
C     explicit boundary conditions of the form Delta(Q) boundaries = 0.
C
      DIMENSION A     (2:JDIM,2:KDIM),
     1          B     (2:JDIM,2:KDIM),
     2          C     (2:JDIM,2:KDIM),
     3          RT    (2:JDIM,2:KDIM),
     4          RTEMP (2:JDIM,1:KDIM+1,NF),
     5          RTEMP2(2:JDIM,1:KDIM+1,NF)
C
C     Storage for the eigenvalues used to calculate A,B and C
C
      DIMENSION EIGENP(2:JDIM,1:KDIM+1,3),
     1          EIGENM(2:JDIM,1:KDIM+1,3),
     2          EIGEN (2:JDIM,1:KDIM+1,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      LOGICAL   VISCOS
C
C
C     Store Q in QI at cell centers
C
      DO 130 L = 1, NQ
         DO 120 K = 1, KDIM+1
            DO 110 J = 2, JDIM
               QI(J,K,L) = Q(J,K,I,L)
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C     Calculate the metrics for the at the cell centers for the eigenvalues 
C     When prevent division by zero, don't set S=<1,1,1> since TANVEC 
C     calculates the tangency vectors by crossing S with <1,1,1>.
C
      DO 150 K = 1, KDIM+1
         DO 140 J = 2, JDIM
            SX           = 0.5E0*( S(J,K  ,I,1,3)*S(J,K  ,I,4,3)
     1                           + S(J,K-1,I,1,3)*S(J,K-1,I,4,3) )
            SY           = 0.5E0*( S(J,K  ,I,2,3)*S(J,K  ,I,4,3)
     1                           + S(J,K-1,I,2,3)*S(J,K-1,I,4,3) )
            SZ           = 0.5E0*( S(J,K  ,I,3,3)*S(J,K  ,I,4,3)
     1                           + S(J,K-1,I,3,3)*S(J,K-1,I,4,3) )
            STEMP(J,K,4) = SQRT (SX*SX + SY*SY + SZ*SZ)
            IF (STEMP(J,K,4) .GT. 0.E0) THEN
               STEMP(J,K,1) = SX / STEMP(J,K,4)
               STEMP(J,K,2) = SY / STEMP(J,K,4)
               STEMP(J,K,3) = SZ / STEMP(J,K,4)
            ELSE
               STEMP(J,K,1) =   SHATX
               STEMP(J,K,2) = - SHATY
               STEMP(J,K,3) =   SHATZ
            ENDIF
  140    CONTINUE
  150 CONTINUE
C
C     Calculate unit tangency vectors
C
      NPTS = (KDIM + 1) * (JDIM - 1)
      CALL TANVEC (NPTS, STEMP, EL, EM)
C
C     Calculate RTEMP = T[-1] * RES
C
      DO 230 L = 1, NF
         DO 220 K = 1, KDIM+1
            DO 210 J = 2, JDIM
               RTEMP2(J,K,L) = RES(J,K,I,L)
  210       CONTINUE
  220    CONTINUE
  230 CONTINUE
C
      CALL GETTKE (ITURB, NPTS, QI, TKE)
      CALL TINVR  (ITURB, NPTS, QI, STEMP, EL, EM, TKE, RTEMP2, RTEMP)
C
C     Calculate eigenvalues 
C
      CALL EIGENI (ITURB, NPTS, QI, STEMP, EIGENP, EIGENM)
C
C     Viscous eigenvalues
C
      IF (VISCOS) THEN
         DO 320 L = 1, NP
            DO 310 K = 1, KDIM+1
               DO 300 J = 2, JDIM
                  PROPSK(J,K,L) = PROPS(J,K,I,L)
  300          CONTINUE
  310       CONTINUE
  320    CONTINUE
C
         DO 340 K = 1, KDIM+1
            DO 330 J = 2, JDIM
               VOLK(J,K) = VOL(J,K,I)
  330       CONTINUE
  340    CONTINUE
C
         CALL EIGENV (ITURB, NPTS, QI, PROPSK, STEMP, VOLK, EIGEN)
      ENDIF
C
C     Setup and invert each system
C
      DO 490 L = 1, NRANK
         IF (L .EQ. 4) THEN
            M = 2
         ELSE IF (L .EQ. 5) THEN
            M = 3
         ELSE
            M = 1
         ENDIF
C
C     Viscous eigenvalue
C
         IF (L .LE. 5) THEN
            MV = 1
         ELSE IF (L .EQ. NQ) THEN
            MV = 3
         ELSE
            MV = 2
         ENDIF
C
C     Fourth-order central scheme
C
         IF (FOURTH) THEN
            DO 410 K = 2, KDIM
               DO 400 J = 2, JDIM
                  DTTERM  = 1.0E0 /
     1                      ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                      + TDPDEL * VOL(J,K,I) / TIMEDT )
                  A (J,K) = - DTTERM*(EIGENP(J,K-1,M) + EIGENM(J,K-1,M))
                  B (J,K) =   1.0E0  
                  C (J,K) =   DTTERM*(EIGENP(J,K+1,M) + EIGENM(J,K+1,M))
                  RT(J,K) =   RTEMP(J,K,L)
  400          CONTINUE
  410       CONTINUE
C
C    Second-order MUSCL scheme
C
         ELSE
            DO 430 K = 2, KDIM
               DO 420 J = 2, JDIM
                  DTTERM  = 1.0E0 /
     1                      ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                      + TDPDEL * VOL(J,K,I) / TIMEDT )
                  A (J,K) = - DTTERM * EIGENP(J,K-1,M)
                  B (J,K) =   1.0E0  + DTTERM *
     1                                 (EIGENP(J,K,M) - EIGENM(J,K,M))
                  C (J,K) =   DTTERM * EIGENM(J,K+1,M)
                  RT(J,K) =   RTEMP(J,K,L)
  420          CONTINUE
  430       CONTINUE
         ENDIF
C
C     Viscous contribution
C
         IF (VISCOS) THEN
            DO 450 K = 2, KDIM
               DO 440 J = 2, JDIM
                  DTTERM  = 1.0E0 /
     1                      ( TDPDTJ / (TDTHET*DTJ(J,K,I))
     2                      + TDPDEL * VOL(J,K,I) / TIMEDT )
                  A (J,K) = A(J,K)-     DTTERM * EIGEN(J,K,MV)
                  B (J,K) = B(J,K)+2.E0*DTTERM * EIGEN(J,K,MV)
                  C (J,K) = C(J,K)-     DTTERM * EIGEN(J,K,MV)
  440          CONTINUE
  450       CONTINUE
         ENDIF
C
C     Invert the tridiagonal systems:
C     A*DQ(K-1) + B*DQ(K) + C*DQ(K+1) = RT
C
         NPTS = KDIM - 1
         NSIM = JDIM - 1
         CALL VSTRI (NPTS, NSIM, A, B, C, RT)
C
C     Store inversion back into RTEMP
C
         DO 470 K = 2, KDIM
            DO 460 J = 2, JDIM
               RTEMP(J,K,L) = RT(J,K)
  460       CONTINUE
  470    CONTINUE
C
C     Next equation
C
  490 CONTINUE
C
C     Calculate RTEMP = T * RTEMP
C
      NPTS = (KDIM + 1) * (JDIM - 1)
      CALL TR     (ITURB, NPTS, QI, STEMP, EL, EM, TKE, RTEMP, RTEMP2)
C
C     Update the residual array
C
      DO 830 L = 1, NF
         DO 820 K = 2, KDIM
            DO 810 J = 2, JDIM
               RES(J,K,I,L) = RTEMP2(J,K,L)
  810       CONTINUE
  820    CONTINUE
  830 CONTINUE
C
C     The reset of non-zero residuals caused by roundoff of transform-solve-
C     transform procedure for the two-dimensional cases is not needed here
C     as this routine is not called if doing a 2D case!
C
C     Finished with the K implicit sweep of diagonalized three factor AF
C
      RETURN
      END





      SUBROUTINE EIGENI (ITURB, NPTS, Q, S, EIGENP, EIGENM)
C
C     Routine to calculate the inviscid eigenvalues
C
C     ITURB          : Turbulence model
C     NPTS           : Size of array to calculate on
C     Q              : Flow variables
C     S              : Metrics - including area term
C     EIGENP,EIGENM  : Eigenvalues (+,-)
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          S     (NPTS,4),
     2          EIGENP(NPTS,3),
     3          EIGENM(NPTS,3)
C
C     Calculate the Eigenvalues
C
      DO 100 I = 1, NPTS
         RHO   = Q(I,1)
         U     = Q(I,2)
         V     = Q(I,3)
         W     = Q(I,4)
         P     = Q(I,5)
C
         A     = SQRT ( GAMMA * P / RHO )
C
         AREA  = S(I,4)
C
         UBAR        = (S(I,1) * U + S(I,2) * V + S(I,3) * W) * AREA
         UBARPA      = (UBAR + A * AREA)
         UBARMA      = (UBAR - A * AREA)
C
         EIGENP(I,1) = 0.5E0 * (UBAR   + ABS(UBAR  ))
         EIGENP(I,2) = 0.5E0 * (UBARPA + ABS(UBARPA))
         EIGENP(I,3) = 0.5E0 * (UBARMA + ABS(UBARMA))
C
         EIGENM(I,1) = 0.5E0 * (UBAR   - ABS(UBAR  ))
         EIGENM(I,2) = 0.5E0 * (UBARPA - ABS(UBARPA))
         EIGENM(I,3) = 0.5E0 * (UBARMA - ABS(UBARMA))
c     write (iout,'(a,i4,6e15.7)') ' EIGENI: i ', i, 
c    1                            (eigenp(i,ll),ll=1,3),
c    2                            (eigenm(i,lk),lk=1,3)
C
  100 CONTINUE
C
      RETURN
      END





      SUBROUTINE EIGENV (ITURB, NPTS, Q, PROPS, S, VOL, EIGEN)
C
C     Routine to calculate the viscous eigenvalues
C
C     VISCOS         : Viscous switch
C     ITURB          : Turbulence model
C     NPTS           : Size of array to calculate on
C     Q              : Flow variables
C     PROPS          : Properties array
C     S              : Metrics - including area term
C     VOL            : Cell volume
C     EIGEN          : Maximum viscous eigenvalues
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          PROPS (NPTS,NP),
     2          S     (NPTS,4),
     3          VOL   (NPTS),
     4          EIGEN (NPTS,3)
C
C     Calculate the Mean Flow Eigenvalue
C
      DO 100 I = 1, NPTS
         RHO      = Q(I,1)
C
         RMU      = PROPS(I,1)
         RMUT     = PROPS(I,2)
C
         AREA     = S(I,4)
C
         EIGEN(I,1) = AREA * AREA / (RHO*VOL(I) + RSMASQ) * FSMACH/RE *
     1            MAX ( 4.E0/3.E0*(RMU+RMUT), GAMMA*(RMU/PR+RMUT/PRT) )
c     write (iout,'(a,i4,e15.7)') ' EIGENV: i ', i, eigen(i)
C
  100 CONTINUE
C
C     Calculate the Additional Eigenvalues
C     II-b. Two-Equation Models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 200 I = 1, NPTS
            RHO        = Q(I,1)
C
            RMU        = PROPS(I,1)
            RMUT       = PROPS(I,2)
C
            AREA       = S(I,4)
C
            EIGEN(I,2) = AREA * AREA / (RHO*VOL(I)+RSMASQ) * FSMACH/RE * 
     1                   (RMU + SIGK  * RMUT) 
            EIGEN(I,3) = AREA * AREA / (RHO*VOL(I)+RSMASQ) * FSMACH/RE *
     1                   (RMU + SIGT2 * RMUT)
C
  200    CONTINUE
C
C     II-c. Reynolds Stress Models
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         T23 = 2.E0 / 3.E0
         DO 210 I = 1, NPTS
            RHO        = Q(I, 1)
            TAUXX      = Q(I, 6)
            TAUYY      = Q(I, 7)
            TAUZZ      = Q(I, 8)
            EPSLN      = Q(I,12)
C
            TKE        = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
            RMU        = PROPS(I,1)
            RMUT       = PROPS(I,2)
C
            AREA       = S(I,4)
C
            TMP        = RHO * TKE / EPSLN *
     1             MAX( (S(I,1)*TAUXX+S(I,2)*TAUYY+S(I,3)*TAUZZ), 0.E0 )
            EIGEN(I,2) = AREA * AREA / (RHO*VOL(I)+RSMASQ) * 
     1                   ( RMU*FSMACH/RE + CSUBS *TMP )
            EIGEN(I,3) = AREA * AREA / (RHO*VOL(I)+RSMASQ) * 
     1                   ( RMU*FSMACH/RE + CEPSLN*TMP )
C
  210    CONTINUE
      ENDIF
C
      RETURN
      END





      SUBROUTINE TINVR  (ITURB, NPTS, Q, S, EL, EM, TKE, R, TIR)
C
C     Routine to calculate T[-1] R.  This routine is done in two steps:
C          I.  Calculate T[-1] R for mean-flow equations
C          II. Calculate T[-1] R for additional equations
C              a. Inviscid, Laminar, Algebraic turb. => None
C              b. Two-Equation turb. model           => k, phi equations
C              c. Reynolds stress model              => tau_ij, epsilon
C
C     ITURB          : Turbulence model
C     NPTS           : Size of array to calculate on
C     Q              : Flow variables
C     S              : Metrics - including area term
C     EL, EM         : Unit tangency vectors
C     TKE            : Tubulent kinetic energy
C     R              : Residual to be multiplied by T[-1]
C     TIR            : T[-1] R
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          S     (NPTS,4),
     2          EL    (NPTS,3),
     3          EM    (NPTS,3),
     4          TKE   (NPTS),
     5          R     (NPTS,NF),
     6          TIR   (NPTS,NF)
C
C     I. Calculate the T[-1] R for the mean-flow equations
C
      DO 100 I = 1, NPTS
         RHO   = Q(I,1)
         U     = Q(I,2)
         V     = Q(I,3)
         W     = Q(I,4)
         P     = Q(I,5)
C
         QSQ   = U*U + V*V + W*W
         ASQ   = GAMMA * P / RHO
         A     = SQRT ( ASQ )
C
         UBAR  = S (I,1) * U + S (I,2) * V + S (I,3) * W
         VBAR  = EL(I,1) * U + EL(I,2) * V + EL(I,3) * W
         WBAR  = EM(I,1) * U + EM(I,2) * V + EM(I,3) * W
C
         CAPG  = GAMM1 / ASQ
C
         TIR(I, 1) = (1.E0 - 0.5E0*CAPG*QSQ)*R(I,1)
     1             + CAPG * (U*R(I,2) + V*R(I,3) + W*R(I,4) - R(I,5))
         TIR(I, 2) = -VBAR*R(I,1) + EL(I,1)*R(I,2) + EL(I,2)*R(I,3) 
     1                            + EL(I,3)*R(I,4)
         TIR(I, 3) = -WBAR*R(I,1) + EM(I,1)*R(I,2) + EM(I,2)*R(I,3) 
     1                            + EM(I,3)*R(I,4)
         TIR(I, 4) = 0.5E0 * ( (0.5E0*CAPG*QSQ-UBAR/A)*R(I,1)
     1                       + (S(I,1)/A-CAPG*U)*R(I,2) 
     2                       + (S(I,2)/A-CAPG*V)*R(I,3)
     3                       + (S(I,3)/A-CAPG*W)*R(I,4) + CAPG*R(I,5) )
         TIR(I, 5) = 0.5E0 * ( (0.5E0*CAPG*QSQ+UBAR/A)*R(I,1)
     1                       - (S(I,1)/A+CAPG*U)*R(I,2) 
     2                       - (S(I,2)/A+CAPG*V)*R(I,3)
     3                       - (S(I,3)/A+CAPG*W)*R(I,4) + CAPG*R(I,5) )
  100 CONTINUE
C
C     II.   Additional equations for T[-1] R
C
C     II-b. Two-Equation model
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 210 I = 1, NPTS
            RHO   = Q(I, 1)
            P     = Q(I, 5)
            EPSLN = Q(I, 7)
C
            ASQ   = GAMMA * P / RHO
            A     = SQRT ( ASQ )
C
            CAPG  = GAMM1 / ASQ
C
            TIR(I, 1) =  TIR(I, 1) +         CAPG * R(I, 6)
            TIR(I, 4) =  TIR(I, 4) - 0.5E0 * CAPG * R(I, 6)
            TIR(I, 5) =  TIR(I, 5) - 0.5E0 * CAPG * R(I, 6)
            TIR(I, 6) = -TKE(I)*R(I, 1) + R(I, 6)
            TIR(I, 7) = -EPSLN *R(I, 1) + R(I, 7)
  210    CONTINUE
C
C     II-c. Reynolds Stress model
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 220 I = 1, NPTS
            RHO   = Q(I, 1)
            P     = Q(I, 5)
            TAUXX = Q(I, 6)
            TAUYY = Q(I, 7)
            TAUZZ = Q(I, 8)
            TAUXY = Q(I, 9)
            TAUXZ = Q(I,10)
            TAUYZ = Q(I,11)
            EPSLN = Q(I,12)
C
            ASQ   = GAMMA * P / RHO
            A     = SQRT ( ASQ )
C
            CAPG  = GAMM1 / ASQ
C
            TRACE     =  R(I,6) + R(I,7) + R(I,8)
            TIR(I, 1) =  TIR(I, 1) + 0.50E0 * CAPG * TRACE
            TIR(I, 4) =  TIR(I, 4) - 0.25E0 * CAPG * TRACE
            TIR(I, 5) =  TIR(I, 5) - 0.25E0 * CAPG * TRACE
            TIR(I, 6) = -TAUXX*R(I, 1) + R(I, 6)
            TIR(I, 7) = -TAUYY*R(I, 1) + R(I, 7)
            TIR(I, 8) = -TAUZZ*R(I, 1) + R(I, 8)
            TIR(I, 9) = -TAUXY*R(I, 1) + R(I, 9)
            TIR(I,10) = -TAUXZ*R(I, 1) + R(I,10)
            TIR(I,11) = -TAUYZ*R(I, 1) + R(I,11)
            TIR(I,12) = -EPSLN*R(I, 1) + R(I,12)
  220    CONTINUE
      ENDIF
C
      RETURN
      END





      SUBROUTINE TR     (ITURB, NPTS, Q, S, EL, EM, TKE, R, TTR)
C
C     Routine to calculate T R.  This routine is done in two steps:
C          I.  Calculate T R for mean-flow equations
C          II. Calculate T R for additional equations
C              a. Inviscid, Laminar, Algebraic turb. => None
C              b. Two-Equation turb. model           => k, phi equations
C              c. Reynolds stress model              => tau_ij, epsilon
C
C     ITURB          : Turbulence model
C     NPTS           : Size of array to calculate on
C     Q              : Flow variables
C     S              : Metrics - including area term
C     EL, EM         : Unit tangency vectors
C     TKE            : Tubulent kinetic energy
C     R              : Residual to be multiplied by T[-1]
C     TTR            : T R (T times R)
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          S     (NPTS,4),
     2          EL    (NPTS,3),
     3          EM    (NPTS,3),
     4          TKE   (NPTS),
     5          R     (NPTS,NF),
     6          TTR   (NPTS,NF)
C
C     I. Calculate the T R for the mean-flow equations
C
      DO 100 I = 1, NPTS
         RHO   = Q(I,1)
         U     = Q(I,2)
         V     = Q(I,3)
         W     = Q(I,4)
         P     = Q(I,5)
C
         H     = GAMMA * P / (RHO * GAMM1)
     1            + 0.5E0 * (U*U + V*V + W*W) + TKE(I)
C
         QSQ   = U*U + V*V + W*W
         ASQ   = GAMMA * P / RHO
         A     = SQRT ( ASQ )
C
         UBAR  = S (I,1) * U + S (I,2) * V + S (I,3) * W
         VBAR  = EL(I,1) * U + EL(I,2) * V + EL(I,3) * W
         WBAR  = EM(I,1) * U + EM(I,2) * V + EM(I,3) * W
C
         CAPG  = GAMM1 / ASQ
C
         TTR(I, 1) = R(I,1) + R(I,4) + R(I,5)
         TTR(I, 2) = U*R(I,1) + EL(I,1)*R(I,2) + EM(I,1)*R(I,3) 
     1               + (U+S(I,1)*A)*R(I,4) + (U-S(I,1)*A)*R(I,5)
         TTR(I, 3) = V*R(I,1) + EL(I,2)*R(I,2) + EM(I,2)*R(I,3) 
     1               + (V+S(I,2)*A)*R(I,4) + (V-S(I,2)*A)*R(I,5)
         TTR(I, 4) = W*R(I,1) + EL(I,3)*R(I,2) + EM(I,3)*R(I,3) 
     1               + (W+S(I,3)*A)*R(I,4) + (W-S(I,3)*A)*R(I,5)
         TTR(I, 5) = (0.5E0*QSQ+TKE(I))*R(I,1) 
     1               + VBAR*R(I,2) + WBAR*R(I,3)
     2               + (H+UBAR*A)*R(I,4) + (H-UBAR*A)*R(I,5)
  100 CONTINUE
C
C     II.   Additional equations for T R
C
C     II-b. Two-Equation model
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 210 I = 1, NPTS
            RHO   = Q(I, 1)
            P     = Q(I, 5)
            EPSLN = Q(I, 7)
C
            ASQ   = GAMMA * P / RHO
            A     = SQRT ( ASQ )
C
            CAPG  = GAMM1 / ASQ
C
            TTR(I, 5) = TTR(I, 5) + R(I, 6)
            TTR(I, 6) = TKE(I)*(R(I, 1) + R(I, 4) + R(I, 5)) + R(I, 6)
            TTR(I, 7) = EPSLN *(R(I, 1) + R(I, 4) + R(I, 5)) + R(I, 7)
  210    CONTINUE
C
C     II-c. Reynolds Stress model
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 220 I = 1, NPTS
            TAUXX = Q(I, 6)
            TAUYY = Q(I, 7)
            TAUZZ = Q(I, 8)
            TAUXY = Q(I, 9)
            TAUXZ = Q(I,10)
            TAUYZ = Q(I,11)
            EPSLN = Q(I,12)
C
            TERM      = R(I, 1) + R(I, 4) + R(I, 5)
            TRACE     = R(I, 6) + R(I, 7) + R(I, 8)
            TTR(I, 5) = TTR(I, 5) + 0.5E0 * TRACE
            TTR(I, 6) = TAUXX*TERM + R(I, 6)
            TTR(I, 7) = TAUYY*TERM + R(I, 7)
            TTR(I, 8) = TAUZZ*TERM + R(I, 8)
            TTR(I, 9) = TAUXY*TERM + R(I, 9)
            TTR(I,10) = TAUXZ*TERM + R(I,10)
            TTR(I,11) = TAUYZ*TERM + R(I,11)
            TTR(I,12) = EPSLN*TERM + R(I,12)
  220    CONTINUE
      ENDIF
C
      RETURN
      END






      SUBROUTINE TANVEC (NPTS, S, EL, EM)
C
C     Routine to calculate the unit tangency vectors lhat (EL) and mhat (EM)
C     such that S dot EL = 0 and S dot EM = 0.
C     Formula adapted from Thomas based on vector identity that
C     A . B X C = - C . B X A so that lhat = nhat X (1, 1, 1)
C     and then mhat = nhat X lhat.
C     This method will fail if nhat = (1, 1, 1) and then a correction will
C     need to be applied.
C
C     NPTS           : Size of array to calculate on
C     S              : Metrics - including area term
C     EL             : Unit tangency vector
C     EM             : Unit tangency vector
C
      include '../header/common.h'
C
      DIMENSION S     (NPTS,4),
     1          EL    (NPTS,3),
     2          EM    (NPTS,3)
C
C     Loop over all points and calculate tangency vectors
C
      DO 100 I = 1, NPTS
C
C     Calculate lhat = nhat X (1, 1, 1)
C
         ELX     = S(I,2) - S(I,3)
         ELY     = S(I,3) - S(I,1)
         ELZ     = S(I,1) - S(I,2)
         ELNORM  = SQRT (ELX*ELX + ELY*ELY + ELZ*ELZ)
c->      if (elnorm .le. 0.e0) write (iout,'(a,i5,3e15.7)')
c->  1       ' lhat zero in TANVEC i,lhat ', i, elx, ely, elz
         EL(I,1) = ELX / ELNORM
         EL(I,2) = ELY / ELNORM
         EL(I,3) = ELZ / ELNORM
C
C     Calculate mhat = nhat X lhat (already normalized)
C
         EM(I,1) =   S(I,2) * EL(I,3) - S(I,3) * EL(I,2)
         EM(I,2) = - S(I,1) * EL(I,3) + S(I,3) * EL(I,1)
         EM(I,3) =   S(I,1) * EL(I,2) - S(I,2) * EL(I,1)
  100 CONTINUE
C
      RETURN
      END
