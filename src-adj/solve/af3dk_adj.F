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
      SUBROUTINE AF3DK_ADJ  (IDIM, JDIM, KDIM, I,
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
      include '../../src-flo/header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          DTJ   (JDIM,KDIM,IDIM)
C
      DIMENSION RES   (0:JDIM+2,0:KDIM+2,0:IDIM+2,NF),
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
C     Calculate RTEMP = T^T * RES
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
      CALL TTrR   (ITURB, NPTS, QI, STEMP, EL, EM, TKE, RTEMP2, RTEMP)
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
C     Second-order MUSCL scheme
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
C     Calculate RTEMP = T^[-T] * RTEMP
C
      NPTS = (KDIM + 1) * (JDIM - 1)
      CALL TTrINVR(ITURB, NPTS, QI, STEMP, EL, EM, TKE, RTEMP, RTEMP2)
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
