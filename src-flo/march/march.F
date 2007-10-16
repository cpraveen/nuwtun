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



      SUBROUTINE LGSFBJ (IDIM, JDIM, KDIM, I, IPSTRT, IPEND,
     1                   RES, A, B, C, RTEMP)
C
C     Routine to do the J implicit forward/backward sweep of LU
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : The current plane
C     IPSTRT         : I Plane (cell center) to start calculations at 
C                      (Minimum Value: 2)
C     IPEND          : I Plane (cell center) to end caclulations at
C                      (Maximum Value: IDIM)
C     RES            : Residual - must be over the entire field for
C                      this time integration scheme
C     A,B,C          : Tridiagonal systems to invert
C                      A*DQ(J-1) + B*DQ(J) + C*DQ(J+1) = RTEMP
C                      Implicit in the J direction.  Vectorization
C                      by inverting multiple systems at one time to
C                      avoid the recursion of the inversion.  The 
C                      vectorization direction is K.  This is
C                      stored in the first index for a unit stride
C                      in the inversion.
C     RTEMP          : Storage for the residual array in the form 
C                      appropriate for the inversion routine
C
      include '../header/common.h'
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF)
C
C     Storage for the implicit system to invert
C     Store the system only in the internal portion of
C     the field and use explicit boundary conditions
C     of the form Delta(Q) boundaries = 0.
C
      DIMENSION A     (2:KDIM,2:JDIM,NRANK,NRANK),
     1          B     (2:KDIM,2:JDIM,NRANK,NRANK),
     2          C     (2:KDIM,2:JDIM,NRANK,NRANK),
     3          RTEMP (2:KDIM,2:JDIM,NF)
C
C     Store Residual array for inversion
C
      DO 490 L = 1, NF
         DO 480 K = 2, KDIM
            DO 470 J = 2, JDIM
               RTEMP(K,J,L) = RES(J,K,I,L)
  470       CONTINUE
  480    CONTINUE
  490 CONTINUE
C
C     Forward/Backward sweeps of the LU decomposition
C
      NPTS = JDIM - 1
      NSIM = KDIM - 1
      CALL VFBTR (NPTS, NSIM, NRANK, A, B, C, RTEMP)
C
C     Update the residual array
C
      DO 880 L = 1, NF
         DO 870 K = 2, KDIM
            DO 860 J = 2, JDIM
               RES(J,K,I,L) = RTEMP(K,J,L)
  860       CONTINUE
  870    CONTINUE
  880 CONTINUE
C
C     Finished with the J implicit forward/backward sweep of LU
C
      RETURN
      END






      SUBROUTINE LGSFBK (IDIM, JDIM, KDIM, I, IPSTRT, IPEND,
     1                   RES, A, B, C, RTEMP)
C
C     Routine to do the K implicit forward/backward sweep of LU
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : The current plane
C     IPSTRT         : I Plane (cell center) to start calculations at 
C                      (Minimum Value: 2)
C     IPEND          : I Plane (cell center) to end caclulations at
C                      (Maximum Value: IDIM)
C     RES            : Residual - must be over the entire field for
C                      this time integration scheme
C     A,B,C          : Tridiagonal systems to invert
C                      A*DQ(K-1) + B*DQ(K) + C*DQ(K+1) = RTEMP
C                      Implicit in the I direction.  Vectorization
C                      by inverting multiple systems at one time to
C                      avoid the recursion of the inversion.  The 
C                      vectorization direction is J.  This is
C                      stored in the first index for a unit stride
C                      in the inversion.
C     RTEMP          : Storage for the residual array in the form 
C                      appropriate for the inversion routine
C
      include '../header/common.h'
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF)
C
C     Storage for the implicit system to invert
C     Store the system only in the internal portion of
C     the field and use explicit boundary conditions
C     of the form Delta(Q) boundaries = 0.
C
      DIMENSION A     (2:JDIM,2:KDIM,NRANK,NRANK),
     1          B     (2:JDIM,2:KDIM,NRANK,NRANK),
     2          C     (2:JDIM,2:KDIM,NRANK,NRANK),
     3          RTEMP (2:JDIM,2:KDIM,NF)
C
C     Store Residual array for inversion
C
      DO 490 L = 1, NF
         DO 480 K = 2, KDIM
            DO 470 J = 2, JDIM
               RTEMP(J,K,L) = RES(J,K,I,L)
  470       CONTINUE
  480    CONTINUE
  490 CONTINUE
C
C     Forward/Backward sweeps of the LU decomposition to invert
C
      NPTS = KDIM - 1
      NSIM = JDIM - 1
      CALL VFBTR (NPTS, NSIM, NRANK, A, B, C, RTEMP)
C
C     Update the residual array
C
      DO 880 L = 1, NF
         DO 870 K = 2, KDIM
            DO 860 J = 2, JDIM
               RES(J,K,I,L) = RTEMP(J,K,L)
  860       CONTINUE
  870    CONTINUE
  880 CONTINUE
C
C     Finished with the K implicit forward/backward sweep of the LU
C
      RETURN
      END






      SUBROUTINE LGSLUJ (IDIM, JDIM, KDIM, I, IPSTRT, IPEND,
     1                   Q, PROPS, NPPTS, S, VOL, DTJ, QP, QM, 
     2                   STEMP, IBCDAT, NBCS, ICUTS, NCUTS, 
     3                   IFLUX, VISCOS, ITURB,
     4                   A, B, C, RM, APT, AMT, ART, AI, WORK)
C
C     Routine to form the LU decomposition of the J implicit sweep of 
C     Line Gauss Seidel with Approximate Factorisation in 
C     the I = constant plane.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : The current plane
C     IPSTRT         : I Plane (cell center) to start calculations at 
C                      (Minimum Value: 2)
C     IPEND          : I Plane (cell center) to end caclulations at
C                      (Maximum Value: IDIM)
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
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     STEMP          : Storage for metrics passed to Jacobian routine
C                      including the area term
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     ICUTS          : Data controlling zone communication
C     NCUTS          : Number of cuts for current block
C     IFLUX          : Flux evaluation scheme
C     VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C     ITURB          : Turbulence model
C     A,B,C          : LU decomposition of the tridiagonal systems to invert
C                      A*DQ(J-1) + B*DQ(J) + C*DQ(J+1) = RTEMP
C                      Implicit in the J direction.  Vectorization
C                      by inverting multiple systems at one time to
C                      avoid the recursion of the inversion.  The 
C                      vectorization direction is K.  This is
C                      stored in the first index for a unit stride
C                      in the inversion.
C     RM             : The I direction contribution to the implicit scheme
C     APT,AMT,ART    : Storage on a plane (JxK) to calculate the
C                      Jacobians
C     AI             : Identity Matrix
C     WORK           : Temporary array of size NPTS*4 for Jacobian calculation
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          DTJ   (JDIM,KDIM,IDIM)
C
      DIMENSION QP    (JDIM,2:KDIM,NQ), 
     1          QM    (JDIM,2:KDIM,NQ),
     2          STEMP (JDIM,2:KDIM,4),
     3          WORK  (JDIM,2:KDIM,4)
C
C     Storage for the implicit system to invert
C     Store the system only in the internal portion of
C     the field and use explicit boundary conditions
C     of the form Delta(Q) boundaries = 0.
C
      DIMENSION A     (2:KDIM,2:JDIM,NRANK,NRANK),
     1          B     (2:KDIM,2:JDIM,NRANK,NRANK),
     2          C     (2:KDIM,2:JDIM,NRANK,NRANK) 
C
      DIMENSION RM    (2:JDIM,2:KDIM,NRANK,NRANK)
C
C     Storage for the Jacobians used to calculate A,B and C
C
      DIMENSION APT   (JDIM,2:KDIM,NRANK,NRANK),
     1          AMT   (JDIM,2:KDIM,NRANK,NRANK),
     2          ART   (JDIM,2:KDIM,NRANK,NRANK)
C
      DIMENSION AI    (NRANK,NRANK)
C
      DIMENSION IBCDAT(10,NBCS), 
     1          ICUTS (21,NCUTS)
C
      LOGICAL   VISCOS
C
C     Calculate QP,QM using first order
C
      DO 130 L = 1, NQ
         DO 120 K = 2, KDIM
            DO 110 J = 1, JDIM
               QP(J,K,L) = Q(J+1,K,I,L)
               QM(J,K,L) = Q(J,K,I,L)
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C     Store the metrics for the Jacobian calculation
C     including the area term
C
      DO 190 L = 1, 4
         DO 180 K = 2, KDIM
            DO 170 J = 1, JDIM
               STEMP(J,K,L) = S(J,K,I,L,2)
  170       CONTINUE
  180    CONTINUE
  190 CONTINUE
C
C     Calculate Jacobians and set up matrices to invert
C
      NPTS = JDIM * (KDIM - 1)
C
C     Roe's scheme
C          APT = dF/dQ(QP)
C          AMT = dF/dQ(QM)
C          ART = T |Lambda| T(-1)
C
      IF (IFLUX .EQ. IFROE) THEN
         CALL DFROE (ITURB, NPTS, QP, QM, STEMP, APT, AMT, ART, WORK)
C
         DO 214 L = 1, NRANK
            DO 213 M = 1, NRANK
               DO 212 K = 2, KDIM
                  DO 211 J = 2, JDIM
                     A(K,J,L,M) = -0.5E0*(AMT(J-1,K,L,M)+ART(J-1,K,L,M))
                     B(K,J,L,M) =  RM(J,K,L,M)
     1                            +0.5E0*( AMT(J,K,L,M)-APT(J-1,K,L,M)
     2                                   + ART(J,K,L,M)+ART(J-1,K,L,M) )
                     C(K,J,L,M) =  0.5E0*(APT(J,K,L,M)-ART(J,K,L,M))
  211             CONTINUE
  212          CONTINUE
  213       CONTINUE
  214    CONTINUE
C
      ELSE
         write (iout,'(a)') ' *** Invalid Flux Function ***'
         stop
      ENDIF
C
C     Implicit Viscous Terms
C          AP = DF(j+1/2)/DQ(j+1)
C          AM = DF(j+1/2)/DQ(j)
C          PROPS(1) is RMU and PROPS(2) is RMUT
C          Use QM as temporary space to store U,V,W,RMUJ,RMUTJ
C          Use QP as temporary space to store T - being
C          careful as T is dimensioned as JDIM+1, therefore
C          store T in QP(1) and leave QP(2) unused in this call.
C
      IF ( VISCOS ) THEN
         CALL DFVISJ (IDIM, JDIM, KDIM, I, ITURB,
     1                Q, PROPS(1,1), PROPS(1,2), S, VOL, IBCDAT, NBCS, 
     2                QP(1,2,1), QM(1,2,1), QM(1,2,2), QM(1,2,3), 
     3                QM(1,2,4), QM(1,2,5), APT, AMT)
C
         DO 324 L = 1, NRANK
            DO 323 M = 1, NRANK
               DO 322 K = 2, KDIM
                  DO 321 J = 2, JDIM
                     A(K,J,L,M) = A(K,J,L,M)
     1                    + AMT(J-1,K,L,M)
                     B(K,J,L,M) = B(K,J,L,M)
     1                    + ( APT(J,K,L,M)-AMT(J-1,K,L,M) )
                     C(K,J,L,M) = C(K,J,L,M) 
     1                    - APT(J,K,L,M)
  321             CONTINUE
  322          CONTINUE
  323       CONTINUE
  324    CONTINUE
      ENDIF
C
C     Form the LU decomposition of the block tridiagonal system:
C     A*DQ(J-1) + B*DQ(J) + C*DQ(J+1)
C
      NPTS = JDIM - 1
      NSIM = KDIM - 1
      CALL VLUTR (NPTS, NSIM, NRANK, A, B, C)
C
C     Finished forming the LU decomposition of the J implicit sweep of 
C     Line Gauss Seidel with AF in the plane.
C
      RETURN
      END






      SUBROUTINE LGSLUK (IDIM, JDIM, KDIM, I, IPSTRT, IPEND,
     1                   Q, PROPS, NPPTS, S, VOL, DTJ, QP, QM, 
     2                   STEMP, IBCDAT, NBCS, ICUTS, NCUTS, 
     3                   IFLUX, VISCOS, ITURB,
     4                   A, B, C, RM, APT, AMT, ART, AI, WORK)
C
C     Routine to form the LU decomposition of the K implicit sweep of 
C     Line Gauss Seidel with Approximate Factorisation in 
C     the I = constant plane.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : The current plane
C     IPSTRT         : I Plane (cell center) to start calculations at 
C                      (Minimum Value: 2)
C     IPEND          : I Plane (cell center) to end caclulations at
C                      (Maximum Value: IDIM)
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
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     STEMP          : Storage for metrics passed to Jacobian routine
C                      including the area term
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     ICUTS          : Data controlling zone communication
C     NCUTS          : Number of cuts for current block
C     IFLUX          : Flux evaluation scheme
C     VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C     ITURB          : Turbulence model
C     A,B,C          : LU decomposition of the tridiagonal systems to invert
C                      A*DQ(K-1) + B*DQ(K) + C*DQ(K+1) = RTEMP
C                      Implicit in the K direction.  Vectorization
C                      by inverting multiple systems at one time to
C                      avoid the recursion of the inversion.  The 
C                      vectorization direction is J.  This is
C                      stored in the first index for a unit stride
C                      in the inversion.
C     RM             : The I direction contribution to the implicit scheme
C     APT,AMT,ART    : Storage on a plane (KxJ) to calculate the
C                      Jacobians
C     AI             : Identity Matrix
C     WORK           : Temporary array of size NPTS*4 for Jacobian calculation
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          DTJ   (JDIM,KDIM,IDIM)
C
      DIMENSION QP    (KDIM,2:JDIM,NQ), 
     1          QM    (KDIM,2:JDIM,NQ),
     2          STEMP (KDIM,2:JDIM,4),
     3          WORK  (KDIM,2:JDIM,4)
C
C     Storage for the implicit system to invert
C     Store the system only in the internal portion of
C     the field and use explicit boundary conditions
C     of the form Delta(Q) boundaries = 0.
C
      DIMENSION A     (2:JDIM,2:KDIM,NRANK,NRANK),
     1          B     (2:JDIM,2:KDIM,NRANK,NRANK),
     2          C     (2:JDIM,2:KDIM,NRANK,NRANK)
C
      DIMENSION RM    (2:JDIM,2:KDIM,NRANK,NRANK)
C
C     Storage for the Jacobians used to calculate A,B and C
C
      DIMENSION APT   (KDIM,2:JDIM,NRANK,NRANK),
     1          AMT   (KDIM,2:JDIM,NRANK,NRANK),
     2          ART   (KDIM,2:JDIM,NRANK,NRANK)
C
      DIMENSION AI    (NRANK,NRANK)
C
      DIMENSION IBCDAT(10,NBCS), 
     1          ICUTS (21,NCUTS)
C
      LOGICAL   VISCOS
C
C     Calculate QP,QM using first order
C
      DO 130 L = 1, NQ
         DO 120 J = 2, JDIM
            DO 110 K = 1, KDIM
               QP(K,J,L) = Q(J,K+1,I,L)
               QM(K,J,L) = Q(J,K,I,L)
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C     Store the metrics for the Jacobian calculation
C     including the area term
C
      DO 190 L = 1, 4
         DO 180 J = 2, JDIM
            DO 170 K = 1, KDIM
               STEMP(K,J,L) = S(J,K,I,L,3)
  170       CONTINUE
  180    CONTINUE
  190 CONTINUE
C
C     Calculate Jacobians and set up matrices to invert
C
      NPTS = KDIM * (JDIM - 1)
C
C     Roe's scheme
C          APT = dF/dQ(QP)
C          AMT = dF/dQ(QM)
C          ART = T |Lambda| T(-1)
C
      IF (IFLUX .EQ. IFROE) THEN
         CALL DFROE (ITURB, NPTS, QP, QM, STEMP, APT, AMT, ART, WORK) 
C
         DO 214 L = 1, NRANK
            DO 213 M = 1, NRANK
               DO 212 K = 2, KDIM
                  DO 211 J = 2, JDIM
                     A(J,K,L,M) = -0.5E0*(AMT(K-1,J,L,M)+ART(K-1,J,L,M))
                     B(J,K,L,M) =  RM(J,K,L,M)
     1                            +0.5E0*( AMT(K,J,L,M)-APT(K-1,J,L,M)
     2                                   + ART(K,J,L,M)+ART(K-1,J,L,M) )
                     C(J,K,L,M) =  0.5E0*(APT(K,J,L,M)-ART(K,J,L,M))
  211             CONTINUE
  212          CONTINUE
  213       CONTINUE
  214    CONTINUE
C
      ELSE
         write (iout,'(a)') ' *** Invalid Flux Function ***'
         stop
      ENDIF
C
C     Implicit Viscous Terms
C          AP = DF(k+1/2)/DQ(k+1)
C          AM = DF(k+1/2)/DQ(k)
C          PROPS(1) is RMU and PROPS(2) is RMUT
C          Use QM as temporary space to store U,V,W,RMUK,RMUTK
C          Use QP as temporary space to store T - being
C          careful as T is dimensioned as KDIM+1, therefore
C          store T in QP(1) and leave QP(2) unused in this call.
C
      IF ( VISCOS ) THEN
         CALL DFVISK (IDIM, JDIM, KDIM, I, ITURB,
     1                Q, PROPS(1,1), PROPS(1,2), S, VOL, IBCDAT, NBCS, 
     2                QP(1,2,1), QM(1,2,1), QM(1,2,2), QM(1,2,3), 
     3                QM(1,2,4), QM(1,2,5), APT, AMT)
C
         DO 324 L = 1, NRANK
            DO 323 M = 1, NRANK
               DO 322 K = 2, KDIM
                  DO 321 J = 2, JDIM
                     A(J,K,L,M) = A(J,K,L,M)
     1                    + AMT(K-1,J,L,M)
                     B(J,K,L,M) = B(J,K,L,M)
     1                    + ( APT(K,J,L,M)-AMT(K-1,J,L,M) )
                     C(J,K,L,M) = C(J,K,L,M) 
     1                    - APT(K,J,L,M)
  321             CONTINUE
  322          CONTINUE
  323       CONTINUE
  324    CONTINUE
      ENDIF
C
C     Form the LU decomposition of the block tridiagonal system:
C     A*DQ(K-1) + B*DQ(K) + C*DQ(K+1)
C
      NPTS = KDIM - 1
      NSIM = JDIM - 1
      CALL VLUTR (NPTS, NSIM, NRANK, A, B, C)
C
C     Finished forming the LU decomposition of the K implicit sweep of 
C     Line Gauss Seidel with AF in the plane.
C
      RETURN
      END






      SUBROUTINE LGSRES (IDIM, JDIM, KDIM, I, IPSTRT, IPEND,
     1                   RES, RM)
C
C     Routine to multiply RM by RES as the step between the two
C     Approximate Factorisation sweeps of implicit Line Gauss Seidel 
C     in the I = constant plane.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : The current plane
C     IPSTRT         : I Plane (cell center) to start calculations at 
C                      (Minimum Value: 2)
C     IPEND          : I Plane (cell center) to end caclulations at
C                      (Maximum Value: IDIM)
C     RES            : Residual - must be over the entire field for
C                      this time integration scheme
C     RM             : The I direction contribution to the implicit scheme
C
      include '../header/common.h'
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF)
C
C     Storage for the I direction contribution to the implicit scheme
C
      DIMENSION RM    (2:JDIM,2:KDIM,NRANK,NRANK)
C
C     Multiply RM by RES
C
C     Rank = 5
C
      IF (NRANK .EQ. 5) THEN
         DO 120 K = 2, KDIM
            DO 110 J = 2, JDIM
               R1 = RES(J,K,I,1)
               R2 = RES(J,K,I,2)
               R3 = RES(J,K,I,3)
               R4 = RES(J,K,I,4)
               R5 = RES(J,K,I,5)
               RES(J,K,I,1) = RM(J,K,1,1)*R1 + RM(J,K,1,2)*R2
     1                      + RM(J,K,1,3)*R3 + RM(J,K,1,4)*R4
     2                      + RM(J,K,1,5)*R5
               RES(J,K,I,2) = RM(J,K,2,1)*R1 + RM(J,K,2,2)*R2
     1                      + RM(J,K,2,3)*R3 + RM(J,K,2,4)*R4
     2                      + RM(J,K,2,5)*R5
               RES(J,K,I,3) = RM(J,K,3,1)*R1 + RM(J,K,3,2)*R2
     1                      + RM(J,K,3,3)*R3 + RM(J,K,3,4)*R4
     2                      + RM(J,K,3,5)*R5
               RES(J,K,I,4) = RM(J,K,4,1)*R1 + RM(J,K,4,2)*R2
     1                      + RM(J,K,4,3)*R3 + RM(J,K,4,4)*R4
     2                      + RM(J,K,4,5)*R5
               RES(J,K,I,5) = RM(J,K,5,1)*R1 + RM(J,K,5,2)*R2
     1                      + RM(J,K,5,3)*R3 + RM(J,K,5,4)*R4
     2                      + RM(J,K,5,5)*R5
  110       CONTINUE
  120    CONTINUE
C
C     Rank = 7
C
      ELSE IF (NRANK .EQ. 7) THEN
         DO 220 K = 2, KDIM
            DO 210 J = 2, JDIM
               R1 = RES(J,K,I,1)
               R2 = RES(J,K,I,2)
               R3 = RES(J,K,I,3)
               R4 = RES(J,K,I,4)
               R5 = RES(J,K,I,5)
               R6 = RES(J,K,I,6)
               R7 = RES(J,K,I,7)
               RES(J,K,I,1) = RM(J,K,1,1)*R1 + RM(J,K,1,2)*R2
     1                      + RM(J,K,1,3)*R3 + RM(J,K,1,4)*R4
     2                      + RM(J,K,1,5)*R5 + RM(J,K,1,6)*R6
     3                      + RM(J,K,1,7)*R7
               RES(J,K,I,2) = RM(J,K,2,1)*R1 + RM(J,K,2,2)*R2
     1                      + RM(J,K,2,3)*R3 + RM(J,K,2,4)*R4
     2                      + RM(J,K,2,5)*R5 + RM(J,K,2,5)*R5
     3                      + RM(J,K,2,7)*R7
               RES(J,K,I,3) = RM(J,K,3,1)*R1 + RM(J,K,3,2)*R2
     1                      + RM(J,K,3,3)*R3 + RM(J,K,3,4)*R4
     2                      + RM(J,K,3,5)*R5 + RM(J,K,3,5)*R5
     3                      + RM(J,K,3,7)*R7
               RES(J,K,I,4) = RM(J,K,4,1)*R1 + RM(J,K,4,2)*R2
     1                      + RM(J,K,4,3)*R3 + RM(J,K,4,4)*R4
     2                      + RM(J,K,4,5)*R5 + RM(J,K,4,5)*R5
     3                      + RM(J,K,4,7)*R7
               RES(J,K,I,5) = RM(J,K,5,1)*R1 + RM(J,K,5,2)*R2
     1                      + RM(J,K,5,3)*R3 + RM(J,K,5,4)*R4
     2                      + RM(J,K,5,5)*R5 + RM(J,K,5,5)*R5
     3                      + RM(J,K,5,7)*R7
               RES(J,K,I,6) = RM(J,K,6,1)*R1 + RM(J,K,6,2)*R2
     1                      + RM(J,K,6,3)*R3 + RM(J,K,6,4)*R4
     2                      + RM(J,K,6,5)*R5 + RM(J,K,6,5)*R5
     3                      + RM(J,K,6,7)*R7
               RES(J,K,I,7) = RM(J,K,7,1)*R1 + RM(J,K,7,2)*R2
     1                      + RM(J,K,7,3)*R3 + RM(J,K,7,4)*R4
     2                      + RM(J,K,7,5)*R5 + RM(J,K,7,5)*R5
     3                      + RM(J,K,7,7)*R7
  210       CONTINUE
  220    CONTINUE
C
C
C     Error
C
      ELSE 
         WRITE (IOUT,1000) NRANK
         STOP
      ENDIF
C
C     Finished with the multiplication of RM and RES.
C
      RETURN
 1000 FORMAT (' ','LGSRES : ERROR-> Unsupported matrix size for ',
     1            'matrix-vector multiply.'
     2       /' ','                 Contact programmer.')
      END






      SUBROUTINE MARCH (IDIM, JDIM, KDIM, MXSECT, MXABC, 
     1                  Q, QN, NQNPTS, PROPS, NPPTS, S, VOL,RC,DTJ,RES, 
     2                  FHAT, DQP, DQM, QP, QM, TAU, F, FM, STEMP,WORKX,
     3                  IBCDAT, NBCS, NXPROF, NPROFL, IOFQBC, QBC, 
     4                  ICUTS, NCUTS, R2NORM, ITRSMX, RESMAX, IRESMX,
     5                  IORDER, LIMIT, RKAP,IFLUX,VISCOS,IMODEL,IFDELQ,
     6                  ABC, RM, RTEMP, AT, AI,
     7                  IMSTRT, IMEND, NITS, NITFO, NITALG, ITUPJ, 
     8                  TOLER, VGNRN, SIGMA, ITER, 
     9                  QMIN, NFIXQ, NRELIZ, NPRLIM, NPRNEG, 
     A                  IFL2DQ, R2DQDT, CL, RCNTR, IERRCD)
C
C     Routine to do Marching.
C     Presently, this routine will march only for one block.
C
C
C IDIM,JDIM,KDIM : Dimensions of current block
C MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                  (used in dimensioning of temporary space used
C                   for calculating q(+,-) and f.
C MXABC          : Amount of memory to use to store the tridiagonal
C                  systems to invert.  Must be at least
C                  max (IDIM*JDIM, JDIM*KDIM) * number of planes to
C                  vectorize inversion on
C Q              : Primitive variables at cell centers
C QN             : Conserved variables at time levels n, n-1, etc.
C NQNPTS         : Number of QN points. This allows avoiding storing QN 
C                  when running time accurate.
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity 
C                  PROPS(3) = Y+
C NPPTS          : Number of points at which PROPS are stored.  This allows
C                  the storage of PROPS for viscous flows and the avoidance
C                  of the additional storage for inviscid flows.
C S              : Metrics
C VOL            : Cell volumes
C RC             : Position vector at cell centers
C DTJ            : Time step divided by cell volume
C RES            : Residual
C FHAT           : Storage for fluxes for 4th order
C DQP,DQM        : Storage for delta q = q(i) - q(i-1)
C QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C TAU            : Storage for Reynolds stresses
C F,FM           : Fluxes at cell interfaces -
C                  if Flux Difference Splitting then F = Flux
C                  if Flux Vector Splitting then F = Fplus, FM = Fminus
C STEMP          : Storage for metrics passed to flux evaluation routine.
C WORKX          : Work array for cross-section (used for damping, etc.)
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C NXPROF         : Maximum number of profile boundary condition segments
C NPROFL         : Number of profile boundary condition segments
C IOFQBC         : Offsets for QBC array
C QBC            : Q array for PROFILE bc
C ICUTS          : Data controlling zone communication
C NCUTS          : Number of cuts for current block
C R2NORM         : Array of residuals for each equation and
C                  L2 norm of all equations
C ITRSMX         : Number of iterations to output maximum residual
C RESMAX         : Maximum residual
C IRESMX         : Location of maximum residual 
C IORDER         : Order of MUSCL interpolation
C                  IORDER = 1 first order
C                         = 2 second order (or 'third' order)
C LIMIT          : Limiter to use for MUSCL interpolation
C RKAP           : Kappa to use in MUSCL interpolation
C IFLUX          : Flux evaluation scheme
C VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C IMODEL         : Turbulence model scheme array
C IFDELQ         : Array controlling calculation of delq
C ABC            : Storage for the tridiagonal system to invert.
C                  This is permanent storage as it stores the LU
C                  decomposition of the two approximate factorization
C                  sweeps which are reused for every NITUPJ iterations.
C RM             : The I direction contribution to the implicit scheme
C RTEMP          : Storage for the residual to be put in the proper
C                  form to use with the inversion routine
C AT             : Storage for the Jacobian calculations
C AI             : Identity Matrix
C IMSTRT         : I plane to begin marching at
C IMEND          : I plane to end marching
C NITS           : Number of iterations to perform on each plane
C NITFO          : Number of iterations to perform first order accurate
C NITALG         : Number of iterations to update turbulent eddy viscosity
C                  with algebraic model while doing higher order turb
C                  model for the initial iterations.
C ITUPJ          : Number of iterations to update the LU decomposition
C                  of the implicit terms
C TOLER          : Specified tolerance to converge the solution on a
C                  plane - RES/RES(iter=1) <= TOLER
C VGNRN          : Logical for Vigneron approximation (TRUE) or not (FALSE)
C SIGMA          : Safety factor for the Vigneron term in the pressure
C ITER           : The iteration counter used back in the main routine
C                  to write out the restart file.
C QMIN           : Control array for minimum value of Q
C NFIXQ          : Number of times reset Q
C NRELIZ         : Number of times violate realizability
C NPRLIM         : Number of times reset production
C NPRNEG         : Number of times invoke positivity preservation for production
C IFL2DQ         : Logical to calculate L2 norm of dQ/dt term for unsteady
C R2DQDT         : L2 norm of dQ/dt term for unsteady
C CL             : Lift coefficient. Used for 2D farfield vortex bc.
C RCNTR          : Center of airfoil. Used as location of point vortex.
C                  Also used as airfoil location at which to calculate moments.
C IERRCD         : Error code
C
      include '../header/common.h'
      include '../header/histry.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QN    (NQNPTS,NQ,NTMLVL),
     2          PROPS (NPPTS,NP),
     3          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     4          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     5          RC    (JDIM+1,KDIM+1,IDIM+1,3),
     6          DTJ   (JDIM,KDIM,IDIM)
C
C     Dimension RES with I = 3 to correspond to the dimensioning
C     as IPSTRT-1:IPEND+1 in routines below this where IPSTRT = IPEND.
C
      DIMENSION RES   (JDIM+1,KDIM+1,3,NF),
     1          FHAT  (JDIM+1,KDIM+1,3,NF),
     2          DQP   (MXSECT,NQ),
     3          DQM   (MXSECT,NQ),
     4          QP    (MXSECT,NQ),
     5          QM    (MXSECT,NQ),
     6          TAU   (MXSECT,6),
     7          F     (MXSECT,NF),
     8          FM    (MXSECT,NF),
     9          STEMP (MXSECT,4),
     A          WORKX (MXSECT,NW)
C
      DIMENSION ABC   (MXABC*NRANK*NRANK,3,2),
     1          RTEMP (MXABC*NF),
     2          RM    (MXABC*NRANK*NRANK),
     3          AT    (MXSECT*NRANK*NRANK,3)
C
      DIMENSION AI    (NRANK,NRANK)
C
      DIMENSION IBCDAT(10,NBCS), 
     1          ICUTS (21,NCUTS)
C
      DIMENSION IOFQBC(NXPROF)
      DIMENSION QBC   (*)
C
      DIMENSION IORDER(3), RKAP(3), LIMIT(3)
C
      DIMENSION IMODEL(NMDL),
     1          IFDELQ(NQ)
C
      LOGICAL   VISCOS(3), INITRS, MRCHNG, VGNRN, NCSTOP, NCCHCK, IFL2DQ
C
      DIMENSION R2NORM(NF+1),
     1          R2DQDT(NF+1)
C
      DIMENSION IRESMX(0:3)
C
      DIMENSION QMIN  (NQ)
      DIMENSION NFIXQ (NQ),
     1          NRELIZ(NQ)
C
      DIMENSION RCNTR (3)
C
C     Set to march and initialize residual to zero
C
      MRCHNG = .TRUE.
      INITRS = .TRUE.
C
      ITURB  = IMODEL(1)
C
C     Test to asure that sufficient room exists to store the entire
C     LU in the plane to reuse.
C
      ISIZE = JDIM * KDIM
      IF (ISIZE .GT. MXABC) THEN
         WRITE (IOUT,1900) ISIZE
         STOP
      ENDIF
 1900 FORMAT (' ','MARCH : ERROR-> Insufficient room in ABC array to'
     1       /' ','        allow for storage of the LU decomposition in'
     2       /' ','        the marching algorithm.'
     3       /' ','        Increase MXABC in MAIN to at least ',I20)
C
C     Save the values of RKAP and LIMIT for the I direction to restore after
C     finished with marching to allow contnuation with global solution.
C     Set RKAP and LIMIT for the I direction to be fully upwind unlimited
C     for marching.
C
      RKAPS    = RKAP(1)
      LIMITS   = LIMIT(1)
      RKAP(1)  = - 1.E0
      LIMIT(1) =   0
C
C     Set variable to abort marching when a plane does not converge.
C
      NCSTOP = .TRUE.
C
C     Make sure that are marching on the available planes (2 to IDIM)
C
      IMSTRT = MAX ( IMSTRT, 2 )
      IMEND  = MIN ( IMEND,  IDIM )
C
C     March plane to plane doing NITS iterations per plane.
C
      DO 900 I = IMSTRT, IMEND
         IPSTRT = I
         IPEND  = I
C
C     Initialize the cross section to the solution from the previous
C     cross section UNLESS are doing a restart.
C     Then call boundary conditions. Use QP as work array
C
         IF (ITTOT .LE. 0) THEN
            WRITE (IOUT,1020) I
 1020 FORMAT (//' ','Initializing cross section I = ',I5)
C
            DO 130 L = 1, NQ
               DO 120 K = 1, KDIM
                  DO 110 J = 1, JDIM
                     Q(J,K,I,L) = Q(J,K,I-1,L)
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
C
            CALL BC (IDIM, JDIM, KDIM, IPSTRT, IPEND,
     1               Q, PROPS, NPPTS, S, VOL, RC, NBCS, IBCDAT, 
     2               NXPROF, NPROFL, IOFQBC, QBC,
     3               VISCOS, IMODEL, QP, MXSECT, CL, RCNTR)
         ENDIF
C
C     End of initialization section
C
         WRITE (IOUT,1000) I
 1000 FORMAT (/' ','Iteration History for Cross Section I = ',I5/)
C
C     Loop for NITS iterations on the current plane
C
         NCCHCK = .FALSE.
         DO 800 ITER = 1, NITS
C                       
C     Check for presence of STOP file; if it exists, write restart and quit
C    
            OPEN (ISTOPFL,FILE='STOP',STATUS='OLD',ERR=701)
            WRITE (IOUT,1490)
            GO TO 999  
 1490 FORMAT (//' ','MARCH : STOP file encountered! Run stopping.'//)
C    
C     STOP does not exist; continue run
C    
  701       CONTINUE  
C
C     Set up switch from first order scheme to higher order
C     This sets first order to be ONLY streamwise for now.
C
            IF (ITER .GT. NITFO) THEN
               IORDER(1) = 2
               IORDER(2) = 2
               IORDER(3) = 2
            ELSE
               IORDER(1) = 1
               IORDER(2) = 2
               IORDER(3) = 2
            ENDIF
C
C     Set up switch to use the algebraic turbulent eddy viscosity in
C     the initial phases of the iteration to improve convergence.
C
            IF (ITER .GT. NITALG) THEN
               MUTALG = .FALSE.
            ELSE
               MUTALG = .TRUE.
            ENDIF
C
C     Update time step.
C     Use DQP as temporary storage for PFACTR.
C  
            IF ( MOD((ITER-1), ITUPDT) .EQ. 0) THEN
               IBLK  =  1
               DTMIN =  999.E0
               DTMAX = -999.E0
               CALL DELTAT (IBLK, IDIM, JDIM, KDIM, IPSTRT, IPEND,
     1                      Q, PROPS, DTJ, VOL, S, DQP,
     2                      VISCOS, ITURB, ITER, CFLCUR,
     3                      DTMIN, IBLKMN, IDTMIN, JDTMIN, KDTMIN,
     4                      DTMAX, IBLKMX, IDTMAX, JDTMAX, KDTMAX)
C-need to output this only when running in CFL specified mode
               WRITE (IOUT,1450) CFLCUR,
     1                           DTMIN, IBLKMN, IDTMIN, JDTMIN, KDTMIN,
     2                           DTMAX, IBLKMX, IDTMAX, JDTMAX, KDTMAX
            ENDIF
 1450 FORMAT (' ','MARCH : CFL =',1PE10.3,'  DTMIN =',1PE10.3,
     1            ' At BLOCK:I,J,K=(',I3,':',I4,',',I4,',',I4,')',
     2       /' ',23X,'  DTMAX =',1PE10.3,
     3            ' At BLOCK:I,J,K=(',I3,':',I4,',',I4,',',I4,')')
C
C     Reset norms to zero
C
            DO 400 L = 1, NF + 1
               R2NORM(L) = 0.E0
 400        CONTINUE
C
            DO 410 L = 1, NQ
               NFIXQ (L) = 0
               NRELIZ(L) = 0
 410        CONTINUE
            NPRLIM = 0
            NPRNEG = 0
            RESMAX = 0.E0
C
C     Set up control parameter for Jacobian update.
C
            IF ( MOD((ITER-1), ITUPJ) .EQ. 0) THEN
               UPDATJ = .TRUE.
            ELSE
               UPDATJ = .FALSE.
            ENDIF
C
C     If running a viscous case, calculate the PROPS array
C     Use QP, QM, F and DQP as temporary storage.
C
            IF (ITURB .GE. ITLMNR) THEN
               CALL PRPRTY (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                      Q, PROPS, S, VOL, 
     2                      IBCDAT, NBCS, IMODEL, IFDELQ,
     3                      STEMP, QP, F, FM, QM, DQP, DQM, WORKX)
            ENDIF
C
C     Call Residual routine for the current plane
C
            CALL RESID (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                  Q, QN, NQNPTS, PROPS, NPPTS, S, VOL, RC,DTJ,RES, 
     2                  FHAT, DQP, DQM, QP, QM, TAU, F, FM, STEMP,WORKX,
     3                  IBCDAT, NBCS, IORDER, LIMIT, RKAP, IFLUX,VISCOS, 
     4                  IMODEL, INITRS, IFDELQ, MRCHNG, VGNRN, SIGMA, 
     5                  NPRLIM, NPRNEG, IFL2DQ, R2DQDT)
C
C     Output if reset turbulence production
C
            IF (NPRNEG .GT. 0) THEN
               WRITE (IOUT,2400) NPRNEG
 2400 FORMAT (' ','MAIN  : from SRCRES: Turb. Prod reset posit',
     1            'ive   at ',I6,' locations')
            ENDIF
C
            IF (NPRLIM .GT. 0) THEN
               WRITE (IOUT,2401) PRDLIM, NPRLIM
 2401 FORMAT (' ','MAIN  : from SRCRES: Turb. Prod limited to ',1PE10.2,
     1            '*diss at ',I6,' locations')
            ENDIF
C
C     Calculate L2 Norm of residual on the current plane and add
C     it to the residual for the entire flowfield
C
            CALL L2NORM (IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   RES, VOL, R2NORM)
C
            DO 405 LL = 1, NF + 1
               R2NORM(LL) = SQRT ( R2NORM(LL) )
  405       CONTINUE
            IT   = ITER + ITTOT
C
C     Look at continuity equation for convergence
C
C--         R2(IT) = R2NORM(NF+1)
            R2(IT) = R2NORM(1)
C
C     Store the first residual for comparison of convergence; If the first
C     residual is zero, then store the second residual.
C
            IF (IT .EQ. 1) THEN
               R2ONE  = R2(IT)
               R2TEST = 1.E0
            ELSE IF (IT .EQ. 2) THEN
               IF (R2ONE .EQ. 0.E0) THEN
                  R2ONE  = R2(IT)
                  R2TEST = 1.E0
               ENDIF
            ELSE
               R2TEST = R2(IT) / R2ONE
            ENDIF
C
C     Calculate the maximum residual
C
            IF ( MOD( ITER, ITRSMX ) .EQ. 0) THEN
               CALL L2MAX  (IBLK, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                      RES, VOL, RESMAX, IRESMX)
               RESMAX = SQRT (RESMAX)
               WRITE (IOUT,2490) ITER, IT, RESMAX, (IRESMX(L),L=0,3)
            ENDIF
 2490 FORMAT (' ','Iter = ',2I6,' Max Res = ',1PE12.5,
     1            ' At BLOCK:I,J,K=(',I3,':',I4,',',I4,',',I4,')')
C
            WRITE (IOUT,2500) ITER, IT, (R2NORM(LL),LL=1,NF+1)
 2500 FORMAT (' ','ITER = ',2I8,' RES = ',13(1PE14.7,2X))
C
C     Set RES = - RES
C     RES is dimensioned as I = 3
C
            DO 413 L = 1, NF
               DO 412 K = 2, KDIM
                  DO 411 J = 2, JDIM
                     RES(J,K,1,L) = - RES(J,K,1,L)
                     RES(J,K,2,L) = - RES(J,K,2,L)
                     RES(J,K,3,L) = - RES(J,K,3,L)
  411             CONTINUE
  412          CONTINUE
  413       CONTINUE
C
C     Update the LU decomposition of the approximate factorization terms
C
 1010 FORMAT (' ','Jacobian LU Update')
            IF (UPDATJ) THEN
               WRITE (IOUT,1010)
C
C     Calculate the I direction contribution to the implicit scheme
C
C
               CALL MARRM (IDIM, JDIM, KDIM, I,
     1                     Q, PROPS, NPPTS, S, VOL, DTJ, 
     2                     MXSECT, DQP, DQM, QP, QM, TAU, F, FM,
     3                     STEMP, WORKX, IBCDAT, NBCS, ICUTS, NCUTS,
     4                     IFLUX, VISCOS(1), IMODEL, IFDELQ,
     5                     RM, AT(1,1), AT(1,2), AT(1,3), AI)
C
C     Form LU decomposition for J implicit sweep
C     K is the vectorization direction.
C
               CALL LGSLUJ (IDIM, JDIM, KDIM, I, IPSTRT, IPEND,
     1                      Q, PROPS, NPPTS, S, VOL, DTJ, QP, QM, 
     2                      STEMP, IBCDAT, NBCS, ICUTS, NCUTS,
     3                      IFLUX, VISCOS(2), ITURB,
     4                      ABC(1,1,1), ABC(1,2,1), ABC(1,3,1), RM,
     5                      AT(1,1), AT(1,2), AT(1,3), AI, DQM)
C
C     Form LU decomposition for K implicit sweep
C     J is the vectorization direction.
C
               IF ( THREED ) THEN
                  CALL LGSLUK (IDIM, JDIM, KDIM, I, IPSTRT, IPEND,
     1                         Q, PROPS, NPPTS, S, VOL, DTJ, QP, QM, 
     2                         STEMP, IBCDAT, NBCS, ICUTS, NCUTS,
     3                         IFLUX, VISCOS(3), ITURB,
     4                         ABC(1,1,2), ABC(1,2,2), ABC(1,3,2), RM,
     5                         AT(1,1), AT(1,2), AT(1,3), AI, DQM)
               ENDIF
            ENDIF
C
C     Forward/Backward sweeps of LU decomposition
C
C     J implicit sweep
C     K is the vectorization direction.
C
            CALL LGSFBJ (IDIM, JDIM, KDIM, I, IPSTRT, IPEND,
     1                   RES, ABC(1,1,1), ABC(1,2,1), ABC(1,3,1), 
     2                   RTEMP)
C
C     K implicit sweep
C     J is the vectorization direction.
C
            IF ( THREED ) THEN
               CALL LGSRES (IDIM, JDIM, KDIM, I, IPSTRT, IPEND,
     1                      RES, RM)
               CALL LGSFBK (IDIM, JDIM, KDIM, I, IPSTRT, IPEND,
     1                      RES, ABC(1,1,2), ABC(1,2,2), ABC(1,3,2), 
     2                      RTEMP)
            ENDIF
C
C     Conservative Variable Update
C
            CALL UPDATE (IDIM, JDIM, KDIM, IPSTRT, IPEND, Q, RES, ITURB,
     1                   QP, QM, QMIN, NFIXQ, NRELIZ, FHAT, IERRCD)
cN.B. Q is it not NaN until after this call!
c     kk = 2
c     do jj = 1, jdim
c        write (iout,'(a,3i3,5e12.4)') 'i,q', i, jj, kk,
c    1         (q(jj,kk,i,ll),ll=1,5)
c     end do
C
C     Check error code and abort run
C
            IF (IERRCD .NE. INO) THEN
               GO TO 999
            ENDIF
C
C     Call Boundary Conditions
C     Use QP as work array
C
            CALL BC (IDIM, JDIM, KDIM, IPSTRT, IPEND,
     1               Q, PROPS, NPPTS, S, VOL, RC, NBCS, IBCDAT, 
     2               NXPROF, NPROFL, IOFQBC, QBC,
     3               VISCOS, IMODEL, QP, MXSECT, CL, RCNTR)
C
C     Check for convergence of present cross section
C
            IF (R2TEST .LE. TOLER) GO TO 805
 800     CONTINUE
C
C     Set up check for non-convergence of cross section
C
          NCCHCK = .TRUE.
 805     CONTINUE
C
C     Write body data
C
         WRITE (IOUT,1100) I
 1100 FORMAT (///' ','  Solution completed at I = ',I5/)
C        DO 897 J = 1, JDIM
C           WRITE (IOUT,1105) J,(Q(J,2,I,L),L=1,NQ)
C1105 FORMAT (' ',I3,5X,6(E14.7,2X))
C897     CONTINUE
C
C     Stop the marching procedure if the present cross section did not converge
C     and the variable NCSTOP is true to enable this check
C
         IF ( NCCHCK ) THEN
            WRITE (IOUT,1200)
 1200 FORMAT (///' ','Marching did not converge on present plane',
     1          /' ','*** MARCHING ABORTED ***')
            IF ( NCSTOP ) GO TO 999
         ENDIF
C
C     Reset ITTOT to zero for the next cross section (this one is converged)
C
         ITTOT = 0
C
  900 CONTINUE
C
C     If all planes converged fine then will get here.  If this occurs,
C     set ITER to zero on return to the main routine so that the convergence
C     history is not saved for a restart of the current plane.
C
      ITER  = 0
C
C     Rewind the restart file that is written to avoid storing residual
C     history
C
      REWIND (IWRRST)
C
C     Jump here if plane doesn't converge and NCSTOP is set or error detected
C     Jump here if STOP file is detected
C
  999 CONTINUE
C
C     Restore RKAP and LIMIT for the I direction to continue global
C     solution.
C
      RKAP(1)  = RKAPS
      LIMIT(1) = LIMITS
C
C     Finished with Marching
C
      RETURN
      END






      SUBROUTINE MARRM (IDIM, JDIM, KDIM, I,
     1                  Q, PROPS, NPPTS, S, VOL, DTJ, 
     2                  MXSECT, DQP, DQM, QP, QM, TAU, F, FM,
     3                  STEMP, WORKX, IBCDAT, NBCS, ICUTS, NCUTS, 
     4                  IFLUX, VISCOS, IMODEL, IFDELQ,
     5                  RM, APT, AMT, ART, AI)
C
C     Routine to calculate the I direction contribution to the 
C     implicit Line Gauss Seidel with Approximate Factorisation 
C     in the I = constant plane for MARCHING.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : The current plane
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
C     MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                      (used in dimensioning of temporary space used
C                       for calculating q(+,-) and f.
C     DQP,DQM        : Storage for delta q = q(j) - q(j-1)
C     QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C     TAU            : Storage for Reynolds stresses
C     F,FM           : Fluxes at cell interfaces - currently temp storage
C     STEMP          : Storage for metrics passed to Jacobian routine
C                      including the area term
C     WORKX          :
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     ICUTS          : Data controlling zone communication
C     NCUTS          : Number of cuts for current block
C     IFLUX          : Flux evaluation scheme
C     VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C     IMODEL         : Turbulence model array
C     IFDELQ         : Array controlling calculation of delq
C     RM             : The I direction contribution to the implicit scheme
C     APT,AMT,ART    : Storage on a plane (JxK) to calculate the
C                      Jacobians
C     AI             : Identity Matrix
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          DTJ   (JDIM,KDIM,IDIM)
C
      DIMENSION DQP   (MXSECT,NQ),
     1          DQM   (MXSECT,NQ),
     2          QP    (2:JDIM,2:KDIM,NQ), 
     3          QM    (2:JDIM,2:KDIM,NQ),
     4          TAU   (2:JDIM,2:KDIM,6),
     5          F     (2:JDIM,2:KDIM,NQ),
     6          FM    (2:JDIM,2:KDIM,NQ),
     7          STEMP (2:JDIM,2:KDIM,4),
     8          WORKX (MXSECT,NW)
C
C     Storage for the I direction contribution to the implicit scheme
C
      DIMENSION RM    (2:JDIM,2:KDIM,NRANK,NRANK)
C
C     Storage for the Jacobians used to calculate A,B and C
C
      DIMENSION APT   (2:JDIM,2:KDIM,NRANK,NRANK),
     1          AMT   (2:JDIM,2:KDIM,NRANK,NRANK),
     2          ART   (2:JDIM,2:KDIM,NRANK,NRANK)
C
      DIMENSION AI    (NRANK,NRANK)
C
      DIMENSION IBCDAT(10,NBCS), 
     1          ICUTS (21,NCUTS)
C
      DIMENSION IMODEL(NMDL),
     1          IFDELQ(NQ)
C
      LOGICAL   VISCOS
C
      ITURB  = IMODEL(1)
C
C     Contribution to RM from I+1/2
C
C     Calculate QM using first order
C
      DO 113 L = 1, NQ
         DO 112 K = 2, KDIM
            DO 111 J = 2, JDIM
               QM(J,K,L) = Q(J,K,I,L)
  111       CONTINUE
  112    CONTINUE
  113 CONTINUE
C
C     Store the metrics for the Jacobian calculation
C     including the area term
C
      DO 163 L = 1, 4
         DO 162 K = 2, KDIM
            DO 161 J = 2, JDIM
               STEMP(J,K,L) = S(J,K,I,L,1)
  161       CONTINUE
  162    CONTINUE
  163 CONTINUE
C
C     Calculate Jacobians and add I+1/2 contribution to RM matrix
C
      NPTS = (JDIM - 1) * (KDIM - 1)
      CALL DFDQ (ITURB, NPTS, QM, STEMP, DQM(1,1), DQM(1,2), AMT)
C
      DO 214 L = 1, NRANK
         DO 213 M = 1, NRANK
            DO 212 K = 2, KDIM
               DO 211 J = 2, JDIM
                  RM(J,K,L,M) = AI(L,M) / DTJ(J,K,I) + AMT(J,K,L,M)
  211          CONTINUE
  212       CONTINUE
  213    CONTINUE
  214 CONTINUE
C
C     Source term contribution to RM matrix
C     Temporary variables passed to the source term evaluation routine:
C         QP     for QC (variables in one-dimensional array)
C         FM     for PROPSC (properties in one-dimensional array - this is 
C                            passed in using FM.  Must always be sure that 
C                            NP <= NQ or else will need to fix this)
C         QM     for DQDX (derivatives of Q in x direction at cell centers)
C         DQP    for DQDY (derivatives of Q in y direction at cell centers)
C         DQM    for DQDZ (derivatives of Q in z direction at cell centers)
C         STEMP  for DEL2K (Laplacian of turbulent kinetic energy)
C         WORKX(1)     for TKE
C         WORKX(2)     for ETILDE
C         WORKX(3)     for XSIW
C         WORKX(4-5)   for FDAMP
C         WORKX(6)     for EPSC
C         WORKX(7)     for PDIL
C         WORKX( 8-13) for SIJ
C         WORKX(14-19) for WIJ
C         F            for SRC
C
      IF (ISRCE) THEN
         IF ( ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW .OR.
     1        ITURB .EQ. ITRS ) THEN
            CALL DSRCDQ (IDIM, JDIM, KDIM, I, IMODEL, 
     1                   Q, PROPS, NPPTS, S, VOL, 
     2                   QP, FM, QM, DQP, DQM, STEMP, 
     3                   WORKX(1,1), WORKX(1,2), WORKX(1,3),WORKX(1, 4), 
     4                   WORKX(1,6), WORKX(1,7), WORKX(1,8),WORKX(1,14), 
     5                   TAU, F, IBCDAT, NBCS, APT)
C 
            DO 314 L = 1, NRANK 
               DO 313 M = 1, NRANK
                  DO 312 K = 2, KDIM
                     DO 311 J = 2, JDIM
                        RM(J,K,L,M) = RM(J,K,L,M) - APT(J,K,L,M)
  311                CONTINUE
  312             CONTINUE
  313          CONTINUE
  314       CONTINUE
         ENDIF
      ENDIF
C
C     Finished with the RM calculation.
C
      RETURN
      END
