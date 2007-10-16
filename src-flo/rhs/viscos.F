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


      SUBROUTINE DELQC (IDIM, JDIM, KDIM, I, IFDELQ, NBOUND,
     1                  Q, PROPS, S, VOL, IBCDAT, NBCS,
     2                  QJ, QK, QPI, QMI,
     3                  QC, PROPSC, DQDX, DQDY, DQDZ)
C
C     Routine to calculate the derivatives of the primitive variables
C     at the cell centers.  It also stores Q, RMU, RMUT in a one-dimensional
C     array stored the same way as the derivative data.
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : Current I Plane
C     IFDELQ         : Array to control calculation of delq
C     NBOUND         : Number of boundaries to include in derivative calculation
C                      = 0 : Calculate only at interior cell centers
C                      = 1 : Calculate at interior cell centers and on boundary
C     Q              : Primitive variables at cell centers
C     PROPS          : Properties at cell centers 
C                      PROPS(1) = molecular viscosity
C                      PROPS(2) = turbulent eddy viscosity
C                      PROPS(3) = Y+
C     S              : Metrics
C     VOL            : Cell volumes
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     QJ,QK,QPI,QMI  : Temporary storage of Q velocities on the J,K,I+,I- faces
C     QC             : Q in one-dimensional array at cell centers
C     PROPSC         : Properties in one-dimensional array
C     DQDX,DQDY,DQDZ : Derivatives of Q at cell centers
C
      include '../header/common.h'
C
      DIMENSION IFDELQ(NQ)
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP), 
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION QJ    (1-NBOUND:JDIM+NBOUND,2-NBOUND:KDIM+NBOUND),
     1          QK    (2-NBOUND:JDIM+NBOUND,1-NBOUND:KDIM+NBOUND),
     2          QPI   (2-NBOUND:JDIM+NBOUND,2-NBOUND:KDIM+NBOUND),
     3          QMI   (2-NBOUND:JDIM+NBOUND,2-NBOUND:KDIM+NBOUND)
C
      DIMENSION QC    (2-NBOUND:JDIM+NBOUND,2-NBOUND:KDIM+NBOUND,NQ),
     1          PROPSC(2-NBOUND:JDIM+NBOUND,2-NBOUND:KDIM+NBOUND,NP),
     2          DQDX  (2-NBOUND:JDIM+NBOUND,2-NBOUND:KDIM+NBOUND,NQ),
     3          DQDY  (2-NBOUND:JDIM+NBOUND,2-NBOUND:KDIM+NBOUND,NQ),
     4          DQDZ  (2-NBOUND:JDIM+NBOUND,2-NBOUND:KDIM+NBOUND,NQ)
C
      DIMENSION IBCDAT(10,NBCS)
C
C     Loop over all of the variables and calculate derivatives at cell centers
C
      DO 600 L = 1, NQ
C
C     Only calculate required derivatives
C
         IF (IFDELQ(L) .EQ. 1) THEN
C
C     Calculate Q at cell interfaces as averages of cell centers in J
C
            DO 220 K = 2-NBOUND, KDIM+NBOUND
               DO 210 J = 1-NBOUND, JDIM+NBOUND
                  QJ(J,K) = 0.5E0 * ( Q(J,K,I,L) + Q(J+1,K,I,L) )
  210          CONTINUE
  220       CONTINUE
C
C     Loop over boundary conditions and fix the wall data
C
            DO 250 IBC = 1, NBCS
               IBCTYP = IBCDAT(1,IBC)
               IDIR   = IBCDAT(2,IBC)
               IF ( (IBCTYP .EQ. IWALL .OR. IBCTYP .EQ. IWALFN)
     1              .AND. IDIR .EQ. 2) THEN
                  INROUT = IBCDAT(3,IBC)
                  ISTRT  = IBCDAT(4,IBC)
                  IEND   = IBCDAT(5,IBC)
                  KSTRT  = IBCDAT(8,IBC)
                  KEND   = IBCDAT(9,IBC)
                  IF (INROUT .EQ. 1) THEN
                     JNTR  = 1
                     JNTR2 = 2
                     JBCPT = 1
                  ELSE
                     JNTR  = JDIM
                     JNTR2 = JDIM - 1
                     JBCPT = JDIM + 1
                  ENDIF
                  IF (I .GE. ISTRT .AND. I .LE. IEND) THEN
                     DO 230 K = KSTRT, KEND
                        QJ(JNTR,K) = Q(JBCPT,K,I,L)
  230                CONTINUE
                  ENDIF
               ENDIF
  250       CONTINUE
C
C     If using the MARCHing algorithm, set the QPI and QMI
C     values to be Q(i) to remove the spatial variation of Q
C     which will not be correct as the downstream values are not updated
C
            IF (ISOLVR .EQ. IMARCH) THEN
               DO 310 K = 2-NBOUND, KDIM+NBOUND
                  DO 305 J = 2-NBOUND, JDIM+NBOUND
                     QPI(J,K) = Q(J,K,I,L) 
                     QMI(J,K) = Q(J,K,I,L) 
  305             CONTINUE
  310          CONTINUE
C
C     Calculate QPI & QMI at cell interfaces as averages of cell centers in I
C
            ELSE
               DO 320 K = 2-NBOUND, KDIM+NBOUND
                  DO 315 J = 2-NBOUND, JDIM+NBOUND
                     QPI(J,K) = 0.5E0 * ( Q(J,K,I,L) + Q(J,K,I+1,L) )
                     QMI(J,K) = 0.5E0 * ( Q(J,K,I,L) + Q(J,K,I-1,L) )
  315             CONTINUE
  320          CONTINUE
            ENDIF
C
C     Loop over boundary conditions and fix the wall data - 
C     Only if I = 2 or I = IDIM
C
            DO 350 IBC = 1, NBCS
               IBCTYP = IBCDAT(1,IBC)
               IDIR   = IBCDAT(2,IBC)
               IF ( (IBCTYP .EQ. IWALL .OR. IBCTYP .EQ. IWALFN)
     1              .AND. IDIR .EQ. 1) THEN
                  INROUT = IBCDAT(3,IBC)
                  JSTRT  = IBCDAT(6,IBC)
                  JEND   = IBCDAT(7,IBC)
                  KSTRT  = IBCDAT(8,IBC)
                  KEND   = IBCDAT(9,IBC)
                  IF (INROUT .EQ. 1 .AND. I .EQ. 2) THEN
                     IBCPT = 1
                     DO 335 K = KSTRT, KEND
                        DO 330 J = JSTRT, JEND
                           QMI(J,K) = Q(J,K,IBCPT,L)
  330                   CONTINUE
  335                CONTINUE
                  ELSE IF (INROUT .NE. 1 .AND. I .EQ. IDIM) THEN
                     IBCPT = IDIM + 1
                     DO 345 K = KSTRT, KEND
                        DO 340 J = JSTRT, JEND
                           QPI(J,K) = Q(J,K,IBCPT,L)
  340                   CONTINUE
  345                CONTINUE
                  ENDIF
               ENDIF
  350       CONTINUE
C
C     Loop over J and K forming the derivatives at the cell centers
C     Calculate first derivative terms using Divergence Theorem
C     Use symmetric differencing of the cross-derivative terms
C     from Chakravarthy et al AIAA 85-0165.
C
            DO 420 K = 2-NBOUND, KDIM+NBOUND
               DO 410 J = 2-NBOUND, JDIM+NBOUND
C
                  DQDX(J,K,L) = 
     1                   ( QPI(J,K)  * S(J,K,I  ,1,1) * S(J,K,I  ,4,1) 
     2                   - QMI(J,K)  * S(J,K,I-1,1,1) * S(J,K,I-1,4,1)
     3                   + QJ(J,K)   * S(J  ,K,I,1,2) * S(J  ,K,I,4,2)  
     4                   - QJ(J-1,K) * S(J-1,K,I,1,2) * S(J-1,K,I,4,2) )
                  DQDY(J,K,L) = 
     1                   ( QPI(J,K)  * S(J,K,I  ,2,1) * S(J,K,I  ,4,1) 
     2                   - QMI(J,K)  * S(J,K,I-1,2,1) * S(J,K,I-1,4,1)
     3                   + QJ(J,K)   * S(J  ,K,I,2,2) * S(J  ,K,I,4,2)  
     4                   - QJ(J-1,K) * S(J-1,K,I,2,2) * S(J-1,K,I,4,2) )
                  DQDZ(J,K,L) = 0.E0
  410          CONTINUE
  420       CONTINUE
C
C     If three-dimensional calculation, then calculate QK at cell 
C     interfaces as averages of cell centers in K
C
            IF (THREED .OR. AXISYM) THEN
               DO 520 J = 2-NBOUND, JDIM+NBOUND
                  DO 510 K = 1-NBOUND, KDIM+NBOUND
                     QK(J,K) = 0.5E0 * ( Q(J,K,I,L) + Q(J,K+1,I,L) )
  510             CONTINUE
  520          CONTINUE
C
C     Loop over boundary conditions and fix the wall data
C
               DO 550 IBC = 1, NBCS
                  IBCTYP = IBCDAT(1,IBC)
                  IDIR   = IBCDAT(2,IBC)
                  IF ((IBCTYP .EQ. IWALL .OR. IBCTYP .EQ. IWALFN)
     1                 .AND. IDIR .EQ. 3) THEN
                     INROUT = IBCDAT(3,IBC)
                     ISTRT  = IBCDAT(4,IBC)
                     IEND   = IBCDAT(5,IBC)
                     JSTRT  = IBCDAT(6,IBC)
                     JEND   = IBCDAT(7,IBC)
                     IF (INROUT .EQ. 1) THEN
                        KNTR  = 1
                        KNTR2 = 2
                        KBCPT = 1
                     ELSE
                        KNTR  = KDIM
                        KNTR2 = KDIM - 1
                        KBCPT = KDIM + 1
                     ENDIF
                     IF (I .GE. ISTRT .AND. I .LE. IEND) THEN
                        DO 530 J = JSTRT, JEND
                           QK(J,KNTR) = Q(J,KBCPT,I,L)
  530                   CONTINUE
                     ENDIF
                  ENDIF
  550          CONTINUE
C
C     Add K contribution to derivative calculation
C
               DO 570 K = 2-NBOUND, KDIM+NBOUND
                  DO 560 J = 2-NBOUND, JDIM+NBOUND
C
C     Metrics and volume
C
                     DQDX(J,K,L) = DQDX(J,K,L) +
     1                   ( QK(J,K)   * S(J,K  ,I,1,3) * S(J,K  ,I,4,3)  
     2                   - QK(J,K-1) * S(J,K-1,I,1,3) * S(J,K-1,I,4,3) )
                     DQDY(J,K,L) = DQDY(J,K,L) +
     1                   ( QK(J,K)   * S(J,K  ,I,2,3) * S(J,K  ,I,4,3)  
     2                   - QK(J,K-1) * S(J,K-1,I,2,3) * S(J,K-1,I,4,3) )
                     DQDZ(J,K,L) = 
     1                   ( QPI(J,K)  * S(J,K,I  ,3,1) * S(J,K,I  ,4,1) 
     2                   - QMI(J,K)  * S(J,K,I-1,3,1) * S(J,K,I-1,4,1)
     3                   + QJ(J,K)   * S(J  ,K,I,3,2) * S(J  ,K,I,4,2)  
     4                   - QJ(J-1,K) * S(J-1,K,I,3,2) * S(J-1,K,I,4,2) 
     5                   + QK(J,K)   * S(J,K  ,I,3,3) * S(J,K  ,I,4,3)  
     6                   - QK(J,K-1) * S(J,K-1,I,3,3) * S(J,K-1,I,4,3) )
  560             CONTINUE
  570          CONTINUE
            ENDIF
C
C     Divide by volume
C
            DO 590 K = 2-NBOUND, KDIM+NBOUND
               DO 580 J = 2-NBOUND, JDIM+NBOUND
                  VOLI        = 1.E0 / (VOL(J,K,I) + RSMASQ)
                  DQDX(J,K,L) = DQDX(J,K,L) * VOLI
                  DQDY(J,K,L) = DQDY(J,K,L) * VOLI
                  DQDZ(J,K,L) = DQDZ(J,K,L) * VOLI
  580          CONTINUE
  590       CONTINUE
C
         ENDIF
  600 CONTINUE
C
C     Put Q into a one-dimensional array
C
      DO 730 L = 1, NQ
         DO 720 K = 2-NBOUND, KDIM+NBOUND
            DO 710 J = 2-NBOUND, JDIM+NBOUND
               QC(J,K,L) = Q(J,K,I,L)
  710       CONTINUE
  720    CONTINUE
  730 CONTINUE
C
C     Put PROPS into a one-dimensional array
C
      DO 830 L = 1, NP
         DO 820 K = 2-NBOUND, KDIM+NBOUND
            DO 810 J = 2-NBOUND, JDIM+NBOUND
               PROPSC(J,K,L) = PROPS(J,K,I,L)
  810       CONTINUE
  820    CONTINUE
  830 CONTINUE
C
C     Finished with derivatives
C
      RETURN
      END





      SUBROUTINE DEL2KC (IDIM, JDIM, KDIM, I, ITURB,
     1                   Q, S, VOL, IBCDAT, NBCS,
     2                   TKEJ, TKEK, TKEPI, TKEMI, 
     3                   DKDX, DKDY, DKDZ,  DEL2K)
C
C Routine to calculate the second derivative of k term required in several 
C turbulence models.  This routine approximates the term as:
C          2 {|del[sqrt(k)]|}^2
C
C IDIM,JDIM,KDIM : Dimensions of current block
C I              : Current I Plane
C ITURB          : Turbulence model
C Q              : Primitive variables at cell centers
C S              : Metrics
C VOL            : Cell volumes
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C TKEJ,TKEK,
C TKEPI,TKEMI    : Temporary storage of K on the J,K,I+,I- faces
C DEL2K          : Approximation to Laplacian of K
C
C This routine is arranged as:
C   1. Evaluate k at faces of computational cell as averages of cell centers
C      a. Two-equation
C      b. Reynolds stress
C   2. Correct k at faces for walls
C   3. Correct streamwise faces for marching scheme
C   4. Evaluate first derivative using Gauss' divergence theorem
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     2          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION TKEJ  (JDIM,2:KDIM),
     1          TKEK  (2:JDIM,KDIM),
     2          TKEPI (2:JDIM,2:KDIM),
     3          TKEMI (2:JDIM,2:KDIM)
C
      DIMENSION DKDX  (2:JDIM,2:KDIM),
     1          DKDY  (2:JDIM,2:KDIM),
     2          DKDZ  (2:JDIM,2:KDIM)
C
      DIMENSION DEL2K (2:JDIM,2:KDIM)
C
      DIMENSION IBCDAT(10,NBCS)
C
C 1. Calculate square root of K at cell interfaces as averages of cell centers
C
C 1-a. Two-Equation model
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
C
C Calculate I interfaces 
C
         DO 110 K = 2, KDIM
            DO 105 J = 2, JDIM
               TKEPI(J,K) = SQRT (0.5E0*(Q(J,K,I,6) + Q(J,K,I+1,6)))
               TKEMI(J,K) = SQRT (0.5E0*(Q(J,K,I,6) + Q(J,K,I-1,6)))
  105       CONTINUE
  110    CONTINUE
C
C Calculate J interfaces
C
         DO 120 K = 2, KDIM
            DO 115 J = 1, JDIM
               TKEJ(J,K) = SQRT (0.5E0*(Q(J,K,I,6) + Q(J+1,K,I,6)))
  115       CONTINUE
  120    CONTINUE
C
C Calculate K interfaces 
C
         IF (THREED .OR. AXISYM) THEN
            DO 130 J = 2, JDIM
               DO 125 K = 1, KDIM
                  TKEK(J,K) = SQRT (0.5E0*(Q(J,K,I,6) + Q(J,K+1,I,6)))
  125          CONTINUE
  130       CONTINUE
         ENDIF
C
C 1-b. Reynolds Stress model (K = (T_xx + T_yy + T_zz) / 2)
C
      ELSE IF (ITURB .EQ. ITRS) THEN
C
C Calculate I interfaces 
C
         DO 160 K = 2, KDIM
            DO 155 J = 2, JDIM
               TKEPI(J,K) = SQRT (0.25E0 * 
     1                ( Q(J,K,I  ,6) + Q(J,K,I  ,7) + Q(J,K,I  ,8)
     2                + Q(J,K,I+1,6) + Q(J,K,I+1,7) + Q(J,K,I+1,8) ) )
               TKEMI(J,K) = SQRT (0.25E0 * 
     1                ( Q(J,K,I  ,6) + Q(J,K,I  ,7) + Q(J,K,I  ,8)
     2                + Q(J,K,I-1,6) + Q(J,K,I-1,7) + Q(J,K,I-1,8) ) )
  155       CONTINUE
  160    CONTINUE
C
C Calculate J interfaces
C
         DO 170 K = 2, KDIM
            DO 165 J = 1, JDIM
               TKEJ(J,K) = SQRT (
     1               0.25E0*( Q(J  ,K,I,6)+Q(J  ,K,I,7)+Q(J  ,K,I,8) 
     2                      + Q(J+1,K,I,6)+Q(J+1,K,I,7)+Q(J+1,K,I,8) ) )
  165       CONTINUE
  170    CONTINUE
C
C Calculate K interfaces 
C
         IF (THREED .OR. AXISYM) THEN
            DO 180 J = 2, JDIM
               DO 175 K = 1, KDIM
                  TKEK(J,K) = SQRT (
     1               0.25E0*( Q(J,K  ,I,6)+Q(J,K  ,I,7)+Q(J,K  ,I,8) 
     2                      + Q(J,K+1,I,6)+Q(J,K+1,I,7)+Q(J,K+1,I,8) ) )
  175          CONTINUE
  180       CONTINUE
         ENDIF
C
C Error - Invalid Turbulence model
C
      ELSE
         WRITE (IOUT,1000) ITURB
         STOP
      ENDIF
C
C 2. Loop over boundary conditions and fix the wall data
C
      DO 250 IBC = 1, NBCS
         IBCTYP = IBCDAT(1,IBC)
         IF (IBCTYP .EQ. IWALL .OR. IBCTYP .EQ. IWALFN) THEN
            IDIR   = IBCDAT(2,IBC)
            INROUT = IBCDAT(3,IBC)
            ISTRT  = IBCDAT(4,IBC)
            IEND   = IBCDAT(5,IBC)
            JSTRT  = IBCDAT(6,IBC)
            JEND   = IBCDAT(7,IBC)
            KSTRT  = IBCDAT(8,IBC)
            KEND   = IBCDAT(9,IBC)
C
C I direction
C
            IF (IDIR .EQ. 1) THEN
               IF (I .EQ. 2 .AND. INROUT .EQ. 1) THEN
                  DO 210 K = KSTRT, KEND
                     DO 205 J = JSTRT, JEND
                        TKEMI(J,K) = 0.E0
  205                CONTINUE
  210             CONTINUE
               ELSE IF (I .EQ. IDIM .AND. INROUT .NE. 1) THEN
                  DO 220 K = KSTRT, KEND
                     DO 215 J = JSTRT, JEND
                        TKEPI(J,K) = 0.E0
  215                CONTINUE
  220             CONTINUE
               ENDIF
C
C J direction
C
            ELSE IF (IDIR .EQ. 2) THEN
               IF (INROUT .EQ. 1) THEN
                  JNTR  = 1
                  JBCPT = 1
               ELSE
                  JNTR  = JDIM
                  JBCPT = JDIM + 1
               ENDIF
               IF (I .GE. ISTRT .AND. I .LE. IEND) THEN
                  DO 230 K = KSTRT, KEND
                     TKEJ(JNTR,K) = 0.E0
  230             CONTINUE
               ENDIF
C
C K direction
C
            ELSE IF (IDIR .EQ. 3) THEN
               IF (INROUT .EQ. 1) THEN
                  KNTR  = 1
                  KBCPT = 1
               ELSE
                  KNTR  = KDIM
                  KBCPT = KDIM + 1
               ENDIF
               IF (I .GE. ISTRT .AND. I .LE. IEND) THEN
                  DO 240 J = JSTRT, JEND
                     TKEK(J,KNTR) = 0.E0
  240             CONTINUE
               ENDIF
            ENDIF
         ENDIF
  250 CONTINUE
C
C 3. If using the MARCHing algorithm, set the TKEPI and TKEMI equal
C    to remove the spatial variation of K which will not be correct as 
C    the downstream values are not updated
C
      IF (ISOLVR .EQ. IMARCH) THEN
         DO 310 K = 2, KDIM
            DO 305 J = 2, JDIM
               TKEPI(J,K) = TKEMI(J,K)
  305       CONTINUE
  310    CONTINUE
      ENDIF
C
C 4. Loop over J and K forming the derivatives at the cell centers
C    Calculate first derivative terms using Divergence Theorem
C    Use symmetric differencing of the cross-derivative terms
C    from Chakravarthy et al AIAA 85-0165.
C
      DO 420 K = 2, KDIM
         DO 410 J = 2, JDIM
C
            DKDX(J,K) =( TKEPI(J,K)  * S(J,K,I  ,1,1) * S(J,K,I  ,4,1) 
     1                 - TKEMI(J,K)  * S(J,K,I-1,1,1) * S(J,K,I-1,4,1)
     2                 + TKEJ(J,K)   * S(J  ,K,I,1,2) * S(J  ,K,I,4,2)  
     3                 - TKEJ(J-1,K) * S(J-1,K,I,1,2) * S(J-1,K,I,4,2) )
            DKDY(J,K) =( TKEPI(J,K)  * S(J,K,I  ,2,1) * S(J,K,I  ,4,1) 
     1                 - TKEMI(J,K)  * S(J,K,I-1,2,1) * S(J,K,I-1,4,1)
     2                 + TKEJ(J,K)   * S(J  ,K,I,2,2) * S(J  ,K,I,4,2)  
     3                 - TKEJ(J-1,K) * S(J-1,K,I,2,2) * S(J-1,K,I,4,2) )
            DKDZ(J,K) = 0.E0
  410    CONTINUE
  420 CONTINUE
C
C     Three-Dimensional Section
C
      IF (THREED .OR. AXISYM) THEN
         DO 440 K = 2, KDIM
            DO 430 J = 2, JDIM
               DKDX(J,K) = DKDX(J,K) +
     1                 ( TKEK(J,K)   * S(J,K  ,I,1,3) * S(J,K  ,I,4,3)  
     2                 - TKEK(J,K-1) * S(J,K-1,I,1,3) * S(J,K-1,I,4,3) )
               DKDY(J,K) = DKDY(J,K) +
     1                 ( TKEK(J,K)   * S(J,K  ,I,2,3) * S(J,K  ,I,4,3)  
     2                 - TKEK(J,K-1) * S(J,K-1,I,2,3) * S(J,K-1,I,4,3) )
               DKDZ(J,K) = 
     1                 ( TKEPI(J,K)  * S(J,K,I  ,3,1) * S(J,K,I  ,4,1) 
     2                 - TKEMI(J,K)  * S(J,K,I-1,3,1) * S(J,K,I-1,4,1)
     3                 + TKEJ(J,K)   * S(J  ,K,I,3,2) * S(J  ,K,I,4,2)  
     4                 - TKEJ(J-1,K) * S(J-1,K,I,3,2) * S(J-1,K,I,4,2) 
     5                 + TKEK(J,K)   * S(J,K  ,I,3,3) * S(J,K  ,I,4,3)  
     6                 - TKEK(J,K-1) * S(J,K-1,I,3,3) * S(J,K-1,I,4,3) )
  430       CONTINUE
  440    CONTINUE
      ENDIF
C
C     Divide by volume and form 2*|del[sqrt(k)]|^2
C     Magnitude of del w is sqrt of sum of squares.
C
      DO 460 K = 2, KDIM
         DO 450 J = 2, JDIM
            VOLI  = 1.E0 / VOL(J,K,I)
            DRKDX = DKDX(J,K) * VOLI
            DRKDY = DKDY(J,K) * VOLI
            DRKDZ = DKDZ(J,K) * VOLI
            DEL2K(J,K) = 2.E0 * (DRKDX*DRKDX + DRKDY*DRKDY + 
     1                           DRKDZ*DRKDZ)
  450    CONTINUE
  460 CONTINUE
C
C     Finished with K derivative
C
      RETURN
 1000 FORMAT (' ','Invalid turbulence model ',I3,' specified in DEL2KC',
     1       /' ','Run Aborting.')
      END





      SUBROUTINE DELQI (IDIM, JDIM, KDIM, IPSTRT, IPEND, K,
     1                  Q, PROPS, S, VOL, IBCDAT, NBCS,
     2                  STEMP, VOLAUX, QI, PROPSI,
     3                  DQDX, DQDY, DQDZ)
C
C     Routine to calculate derivatives of the primitive variables at
C     cell interfaces and to calculate interface averages for the 
C     I direction.
C     Thin Layer:         Calculate DQDX,DQDY,DQDZ as metric times DQDXSI
C     Full Navier-Stokes: Use divergence theorem to calculate derivatives
C
C     This routine is not presently suited for the multiblock
C     approach.  It is close to being ready but is not.  The problem
C     lies that the metrics for one cell into the adjacent block
C     are needed to be completely general.  These are currently not available.
C     The routine metbc needs to be completed.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     IPSTRT,IPEND   : I plane to start, I plane to end on
C     K              : Current K Plane
C     Q              : Primitive variables at cell centers
C     PROPS          : Properties at cell centers
C                      PROPS(1) = RMU
C                      PROPS(2) = RMUT
C                      PROPS(3) = Y+
C     S              : Metrics
C     VOL            : Cell volumes
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     STEMP          : Temporary storage for metrics 
C     QI             : Primitive variables at cell interfaces
C     PROPSI         : Properties at cell interfaces
C     DQDX,DQDY...   : Derivatives of primitive variables at interfaces.
C                      This routine calculates either thin-layer or full NS
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP), 
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION STEMP (IPSTRT-1:IPEND+1,JDIM,3),
     1          VOLAUX(IPSTRT-1:IPEND,2:JDIM)
C
      DIMENSION QI    (IPSTRT-1:IPEND,2:JDIM,NQ), 
     1          PROPSI(IPSTRT-1:IPEND,2:JDIM,NP)
C
      DIMENSION DQDX  (IPSTRT-1:IPEND,2:JDIM,NQ), 
     1          DQDY  (IPSTRT-1:IPEND,2:JDIM,NQ), 
     2          DQDZ  (IPSTRT-1:IPEND,2:JDIM,NQ)
C
      DIMENSION IBCDAT(10,NBCS)
C
C     Loop over J forming derivative terms along the I direction
C     on the I = IPSTRT-1, IPEND faces.  
C     The storage of SJ reduces the overhead of recalculating the average
C     metric terms required for the derivative calculation.
C
C     VOLAUX is the inverse of the average auxilliary cell volume
C
      DO 120 J = 2, JDIM
         DO 110 I = IPSTRT-1, IPEND
            VOLAUX(I,J) = 2.0E0 / ( VOL(J,K,I) + VOL(J,K,I+1) )
  110    CONTINUE
  120 CONTINUE
C
C     For WALL FUNCTIONS ONLY: Modify the auxilliary cell volume at the
C     wall function cell interface to reflect the wall function grid.
C
      IF ( IFWALF ) THEN
C
C     Loop over boundary conditions and fix the wall function volume
C
         DO 125 IBC = 1, NBCS
            IBCTYP = IBCDAT(1,IBC)
            IDIR   = IBCDAT(2,IBC)
            IF (IBCTYP .EQ. IWALFN .AND. IDIR .EQ. 1) THEN
               INROUT = IBCDAT(3,IBC)
               JSTRT  = IBCDAT(6,IBC)
               JEND   = IBCDAT(7,IBC)
               KSTRT  = IBCDAT(8,IBC)
               KEND   = IBCDAT(9,IBC)
               IF (INROUT .EQ. 1) THEN
                  IINT   = 2
                  IBCN2  = 3
               ELSE
                  IINT   = IDIM - 1
                  IBCN2  = IDIM - 2
               ENDIF
               IF (IINT .GE. IPSTRT-1 .AND. IINT .LE. IPEND .AND.
     1                K .GE. KSTRT    .AND.    K .LE. KEND) THEN
                  DO 122 J = JSTRT, JEND
                     VOLAUX(IINT,J) = 1.0E0 / VOL(J,K,IBCN2)
  122             CONTINUE
               ENDIF
            ENDIF
  125    CONTINUE
C
      ENDIF
C
C     Thin Layer Derivatives
C
      IF ( THNLYR ) THEN
         DO 150 L = 1, NQ
            DO 140 J = 2, JDIM
               DO 130 I = IPSTRT-1, IPEND
                  DXSI        = (Q(J,K,I+1,L) - Q(J,K,I,L)) 
     1                          * VOLAUX(I,J) * S(J,K,I,4,1)
                  DQDX(I,J,L) = S(J,K,I,1,1) * DXSI
                  DQDY(I,J,L) = S(J,K,I,2,1) * DXSI
                  DQDZ(I,J,L) = S(J,K,I,3,1) * DXSI
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
      ELSE
C
C     Full Navier-Stokes terms
C     Calculate first derivative terms using Divergence Theorem
C     Use symmetric differencing of the cross-derivative terms
C     from Chakravarthy et al AIAA 85-0165.
C
C     Average metrics in I direction 
C
         DO 220 J = 2, JDIM
            DO 210 I = IPSTRT-1, IPEND+1
               STEMP(I,J,1) = 0.5E0 * ( S(J,K,I-1,1,1)*S(J,K,I-1,4,1)
     1                                + S(J,K,I,1,1)*S(J,K,I,4,1) )
               STEMP(I,J,2) = 0.5E0 * ( S(J,K,I-1,2,1)*S(J,K,I-1,4,1)
     1                                + S(J,K,I,2,1)*S(J,K,I,4,1) )
               STEMP(I,J,3) = 0.5E0 * ( S(J,K,I-1,3,1)*S(J,K,I-1,4,1)
     1                                + S(J,K,I,3,1)*S(J,K,I,4,1) )
  210       CONTINUE
  220    CONTINUE
C
C     Calculate first derivative terms using Divergence Theorem
C
         DO 250 L = 1, NQ
            DO 240 J = 2, JDIM
               DO 230 I = IPSTRT-1, IPEND
                  DQDX(I,J,L) = VOLAUX(I,J)*(Q(J,K,I+1,L)*STEMP(I+1,J,1)
     1                                      -Q(J,K,I,L)*STEMP(I,J,1) )
                  DQDY(I,J,L) = VOLAUX(I,J)*(Q(J,K,I+1,L)*STEMP(I+1,J,2)
     1                                      -Q(J,K,I,L)*STEMP(I,J,2) )
                  DQDZ(I,J,L) = VOLAUX(I,J)*(Q(J,K,I+1,L)*STEMP(I+1,J,3)
     1                                      -Q(J,K,I,L)*STEMP(I,J,3) )
  230          CONTINUE
  240       CONTINUE
  250    CONTINUE
C
C     Average metric terms in J direction
C
         DO 320 J = 1, JDIM
            DO 310 I = IPSTRT-1, IPEND
               STEMP(I,J,1) = 0.5E0 * ( S(J,K,I,1,2)*S(J,K,I,4,2)
     1                                + S(J,K,I+1,1,2)*S(J,K,I+1,4,2) )
               STEMP(I,J,2) = 0.5E0 * ( S(J,K,I,2,2)*S(J,K,I,4,2)
     1                                + S(J,K,I+1,2,2)*S(J,K,I+1,4,2) )
               STEMP(I,J,3) = 0.5E0 * ( S(J,K,I,3,2)*S(J,K,I,4,2)
     1                                + S(J,K,I+1,3,2)*S(J,K,I+1,4,2) )
  310       CONTINUE
  320    CONTINUE
C
C     Add the derivative data from the J direction
C
         DO 350 L = 1, NQ
            DO 340 J = 2, JDIM
               DO 330 I = IPSTRT-1, IPEND
                  QPJ = 0.25E0 * ( Q(J,K,I,L) + Q(J+1,K,I+1,L)
     1                           + Q(J+1,K,I,L) + Q(J,K,I+1,L) )
                  QMJ = 0.25E0 * ( Q(J,K,I,L) + Q(J-1,K,I+1,L)
     1                           + Q(J-1,K,I,L) + Q(J,K,I+1,L) )
C
                  DQDX(I,J,L) = DQDX(I,J,L)
     1                        + VOLAUX(I,J) * ( QPJ*STEMP(I,J,1) 
     2                                        - QMJ*STEMP(I,J-1,1) )
                  DQDY(I,J,L) = DQDY(I,J,L)
     1                        + VOLAUX(I,J) * ( QPJ*STEMP(I,J,2) 
     2                                        - QMJ*STEMP(I,J-1,2) )
                  DQDZ(I,J,L) = DQDZ(I,J,L)
     1                        + VOLAUX(I,J) * ( QPJ*STEMP(I,J,3) 
     2                                        - QMJ*STEMP(I,J-1,3) )
  330          CONTINUE
  340       CONTINUE
  350    CONTINUE
C
C     Add the contribution from the K+ term
C
C     Average metric terms in K+ direction
C
         IF (THREED .OR. AXISYM) THEN
            DO 420 J = 2, JDIM
               DO 410 I = IPSTRT-1, IPEND
                  STEMP(I,J,1) = 0.5E0*( S(J,K,I,1,3)*S(J,K,I,4,3)
     1                                 + S(J,K,I+1,1,3)*S(J,K,I+1,4,3) )
                  STEMP(I,J,2) = 0.5E0*( S(J,K,I,2,3)*S(J,K,I,4,3)
     1                                 + S(J,K,I+1,2,3)*S(J,K,I+1,4,3) )
                  STEMP(I,J,3) = 0.5E0*( S(J,K,I,3,3)*S(J,K,I,4,3)
     1                                 + S(J,K,I+1,3,3)*S(J,K,I+1,4,3) )
  410          CONTINUE
  420       CONTINUE
C
C     Add the derivative data from the K+ term
C
            DO 450 L = 1, NQ
               DO 440 J = 2, JDIM
                  DO 430 I = IPSTRT-1, IPEND
                     QPK = 0.25E0 * ( Q(J,K,I,L) + Q(J,K+1,I+1,L)
     1                              + Q(J,K+1,I,L) + Q(J,K,I+1,L) )
C
                     DQDX(I,J,L) = DQDX(I,J,L)
     1                           + VOLAUX(I,J) * ( QPK*STEMP(I,J,1) )
                     DQDY(I,J,L) = DQDY(I,J,L)
     1                           + VOLAUX(I,J) * ( QPK*STEMP(I,J,2) )
                     DQDZ(I,J,L) = DQDZ(I,J,L)
     1                           + VOLAUX(I,J) * ( QPK*STEMP(I,J,3) )
  430             CONTINUE
  440          CONTINUE
  450       CONTINUE
C
C     Add the contribution from the K- term
C
C     Average metric terms in K- direction
C
            DO 520 J = 2, JDIM
               DO 510 I = IPSTRT-1, IPEND
                  STEMP(I,J,1) = 0.5E0 * ( S(J,K-1,I,1,3)*S(J,K-1,I,4,3)
     1                             + S(J,K-1,I+1,1,3)*S(J,K-1,I+1,4,3) )
                  STEMP(I,J,2) = 0.5E0 * ( S(J,K-1,I,2,3)*S(J,K-1,I,4,3)
     1                             + S(J,K-1,I+1,2,3)*S(J,K-1,I+1,4,3) )
                  STEMP(I,J,3) = 0.5E0 * ( S(J,K-1,I,3,3)*S(J,K-1,I,4,3)
     1                             + S(J,K-1,I+1,3,3)*S(J,K-1,I+1,4,3) )
  510          CONTINUE
  520       CONTINUE
C
C     Add the derivative data from the K- term
C
            DO 550 L = 1, NQ
               DO 540 J = 2, JDIM
                  DO 530 I = IPSTRT-1, IPEND
                     QMK = 0.25E0 * ( Q(J,K,I,L) + Q(J,K-1,I+1,L)
     1                              + Q(J,K,I+1,L) + Q(J,K-1,I,L) )
C
                     DQDX(I,J,L) = DQDX(I,J,L)
     1                           - VOLAUX(I,J) * ( QMK*STEMP(I,J,1) )
                     DQDY(I,J,L) = DQDY(I,J,L)
     1                           - VOLAUX(I,J) * ( QMK*STEMP(I,J,2) )
                     DQDZ(I,J,L) = DQDZ(I,J,L)
     1                           - VOLAUX(I,J) * ( QMK*STEMP(I,J,3) )
  530             CONTINUE
  540          CONTINUE
  550       CONTINUE
         ENDIF
      ENDIF
C
C     Calculate Q at cell interfaces as averages of cell centers
C
      DO 630 L = 1, NQ
         DO 620 J = 2, JDIM
            DO 610 I = IPSTRT-1, IPEND
               QI(I,J,L) = 0.5E0 * ( Q(J,K,I,L)  + Q(J,K,I+1,L) )
  610       CONTINUE
  620    CONTINUE
  630 CONTINUE
C
C     Calculate PROPSI at cell interfaces as averages of cell centers
C
      DO 660 L = 1, NP
         DO 650 J = 2, JDIM
            DO 640 I = IPSTRT-1, IPEND
               PROPSI(I,J,L)  = 0.5E0*(PROPS(J,K,I,L)+PROPS(J,K,I+1,L))
  640       CONTINUE
  650    CONTINUE
  660 CONTINUE
C
C     Loop over boundary conditions and fix the wall data
C
      DO 750 IBC = 1, NBCS
         IBCTYP = IBCDAT(1,IBC)
         IDIR   = IBCDAT(2,IBC)
         IF ( (IBCTYP .EQ. IWALL .OR. IBCTYP .EQ. IWALFN)
     1        .AND. IDIR .EQ. 1) THEN
            INROUT = IBCDAT(3,IBC)
            JSTRT  = IBCDAT(6,IBC)
            JEND   = IBCDAT(7,IBC)
            KSTRT  = IBCDAT(8,IBC)
            KEND   = IBCDAT(9,IBC)
            IF (INROUT .EQ. 1) THEN
               IBCPT = 1
               INTR  = 1
            ELSE
               IBCPT = IDIM + 1
               INTR  = IDIM
            ENDIF
            IF (IBCPT .GE. IPSTRT-1 .AND. IBCPT .LE. IPEND .AND.
     1              K .GE. KSTRT    .AND.     K .LE. KEND) THEN
               DO 720 L = 1, NQ
                  DO 710 J = JSTRT, JEND
                     QI(INTR,J,L) = Q(J,K,IBCPT,L)
  710             CONTINUE
  720          CONTINUE
C
C     Also do properties
C
               DO 740 L = 1, NP
                  DO 730 J = JSTRT, JEND
                     PROPSI(INTR,J,L) = PROPS(J,K,IBCPT,L)
  730             CONTINUE
  740          CONTINUE
            ENDIF
         ENDIF
  750 CONTINUE
C
C     Finished calculating derivatives in the I direction
C
      RETURN
      END





      SUBROUTINE DELQJ (IDIM, JDIM, KDIM, I,
     1                  Q, PROPS, S, VOL, IBCDAT, NBCS,
     2                  STEMP, VOLAUX, QI, PROPSI,
     3                  DQDX, DQDY, DQDZ)
C
C     Routine to calculate derivatives of the primitive variables at
C     cell interfaces and to calculate interface averages for the 
C     J direction.
C     Thin Layer:         Calculate DQDX,DQDY,DQDZ as metric times DQDETA
C     Full Navier-Stokes: Use divergence theorem to calculate derivatives
C
C     This routine is not presently suited for the multiblock
C     approach.  It is close to being ready but is not.  The problem
C     lies that the metrics for one cell into the adjacent block
C     are needed to be completely general.  These are currently not available.
C     The routine metbc needs to be completed.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : Current I Plane
C     Q              : Primitive variables at cell centers
C     PROPS          : Properties at cell centers
C                      PROPS(1) = RMU
C                      PROPS(2) = RMUT
C                      PROPS(3) = Y+
C     S              : Metrics
C     VOL            : Cell volumes
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     STEMP          : Temporary storage for metrics 
C     QI             : Primitive variables at cell interfaces
C     PROPSI         : Properties at cell interfaces
C     DQDX,DQDY...   : Derivatives of primitive variables at interfaces.
C                      This routine calculates either thin-layer or full NS
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP), 
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION STEMP (JDIM+1,KDIM,3),
     1          VOLAUX(JDIM,2:KDIM)
C
      DIMENSION QI    (JDIM,2:KDIM,NQ), 
     1          PROPSI(JDIM,2:KDIM,NP)
C
      DIMENSION DQDX  (JDIM,2:KDIM,NQ), 
     1          DQDY  (JDIM,2:KDIM,NQ), 
     2          DQDZ  (JDIM,2:KDIM,NQ)
C
      DIMENSION IBCDAT(10,NBCS)
C
C     Loop over K forming derivative terms along the J direction
C     on the J = 1, JDIM faces.  
C     The storage of SJ reduces the overhead of recalculating the average
C     metric terms required for the derivative calculation.
C
C     VOLAUX is the inverse of the average auxilliary cell volume
C
      DO 120 K = 2, KDIM
         DO 110 J = 1, JDIM
            VOLAUX(J,K) = 2.0E0 / ( VOL(J,K,I) + VOL(J+1,K,I) )
  110    CONTINUE
  120 CONTINUE
C
C     For WALL FUNCTIONS ONLY: Modify the auxilliary cell volume at the
C     wall function cell interface to reflect the wall function grid.
C
      IF ( IFWALF ) THEN
C
C     Loop over boundary conditions and fix the wall function volume
C
         DO 125 IBC = 1, NBCS
            IBCTYP = IBCDAT(1,IBC)
            IDIR   = IBCDAT(2,IBC)
            IF (IBCTYP .EQ. IWALFN .AND. IDIR .EQ. 2) THEN
               INROUT = IBCDAT(3,IBC)
               ISTRT  = IBCDAT(4,IBC)
               IEND   = IBCDAT(5,IBC)
               KSTRT  = IBCDAT(8,IBC)
               KEND   = IBCDAT(9,IBC)
               IF (INROUT .EQ. 1) THEN
                  JINT   = 2
                  JBCN2  = 3
               ELSE
                  JINT   = JDIM - 1
                  JBCN2  = JDIM - 2
               ENDIF
               IF (I .GE. ISTRT .AND. I .LE. IEND) THEN
                  DO 122 K = KSTRT, KEND
                     VOLAUX(JINT,K) = 1.0E0 / VOL(JBCN2,K,I)
  122             CONTINUE
               ENDIF
            ENDIF
  125    CONTINUE
C
      ENDIF
C
C     End of corrections for wall functions
C
C     Thin Layer Derivatives
C
      IF ( THNLYR ) THEN
         DO 150 L = 1, NQ
            DO 140 K = 2, KDIM
               DO 130 J = 1, JDIM
                  DETA        = (Q(J+1,K,I,L) - Q(J,K,I,L))
     1                          * VOLAUX(J,K) * S(J,K,I,4,2)
                  DQDX(J,K,L) = S(J,K,I,1,2) * DETA
                  DQDY(J,K,L) = S(J,K,I,2,2) * DETA
                  DQDZ(J,K,L) = S(J,K,I,3,2) * DETA
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
      ELSE
C
C     Full Navier-Stokes terms
C     Calculate first derivative terms using Divergence Theorem
C     Use symmetric differencing of the cross-derivative terms
C     from Chakravarthy et al AIAA 85-0165.
C
C     Average metrics in J direction 
C
         DO 220 K = 2, KDIM
            DO 210 J = 1, JDIM+1
               STEMP(J,K,1) = 0.5E0 * ( S(J-1,K,I,1,2)*S(J-1,K,I,4,2)
     1                                + S(J,K,I,1,2)*S(J,K,I,4,2) )
               STEMP(J,K,2) = 0.5E0 * ( S(J-1,K,I,2,2)*S(J-1,K,I,4,2)
     1                                + S(J,K,I,2,2)*S(J,K,I,4,2) )
               STEMP(J,K,3) = 0.5E0 * ( S(J-1,K,I,3,2)*S(J-1,K,I,4,2)
     1                                + S(J,K,I,3,2)*S(J,K,I,4,2) )
  210       CONTINUE
  220    CONTINUE
C
C     Calculate first derivative terms using Divergence Theorem
C
         DO 250 L = 1, NQ
            DO 240 K = 2, KDIM
               DO 230 J = 1, JDIM
                  DQDX(J,K,L) = VOLAUX(J,K)*(Q(J+1,K,I,L)*STEMP(J+1,K,1) 
     1                                      -Q(J,K,I,L)*STEMP(J,K,1) )
                  DQDY(J,K,L) = VOLAUX(J,K)*(Q(J+1,K,I,L)*STEMP(J+1,K,2) 
     1                                      -Q(J,K,I,L)*STEMP(J,K,2) )
                  DQDZ(J,K,L) = VOLAUX(J,K)*(Q(J+1,K,I,L)*STEMP(J+1,K,3) 
     1                                      -Q(J,K,I,L)*STEMP(J,K,3) )
  230          CONTINUE
  240       CONTINUE
  250    CONTINUE
C
C     Average metric terms in K direction
C
         IF (THREED .OR. AXISYM) THEN
            DO 320 K = 1, KDIM
               DO 310 J = 1, JDIM
                  STEMP(J,K,1) = 0.5E0*( S(J,K,I,1,3)*S(J,K,I,4,3)
     1                                 + S(J+1,K,I,1,3)*S(J+1,K,I,4,3) )
                  STEMP(J,K,2) = 0.5E0*( S(J,K,I,2,3)*S(J,K,I,4,3)
     1                                 + S(J+1,K,I,2,3)*S(J+1,K,I,4,3) )
                  STEMP(J,K,3) = 0.5E0*( S(J,K,I,3,3)*S(J,K,I,4,3)
     1                                 + S(J+1,K,I,3,3)*S(J+1,K,I,4,3) )
  310          CONTINUE
  320       CONTINUE
C
C     Add the derivative data from the K direction
C
            DO 350 L = 1, NQ
               DO 340 K = 2, KDIM
                  DO 330 J = 1, JDIM
                     QPK = 0.25E0 * ( Q(J,K,I,L) + Q(J+1,K+1,I,L)
     1                              + Q(J+1,K,I,L) + Q(J,K+1,I,L) )
                     QMK = 0.25E0 * ( Q(J,K,I,L) + Q(J+1,K-1,I,L)
     1                              + Q(J+1,K,I,L) + Q(J,K-1,I,L) )
C
                     DQDX(J,K,L) = DQDX(J,K,L)
     1                           + VOLAUX(J,K) * ( QPK*STEMP(J,K,1) 
     2                                           - QMK*STEMP(J,K-1,1) )
                     DQDY(J,K,L) = DQDY(J,K,L)
     1                           + VOLAUX(J,K) * ( QPK*STEMP(J,K,2) 
     2                                           - QMK*STEMP(J,K-1,2) )
                     DQDZ(J,K,L) = DQDZ(J,K,L)
     1                           + VOLAUX(J,K) * ( QPK*STEMP(J,K,3) 
     2                                           - QMK*STEMP(J,K-1,3) )
  330             CONTINUE
  340          CONTINUE
  350       CONTINUE
         ENDIF
C
C     Add the contribution from the I+ term
C
C     Average metric terms in I+ direction
C
         DO 420 K = 2, KDIM
            DO 410 J = 1, JDIM
               STEMP(J,K,1) = 0.5E0 * ( S(J,K,I,1,1)*S(J,K,I,4,1)
     1                                + S(J+1,K,I,1,1)*S(J+1,K,I,4,1) )
               STEMP(J,K,2) = 0.5E0 * ( S(J,K,I,2,1)*S(J,K,I,4,1)
     1                                + S(J+1,K,I,2,1)*S(J+1,K,I,4,1) )
               STEMP(J,K,3) = 0.5E0 * ( S(J,K,I,3,1)*S(J,K,I,4,1)
     1                                + S(J+1,K,I,3,1)*S(J+1,K,I,4,1) )
  410       CONTINUE
  420    CONTINUE
C
C     Add the derivative data from the I+ term
C
         DO 450 L = 1, NQ
            DO 440 K = 2, KDIM
               DO 430 J = 1, JDIM
                  QPI = 0.25E0 * ( Q(J,K,I,L) + Q(J+1,K,I+1,L)
     1                           + Q(J+1,K,I,L) + Q(J,K,I+1,L) )
C
                  DQDX(J,K,L) = DQDX(J,K,L)
     1                        + VOLAUX(J,K) * ( QPI*STEMP(J,K,1) )
                  DQDY(J,K,L) = DQDY(J,K,L)
     1                        + VOLAUX(J,K) * ( QPI*STEMP(J,K,2) )
                  DQDZ(J,K,L) = DQDZ(J,K,L)
     1                        + VOLAUX(J,K) * ( QPI*STEMP(J,K,3) )
  430          CONTINUE
  440       CONTINUE
  450    CONTINUE
C
C     Add the contribution from the I- term
C
C     Average metric terms in I- direction
C
         DO 520 K = 2, KDIM
            DO 510 J = 1, JDIM
               STEMP(J,K,1) = 0.5E0 * ( S(J,K,I-1,1,1)*S(J,K,I-1,4,1)
     1                          + S(J+1,K,I-1,1,1)*S(J+1,K,I-1,4,1) )
               STEMP(J,K,2) = 0.5E0 * ( S(J,K,I-1,2,1)*S(J,K,I-1,4,1)
     1                          + S(J+1,K,I-1,2,1)*S(J+1,K,I-1,4,1) )
               STEMP(J,K,3) = 0.5E0 * ( S(J,K,I-1,3,1)*S(J,K,I-1,4,1)
     1                          + S(J+1,K,I-1,3,1)*S(J+1,K,I-1,4,1) )
  510       CONTINUE
  520    CONTINUE
C
C     Add the derivative data from the I- term
C
         DO 550 L = 1, NQ
            DO 540 K = 2, KDIM
               DO 530 J = 1, JDIM
                  QMI = 0.25E0 * ( Q(J,K,I,L) + Q(J+1,K,I-1,L)
     1                           + Q(J+1,K,I,L) + Q(J,K,I-1,L) )
C
                  DQDX(J,K,L) = DQDX(J,K,L)
     1                        - VOLAUX(J,K) * ( QMI*STEMP(J,K,1) )
                  DQDY(J,K,L) = DQDY(J,K,L)
     1                        - VOLAUX(J,K) * ( QMI*STEMP(J,K,2) )
                  DQDZ(J,K,L) = DQDZ(J,K,L)
     1                        - VOLAUX(J,K) * ( QMI*STEMP(J,K,3) )
  530          CONTINUE
  540       CONTINUE
  550    CONTINUE
      ENDIF
C
C     Calculate Q at cell interfaces as averages of cell centers
C
      DO 630 L = 1, NQ
         DO 620 K = 2, KDIM
            DO 610 J = 1, JDIM
               QI(J,K,L) = 0.5E0 * ( Q(J,K,I,L)  + Q(J+1,K,I,L) )
  610       CONTINUE
  620    CONTINUE
  630 CONTINUE
C
C     Calculate PROPSI at cell interfaces as averages of cell centers
C
      DO 660 L = 1, NP
         DO 650 K = 2, KDIM
            DO 640 J = 1, JDIM
               PROPSI(J,K,L)  = 0.5E0*(PROPS(J,K,I,L)+PROPS(J+1,K,I,L))
  640       CONTINUE
  650    CONTINUE
  660 CONTINUE
C
C     Loop over boundary conditions and fix the wall data
C
      DO 750 IBC = 1, NBCS
         IBCTYP = IBCDAT(1,IBC)
         IDIR   = IBCDAT(2,IBC)
         IF ( (IBCTYP .EQ. IWALL .OR. IBCTYP .EQ. IWALFN)
     1        .AND. IDIR .EQ. 2) THEN
            INROUT = IBCDAT(3,IBC)
            ISTRT  = IBCDAT(4,IBC)
            IEND   = IBCDAT(5,IBC)
            KSTRT  = IBCDAT(8,IBC)
            KEND   = IBCDAT(9,IBC)
            IF (INROUT .EQ. 1) THEN
               JBCPT = 1
               JNTR  = 1
            ELSE
               JBCPT = JDIM + 1
               JNTR  = JDIM
            ENDIF
            IF (I .GE. ISTRT .AND. I .LE. IEND) THEN
               DO 720 L = 1, NQ
                  DO 710 K = KSTRT, KEND
                     QI(JNTR,K,L) = Q(JBCPT,K,I,L)
  710             CONTINUE
  720          CONTINUE
C
C     Also do properties
C
               DO 740 L = 1, NP
                  DO 730 K = KSTRT, KEND
                     PROPSI(JNTR,K,L) = PROPS(JBCPT,K,I,L)
  730             CONTINUE
  740          CONTINUE
            ENDIF
         ENDIF
  750 CONTINUE
C
C     Finished calculating derivatives in the J direction
C
      RETURN
      END








      SUBROUTINE DELQK (IDIM, JDIM, KDIM, I,
     1                  Q, PROPS, S, VOL, IBCDAT, NBCS,
     2                  STEMP, VOLAUX, QI, PROPSI,
     3                  DQDX, DQDY, DQDZ)
C
C     Routine to calculate derivatives of the primitive variables at
C     cell interfaces and to calculate interface averages of primitive
C     variables in the K direction.
C     Thin Layer:         Calculate DQDX,DQDY,DQDZ as metric times DQDZETA
C     Full Navier-Stokes: Use divergence theorem to calculate derivatives
C
C     This routine is not presently suited for the multiblock
C     approach.  It is close to being ready but is not.  The problem
C     lies that the metrics for one cell into the adjacent block
C     are needed to be completely general.  These are currently not available.
C     The routine metbc needs to be completed.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : Current I Plane
C     Q              : Primitive variables at cell centers
C     PROPS          : Properties at cell centers
C                      PROPS(1) = RMU
C                      PROPS(2) = RMUT
C                      PROPS(3) = Y+
C     S              : Metrics
C     VOL            : Cell volumes
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     STEMP          : Temporary storage for metrics 
C     QI             : Primitive variables at cell interfaces
C     PROPSI         : Properties at cell interfaces
C     DQDX,DQDY...   : Derivatives of primitive variables at interfaces.
C                      This routine calculates either thin-layer or full NS
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP), 
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION STEMP (KDIM+1,JDIM,3),
     1          VOLAUX(KDIM,2:JDIM)
C
      DIMENSION QI    (KDIM,2:JDIM,NQ), 
     1          PROPSI(KDIM,2:JDIM,NP)
C
      DIMENSION DQDX  (KDIM,2:JDIM,NQ), 
     1          DQDY  (KDIM,2:JDIM,NQ), 
     2          DQDZ  (KDIM,2:JDIM,NQ)
C
      DIMENSION IBCDAT(10,NBCS)
C
C     Loop over J forming derivative terms along the K direction
C     on the K = 1, KDIM faces.  
C     The storage of SK reduces the overhead of recalculating the average
C     metric terms required for the derivative calculation.
C
C     VOLAUX is the inverse of the average auxilliary cell volume
C
      DO 120 J = 2, JDIM
         DO 110 K = 1, KDIM
            VOLAUX(K,J) = 2.0E0 / ( VOL(J,K,I) + VOL(J,K+1,I) )
  110    CONTINUE
  120 CONTINUE
C
C     For WALL FUNCTIONS ONLY: Modify the auxilliary cell volume at the
C     wall function cell interface to reflect the wall function grid.
C
      IF ( IFWALF ) THEN
C
C     Loop over boundary conditions and fix the wall function volume
C
         DO 125 IBC = 1, NBCS
            IBCTYP = IBCDAT(1,IBC)
            IDIR   = IBCDAT(2,IBC)
            IF (IBCTYP .EQ. IWALFN .AND. IDIR .EQ. 3) THEN
               INROUT = IBCDAT(3,IBC)
               ISTRT  = IBCDAT(4,IBC)
               IEND   = IBCDAT(5,IBC)
               JSTRT  = IBCDAT(6,IBC)
               JEND   = IBCDAT(7,IBC)
               IF (INROUT .EQ. 1) THEN
                  KINT   = 2
                  KBCN2  = 3
               ELSE
                  KINT   = KDIM - 1
                  KBCN2  = KDIM - 2
               ENDIF
               IF (I .GE. ISTRT .AND. I .LE. IEND) THEN
                  DO 122 J = JSTRT, JEND
                     VOLAUX(KINT,J) = 1.0E0 / VOL(J,KBCN2,I)
  122             CONTINUE
               ENDIF
            ENDIF
  125    CONTINUE
C
      ENDIF
C
C     Thin Layer Derivatives
C
      IF ( THNLYR ) THEN
         DO 150 L = 1, NQ
            DO 140 J = 2, JDIM
               DO 130 K = 1, KDIM
                  DZETA       = (Q(J,K+1,I,L) - Q(J,K,I,L))
     1                          * VOLAUX(K,J) * S(J,K,I,4,3)
                  DQDX(K,J,L) = S(J,K,I,1,3) * DZETA
                  DQDY(K,J,L) = S(J,K,I,2,3) * DZETA
                  DQDZ(K,J,L) = S(J,K,I,3,3) * DZETA
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
      ELSE
C
C     Full Navier-Stokes terms
C     Calculate first derivative terms using Divergence Theorem
C     Use symmetric differencing of the cross-derivative terms
C     from Chakravarthy et al AIAA 85-0165.
C
C     Average metrics in K direction 
C
         DO 220 J = 2, JDIM
            DO 210 K = 1, KDIM+1
               STEMP(K,J,1) = 0.5E0 * ( S(J,K-1,I,1,3)*S(J,K-1,I,4,3)
     1                                + S(J,K,I,1,3)*S(J,K,I,4,3) )
               STEMP(K,J,2) = 0.5E0 * ( S(J,K-1,I,2,3)*S(J,K-1,I,4,3)
     1                                + S(J,K,I,2,3)*S(J,K,I,4,3) )
               STEMP(K,J,3) = 0.5E0 * ( S(J,K-1,I,3,3)*S(J,K-1,I,4,3)
     1                                + S(J,K,I,3,3)*S(J,K,I,4,3) )
  210       CONTINUE
  220    CONTINUE
C
C     Calculate first derivative terms using Divergence Theorem
C
         DO 250 L = 1, NQ
            DO 240 J = 2, JDIM
               DO 230 K = 1, KDIM
                  DQDX(K,J,L) = VOLAUX(K,J)*(Q(J,K+1,I,L)*STEMP(K+1,J,1) 
     1                                      -Q(J,K,I,L)*STEMP(K,J,1) )
                  DQDY(K,J,L) = VOLAUX(K,J)*(Q(J,K+1,I,L)*STEMP(K+1,J,2) 
     1                                      -Q(J,K,I,L)*STEMP(K,J,2) )
                  DQDZ(K,J,L) = VOLAUX(K,J)*(Q(J,K+1,I,L)*STEMP(K+1,J,3) 
     1                                      -Q(J,K,I,L)*STEMP(K,J,3) )
  230          CONTINUE
  240       CONTINUE
  250    CONTINUE
C
C     Average metric terms in J direction
C
         DO 320 J = 1, JDIM
            DO 310 K = 1, KDIM
               STEMP(K,J,1) = 0.5E0 * ( S(J,K,I,1,2)*S(J,K,I,4,2)
     1                                + S(J,K+1,I,1,2)*S(J,K+1,I,4,2) )
               STEMP(K,J,2) = 0.5E0 * ( S(J,K,I,2,2)*S(J,K,I,4,2)
     1                                + S(J,K+1,I,2,2)*S(J,K+1,I,4,2) )
               STEMP(K,J,3) = 0.5E0 * ( S(J,K,I,3,2)*S(J,K,I,4,2)
     1                                + S(J,K+1,I,3,2)*S(J,K+1,I,4,2) )
  310       CONTINUE
  320    CONTINUE
C
C     Add the derivative data from the J direction
C
         DO 350 L = 1, NQ
            DO 340 J = 2, JDIM
               DO 330 K = 1, KDIM
                  QPJ = 0.25E0 * ( Q(J,K,I,L) + Q(J+1,K+1,I,L)
     1                           + Q(J+1,K,I,L) + Q(J,K+1,I,L) )
                  QMJ = 0.25E0 * ( Q(J,K,I,L) + Q(J-1,K+1,I,L)
     1                           + Q(J-1,K,I,L) + Q(J,K+1,I,L) )
C
                  DQDX(K,J,L) = DQDX(K,J,L)
     1                        + VOLAUX(K,J) * ( QPJ*STEMP(K,J,1) 
     2                                        - QMJ*STEMP(K,J-1,1) )
                  DQDY(K,J,L) = DQDY(K,J,L)
     1                        + VOLAUX(K,J) * ( QPJ*STEMP(K,J,2) 
     2                                        - QMJ*STEMP(K,J-1,2) )
                  DQDZ(K,J,L) = DQDZ(K,J,L)
     1                        + VOLAUX(K,J) * ( QPJ*STEMP(K,J,3) 
     2                                        - QMJ*STEMP(K,J-1,3) )
  330          CONTINUE
  340       CONTINUE
  350    CONTINUE
C
C     Add the contribution from the I+ term
C
C     Average metric terms in I+ direction
C
         DO 420 J = 2, JDIM
            DO 410 K = 1, KDIM
               STEMP(K,J,1) = 0.5E0 * ( S(J,K,I,1,1)*S(J,K,I,4,1) 
     1                                + S(J,K+1,I,1,1)*S(J,K+1,I,4,1) )
               STEMP(K,J,2) = 0.5E0 * ( S(J,K,I,2,1)*S(J,K,I,4,1)
     1                                + S(J,K+1,I,2,1)*S(J,K+1,I,4,1) )
               STEMP(K,J,3) = 0.5E0 * ( S(J,K,I,3,1)*S(J,K,I,4,1)
     1                                + S(J,K+1,I,3,1)*S(J,K+1,I,4,1) )
  410       CONTINUE
  420    CONTINUE
C
C     Add the derivative data from the I+ term
C
         DO 450 L = 1, NQ
            DO 440 J = 2, JDIM
               DO 430 K = 1, KDIM
                  QPI = 0.25E0 * ( Q(J,K,I,L) + Q(J,K+1,I+1,L)
     1                           + Q(J,K+1,I,L) + Q(J,K,I+1,L) )
C
                  DQDX(K,J,L) = DQDX(K,J,L)
     1                        + VOLAUX(K,J) * ( QPI*STEMP(K,J,1) )
                  DQDY(K,J,L) = DQDY(K,J,L)
     1                        + VOLAUX(K,J) * ( QPI*STEMP(K,J,2) )
                  DQDZ(K,J,L) = DQDZ(K,J,L)
     1                        + VOLAUX(K,J) * ( QPI*STEMP(K,J,3) )
  430          CONTINUE
  440       CONTINUE
  450    CONTINUE
C
C     Add the contribution from the I- term
C
C     Average metric terms in I- direction
C
         DO 520 J = 2, JDIM
            DO 510 K = 1, KDIM
               STEMP(K,J,1) = 0.5E0*( S(J,K,I-1,1,1)*S(J,K,I-1,4,1) 
     1                              + S(J,K+1,I-1,1,1)*S(J,K+1,I-1,4,1))
               STEMP(K,J,2) = 0.5E0*( S(J,K,I-1,2,1)*S(J,K,I-1,4,1)
     1                              + S(J,K+1,I-1,2,1)*S(J,K+1,I-1,4,1))
               STEMP(K,J,3) = 0.5E0*( S(J,K,I-1,3,1)*S(J,K,I-1,4,1)
     1                              + S(J,K+1,I-1,3,1)*S(J,K+1,I-1,4,1))
  510       CONTINUE
  520    CONTINUE
C
C     Add the derivative data from the I- term
C
         DO 550 L = 1, NQ
            DO 540 J = 2, JDIM
               DO 530 K = 1, KDIM
                  QMI = 0.25E0 * ( Q(J,K,I,L) + Q(J,K+1,I-1,L)
     1                           + Q(J,K+1,I,L) + Q(J,K,I-1,L) )
C
                  DQDX(K,J,L) = DQDX(K,J,L)
     1                        - VOLAUX(K,J) * ( QMI*STEMP(K,J,1) )
                  DQDY(K,J,L) = DQDY(K,J,L)
     1                        - VOLAUX(K,J) * ( QMI*STEMP(K,J,2) )
                  DQDZ(K,J,L) = DQDZ(K,J,L)
     1                        - VOLAUX(K,J) * ( QMI*STEMP(K,J,3) )
  530          CONTINUE
  540       CONTINUE
  550    CONTINUE
      ENDIF
C
C     Calculate Q at cell interfaces as averages of cell centers
C
      DO 630 L = 1, NQ
         DO 620 J = 2, JDIM
            DO 610 K = 1, KDIM
               QI(K,J,L) = 0.5E0 * ( Q(J,K,I,L)  + Q(J,K+1,I,L) )
  610       CONTINUE
  620    CONTINUE
  630 CONTINUE
C
C     Calculate PROPSI at cell interfaces as averages of cell centers
C
      DO 660 L = 1, NP
         DO 650 J = 2, JDIM
            DO 640 K = 1, KDIM
               PROPSI(K,J,L) = 0.5E0*(PROPS(J,K,I,L)+PROPS(J,K+1,I,L))
  640       CONTINUE
  650    CONTINUE
  660 CONTINUE
C
C     Loop over boundary conditions and fix the wall data
C
      DO 750 IBC = 1, NBCS
         IBCTYP = IBCDAT(1,IBC)
         IDIR   = IBCDAT(2,IBC)
         IF ( (IBCTYP .EQ. IWALL .OR. IBCTYP .EQ. IWALFN)
     1        .AND. IDIR .EQ. 3) THEN
            INROUT = IBCDAT(3,IBC)
            ISTRT  = IBCDAT(4,IBC)
            IEND   = IBCDAT(5,IBC)
            JSTRT  = IBCDAT(6,IBC)
            JEND   = IBCDAT(7,IBC)
            IF (INROUT .EQ. 1) THEN
               KBCPT = 1
               KNTR  = 1
            ELSE
               KBCPT = KDIM + 1
               KNTR  = KDIM
            ENDIF
            IF (I .GE. ISTRT .AND. I .LE. IEND) THEN
               DO 720 L = 1, NQ
                  DO 710 J = JSTRT, JEND
                     QI(KNTR,J,L) = Q(J,KBCPT,I,L)
  710             CONTINUE
  720          CONTINUE
C
C     Also do properties
C
               DO 740 L = 1, NP
                  DO 730 J = JSTRT, JEND
                     PROPSI(KNTR,J,L) = PROPS(J,KBCPT,I,L)
  730             CONTINUE
  740          CONTINUE
            ENDIF
         ENDIF
  750 CONTINUE
C
C     Finished calculating derivatives for K interfaces
C
      RETURN
      END








      SUBROUTINE FVIS (IMODEL, NPTS, QI, PROPSI, DQDX, DQDY, DQDZ, 
     1                 S, TAU, F)
C
C     Routine to calculate the Navier-Stokes viscous fluxes.  This routine 
C     is for laminar flow or mean equations in Favre averaged form coupled 
C     with a turbulence model.  Currently, algebraic, two-equation and
C     Reynold's stress models are supported by this routine.
C
C     The flow variables, properties and derivatives of the flow
C     variables are provided to this routine from which it calculates
C     the viscous fluxes.
C
C     This routine is organized as follows:
C          1. Calculate the Reynolds Stresses at the Interfaces
C             (this uses a model or the actual Reynolds stresses)
C          2. Calculate the Viscous Fluxes for the Mean Equations
C          3. Calculate Additional Viscous Fluxes
C             a. Laminar, algebraic model => No additional fluxes
C             b. Two-Equation model       => Fluxes for K and T2 equations
C             c. Reynolds Stress model    => Fluxes for Tau_ij and T2 eqns,
C                                            Molecular and Turbulent diffusion
C
C     IMODEL         : Turbulence model array
C     NPTS           : Number of interfaces to calculate viscous fluxes
C     QI             : Primitive variables at interfaces
C     PROPSI         : Properties at the interfaces
C                      PROPSI(1) = RMU
C                      PROPSI(2) = RMUT
C                      PROPSI(3) = Y+
C     DQDX           : X-Derivative of primitive variables at the interface
C     DQDY           : Y-Derivative of primitive variables at the interface
C     DQDZ           : Z-Derivative of primitive variables at the interface
C     S              : Metrics
C     TAU            : Reynold's stresses
C     F              : Viscous fluxes at cell interfaces
C
      include '../header/common.h'
C
      DIMENSION IMODEL(NMDL)
C
      DIMENSION QI    (NPTS,NQ),
     1          PROPSI(NPTS,NP),
     2          DQDX  (NPTS,NQ),
     3          DQDY  (NPTS,NQ),
     4          DQDZ  (NPTS,NQ),
     5          S     (NPTS,3),
     6          TAU   (NPTS,6),
     7          F     (NPTS,NF)
C
C     Two-Thirds
C
      T23 = 2.E0 / 3.E0
C
      ITURB = IMODEL(1)
C
C     1. Calculate the Reynold's Stresses at the Interfaces.  
C                                      ~
C     The array TAU contains {- rho (u" u")}.
C                                     i  j
C     For laminar flow:            TAU_ij => zero.
C         algebraic models:        TAU_ij => modeled with Boussinesq approx.
C         two-equation models:     TAU_ij => modeled with Boussinesq approx.
C         Reynold's stress models: TAU_ij => computed from PDE.
C     Pass in PROPS and DUDX,DUDY,DUDZ for Boussinesq approximation.
C     Pass in QI for average interface data.  Pass in derivative data DQDX,...
C     The Reynold's stresses are stored as:
C         TAU(1) = TAUXX, TAU(2) = TAUYY, TAU(3) = TAUZZ,
C         TAU(4) = TAUXY, TAU(5) = TAUXZ, TAU(6) = TAUYZ
C
      CALL TAUIJ (IMODEL, NPTS, QI, PROPSI, DQDX, DQDY, DQDZ, TAU)
C
C     2. Calculate the Viscous Fluxes for the Mean Flow Equations
C
      DO 100 I = 1, NPTS
C
         RHO   = QI(I,1)
         U     = QI(I,2)
         V     = QI(I,3)
         W     = QI(I,4)
         P     = QI(I,5)
C
         T     = GAMMA * P / RHO
C
         RMU   = PROPSI(I,1)
         RMUT  = PROPSI(I,2)
C
C     Calculate Strain Rates and Shear Stress Terms
C     Use molecular viscosity (RMU) averaged from cell centers to 
C     get the viscosity at the cell interface.
C     T23 is two-thirds.
C
         SXX   = T23 * ( 2.E0*DQDX(I,2) - DQDY(I,3) - DQDZ(I,4) )
         SYY   = T23 * ( 2.E0*DQDY(I,3) - DQDX(I,2) - DQDZ(I,4) )
         SZZ   = T23 * ( 2.E0*DQDZ(I,4) - DQDX(I,2) - DQDY(I,3) )
         SXY   = ( DQDY(I,2) + DQDX(I,3) )
         SXZ   = ( DQDX(I,4) + DQDZ(I,2) )
         SYZ   = ( DQDY(I,4) + DQDZ(I,3) )
C
         VISC  = RMU * FSMACH / RE
         SIGXX = VISC * SXX + TAU(I,1) 
         SIGYY = VISC * SYY + TAU(I,2)
         SIGZZ = VISC * SZZ + TAU(I,3)
         SIGXY = VISC * SXY + TAU(I,4)
         SIGXZ = VISC * SXZ + TAU(I,5)
         SIGYZ = VISC * SYZ + TAU(I,6)
C
C     Calculate Heat Flux Terms.
C     Turbulent heat flux uses the gradient transport assumption here.
C     DTDx_i can be calculated from derivatives of primitive variables as
C          DTDx_i = T [DPDx_i / P - DrhoDx_i / rho]
C
         DTDX  = T * (DQDX(I,5) / P - DQDX(I,1) / RHO)
         DTDY  = T * (DQDY(I,5) / P - DQDY(I,1) / RHO)
         DTDZ  = T * (DQDZ(I,5) / P - DQDZ(I,1) / RHO)
C
         CNDCT =   (RMU/PR + RMUT/PRT) * FSMACH / (RE * GAMM1)
         QDOTX = - CNDCT * DTDX
         QDOTY = - CNDCT * DTDY
         QDOTZ = - CNDCT * DTDZ
C
c--->    CSTAR = SIGRHO * RMUT / (RHO * RHO) * FSMACH / RE
c turn this term off as it is inconsistent to have it here and not in the
c tke or tau_ij equations.  this causes problems for the code - so zero it out.
         CSTAR = 0.e0
C
C     Calculate F(viscous)
C
         F(I,1)  = 0.E0
         F(I,2)  = S(I,1) * SIGXX + S(I,2) * SIGXY + S(I,3) * SIGXZ
         F(I,3)  = S(I,1) * SIGXY + S(I,2) * SIGYY + S(I,3) * SIGYZ
         F(I,4)  = S(I,1) * SIGXZ + S(I,2) * SIGYZ + S(I,3) * SIGZZ
         F(I,5)  = S(I,1) * ( U*SIGXX + V*SIGXY + W*SIGXZ - QDOTX 
     1              +CSTAR*(SXX*DQDX(I,1)+SXY*DQDY(I,1)+SXZ*DQDZ(I,1)) )
     2           + S(I,2) * ( U*SIGXY + V*SIGYY + W*SIGYZ - QDOTY 
     3              +CSTAR*(SXY*DQDX(I,1)+SYY*DQDY(I,1)+SYZ*DQDZ(I,1)) )
     4           + S(I,3) * ( U*SIGXZ + V*SIGYZ + W*SIGZZ - QDOTZ 
     5              +CSTAR*(SXZ*DQDX(I,1)+SYZ*DQDY(I,1)+SZZ*DQDZ(I,1)) )
  100 CONTINUE
C
C     3. Calculate Additional Viscous Fluxes
C
C     3-b. Calculate the viscous fluxes for two-equation turbulence models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 200 I = 1, NPTS
C
            RHO   = QI(I,1)
            P     = QI(I,5)
C
            T     = GAMMA * P / RHO
C
            RMU   = PROPSI(I,1)
            RMUT  = PROPSI(I,2)
C
C     DTDx_i can be calculated from derivatives of primitive variables as
C          DTDx_i = T [DPDx_i / P - DrhoDx_i / rho]
C
            DTDX  = T * (DQDX(I,5) / P - DQDX(I,1) / RHO)
            DTDY  = T * (DQDY(I,5) / P - DQDY(I,1) / RHO)
            DTDZ  = T * (DQDZ(I,5) / P - DQDZ(I,1) / RHO)
C
C     Calculate Pressure Diffusion Correlation terms
C
            VISCPD = RMUT * P / RHO * FSMACH / RE
c           VISCPD = (RMU + RMUT) * P / RHO * FSMACH / RE
            PDIFX  = VISCPD * (CPDFRO*DQDX(I,1)/RHO + CPDFT*DTDX/T)
            PDIFY  = VISCPD * (CPDFRO*DQDY(I,1)/RHO + CPDFT*DTDY/T)
            PDIFZ  = VISCPD * (CPDFRO*DQDZ(I,1)/RHO + CPDFT*DTDZ/T)
C
C     Coefficients for DKDX and DT2DX diffusion terms
C
            VISCK  = (RMU + SIGK  * RMUT) * FSMACH / RE
            VISCT2 = (RMU + SIGT2 * RMUT) * FSMACH / RE
            F(I,5) = F(I,5) + S(I,1) * ( VISCK*DQDX(I,6) )
     1                      + S(I,2) * ( VISCK*DQDY(I,6) )
     2                      + S(I,3) * ( VISCK*DQDZ(I,6) )
            F(I,6) = S(I,1) * ( VISCK*DQDX(I,6) + PDIFX )
     1             + S(I,2) * ( VISCK*DQDY(I,6) + PDIFY )
     2             + S(I,3) * ( VISCK*DQDZ(I,6) + PDIFZ )
            F(I,7) = S(I,1) * ( VISCT2*DQDX(I,7) )
     1             + S(I,2) * ( VISCT2*DQDY(I,7) )
     2             + S(I,3) * ( VISCT2*DQDZ(I,7) )
  200    CONTINUE
C
C     3-c. Calculate viscous fluxes for Reynolds stress models
C          i.   Zero out the fluxes F(6-12)
C          ii.  Call routines to calculate 
C               1. Molecular diffusion 
C               2. Diffusion model for Epsilon_ij_wall from So et al. 1995
C               3. Turbulent diffusion (triple correlation)
C          iii. Calculate additional fluxes including the contraction (trace) 
C               of the molecular and turbulent diffusions for Reynolds 
C               stresses in the energy equation
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         IPSMDL = IMODEL(6)
         IEPSIJ = IMODEL(7)
         IDMMDL = IMODEL(8)
         IDTMDL = IMODEL(9)
C
C     3-c.i   Zero fluxes F(6-12)
C
         DO 300 I = 1, NPTS
            F(I, 6) = 0.E0
            F(I, 7) = 0.E0
            F(I, 8) = 0.E0
            F(I, 9) = 0.E0
            F(I,10) = 0.E0
            F(I,11) = 0.E0
            F(I,12) = 0.E0
  300    CONTINUE
C
C     3-c.ii.1  Add contribution due to Molecular diffusion model
C
         CALL RSMDIF (IDMMDL, NPTS, QI, PROPSI, DQDX, DQDY, DQDZ, S, 
     1                TAU, F)
C
C     3-c.ii.2  Add contribution for Epsilon_ij_wall for So et al. 1995 model
C
         IF (IEPSIJ .EQ. IESO95) THEN
            CALL RSDFEW (IEPSIJ, NPTS, QI, PROPSI, DQDX, DQDY, DQDZ, S, 
     1                   TAU, F)
         ENDIF
C
C     3-c.ii.3  Add contribution due to Turbulent diffusion model
C
         CALL RSTDIF (IDTMDL, NPTS, QI, PROPSI, DQDX, DQDY, DQDZ, S, 
     1                TAU, F)
C
C     3-c.iii Add additional fluxes include contraction of diffusion
C             fluxes to Energy equation
C
         DO 350 I = 1, NPTS
            RHO    = QI(I,1)
            P      = QI(I,5)
C
            T      = GAMMA * P / RHO
C
            TKE    = 0.5E0 * (QI(I,6) + QI(I,7) + QI(I,8))
            EPSLN  = QI(I,12)
C
            RMU    = PROPSI(I,1)
            RMUT   = PROPSI(I,2)
C
C     DTDx_i can be calculated from derivatives of primitive variables as
C          DTDx_i = T [DPDx_i / P - DrhoDx_i / rho]
C
            DTDX  = T * (DQDX(I,5) / P - DQDX(I,1) / RHO)
            DTDY  = T * (DQDY(I,5) / P - DQDY(I,1) / RHO)
            DTDZ  = T * (DQDZ(I,5) / P - DQDZ(I,1) / RHO)
C
C     Calculate Pressure Diffusion Correlation terms
C
            VISCPD = RMUT * P / RHO * FSMACH / RE
c           VISCPD = (RMU + RMUT) * P / RHO * FSMACH / RE
            PDIFX  = VISCPD * (CPDFRO*DQDX(I,1)/RHO + CPDFT*DTDX/T)
            PDIFY  = VISCPD * (CPDFRO*DQDY(I,1)/RHO + CPDFT*DTDY/T)
            PDIFZ  = VISCPD * (CPDFRO*DQDZ(I,1)/RHO + CPDFT*DTDZ/T)
C
C     Coefficients for diffusion terms in epsilon equation
C
            VISCT2 = RMU * FSMACH / RE
            CRKE   = CEPSLN * RHO * TKE / EPSLN
C
C     Calculate F(viscous)
C
           TRACE   = 0.5E0 * (F(I, 6) + F(I, 7) + F(I, 8))
           F(I, 5) = F(I, 5) + TRACE
           F(I, 6) = F(I, 6) + S(I,1)*PDIFX+S(I,2)*PDIFY+S(I,3)*PDIFZ
           F(I, 7) = F(I, 7) + S(I,1)*PDIFX+S(I,2)*PDIFY+S(I,3)*PDIFZ
           F(I, 8) = F(I, 8) + S(I,1)*PDIFX+S(I,2)*PDIFY+S(I,3)*PDIFZ
           F(I, 9) = F(I, 9) + S(I,1)*PDIFX+S(I,2)*PDIFY+S(I,3)*PDIFZ
           F(I,10) = F(I,10) + S(I,1)*PDIFX+S(I,2)*PDIFY+S(I,3)*PDIFZ
           F(I,11) = F(I,11) + S(I,1)*PDIFX+S(I,2)*PDIFY+S(I,3)*PDIFZ
           F(I,12) = F(I,12) 
     1             + S(I,1)*( VISCT2*DQDX(I,12)+
     2                        CRKE  *(QI(I, 6)*DQDX(I,12)+
     3                                QI(I, 9)*DQDY(I,12)+
     4                                QI(I,10)*DQDZ(I,12)) )
     5             + S(I,2)*( VISCT2*DQDY(I,12)+
     6                        CRKE  *(QI(I, 9)*DQDX(I,12)+
     7                                QI(I, 7)*DQDY(I,12)+
     8                                QI(I,11)*DQDZ(I,12)) )
     9             + S(I,3)*( VISCT2*DQDZ(I,12)+
     A                        CRKE  *(QI(I,10)*DQDX(I,12)+
     B                                QI(I,11)*DQDY(I,12)+
     C                                QI(I, 8)*DQDZ(I,12)) )
  350    CONTINUE
      ENDIF
C
C     Finished with Viscous Terms
C
      RETURN
      END






      SUBROUTINE TAUIJ (IMODEL, NPTS, Q, PROPS, DQDX, DQDY, DQDZ, TAU)
C
C     Routine to calculate the Reynold's stress terms.
C                       ~
C     TAU   = {- rho (u"  u" )} => be careful of the density term.
C        ij            i   j
C     For laminar flow:            TAU_ij => zero.
C         algebraic models:        TAU_ij => modeled with Boussinesq approx.
C         two-equation models:     TAU_ij => modeled with Boussinesq approx.
C         Reynold's stress models: TAU_ij => computed from PDE.
C     Pass in PROPS and DQDX,DQDY,DQDZ for Boussinesq approximation.  
C     Pass in Q for primtive variable data at the same location as 
C     DQDX and PROPS.
C         DQDn(2) = DuDn; DQDn(3) = DvDn; DQDn(4) = DwDn, where n = (X,Y,Z).
C     The Reynold's stresses are stored as:
C         TAU(1) = TAUXX, TAU(2) = TAUYY, TAU(3) = TAUZZ,
C         TAU(4) = TAUXY, TAU(5) = TAUXZ, TAU(6) = TAUYZ
C
C     This routine is structured as follows:
C          1. Laminar Flow          => Tau_ij = 0
C          2. Algebraic Model       => Boussinesq assumption
C          3. Two-Equation Model
C             a. Boussinesq assumption
C             b. Algebraic Stress model
C          4. Reynolds Stress Model => Solved directly
C          5. Large-Eddy Simulation => Boussinesq assumption
C
C     IMODEL         : Turbulence model array
C     NPTS           : Number of points
C     Q              : Primitive variables
C     PROPS          : Properties at same location as Q and DQ
C                      PROPS(1) = RMU
C                      PROPS(2) = RMUT
C                      PROPS(3) = Y+
C     DQDX           : X-Derivative of primitive variables
C     DQDY           : Y-Derivative of primitive variables
C     DQDZ           : Z-Derivative of primitive variables
C     TAU            : Reynold's stresses
C
      include '../header/common.h'
C
      DIMENSION IMODEL(NMDL)
C
      DIMENSION Q    (NPTS,NQ),
     1          PROPS(NPTS,NP),
     2          DQDX (NPTS,NQ),
     3          DQDY (NPTS,NQ),
     4          DQDZ (NPTS,NQ),
     5          TAU  (NPTS,6)
C
C     One-Third
C     Two-Thirds
C
      T13 = 1.E0 / 3.E0
      T23 = 2.E0 / 3.E0
C
      ITURB = IMODEL( 1)
      IASM  = IMODEL(10)
C
C     Loop over the points calculating Reynold's stresses
C
C     1. Laminar Flow.
C        TAU = 0.
C
      IF (ITURB .EQ. ITLMNR) THEN
         DO 100 I = 1, NPTS
            TAU(I,1) = 0.E0
            TAU(I,2) = 0.E0
            TAU(I,3) = 0.E0
            TAU(I,4) = 0.E0
            TAU(I,5) = 0.E0
            TAU(I,6) = 0.E0
  100    CONTINUE
C
C     2. Algebraic Turbulence Model.
C        Boussinesq assumption (no rho*k term).
C
      ELSE IF (ITURB .EQ. ITBLMX) THEN
         DO 200 I = 1, NPTS
C
            RMUT     = PROPS(I,2)
C
            SXX      = T23 * ( 2.E0*DQDX(I,2) - DQDY(I,3) - DQDZ(I,4) )
            SYY      = T23 * ( 2.E0*DQDY(I,3) - DQDX(I,2) - DQDZ(I,4) )
            SZZ      = T23 * ( 2.E0*DQDZ(I,4) - DQDX(I,2) - DQDY(I,3) )
            SXY      = ( DQDY(I,2) + DQDX(I,3) )
            SXZ      = ( DQDX(I,4) + DQDZ(I,2) )
            SYZ      = ( DQDY(I,4) + DQDZ(I,3) )
C
            VISC     = RMUT * FSMACH / RE
            TAU(I,1) = VISC * SXX
            TAU(I,2) = VISC * SYY
            TAU(I,3) = VISC * SZZ
            TAU(I,4) = VISC * SXY
            TAU(I,5) = VISC * SXZ
            TAU(I,6) = VISC * SYZ
  200    CONTINUE
C
C     3. Two-Equation Turbulence Model.
C
C     3-a. Boussinesq assumption.
C     3-b. Gatski's algebraic Reynolds stress model
C
      ELSE IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 300 I = 1, NPTS
C
            RHOK     = Q(I,1) * Q(I,6)
            RMUT     = PROPS(I,2)
C
            SXX      = T23 * ( 2.E0*DQDX(I,2) - DQDY(I,3) - DQDZ(I,4) )
            SYY      = T23 * ( 2.E0*DQDY(I,3) - DQDX(I,2) - DQDZ(I,4) )
            SZZ      = T23 * ( 2.E0*DQDZ(I,4) - DQDX(I,2) - DQDY(I,3) )
            SXY      = ( DQDY(I,2) + DQDX(I,3) )
            SXZ      = ( DQDX(I,4) + DQDZ(I,2) )
            SYZ      = ( DQDY(I,4) + DQDZ(I,3) )
C
            VISC     = RMUT * FSMACH / RE
            TAU(I,1) = VISC * SXX - T23 * RHOK
            TAU(I,2) = VISC * SYY - T23 * RHOK
            TAU(I,3) = VISC * SZZ - T23 * RHOK
            TAU(I,4) = VISC * SXY
            TAU(I,5) = VISC * SXZ
            TAU(I,6) = VISC * SYZ
  300    CONTINUE
C
         IF (IASM .EQ. IASMGS) THEN
            DO 320 LOOP = 1, 1
               DO 310 I = 1, NPTS
C
                  RHO      = Q(I,1)
                  TKE      = Q(I,6)
                  RHOK     = RHO * TKE
C
C     Constants from pressure-strain model
C---- Begin SSG constants -------------------------------------------------
                  CC1      = 6.80E0
                  CC2      = 0.36E0
                  CC3      = 1.25E0
                  CC4      = 0.40E0
                  CC5      = 1.88E0
C---- End SSG constants ---------------------------------------------------
C
                  IF (ITURB .EQ. ITKW) THEN
                     TIMSCL   = 1.E0 / (BSTRKW * Q(I,7) + RSMALL)
                  ELSE
                     TIMSCL   = TKE  / (Q(I,7) + RSMALL)
                  ENDIF
C
                  GSMALL   = 1.E0 / (.5E0*CC1 + CC5 - 1.E0)
C
                  RALF1    = 0.50E0 * (4.E0/3.E0-CC2) * GSMALL
                  RALF2    = 0.25E0 * (2.E0 - CC3)**2 * GSMALL * GSMALL
                  RALF3    = 0.25E0 * (2.E0 - CC4)**2 * GSMALL * GSMALL
                  RALF4    = 0.50E0 * (2.E0 - CC4)    * GSMALL
                  RALF5    =          (2.E0 - CC3)    * GSMALL
C
C     S_ij
C           - note: sxx,etc. terms appearing elsewhere in this routine are
C                  NOT S_ij but a leftover from the original viscos routines
C
                  SXX      = DQDX(I,2)
                  SYY      = DQDY(I,3)
                  SZZ      = DQDZ(I,4)
                  SXY      = 0.5E0 * (DQDY(I,2) + DQDX(I,3))
                  SXZ      = 0.5E0 * (DQDZ(I,2) + DQDX(I,4))
                  SYZ      = 0.5E0 * (DQDZ(I,3) + DQDY(I,4))
C
C     S_kk = (S_xx + S_yy + S_zz) / 3
C
                  SKK      = (SXX + SYY + SZZ) / 3.E0
C
C     S_kl * S_lk / 3
C
                  SKLSLK   = ( SXX*SXX + SYY*SYY + SZZ*SZZ
     1                 + 2.E0*(SXY*SXY + SXZ*SXZ + SYZ*SYZ) )/3.E0
C
C     W_ij : for a non-rotating frame
C
                  WYX      = 0.5E0 * (DQDX(I,3) - DQDY(I,2))
                  WZX      = 0.5E0 * (DQDX(I,4) - DQDZ(I,2))
                  WZY      = 0.5E0 * (DQDY(I,4) - DQDZ(I,3))
                  WXY      = 0.5E0 * (DQDY(I,2) - DQDX(I,3))
                  WXZ      = 0.5E0 * (DQDZ(I,2) - DQDX(I,4))
                  WYZ      = 0.5E0 * (DQDZ(I,3) - DQDY(I,4))
C
                  ETASQ    = RALF2 * TIMSCL * TIMSCL *
     1                           ( SXX*SXX + SYY*SYY + SZZ*SZZ
     2                     + 2.E0*(SXY*SXY + SXZ*SXZ + SYZ*SYZ) )
                  ZETSQ    = RALF3 * TIMSCL * TIMSCL *
     1                       2.E0*(WXY*WXY + WXZ*WXZ + WYZ*WYZ) 
C
                  ETA6     = ETASQ * ETASQ * ETASQ
                  ZET6     = ZETSQ * ZETSQ * ZETSQ
C
                  CMUS     = RALF1 * ( 3.E0 * (1.E0 + ETASQ) 
     1                               + 0.2E0 * (ETA6 + ZET6) ) /
     2                  ( 3.E0 + ETASQ + 6.E0*ETASQ*ZETSQ + 6.E0*ZETSQ 
     3                  + ETA6 + ZET6 )
                  CMUSS    = RALF1 * ( 3.E0 * (1.E0 + ETASQ) ) /
     1                  ( 3.E0 + ETASQ + 6.E0*ETASQ*ZETSQ + 6.E0*ZETSQ 
     2                  + ETA6 + ZET6 )
C
                  RMUTS    = RHO * CMUS  * TKE * TIMSCL
                  RMUTSS   = RHO * CMUSS * TKE * TIMSCL
C
C Store C_mu^* (CMUS) to pass into calculation of Mu_T
C
                  PROPS(I,ICMUST) = CMUS
C
                  TAU(I,1) = - T23 * RHOK + 2.E0 * RMUTS * (SXX-SKK)
     1                       + 2.E0 * RMUTSS * TIMSCL * 
     2                   ( RALF4*2.E0*(        SXY*WYX+SXZ*WZX) 
     3                   - RALF5*     (SXX*SXX+SXY*SXY+SXZ*SXZ-SKLSLK) )
                  TAU(I,2) = - T23 * RHOK + 2.E0 * RMUTS * (SYY-SKK)
     1                       + 2.E0 * RMUTSS * TIMSCL * 
     2                   ( RALF4*2.E0*(SXY*WXY        +SYZ*WZY) 
     3                   - RALF5*     (SXY*SXY+SYY*SYY+SYZ*SYZ-SKLSLK) )
                  TAU(I,3) = - T23 * RHOK + 2.E0 * RMUTS * (SZZ-SKK)
     1                       + 2.E0 * RMUTSS * TIMSCL * 
     2                   ( RALF4*2.E0*(SXZ*WXZ+SYZ*WYZ        ) 
     3                   - RALF5*     (SXZ*SXZ+SYZ*SYZ+SZZ*SZZ-SKLSLK) )
                  TAU(I,4) =                2.E0 * RMUTS * (SXY)
     1                       + 2.E0 * RMUTSS * TIMSCL * 
     2                   ( RALF4*     (SXX*WXY          +SXZ*WZY
     3                                          +SYY*WYX+SYZ*WZX) 
     4                   - RALF5*     (SXX*SXY+SXY*SYY+SXZ*SYZ) )
                  TAU(I,5) =                2.E0 * RMUTS * (SXZ)
     1                       + 2.E0 * RMUTSS * TIMSCL * 
     2                   ( RALF4*     (SXX*WXZ+SXY*WYZ          
     3                                        +SYZ*WYX+SZZ*WZX) 
     4                   - RALF5*     (SXX*SXZ+SXY*SYZ+SXZ*SZZ) )
                  TAU(I,6) =                2.E0 * RMUTS * (SYZ)
     1                       + 2.E0 * RMUTSS * TIMSCL * 
     2                   ( RALF4*     (SXY*WXZ+SYY*WYZ          
     3                                +SXZ*WXY        +SZZ*WZY) 
     4                   - RALF5*     (SXY*SXZ+SYY*SYZ+SYZ*SZZ) )
C
C     Enforce positivity of normal stresses 
C     N.B. TAU = - rho tau_ij
C
                  TAU(I,1) = MIN (TAU(I,1), 0.E0)
                  TAU(I,2) = MIN (TAU(I,2), 0.E0)
                  TAU(I,3) = MIN (TAU(I,3), 0.E0)
  310          CONTINUE
  320       CONTINUE
         ENDIF
C
C     4. Reynold's Stress Model.
C        Solved directly.
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 500 I = 1, NPTS
C
            RHO      = Q(I,1)
C
            TAU(I,1) = - RHO * Q(I, 6)
            TAU(I,2) = - RHO * Q(I, 7)
            TAU(I,3) = - RHO * Q(I, 8)
            TAU(I,4) = - RHO * Q(I, 9)
            TAU(I,5) = - RHO * Q(I,10)
            TAU(I,6) = - RHO * Q(I,11)
  500    CONTINUE
C
C     5. Large-Eddy Simulation
C        Boussinesq assumption (no rho*k term).
C
      ELSE IF (ITURB .EQ. ITLES) THEN
         DO 600 I = 1, NPTS
C
            RMUT     = PROPS(I,2)
C
            SXX      = T23 * ( 2.E0*DQDX(I,2) - DQDY(I,3) - DQDZ(I,4) )
            SYY      = T23 * ( 2.E0*DQDY(I,3) - DQDX(I,2) - DQDZ(I,4) )
            SZZ      = T23 * ( 2.E0*DQDZ(I,4) - DQDX(I,2) - DQDY(I,3) )
            SXY      = ( DQDY(I,2) + DQDX(I,3) )
            SXZ      = ( DQDX(I,4) + DQDZ(I,2) )
            SYZ      = ( DQDY(I,4) + DQDZ(I,3) )
C
            VISC     = RMUT * FSMACH / RE
            TAU(I,1) = VISC * SXX
            TAU(I,2) = VISC * SYY
            TAU(I,3) = VISC * SZZ
            TAU(I,4) = VISC * SXY
            TAU(I,5) = VISC * SXZ
            TAU(I,6) = VISC * SYZ
  600    CONTINUE
      ELSE
         WRITE (IOUT,1000) ITURB
         STOP
      ENDIF
C
C     Finished with Reynold's stresses
C
      RETURN
 1000 FORMAT (' ','Invalid turbulence model specified in TAUIJ => ',I5,
     1       /' ','Run Aborting.')
      END







      SUBROUTINE NEWTPN (NORDR, P, ROOT)
C
C     Routine to perform a Newton iteration to solve for the root of
C     a polynomial of degree N.
C
C     NORDR : Degree (order) of polynomial to solve (>= 1).
C             This is set < 0 if iterations do not converge.
C     P     : Coefficients of polynomial to solve
C             F(x) = P(0) + P(1)*x + P(2)*x^2 + P(3)*X^3 + ...
C     ROOT  : The solution for the polynomial
C             ROOT = guess of root on input
C             ROOT = calculated root upon return
C
      include '../header/common.h'
C
      DIMENSION P     (0:NORDR)
C
C     Check to verify that NORDR > 0
C
      IF (NORDR .LE. 0) THEN
         WRITE (IOUT,1000) NORDR
         STOP
      ENDIF
C
C     Set maximum number of iterations and tolerance for convergence
C
      MAXIT = 50
c-->  TOLER = 1.0E-8
      TOLER = 1.0E-6
C
C     Perform MAXIT iterations of Newton iteration to calculate root
C
      DO 900 ITER = 1, MAXIT
         F  = 0.E0
         FP = 0.E0
         DO 100 I = NORDR, 1, -1
            F  = (F + P(I)) * ROOT
            FP = FLOAT(I) * P(I) + FP * ROOT
  100    CONTINUE
         F = F + P(0)
C
C     Update the root guess, prevent divide by zero
C
         IF (FP .EQ. 0.E0) THEN
            DROOT = - 0.1E0 * F
         ELSE
            DROOT = - F / FP
         ENDIF
         ROOT  = ROOT + DROOT
C
c        write (*,'(a,i5,4e15.7)') ' iter,f,fp,dr,root ',
c    1                               iter,f,fp,droot,root
C
C     Check for convergence of ROOT.
C
c-->     IF (ABS(DROOT) .LE. TOLER) GO TO 990
         IF (ABS(DROOT/ROOT) .LE. TOLER) GO TO 990
  900 CONTINUE
C
C     Newton iteration did not converge.  Set NORDR negative and return.
C
      NORDR = -1
  990 RETURN
 1000 FORMAT (' ','RUN ABORTING=> Invalid polynomial in NEWTPN ',I5)
      END
