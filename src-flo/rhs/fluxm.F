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
      SUBROUTINE FLUXM (IDIM, JDIM, KDIM, IPSTRT, IPEND, MXSECT,
     1                  Q, PROPS, NPPTS, S, VOL, RES, 
     2                  DQP, DQM, QP, QM, FP, FM, STEMP,
     3                  IBCDAT, NBCS, 
     4                  IORDER, LIMIT, VGNRN, SIGMA, IMODEL)
C
C     Routine to calculate the contribution to the residual from the
C     i direction for marching.  Uses fully upwind either first or
C     second order accurate in i.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     IPSTRT         : I Plane (cell center) to start calculations at 
C                      (Minimum Value: 2)
C     IPEND          : I Plane (cell center) to end caclulations at
C                      (Maximum Value: IDIM)
C     MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                      (used in dimensioning temporary space used for 
C                       calculating q(+,-) and f.
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
C     RES            : Residual
C     DQP,DQM        : Used as temporary storage for FLUX presently
C     QP,QM          : Q at a cell interface for MUSCL scheme
C                      QP = Q(i+1/2), QM = Q(i-1/2)
C     FP,FM          : Fluxes at cell interfaces
C                      FP = F(i+1/2), FM = F(i-1/2)
C     STEMP          : Storage for metrics passed to flux evaluation routine.
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     IORDER         : Order of MUSCL interpolation for I direction 
C                      IORDER = 1 first order
C                             = 2 second order (or 'third' order)
C     LIMIT          : Limiter to use for MUSCL interpolation
C                      If LIMIT > 0 then do catastrophic limiting in the
C                      I direction for marching
C     VGNRN          : Logical for Vigneron approach in marching flux 
C                      evaluation
C     SIGMA          : Safety factor for Vigneron technique when marching
C     IMODEL         : Turbulence model array
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          DQP   (2:JDIM,2:KDIM,NQ), 
     2          DQM   (2:JDIM,2:KDIM,NQ), 
     3          QP    (2:JDIM,2:KDIM,NQ), 
     4          QM    (2:JDIM,2:KDIM,NQ), 
     5          FP    (2:JDIM,2:KDIM,NF), 
     6          FM    (2:JDIM,2:KDIM,NF), 
     7          STEMP (2:JDIM,2:KDIM,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IMODEL(NMDL)
C
      LOGICAL   VGNRN
C
      ITURB = IMODEL(1)
C
C     ICAT counts the number of times that the catastrophic limiter is on
C
      ICAT = 0
C
C     Calculate everything at the cross section I = IPSTRT = IPEND
C
      IF (IPSTRT .NE. IPEND) THEN
         WRITE (IOUT,1000) IPSTRT, IPEND
         STOP
      ENDIF
C
      I = IPSTRT
C
C     Inviscid Residuals
C
      IF (IORDER .EQ. 1) THEN
         DO 60 L = 1, NQ
            DO 50 K = 2, KDIM
               DO 40 J = 2, JDIM
                  QP(J,K,L) = Q(J,K,I,L) 
                  QM(J,K,L) = Q(J,K,I-1,L) 
   40          CONTINUE
   50       CONTINUE
   60    CONTINUE
      ELSE
         DO 90 L = 1, NQ
            DO 80 K = 2, KDIM
               DO 70 J = 2, JDIM
                  QP(J,K,L) = 1.5E0*Q(J,K,I,L)   - 0.5E0*Q(J,K,I-1,L)
                  QM(J,K,L) = 1.5E0*Q(J,K,I-1,L) - 0.5E0*Q(J,K,I-2,L)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
C
C     Catastrophic limiting for higher order - check density and pressure
C
         IF ( LIMIT .GT. 0 ) THEN
            EPS  = 0.E0
            ICAT = 0
            DO 95 K = 2, KDIM
               DO 94 J = 2, JDIM
                  IF (QP(J,K,1) .LE. EPS .OR. QP(J,K,5) .LE. EPS) THEN
                     ICAT = ICAT + 1
                     DO 91 L = 1, NQ
                        QP(J,K,L) = Q(J,K,I,L)
   91                CONTINUE
                  ENDIF
                  IF (QM(J,K,1) .LE. EPS .OR. QM(J,K,5) .LE. EPS) THEN
                     ICAT = ICAT + 1
                     DO 92 L = 1, NQ
                        QM(J,K,L) = Q(J,K,I-1,L)
   92                CONTINUE
                  ENDIF
   94          CONTINUE
   95       CONTINUE
            IF (ICAT .GT. 0) THEN
               WRITE (IOUT,1200) ICAT
            ENDIF
         ENDIF
      ENDIF
C
C     Calculate the fluxes F(i+1/2) from q(plus). Pass only
C     the appropriate metrics from the current direction.
C
      DO 103 L = 1, 3
         DO 102 K = 2, KDIM
            DO 101 J = 2, JDIM
               STEMP(J,K,L) = S(J,K,I,L,1)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
C
C     Call Flux evaluation routine.  For temporary storage use:
C        DQP => TKE
C        DQM => H
C
      NPTS = (JDIM - 1) * (KDIM - 1)
      CALL FLUX (ITURB, NPTS, QP, STEMP, DQP, DQM, FP)
C
C     Calculate the fluxes F(i-1/2) from q(minus). Pass only
C     the appropriate metrics from the current direction.
C
      DO 113 L = 1, 3
         DO 112 K = 2, KDIM
            DO 111 J = 2, JDIM
               STEMP(J,K,L) = S(J,K,I-1,L,1)
  111       CONTINUE
  112    CONTINUE
  113 CONTINUE
C
C     Call Flux evaluation routine.  For temporary storage use:
C        DQP => TKE
C        DQM => H
C
      NPTS = (JDIM - 1) * (KDIM - 1)
      CALL FLUX (ITURB, NPTS, QM, STEMP, DQP, DQM, FM)
C
C     Add the flux contributions to the residual.
C
      DO 180 L = 1, NF
         DO 160 K = 2, KDIM
            DO 140 J = 2, JDIM
               RES(J,K,I,L) = RES(J,K,I,L) + 
     1                        ( FP(J,K,L) * S(J,K,I,4,1) -
     2                          FM(J,K,L) * S(J,K,I-1,4,1) )
  140       CONTINUE
  160    CONTINUE
  180 CONTINUE
C
C     Source term for the Vigneron technique
C
      IF (VGNRN) THEN
         DO 290 K = 2, KDIM
            DO 280 J = 2, JDIM
               UBAR  = 0.5E0*(S(J,K,I,1,1)+S(J,K,I-1,1,1))*Q(J,K,I,2)
     1               + 0.5E0*(S(J,K,I,2,1)+S(J,K,I-1,2,1))*Q(J,K,I,3)
     2               + 0.5E0*(S(J,K,I,3,1)+S(J,K,I-1,3,1))*Q(J,K,I,4)
               AMSQ  = UBAR*UBAR / (GAMMA * Q(J,K,I,5) / Q(J,K,I,1))
               OMEGA = MIN(1.E0, SIGMA*GAMMA*AMSQ/(1.E0+GAMM1*AMSQ) )
C
               PPI   = QP(J,K,5)
               PMI   = QM(J,K,5)
               PPJ   = 0.5E0*(Q(J,K,I,5) + Q(J+1,K,I,5))
               PMJ   = 0.5E0*(Q(J,K,I,5) + Q(J-1,K,I,5))
               PPK   = 0.5E0*(Q(J,K,I,5) + Q(J,K+1,I,5))
               PMK   = 0.5E0*(Q(J,K,I,5) + Q(J,K-1,I,5))
C
               UPI   = S(J,K,I,1,1)   * QP(J,K,2) +
     1                 S(J,K,I,2,1)   * QP(J,K,3) +
     2                 S(J,K,I,3,1)   * QP(J,K,4)
               UMI   = S(J,K,I-1,1,1) * QM(J,K,2) + 
     1                 S(J,K,I-1,2,1) * QM(J,K,3) +
     2                 S(J,K,I-1,3,1) * QM(J,K,4)
C
               RES(J,K,I,2) = RES(J,K,I,2) + (OMEGA - 1.E0) *
     1                        ( PPI*S(J,K,I,1,1)*S(J,K,I,4,1)
     2                        - PMI*S(J,K,I-1,1,1)*S(J,K,I-1,4,1)
     3                        + PPJ*S(J,K,I,1,2)*S(J,K,I,4,2)
     4                        - PMJ*S(J-1,K,I,1,2)*S(J-1,K,I,4,2)
     5                        + PPK*S(J,K,I,1,3)*S(J,K,I,4,3)
     6                        - PMK*S(J,K-1,I,1,3)*S(J,K-1,I,4,3) )
c              RES(J,K,I,3) = RES(J,K,I,3) + (OMEGA - 1.E0) *
c    1                        ( PPI*S(J,K,I,2,1)*S(J,K,I,4,1)
c    2                        - PMI*S(J,K,I-1,2,1)*S(J,K,I-1,4,1)
c    3                        + PPJ*S(J,K,I,2,2)*S(J,K,I,4,2)
c    4                        - PMJ*S(J-1,K,I,2,2)*S(J-1,K,I,4,2)
c    5                        + PPK*S(J,K,I,2,3)*S(J,K,I,4,3)
c    6                        - PMK*S(J,K-1,I,2,3)*S(J,K-1,I,4,3) )
c              RES(J,K,I,4) = RES(J,K,I,4) + (OMEGA - 1.E0) *
c    1                        ( PPI*S(J,K,I,3,1)*S(J,K,I,4,1)
c    2                        - PMI*S(J,K,I-1,3,1)*S(J,K,I-1,4,1)
c    3                        + PPJ*S(J,K,I,3,2)*S(J,K,I,4,2)
c    4                        - PMJ*S(J-1,K,I,3,2)*S(J-1,K,I,4,2)
c    5                        + PPK*S(J,K,I,3,3)*S(J,K,I,4,3)
c    6                        - PMK*S(J,K-1,I,3,3)*S(J,K-1,I,4,3) )
  280       CONTINUE
  290    CONTINUE
      ENDIF
C
C     Finished with i direction residual computation for marching.
C
      RETURN
C
 1000 FORMAT (' ','Routine FLUXM - invalid IPSTRT not equal IPEND ',2I6,
     1       /' ','RUN STOPPING')
 1200 FORMAT (' ',' Catastrophic limiting for I direction ',I8,' times')
C
      END
