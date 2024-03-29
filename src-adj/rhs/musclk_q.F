      SUBROUTINE MUSCLK_Q (IDIM, JDIM, KDIM, I,
     1                   Q, QB, DQP, DQPB, DQM, DQMB, QP, QPB, QM, QMB,
     2                   IBCDAT, NBCS, IORDER, LIMIT, RKAP)
C
C Routine to do MUSCL interpolation in the K direction. This routine 
C interpolates for Q(+,-) on the K faces for k = 1, kdim by j = 2, jdim 
C on the current I plane. This allows for greater vectorization and also 
C minimizes intermediate storage to a cross section.  
C
C This routine is organized as follows:
C     0. Set up beginning/ending indices for first and second order sections
C     1. Calculate Q(+,-)
C        a. First order
C        b. Second order
C           i.   Calculate delta(Q)
C           ii.  Limit delta(Q)
C           iii. Construct Q(+,-)
C     2. Correct wall values of Q(+,-)
C
C IDIM,JDIM,KDIM : Dimensions of current block
C I              : I Plane to do interpolation on
C Q              : Primitive variables at cell centers
C DQP,DQM        : Storage for delta q = q(k) - q(k-1) and
C                  for the higher order interpolant to be added
C                  to the first order scheme to create QP,QM
C QP,QM          : Q(+),Q(-) at a cell interface for MUSCL scheme
C IBCDAT         : Data controlling boundary conditions:
C                  IBCDAT(1,ibc) = Boundary Condition Type
C                         2      = Direction (1=i, 2=j, 3=k)
C                         3      = k(1 or kdim)
C                         4      = istrt
C                         5      = iend
C                         6      = jstrt
C                         7      = jend
C                         8      = kstrt
C                         9      = kend
C                         10     = profile identifier
C NBCS           : Number of boundary conditions for current block
C IORDER         : Order of MUSCL interpolation for I direction 
C                  IORDER = 1 first order
C                         = 2 second order (or 'third' order)
C LIMIT          : Limiter to use for MUSCL interpolation
C RKAP           : Kappa to use in MUSCL interpolation
C
      include '../../src-flo/header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION DQP   (0:KDIM+1,2:JDIM,NQ),
     1          DQM   (0:KDIM+1,2:JDIM,NQ),
     2          QP    (KDIM,2:JDIM,NQ), 
     3          QM    (KDIM,2:JDIM,NQ) 
C
      DIMENSION QB    (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          DQPB  (0:KDIM+1,2:JDIM,NQ),
     2          DQMB  (0:KDIM+1,2:JDIM,NQ),
     3          QPB   (KDIM,2:JDIM,NQ), 
     4          QMB   (KDIM,2:JDIM,NQ) 
C
      DIMENSION IBCDAT(10,NBCS)
C
C Local Storage
C
      LOGICAL   DOFRST,
     1          DOSCND
C
C 0.  Set up beginning/ending indices for first and second order sections
C
      LFOBEG = 0
      LFOEND = -1
      LSOBEG = 1
      LSOEND = NQ
      DOFRST = .FALSE.
      DOSCND = .TRUE.
C
      IF (IORDER .EQ. 1) THEN
         LFOBEG = 1
         LFOEND = NQ
         LSOBEG = 0
         LSOEND = -1
         DOFRST = .TRUE.
         DOSCND = .FALSE.
      ELSE IF (FOTURB) THEN
         LFOBEG = 6
         LFOEND = NQ
         LSOBEG = 1
         LSOEND = 5
         DOFRST = .TRUE.
         DOSCND = .TRUE.
      END IF
C
C 1.  Calculate Q(+,-)
C     Loop over j,k formulating q(extrapolated). K is stored first so that 
C     successive values of k (k = 1,2,3...) are stored together for limiting.
C
C 1-a. First Order
C
C
C 1-b. Second Order or Third Order
C 1-b-i. Calculate delta(Q)
C        Calculate delta(q(k)) = q(k+1) - q(k) to formulate q(extrapolated) 
C        using limiting.
C
      IF (DOSCND) THEN
         DO 230 L = LSOBEG, LSOEND
            DO 220 J = 2, JDIM
               DO 200 K = 0, KDIM + 1
                  DQP(K,J,L) = Q(J,K+1,I,L) - Q(J,K,I,L)
  200          CONTINUE
  220       CONTINUE
  230    CONTINUE
C
C     Set delta(q) to 2*delta(q) to adjust for the wall boundary condition 
C     being on the wall rather than at a cell center.
C     Set delta(q) = 0 on wall for wall function to adjust for grid spacing.
C
         DO 300 IBC = 1, NBCS
            IBCTYP = IBCDAT(1,IBC)
            IDIR   = IBCDAT(2,IBC)
            IF ( IDIR .EQ. 3 .AND.
     1           ( IBCTYP .EQ. ITAN  .OR. IBCTYP .EQ. IWALL .OR.
     2             IBCTYP .EQ. IWALFN ) ) THEN
               INROUT = IBCDAT(3,IBC)
               ISTRT  = IBCDAT(4,IBC)
               IEND   = IBCDAT(5,IBC)
               JSTRT  = IBCDAT(6,IBC)
               JEND   = IBCDAT(7,IBC)
               IF (INROUT .EQ. 1) THEN
                  KBCPT = 1
               ELSE
                  KBCPT = KDIM
               ENDIF
               IF (IBCTYP .EQ. IWALFN) THEN
                  FACTR = 0.E0
               ELSE
                  FACTR  = 2.E0
               ENDIF
               IF ( I .GE. ISTRT .AND. I .LE. IEND ) THEN
                  DO 290 L = LSOBEG, LSOEND
                     DO 280 J = JSTRT, JEND
                        DQP(KBCPT,J,L) = FACTR * DQP(KBCPT,J,L)
  280                CONTINUE
  290             CONTINUE
               ENDIF
            ENDIF
  300    CONTINUE
C
C 1-b-ii.  Limit delta(Q)
C     Limit from k = 0 to k = kdim+1 by j = 2 to j = jdim times.
C     DQP is used to pass in the values of delta q.
C     The values of the limited deltas are stored in DQP and DQM.
C     These values not valid at the endpoints of each k line, but are
C     calculated here to maximize vectorization on a kxj plane.
C
C
C 1-b-iii. Construct Q(+,-)
C     Calculate q(extrapolated,k+1/2) from q and the limited values of delta(q)
C
      ENDIF
C
C 2.  Correct wall values of Q(+,-)
C     Fix boundary values of q(extrapolated) due to boundary conditions
C     for either the first order or higher order scheme.
C     WALL boundary condition.  Reset QP(wall)=QM(wall)=Q(wall)
C
C BEGIN BACKWARD SWEEP
C
C 2.  Correct wall values of Q(+,-)
C     Fix boundary values of q(extrapolated) due to boundary conditions
C     for either the first order or higher order scheme.
C     WALL boundary condition.  Reset QP(wall)=QM(wall)=Q(wall)
C
      DO 1800 IBC = 1, NBCS
         IBCTYP = IBCDAT(1,IBC)
         IDIR   = IBCDAT(2,IBC)
         IF ( IDIR .EQ. 3 .AND.
     1        ( IBCTYP .EQ. ITAN  .OR. IBCTYP .EQ. IWALL .OR.
     2          IBCTYP .EQ. IWALFN ) ) THEN
            INROUT = IBCDAT(3,IBC)
            ISTRT  = IBCDAT(4,IBC)
            IEND   = IBCDAT(5,IBC)
            JSTRT  = IBCDAT(6,IBC)
            JEND   = IBCDAT(7,IBC)
            IF (INROUT .EQ. 1) THEN
               KQPM  = 1
               KBCPT = 1
            ELSE
               KQPM  = KDIM
               KBCPT = KDIM + 1
            ENDIF
            IF ( I .GE. ISTRT .AND. I .LE. IEND ) THEN
               DO 1790 L = 1, NQ
                  DO 1780 J = JSTRT, JEND
                     QB(J,KBCPT,I,L) = QB(J,KBCPT,I,L) + QPB(KQPM,J,L)
     1                                                 + QMB(KQPM,J,L)
                     QPB(KQPM,J,L)   = 0.0E0
                     QMB(KQPM,J,L)   = 0.0E0
 1780             CONTINUE
 1790          CONTINUE
            ENDIF
         ENDIF
 1800 CONTINUE
C
C 1-b. Second Order or Third Order
C
      IF (DOSCND) THEN
C
C 1-b-iii. Construct Q(+,-)
C     Calculate q(extrapolated,k+1/2) from q and the limited values of delta(q)
C
         DO 1380 L = LSOBEG, LSOEND
            DO 1370 J = 2, JDIM
               DO 1360 K = 1, KDIM
                  QB(J,K,I,L) = QB(J,K,I,L) + QMB(K,J,L)
                  DQMB(K,J,L) = DQMB(K,J,L) + QMB(K,J,L)
                  QMB(K,J,L)  = 0.0E0
                  QB(J,K+1,I,L) = QB(J,K+1,I,L) + QPB(K,J,L)
                  DQPB(K,J,L)   = DQPB(K,J,L)   + QPB(K,J,L)
                  QPB(K,J,L)    = 0.0E0
 1360          CONTINUE
 1370       CONTINUE
 1380    CONTINUE
C
C 1-b-ii.  Limit delta(Q)
C     Limit from k = 0 to k = kdim+1 by j = 2 to j = jdim times.
C     DQP is used to pass in the values of delta q.
C     The values of the limited deltas are stored in DQP and DQM.
C     These values not valid at the endpoints of each k line, but are
C     calculated here to maximize vectorization on a kxj plane.
C
         NPTS = ( KDIM + 2 ) * ( JDIM - 1 )
         CALL LIMITR_Q (NPTS, DQP, DQPB, DQM, DQMB, LIMIT, RKAP, 
     1                  LSOBEG, LSOEND, IERROR)
C
C     Set delta(q) to 2*delta(q) to adjust for the wall boundary condition 
C     being on the wall rather than at a cell center.
C     Set delta(q) = 0 on wall for wall function to adjust for grid spacing.
C
         DO 1300 IBC = 1, NBCS
            IBCTYP = IBCDAT(1,IBC)
            IDIR   = IBCDAT(2,IBC)
            IF ( IDIR .EQ. 3 .AND.
     1           ( IBCTYP .EQ. ITAN  .OR. IBCTYP .EQ. IWALL .OR.
     2             IBCTYP .EQ. IWALFN ) ) THEN
               INROUT = IBCDAT(3,IBC)
               ISTRT  = IBCDAT(4,IBC)
               IEND   = IBCDAT(5,IBC)
               JSTRT  = IBCDAT(6,IBC)
               JEND   = IBCDAT(7,IBC)
               IF (INROUT .EQ. 1) THEN
                  KBCPT = 1
               ELSE
                  KBCPT = KDIM
               ENDIF
               IF (IBCTYP .EQ. IWALFN) THEN
                  FACTR = 0.E0
               ELSE
                  FACTR  = 2.E0
               ENDIF
               IF ( I .GE. ISTRT .AND. I .LE. IEND ) THEN
                  DO 1290 L = LSOBEG, LSOEND
                     DO 1280 J = JSTRT, JEND
                        DQPB(KBCPT,J,L) = FACTR * DQPB(KBCPT,J,L)
 1280                CONTINUE
 1290             CONTINUE
               ENDIF
            ENDIF
 1300    CONTINUE
C
C 1-b-i. Calculate delta(Q)
C        Calculate delta(q(k)) = q(k+1) - q(k) to formulate q(extrapolated) 
C        using limiting.
C
         DO 1230 L = LSOBEG, LSOEND
            DO 1220 J = 2, JDIM
               DO 1200 K = 0, KDIM + 1
                  QB(J,K+1,I,L) = QB(J,K+1,I,L) + DQPB(K,J,L)
                  QB(J,K  ,I,L) = QB(J,K  ,I,L) - DQPB(K,J,L)
                  DQPB(K,J,L)   = 0.0E0
 1200          CONTINUE
 1220       CONTINUE
 1230    CONTINUE
      ENDIF
C
C 1-a. First Order
C
      IF (DOFRST) THEN
         DO 1130 L = LFOBEG, LFOEND
            DO 1120 J = 2, JDIM
               DO 1110 K = 1, KDIM
                  QB(J,K,I,L)   = QB(J,K,I,L) + QMB(K,J,L)
                  QMB(K,J,L)    = 0.0E0
                  QB(J,K+1,I,L) = QB(J,K+1,I,L) + QPB(K,J,L)
                  QPB(K,J,L)    = 0.0E0
 1110          CONTINUE
 1120       CONTINUE
 1130    CONTINUE
      ENDIF
C
C Finished with Q(+/-) calculation for MUSCL scheme in k.
C
      RETURN
      END
