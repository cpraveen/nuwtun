      SUBROUTINE LIMITR_Q (NPTS, DQP, DQPB, DQM, DQMB, LIMIT, RKAP, 
     1                     LBEG, LEND, IERROR)
C
C The second order limited correction to the first order interpolation
C is calculated by this routine and stored in DQP for the q(+) term and 
C DQM for the q(-) term.  Delta(q(i)) = q(i+1) - q(i) is stored in DQP(i) 
C on entry.
C
C N.B. The values of DQP(NPTS) and DQM(1) are not defined properly on exit 
C as this leads to better vectorization and they are unnecessary.
C
C This routine is organized as follows:
C     0. Set constants for the interpolation
C     1. Venkatakrishnan's limiter
C     2. MinMod limiter
C     3. Smooth limiter
C     4. Unlimited scheme
C     5. Error
C
C NPTS           : Length of array to limit successive values of delta(q)
C DQP,DQM        : Storage for delta q = q(i) - q(i-1) on entry.
C                  Storage of the higher order correction to the 
C                  first order scheme on exit
C LIMIT          : Limiter to use for MUSCL interpolation
C                  ILNONE [0] => Unlimited Scheme
C                  ILSMTH [1] => Smooth Limiter
C                  ILMNMD [2] => Minmod Limiter
C                  ILVNKT [3] => Smooth limiter (fixed at kappa = 1/3) 
C                           from Venkatakrishnan, AIAA J. Vol. 29 No. 7
C                           also Thomas, J. L., 5th Copper Mtn. Multigrid
C RKAP           : Kappa to use in MUSCL interpolation
C LBEG,LEND      : Starting,stopping equation index for limiting
C IERROR         : Error code
C
      include '../../src-flo/header/common.h'
C
      DIMENSION DQP   (NPTS,NQ),
     1          DQM   (NPTS,NQ),
     2          DQPB  (NPTS,NQ),
     3          DQMB  (NPTS,NQ)
C
C     Initialize error
C
      IERROR = INO
C
C 0.  Set constants for the interpolation.
C
      OMKAP = 1.E0 - RKAP
      OPKAP = 1.E0 + RKAP
C
C 1.  Venkat's Smooth Limiter
C
      IF (LIMIT .EQ. ILVNKT) THEN
         DO L = LEND, LBEG, -1
            DO I = NPTS, 2, -1
               CALL SMILIM_BQ(DQP(I-1,L), DQPB(I-1,L),
     1                        DQP(I,L),   DQPB(I,L),
     2                        DQM(I,L),   DQMB(I,L))
            END DO
         END DO
C
C 2.  MinMod Limiter
C
      ELSE IF (LIMIT .EQ. ILMNMD) THEN
         B = (3.E0 - RKAP) / (1.E0 - RKAP)
         DO L = LEND, LBEG, -1
            DO I = NPTS, 2, -1
               CALL RMMLIM_BQ(DQP(I-1,L), DQPB(I-1,L),
     1                        DQP(I,L),   DQPB(I,L),
     2                        DQM(I,L),   DQMB(I,L),
     3                        OMKAP, OPKAP, B)
            END DO
         END DO
C
C 3.  Smooth Limiter
C
      ELSE IF (LIMIT .EQ. ILSMTH) THEN
         DO L = LEND, LBEG, -1
            DO I = NPTS, 2, -1
               CALL SMTLIM_BQ(DQP(I-1,L), DQPB(I-1,L), 
     1                        DQP(I,L),   DQPB(I,L), 
     2                        DQM(I,L),   DQMB(I,L), RKAP)
            END DO
         END DO
C
C 4.  Unlimited Scheme
C
      ELSE IF (LIMIT .EQ. ILNONE) THEN
         DO L = LEND, LBEG, -1
            DO I = NPTS, 2, -1
               CALL RNOLIM_BQ(DQP(I-1,L), DQPB(I-1,L), 
     1                        DQP(I,L),   DQPB(I,L), 
     2                        DQM(I,L),   DQMB(I,L), OMKAP, OPKAP)
            END DO
         END DO
C
C 5.  Error
C
      ELSE
         WRITE (IOUT,1000) LIMIT
         IERROR = IYES
         STOP
      ENDIF
C
C Finished with limiting routine.
C
      RETURN
 1000 FORMAT(' ','LIMITR_Q  : Invalid value of LIMIT = ',I10,
     1      /' ','Run Aborting.')
      END
