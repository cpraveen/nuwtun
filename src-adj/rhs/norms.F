      SUBROUTINE L2MAX1 (IBLK, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   RES, VOL, RESMAX, IRESMX)
C
C Routine to calculate the maximum residual (of mean equations).  Actually,
C calculate the square of the residual and take the square root before output.
C RES as passed in is the sum of the fluxes and source term.  Need to 
C divide this by the volume for the residual output.
C
C IBLK           : Current block
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at 
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C RES            : Residual
C VOL            : Cell volume
C R2NORM         : L2 norm of the residual of each equation
C                  and of all the equations
C RESMAX         : Maximum residual
C IRESMX         : Location of maximum residual
C
      include '../../src-flo/header/common.h'
C
      DIMENSION RES   (0:JDIM+2,0:KDIM+2,0:IDIM+2,NF),
     1          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION IRESMX(0:3)
C
C 1.  Loop to find the maximum residual.
C     Do not initialize the maximum residual to 0.  This is done in the main
C     routine so that the maximum is over the entire flowfield (i.e.,
C     including all zones) or by lines for marching.
C
      DO 130 I = 2, IDIM
         DO 120 K = 2, KDIM
            DO 110 J = 2, JDIM
               R2  = ( RES(J,K,I,1) * RES(J,K,I,1) 
     1               + RES(J,K,I,2) * RES(J,K,I,2) 
     2               + RES(J,K,I,3) * RES(J,K,I,3) 
     3               + RES(J,K,I,4) * RES(J,K,I,4) 
     4               + RES(J,K,I,5) * RES(J,K,I,5) ) / 
     5               ( VOL(J,K,I)   * VOL(J,K,I) )
               IF (R2 .GT. RESMAX) THEN
                  RESMAX    = R2
                  IRESMX(0) = IBLK
                  IRESMX(1) = I
                  IRESMX(2) = J
                  IRESMX(3) = K
               ENDIF
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C Finished with maximum residual on this block
C
      RETURN
      END





      SUBROUTINE L2NRM1 (IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   RES, VOL, R2NORM)
C
C     Routine to calculate the various L2 Norms of the Residual.
C     RES as passed in is the sum of the fluxes and source term.  Need
C     to divide this by the volume for the L2NORM of residual output.
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     IPSTRT         : I Plane (cell center) to start calculations at 
C                      (Minimum Value: 2)
C     IPEND          : I Plane (cell center) to end caclulations at
C                      (Maximum Value: IDIM)
C     RES            : Residual
C     VOL            : Cell volume
C     R2NORM         : L2 norm of the residual of each equation
C                      and of all the equations
C
      include '../../src-flo/header/common.h'
C
      DIMENSION RES   (0:JDIM+2,0:KDIM+2,0:IDIM+2,NF),
     1          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION R2NORM(NF+1)
C
C     Loop to form the L2 Norms
C     Do not initialize the norms to 0.  This is done in the main
C     routine so that the norms are over the entire flowfield (i.e.,
C     including all zones) or by lines for marching.
C
      DO 400 L = 1, NF
         DO 300 I = 2, IDIM
            DO 200 K = 2, KDIM
               DO 100 J = 2, JDIM
                  R2NORM(L) = R2NORM(L) + RES(J,K,I,L)*RES(J,K,I,L) /
     1                                     (VOL(J,K,I)*VOL(J,K,I))
  100          CONTINUE
  200       CONTINUE
  300    CONTINUE
  400 CONTINUE
C
C     Form the L2 Norm of the entire residual
C
      R2NORM(NF+1) = 0.E0
      DO 500 L = 1, NF
         R2NORM(NF+1) = R2NORM(NF+1) + R2NORM(L)
  500 CONTINUE
C
C     Finished with L2 Norm
C
      RETURN
      END
