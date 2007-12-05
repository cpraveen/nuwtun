      SUBROUTINE METBC (IDIM, JDIM, KDIM, S, VOL, NBCS, IBCDAT)
C
C     Routine to do boundary conditions on the metric terms.  This
C     routine sets up metric terms outside of the domain that is 
C     necessary for the inviscid terms so that the viscous terms
C     may more easily be evaluated.
C
C     Presently set up only for IWALL, WALFN, ITAN  and IFAR for THIN LAYER
C
C
      include '../header/common.h'
C
      DIMENSION S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     1          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IJKDIM(3)
C
C     Store IDIM, JDIM, KDIM in IJKDIM for looping
C
      IJKDIM(1) = IDIM
      IJKDIM(2) = JDIM
      IJKDIM(3) = KDIM
C
      DO 900 IBC = 1, NBCS
C
C     Store boundary condition data in mnemonic names
C
         IBCTYP = IBCDAT(1,IBC)
         IDIR   = IBCDAT(2,IBC)
         INROUT = IBCDAT(3,IBC)
         ISTRT  = IBCDAT(4,IBC)
         IEND   = IBCDAT(5,IBC)
         JSTRT  = IBCDAT(6,IBC)
         JEND   = IBCDAT(7,IBC)
         KSTRT  = IBCDAT(8,IBC)
         KEND   = IBCDAT(9,IBC)
C
C   a. FAR, FIX, EXT, JET, PROFILE, SUBSONIC PROFILE BCs
C
         IF (IBCTYP .EQ. IFAR   .OR. IBCTYP .EQ. IEXT .OR.
     1       IBCTYP .EQ. IFIX   .OR. IBCTYP .EQ. IJET .OR.
     2       IBCTYP .EQ. IPROFL .OR. IBCTYP .EQ. IPROSB) THEN
C
C     Set needed variables depending on whether the boundary is
C     the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C          NMET  is the location of the metrics on the boundary
C          NINRM is the location of the metrics interior to the boundary
C          NVOL  is the location of the volume on the boundary
C          NINRV is the location of the volume interior to the boundary
C
            IF (INROUT .EQ. 1) THEN
               NMET   = 0
               NINRM  = 1
               NVOL   = 1
               NINRV  = 2
            ELSE
               NMET   = IJKDIM(IDIR) + 1
               NINRM  = IJKDIM(IDIR) 
               NVOL   = IJKDIM(IDIR) + 1
               NINRV  = IJKDIM(IDIR) 
            ENDIF
C
C     I direction
C
            IF (IDIR .EQ. 1) THEN
               DO 120 K = KSTRT, KEND
                  DO 110 J = JSTRT, JEND
                     VOL(J,K,NVOL) = VOL(J,K,NINRV)
  110             CONTINUE
  120          CONTINUE
C
C     J Direction
C
            ELSE IF (IDIR .EQ. 2) THEN
               DO 140 I = ISTRT, IEND
                  DO 130 K = KSTRT, KEND
                     VOL(NVOL,K,I) = VOL(NINRV,K,I)
  130             CONTINUE
  140          CONTINUE
C
C     K Direction
C
            ELSE IF (IDIR .EQ. 3) THEN
               DO 160 I = ISTRT, IEND
                  DO 150 J = JSTRT, JEND
                     VOL(J,NVOL,I) = VOL(J,NINRV,I)
  150             CONTINUE
  160          CONTINUE
            ENDIF
C
C   b. PERIODIC BCs and Periodic HOMOGENEOUS SHEAR BC
C
         ELSE IF (IBCTYP .EQ. IPERD .OR. IBCTYP .EQ. IHSHR) THEN
C
C     Set needed variables depending on whether the boundary is
C     the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C          NMET  is the location of the metrics on the boundary
C          NPRDM is the location of the periodic metrics 
C          NVOL  is the location of the volume on the boundary
C          NPRDV is the location of the periodic volume 
C     N.B. Metrics normal to the BC direction must be set using periodic
C          condition with ISTRT-1,JSTRT-1,KSTRT-1 to account for the 
C          numbering on a cell.
C
            IF (INROUT .EQ. 1) THEN
               NMET   = 0
               NPRDM  = IJKDIM(IDIR) 
               NVOL   = 1
               NPRDV  = IJKDIM(IDIR) 
            ELSE
               NMET   = IJKDIM(IDIR) + 1
               NPRDM  = 1
               NVOL   = IJKDIM(IDIR) + 1
               NPRDV  = 2
            ENDIF
C
C     I direction
C
            IF (IDIR .EQ. 1) THEN
               DO 210 K = KSTRT, KEND
                  DO 205 J = JSTRT, JEND
                     VOL(J,K,NVOL)     = VOL(J,K,NPRDV)
C
                     S  (J,K,NMET,1,1) = S  (J,K,NPRDM,1,1)
                     S  (J,K,NMET,2,1) = S  (J,K,NPRDM,2,1)
                     S  (J,K,NMET,3,1) = S  (J,K,NPRDM,3,1)
                     S  (J,K,NMET,4,1) = S  (J,K,NPRDM,4,1)
  205             CONTINUE
  210          CONTINUE
C
               DO 220 K = KSTRT-1, KEND
                  DO 215 J = JSTRT-1, JEND
                     S  (J,K,NVOL,1,2) = S  (J,K,NPRDV,1,2)
                     S  (J,K,NVOL,2,2) = S  (J,K,NPRDV,2,2)
                     S  (J,K,NVOL,3,2) = S  (J,K,NPRDV,3,2)
                     S  (J,K,NVOL,4,2) = S  (J,K,NPRDV,4,2)
                     S  (J,K,NVOL,1,3) = S  (J,K,NPRDV,1,3)
                     S  (J,K,NVOL,2,3) = S  (J,K,NPRDV,2,3)
                     S  (J,K,NVOL,3,3) = S  (J,K,NPRDV,3,3)
                     S  (J,K,NVOL,4,3) = S  (J,K,NPRDV,4,3)
  215             CONTINUE
  220          CONTINUE
C
C     J Direction
C
            ELSE IF (IDIR .EQ. 2) THEN
               DO 230 I = ISTRT, IEND
                  DO 225 K = KSTRT, KEND
                     VOL(NVOL,K,I)     = VOL(NPRDV,K,I)
C
                     S  (NMET,K,I,1,2) = S  (NPRDM,K,I,1,2)
                     S  (NMET,K,I,2,2) = S  (NPRDM,K,I,2,2)
                     S  (NMET,K,I,3,2) = S  (NPRDM,K,I,3,2)
                     S  (NMET,K,I,4,2) = S  (NPRDM,K,I,4,2)
  225             CONTINUE
  230          CONTINUE
C
               DO 240 I = ISTRT-1, IEND
                  DO 235 K = KSTRT-1, KEND
                     S  (NVOL,K,I,1,1) = S  (NPRDV,K,I,1,1)
                     S  (NVOL,K,I,2,1) = S  (NPRDV,K,I,2,1)
                     S  (NVOL,K,I,3,1) = S  (NPRDV,K,I,3,1)
                     S  (NVOL,K,I,4,1) = S  (NPRDV,K,I,4,1)
                     S  (NVOL,K,I,1,3) = S  (NPRDV,K,I,1,3)
                     S  (NVOL,K,I,2,3) = S  (NPRDV,K,I,2,3)
                     S  (NVOL,K,I,3,3) = S  (NPRDV,K,I,3,3)
                     S  (NVOL,K,I,4,3) = S  (NPRDV,K,I,4,3)
  235             CONTINUE
  240          CONTINUE
C
C     K Direction
C
            ELSE IF (IDIR .EQ. 3) THEN
               DO 250 I = ISTRT, IEND
                  DO 245 J = JSTRT, JEND
                     VOL(J,NVOL,I)     = VOL(J,NPRDV,I)
C
                     S  (J,NMET,I,1,3) = S  (J,NPRDM,I,1,3)
                     S  (J,NMET,I,2,3) = S  (J,NPRDM,I,2,3)
                     S  (J,NMET,I,3,3) = S  (J,NPRDM,I,3,3)
                     S  (J,NMET,I,4,3) = S  (J,NPRDM,I,4,3)
  245             CONTINUE
  250          CONTINUE
C
               DO 260 I = ISTRT, IEND
                  DO 255 J = JSTRT, JEND
                     S  (J,NVOL,I,1,1) = S  (J,NPRDV,I,1,1)
                     S  (J,NVOL,I,2,1) = S  (J,NPRDV,I,2,1)
                     S  (J,NVOL,I,3,1) = S  (J,NPRDV,I,3,1)
                     S  (J,NVOL,I,4,1) = S  (J,NPRDV,I,4,1)
                     S  (J,NVOL,I,1,2) = S  (J,NPRDV,I,1,2)
                     S  (J,NVOL,I,2,2) = S  (J,NPRDV,I,2,2)
                     S  (J,NVOL,I,3,2) = S  (J,NPRDV,I,3,2)
                     S  (J,NVOL,I,4,2) = S  (J,NPRDV,I,4,2)
  255             CONTINUE
  260          CONTINUE
            ENDIF
         ENDIF
  900 CONTINUE
C
C     Finished with metric boundary condition routine
C
      RETURN
 1000 FORMAT (' ','*** WARNING ***',
     1       /' ','Unsupported boundary condition type in METBC ',I6,
     2       /' ','Substituting extrapolation of interior metrics and ',
     3            'volumes and continuing.')
      END
