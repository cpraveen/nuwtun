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

      SUBROUTINE COLQ2H (IDIMF, JDIMF, KDIMF, QF, VOLF,
     1                   IDIMC, JDIMC, KDIMC, QC, VOLC,
     2                   QS, ITURB, TKE)
C
C Routine to collect Q from h grid to 2h grid.  Collection is volume 
C weighted summation of conserved variables:
C      QC = sum(QF * VOLF) / VOLC
C
C This routine is organized as follows:
C     1. Convert QF from primitive to conserved variables
C     2. Loop over coarse I and collect in plane to minimize temporary storage
C        a. Sum in the I direction
C        b. If three-dimensional case, sum in K direction (not for 2D or axi)
C        c. Sum in the J direction and calculate QC = sum(QF * VOLF) / VOLC
C     3. Convert QF and from conserved to primitive variables
C
C IDIMF,JDIMF,KDIMF : Fine   grid dimensions
C QF                : Fine   grid variables at cell centers
C VOLF              : Fine   grid volumes
C IDIMC,JDIMC,KDIMC : Coarse grid dimensions
C QC                : Coarse grid variables at cell centers
C VOLC              : Coarse grid volumes
C QS                : Temporary space for forming summation(QF * VOLF)
C ITURB             : Turbulence model
C TKE               : Temporary storage of TKE to convert from primitive to
C                     conserved variables
C
      include '../header/common.h'
C
      INTEGER   IDIMF, JDIMF, KDIMF,
     1          IDIMC, JDIMC, KDIMC
C
      DIMENSION QF  (0:JDIMF+2,0:KDIMF+2,0:IDIMF+2,NQ),
     1          VOLF(0:JDIMF+1,0:KDIMF+1,0:IDIMF+1),
     2          QC  (0:JDIMC+2,0:KDIMC+2,0:IDIMC+2,NQ),
     3          VOLC(0:JDIMC+1,0:KDIMC+1,0:IDIMC+1)
C
      DIMENSION QS  (0:JDIMF+2,0:KDIMF+2,NQ),
     1          TKE (0:JDIMF+2,0:KDIMF+2,0:IDIMF+2)
C
C Local variables
C
      INTEGER   I , J , K , L,
     1          IC, JC, KC
C
C 1.   Convert QF from primitive to conserved variables
C
      NPTSF = (IDIMF + 3) * (JDIMF + 3) * (KDIMF + 3)
      CALL CONP2C (ITURB, NPTSF, QF, TKE)
      NPTSC = (IDIMC + 3) * (JDIMC + 3) * (KDIMC + 3)
      CALL CONP2C (ITURB, NPTSC, QC, TKE)
C
C 2.   Loop over coarse I and collect in plane to minimize temporary storage
C
      DO 290 IC = 2, IDIMC
         I = 2 * (IC - 1)
C
C 2-a. Sum in the I direction
C
         DO 215 L = 1, NQ
            DO 210 K = 2, KDIMF
               DO 205 J = 2, JDIMF
                  QS(J,K,L) = QF(J,K,I  ,L) * VOLF(J,K,I  ) +
     1                        QF(J,K,I+1,L) * VOLF(J,K,I+1)
  205          CONTINUE
  210       CONTINUE
  215    CONTINUE
C
C 2-b. If three-dimensional case, sum in K direction (not for 2D or axi)
C      The volume multiplication is handled in the I direction - now just sum.
C
         IF (THREED) THEN
            DO 240 L = 1, NQ
               DO 235 KC = 2, KDIMC
                  K = 2 * (KC - 1)
                  DO 230 J = 2, JDIMF
                     QS(J,KC,L) = QS(J,K,L) + QS(J,K+1,L)
  230             CONTINUE
  235          CONTINUE
  240       CONTINUE
         ENDIF
C
C 2-c. Sum in the J direction and calculate QC = sum(QF * VOLF) / sum(VOLF)
C
         DO 280 L = 1, NQ
            DO 275 KC = 2, KDIMC
               DO 270 JC = 2, JDIMC
                  J = 2 * (JC - 1)
                  QC(JC,KC,IC,L) = (QS(J,KC,L) + QS(J+1,KC,L))
     1                             / VOLC(JC,KC,IC)
  270          CONTINUE
  275       CONTINUE
  280    CONTINUE
C
  290 CONTINUE
C
C 3.   Convert QF and from conserved to primitive variables
C
      CALL CONC2P (ITURB, NPTSF, QF, TKE)
      NPTSC = (IDIMC + 3) * (JDIMC + 3) * (KDIMC + 3)
      CALL CONC2P (ITURB, NPTSC, QC, TKE)
C
C     Finished with collection of Q
C
      RETURN
      END




      SUBROUTINE COLR2H (IDIMF, JDIMF, KDIMF, RESF, 
     1                   IDIMC, JDIMC, KDIMC, RESCRS, RESS)
C
C Routine to collect the residual from h grid to 2h grid to form the coarse
C grid forcing function.  Coarse grid forcing function is:
C      RESCRS = sum(RESF) - RESC
C
C This routine is organized as follows:
C     Loop over coarse I and collect in a plane to minimize temporary storage
C     1. Sum in the I direction
C     2. If three-dimensional case, sum in K direction (not for 2D or axi)
C     3. Sum in the J direction and form RESCRS = sum(RESF) - RESC
C
C IDIMF,JDIMF,KDIMF : Fine   grid dimensions
C RESF              : Fine   grid residual
C IDIMC,JDIMC,KDIMC : Coarse grid dimensions
C RESCRS            : Input  : RES on coarse grid
C                     Returns: coarse grid forcing function
C RESS              : Temporary space for forming summation(RESF)
C
      include '../header/common.h'
C
      INTEGER   IDIMF, JDIMF, KDIMF,
     1          IDIMC, JDIMC, KDIMC
C
      DIMENSION RESF  (JDIMF+1,KDIMF+1,IDIMF+1,NF),
     1          RESCRS(JDIMC+1,KDIMC+1,IDIMC+1,NF)
C
      DIMENSION RESS  (JDIMF+1,KDIMF+1,NF)
C
C Local variables
C
      INTEGER   I , J , K , L,
     1          IC, JC, KC
C
C Loop over coarse I and collect in a plane to minimize temporary storage
C
      DO 900 IC = 2, IDIMC
         I = 2 * (IC - 1)
C
C 1.  Sum in the I direction
C
         DO 130 L = 1, NF
            DO 120 K = 2, KDIMF
               DO 110 J = 2, JDIMF
                  RESS(J,K,L) = RESF(J,K,I,L) + RESF(J,K,I+1,L)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
C
C 2.  If three-dimensional case, sum in K direction (not for 2D or axi)
C
         IF (THREED) THEN
            DO 230 L = 1, NF
               DO 230 KC = 2, KDIMC
                  K = 2 * (KC - 1)
                  DO 210 J = 2, JDIMF
                     RESS(J,KC,L) = RESS(J,K,L) + RESS(J,K+1,L)
  210             CONTINUE
  220          CONTINUE
  230       CONTINUE
         ENDIF
C
C 3.  Sum in the J direction and form RESCRS = sum(RESF) - RESC
C
         DO 330 L = 1, NF
            DO 320 KC = 2, KDIMC
               DO 310 JC = 2, JDIMC
                  J = 2 * (JC - 1)
                  RESCRS(JC,KC,IC,L) = RESS  (J  ,KC,L) 
     1                               + RESS  (J+1,KC,L)
c**  2                               - RESCRS(JC ,KC,IC,L)
  310          CONTINUE
  320       CONTINUE
  330    CONTINUE
C
  900 CONTINUE
C
C     Finished with calculating coarse grid forcing function
C
      RETURN
      END




      SUBROUTINE COLV2H (IDIMF, JDIMF, KDIMF, VOLF,
     1                   IDIMC, JDIMC, KDIMC, VOLC, VOLP)
C
C Routine to collect VOL from h grid to 2h grid. N.B. It is important that
C the coarse grid volumes are calculated as the summation of the fine grid 
C volumes so that the volume weighted collection generates a smooth and
C accurate coarse grid representation of the fine grid solution.
C
C This routine is organized as follows:
C     1. Loop over coarse I and collect in plane to minimize temporary storage
C        a. Sum in the I direction
C        b. If three-dimensional case, sum in K direction (not for 2D or axi)
C        c. Sum in the J direction and calculate VOLC = sum(VOLF)
C
C IDIMF,JDIMF,KDIMF : Fine   grid dimensions
C VOLF              : Fine   grid volumes
C IDIMC,JDIMC,KDIMC : Coarse grid dimensions
C VOLC              : Coarse grid volumes
C VOLP              : Temporary space for forming summation(VOLF) in a plane
C
      include '../header/common.h'
C
      INTEGER   IDIMF, JDIMF, KDIMF,
     1          IDIMC, JDIMC, KDIMC
C
      DIMENSION VOLF(0:JDIMF+1,0:KDIMF+1,0:IDIMF+1),
     1          VOLC(0:JDIMC+1,0:KDIMC+1,0:IDIMC+1)
C
      DIMENSION VOLP(0:JDIMF+2,0:KDIMF+2)
C
C Local variables
C
      INTEGER   I , J , K ,
     1          IC, JC, KC
C
C 1.   Loop over coarse I and collect in plane to minimize temporary storage
C
      DO 190 IC = 2, IDIMC
         I = 2 * (IC - 1)
C
C 1-a. Sum in the I direction
C
         DO 110 K = 2, KDIMF
            DO 105 J = 2, JDIMF
               VOLP(J,K) = VOLF(J,K,I) + VOLF(J,K,I+1)
  105       CONTINUE
  110    CONTINUE
C
C 1-b. If three-dimensional case, sum in K direction (not for 2D or axi)
C      The volume multiplication is handled in the I direction - now just sum.
C
         IF (THREED) THEN
            DO 135 KC = 2, KDIMC
               K = 2 * (KC - 1)
               DO 130 J = 2, JDIMF
                  VOLP(J,KC) = VOLP(J,K) + VOLP(J,K+1)
  130          CONTINUE
  135       CONTINUE
         ENDIF
C
C 1-c. Sum in the J direction and calculate VOLC = sum(VOLF)
C
         DO 175 KC = 2, KDIMC
            DO 170 JC = 2, JDIMC
               J = 2 * (JC - 1)
               VOLC(JC,KC,IC) = VOLP(J,KC) + VOLP(J+1,KC)
  170       CONTINUE
  175    CONTINUE
C
  190 CONTINUE
C
C     Finished calculating coarse grid volumes
C
      RETURN
      END





      SUBROUTINE CRSGRD (IDIMF, JDIMF, KDIMF, RF,
     1                   IDIMC, JDIMC, KDIMC, RC)
C
C Routine to form the coarser grid (IDIMC ... RC) by eliminating alternate 
C grid points in each direction of the fine grid (IDIMF ... RF)
C
C This routine is organized as follows:
C     1. Loop over the coarse grid, calculate fine grid location, set RC=RF
C
C IDIMF,JDIMF,KDIMF : Fine   grid dimensions
C RF                : Fine   grid postion vector
C IDIMC,JDIMC,KDIMC : Coarse grid dimensions
C RC                : Coarse grid postion vector
C
      include '../header/common.h'
C
      DIMENSION RF    (0:JDIMF+1,0:KDIMF+1,0:IDIMF+1,3),
     1          RC    (0:JDIMC+1,0:KDIMC+1,0:IDIMC+1,3)
C
C 1.  Loop over the coarse mesh and set each grid point to be every other 
C     grid point in the fine mesh.
C
      DO 130 IC = 1, IDIMC
         IF = 2 * IC - 1
         DO 120 KC = 1, KDIMC
            KF = 2 * KC - 1
            DO 110 JC = 1, JDIMC
               JF = 2 * JC - 1
               RC(JC,KC,IC,1) = RF(JF,KF,IF,1)
               RC(JC,KC,IC,2) = RF(JF,KF,IF,2)
               RC(JC,KC,IC,3) = RF(JF,KF,IF,3)
 110        CONTINUE
 120     CONTINUE
 130  CONTINUE
C
C Finished forming the coarser grid
C
      RETURN
      END




      SUBROUTINE CRSPRF (NPTSF, QBCF, IBCDTF,
     1                   NPTSC, QBCC, IBCDTC, IERRCD)
C
C Routine to calculate the coarse grid profile data from the fine grid data.
C
C NPTSF             : Number of points in the QBC array on fine grid
C QBCF              : Primitive variables on fine grid
C IBCDTF            : Data controlling boundary conditions on fine grid
C NPTSC             : Number of points in the QBC array on coarse grid
C QBCC              : Primitive variables on coarse grid
C IBCDTC            : Data controlling boundary conditions on coarse grid:
C                     IBCDAT( 1,ibc) = Boundary Condition Type
C                             2      = Direction (1=i, 2=j, 3=k)
C                             3      = i(1 or idim)
C                             4      = istrt
C                             5      = iend
C                             6      = jstrt
C                             7      = jend
C                             8      = kstrt
C                             9      = kend
C                            10      = profile identifier
C IERRCD         : Error Code
C                  IERRCD = INO  if NO storage error or error reading file
C                  IERRCD = IYES if    storage error or error reading file
C
      include '../header/common.h'
C
      DIMENSION QBCF  (NPTSF,NQ),
     1          QBCC  (NPTSC,NQ)
C
      DIMENSION IBCDTF(10),
     1          IBCDTC(10)
C
C     Store boundary condition data in mnemonic names
C
      IBCTYP = IBCDTF(1)
      IDIR   = IBCDTF(2)
      INROUT = IBCDTF(3)
      ISTRTF = IBCDTF(4)
      IENDF  = IBCDTF(5)
      JSTRTF = IBCDTF(6)
      JENDF  = IBCDTF(7)
      KSTRTF = IBCDTF(8)
      KENDF  = IBCDTF(9)
C
      ISTRTC = IBCDTC(4)
      IENDC  = IBCDTC(5)
      JSTRTC = IBCDTC(6)
      JENDC  = IBCDTC(7)
      KSTRTC = IBCDTC(8)
      KENDC  = IBCDTC(9)
C
C     I Direction
C
      IF (IDIR .EQ. 1) THEN
         DO 120 KC = KSTRTC, KENDC
            KF = 2 * (KC - 1)
            DO 110 JC = JSTRTC, JENDC
               JF = 2 * (JC - 1)
               JKBCPC = (JC-JSTRTC+1) + (KC-KSTRTC)*(JENDC-JSTRTC+1)
               JKBCPF = (JF-JSTRTF+1) + (KF-KSTRTF)*(JENDF-JSTRTF+1)
               IF (THREED) THEN
                  DO 100 L = 1, NQ
                     QBCC(JKBCPC,L) = 0.25E0 * 
     1                             ( QBCF(JKBCPF,L) 
     2                             + QBCF(JKBCPF+1,L)     
     3                             + QBCF(JKBCPF+(JENDF-JSTRTF+1),L) 
     4                             + QBCF(JKBCPF+(JENDF-JSTRTF+1)+1,L) )
  100             CONTINUE
               ELSE
                  DO 105 L = 1, NQ
                     QBCC(JKBCPC,L) = 0.50E0 * (QBCF(JKBCPF,L) + 
     1                                          QBCF(JKBCPF+1,L))
  105             CONTINUE
               ENDIF
  110       CONTINUE
  120    CONTINUE
C
C     J Direction
C
      ELSE IF (IDIR .EQ. 2) THEN
         DO 220 IC = ISTRTC, IENDC
            IF = 2 * (IC - 1)
            DO 210 KC = KSTRTC, KENDC
               KF = 2 * (KC - 1)
               IKBCPC = (KC-KSTRTC+1) + (IC-ISTRTC)*(KENDC-KSTRTC+1)
               IKBCPF = (KF-KSTRTF+1) + (IF-ISTRTF)*(KENDF-KSTRTF+1)
               IF (THREED) THEN
                  DO 200 L = 1, NQ
                     QBCC(IKBCPC,L) = 0.25E0 * 
     1                             ( QBCF(IKBCPF,L) 
     2                             + QBCF(IKBCPF+1,L)     
     3                             + QBCF(IKBCPF+(KENDF-KSTRTF+1),L) 
     4                             + QBCF(IKBCPF+(KENDF-KSTRTF+1)+1,L) )
  200             CONTINUE
               ELSE
                  DO 205 L = 1, NQ
                     QBCC(IKBCPC,L) = 0.50E0 * (QBCF(IKBCPF,L) + 
     1                                          QBCF(IKBCPF+1,L))
  205             CONTINUE
               ENDIF
  210       CONTINUE
  220    CONTINUE
C
C     K Direction
C
      ELSE IF (IDIR .EQ. 3) THEN
         DO 320 IC = ISTRTC, IENDC
            IF = 2 * (IC - 1)
            DO 310 JC = JSTRTC, JENDC
               JF = 2 * (JC - 1)
               IJBCPC = (JC-JSTRTC+1) + (IC-ISTRTC)*(JENDC-JSTRTC+1)
               IJBCPF = (JF-JSTRTF+1) + (IF-ISTRTF)*(JENDF-JSTRTF+1)
               IF (THREED) THEN
                  DO 300 L = 1, NQ
                     QBCC(IJBCPC,L) = 0.25E0 * 
     1                             ( QBCF(IJBCPF,L) 
     2                             + QBCF(IJBCPF+1,L)     
     3                             + QBCF(IJBCPF+(JENDF-JSTRTF+1),L) 
     4                             + QBCF(IJBCPF+(JENDF-JSTRTF+1)+1,L) )
  300             CONTINUE
               ELSE
                  DO 305 L = 1, NQ
                     QBCC(IJBCPC,L) = 0.50E0 * (QBCF(IJBCPF,L) + 
     1                                          QBCF(IJBCPF+1,L))
  305             CONTINUE
               ENDIF
  310       CONTINUE
  320    CONTINUE
      ELSE
         IERRCD = IYES
         WRITE (IOUT,1200) IDIR
         RETURN
      ENDIF
C
C     Finished reading profile boundary data
C
      RETURN
C
 1200 FORMAT (' ','CRSPRF: ERROR-> Invalid direction specified ',I4,
     1       /' ', 8X,'RUN ABORTING!')
C
      END





      SUBROUTINE INT2HQ (IDIMC, JDIMC, KDIMC, QC,
     1                   IDIMF, JDIMF, KDIMF, QF)
C
C Routine to interpolate (Prolong) from the coarse grid (IDIMC ... QC)
C to the fine grid (IDIMF ... QF).  Interpolate only the interior, use
C boundary condition routine to set boundaries.
C
C This routine is organized as follows:
C     1. Interpolate in the J direction
C     2. Interpolate in the I direction
C     3. If 3D, interpolate in the K direction
C
C IDIMC,JDIMC,KDIMC : Coarse grid dimensions
C QC                : Coarse grid primitive variables at cell centers
C IDIMF,JDIMF,KDIMF : Fine   grid dimensions
C QF                : Fine   grid primitive variables at cell centers
C
      include '../header/common.h'
C
      DIMENSION QC    (0:JDIMC+2,0:KDIMC+2,0:IDIMC+2,NQ),
     1          QF    (0:JDIMF+2,0:KDIMF+2,0:IDIMF+2,NQ)
C
C 1.  Interpolate in the J direction.
C
      DO 140 L = 1, NQ
         DO 130 J = 1, JDIMC
            JF = 2 * J - 1
            DO 120 I = 1, IDIMC+1
               DO 110 K = 1, KDIMC+1
                  QF(JF  ,K,I,L) = .75E0*QC(J,K,I,L)+.25E0*QC(J+1,K,I,L)
                  QF(JF+1,K,I,L) = .25E0*QC(J,K,I,L)+.75E0*QC(J+1,K,I,L)
 110           CONTINUE
 120        CONTINUE
 130     CONTINUE
 140  CONTINUE
C
C 2.  Interpolate in the I direction.
C     Decrement I and interpolate information stored beginning in the IDIMC 
C     location of QF into the IDIMF location of QF, repeat with IDIM-1, etc.
C     Use temporary variables to avoid overwriting data before it is used.
C
      DO 240 L = 1, NQ
         DO 230 I = IDIMC, 1, -1
            IF = 2 * I - 1
            DO 220 K = 1, KDIMC+1
               DO 210 J = 1, JDIMF+1
                  QI             = QF(J,K,I  ,L)
                  QIP1           = QF(J,K,I+1,L)
                  QF(J,K,IF  ,L) = 0.75E0 * QI + 0.25E0 * QIP1
                  QF(J,K,IF+1,L) = 0.25E0 * QI + 0.75E0 * QIP1
 210           CONTINUE
 220        CONTINUE
 230     CONTINUE
 240  CONTINUE
C
C 3.  If 3D, interpolate in the K direction.
C     Decrement K and interpolate information stored beginning in the KDIMC 
C     location of QF into the KDIMF location of QF, repeat with KDIM-1, etc.
C     Use temporary variables to avoid overwriting data before it is used.
C
      IF (THREED) THEN
         DO 340 L = 1, NQ
            DO 330 K = KDIMC, 1, -1
               KF = 2 * K - 1
               DO 320 I = 1, IDIMF+1
                  DO 310 J = 1, JDIMF+1
                     QK             = QF(J,K  ,I,L)
                     QKP1           = QF(J,K+1,I,L)
                     QF(J,KF  ,I,L) = 0.75E0 * QK + 0.25E0 * QKP1
                     QF(J,KF+1,I,L) = 0.25E0 * QK + 0.75E0 * QKP1
 310              CONTINUE
 320           CONTINUE
 330        CONTINUE
 340     CONTINUE
      ENDIF
C
C Finished interpolating from coarse to fine grid
C
      RETURN
      END





      SUBROUTINE PROLNG (IDIMC, JDIMC, KDIMC, QC, QCRS, NBCS, IBCDTC,
     1                   IDIMF, JDIMF, KDIMF, QF, DQ, IBCDTF,
     2                   ITURB, QMIN, NFIXQ, NRELIZ,
     3                   MXSECT, ABC, RTEMP)
C
C Routine to Prolong from the coarse grid (IDIMC ... QC) to the
C fine grid (IDIMF ... QF).
C
C N.B. This routine must maintain the values of QC on the edges for 
C      additional cycles.
C
C This routine is organized as follows:
C     1. Transfer boundary conditions for coarse grids to ghost cells
C     2. Calculate DQ = Q^(n+1) - Q^(n) on coarse grid
C     3. Transfer boundary conditions for QC back to surface to prepare for 
C        next cycle.
C     4. Apply boundary conditions to delta q before prolonging
C     5. Interpolate DQ from coarse grid to fine grid
C     6. Smooth DQ
C     7. Constrain DQ
C     8. Update fine grid solution as Q^(n+1) = Q^(n) + DQ
C     9. Check for positivity of variables
C    10. Enforce realizability of Reynolds stresses for Reynolds stress models
C
C IDIMC,JDIMC,KDIMC : Coarse grid dimensions
C QC                : Coarse grid primitive variables at cell centers
C QCRS              : Restricted fine grid primitive variables
C NBCS              : Number of boundary conditions
C IBCDTC            : Data controlling boundary conditions for coarse grid
C IDIMF,JDIMF,KDIMF : Fine   grid dimensions
C QF                : Fine   grid primitive variables at cell centers
C DQ                : Storage for interpolated delta(Q)
C IBCDTF            : Data controlling boundary conditions for fine grid
C ITURB             : Turbulence model
C QMIN              : Control array for minimum value of Q
C NFIXQ             : Number of times reset Q
C NRELIZ            : Number of times violate realizability
C MXSECT            : Maximum size of any (ixj, jxk or ixk)
C                     (used in dimensioning of temporary space for solver)
C ABC               : Storage for coefficient matrices for tridiagonal system
C RTEMP             : Storage for right hand side of implicit system
C
      include '../header/common.h'
C
      DIMENSION QC    (0:JDIMC+2,0:KDIMC+2,0:IDIMC+2,NQ),
     1          QCRS  (0:JDIMC+2,0:KDIMC+2,0:IDIMC+2,NQ),
     2          QF    (0:JDIMF+2,0:KDIMF+2,0:IDIMF+2,NQ),
     3          DQ    (0:JDIMF+2,0:KDIMF+2,0:IDIMF+2,NQ)
C
      DIMENSION IBCDTC(9,NBCS),
     1          IBCDTF(9,NBCS)
C
      DIMENSION QMIN  (NQ)
      DIMENSION NFIXQ (NQ),
     1          NRELIZ(NQ)
C
      DIMENSION ABC   (MXSECT,3),
     1          RTEMP (MXSECT)
C
C 1.  Transfer boundary conditions for coarse grids to ghost cells
C
      CALL QBCMOD (IDIMC, JDIMC, KDIMC, QC  , NBCS, IBCDTC)
      CALL QBCMOD (IDIMC, JDIMC, KDIMC, QCRS, NBCS, IBCDTC)
C
C 2.  Calculate DQ = Q^(n+1) - Q^(n) on coarse grid
C     Use QCRS to store DQ to maintain QC's values in the array edges for
C     additional cycles.
C
      DO 240 L = 1, NQ
         DO 230 I = 0, IDIMC+2
            DO 220 K = 0, KDIMC+2
               DO 210 J = 0, JDIMC+2
                  QCRS(J,K,I,L) = QC(J,K,I,L) - QCRS(J,K,I,L)
  210          CONTINUE
  220       CONTINUE
  230    CONTINUE
  240 CONTINUE
C
C 3.  Transfer boundary conditions for QC back to surface to prepare for 
C     next cycle.
C
      CALL QBCFIX (IDIMC, JDIMC, KDIMC, QC  , NBCS, IBCDTC)
C
C 4.  Apply boundary conditions to delta q before prolonging
C
      CALL DQBC   (IDIMC, JDIMC, KDIMC, QCRS, NBCS, IBCDTC)
C
C 5.  Interpolate DQ from coarse grid to fine grid
C
      CALL INT2HQ (IDIMC, JDIMC, KDIMC, QCRS,
     1             IDIMF, JDIMF, KDIMF, DQ)
C
C 6.  Smooth DQ
C
      CALL SMOOTH (IDIMF, JDIMF, KDIMF, DQ, 
     1             MXSECT, ABC(1,1), ABC(1,2), ABC(1,3), RTEMP)
C
C 7.  Constrain DQ
C
      NPTS = (IDIMF + 3) * (JDIMF + 3) * (KDIMF + 3)
c     CALL DQLIMT (NPTS, QF, DQ)
      CALL xDQLIMT (idimf, jdimf, kdimf, QF, DQ)
C
C 8.  Update fine grid solution as Q^(n+1) = Q^(n) + DQ
C
      DO 840 L = 1, NQ
         DO 830 I = 2, IDIMF
            DO 820 K = 2, KDIMF
               DO 810 J = 2, JDIMF
                  QF(J,K,I,L) = QF(J,K,I,L) + DQ(J,K,I,L)
  810          CONTINUE
  820       CONTINUE
  830    CONTINUE
  840 CONTINUE
C
C 9.  Check for positivity of flow variables
C
      IPSTRT = 1
      IPEND  = IDIMF
      IF (IFIXQ .NE. INO) THEN
         CALL FIXQ   (ITURB, IDIMF, JDIMF, KDIMF, IPSTRT, IPEND, QF,
     1                QMIN,  NFIXQ, IERRCD)
      ENDIF
C
C 10. Enforce realizability of Reynolds stresses for Reynolds stress models
C
      IF (ITURB .EQ. ITRS) THEN
         CALL REALIZ (ITURB, IDIMF, JDIMF, KDIMF, IPSTRT, IPEND, QF, 
     1                NRELIZ)
      ENDIF
C
C     Finished prolongation
C
      RETURN
      END



      SUBROUTINE DQBC   (IDIM, JDIM, KDIM, DQ, NBCS, IBCDAT)
C
C Routine to calculate delta Q on the boundary data for use in the multigrid
C prolongation operator.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Loop over all boundary conditions
C        If boundary condition then
C        a. Set delta(Q) boundary = 0
C
C IDIM,JDIM,KDIM : Dimensions of current block
C DQ             : Delta Q
C NBCS           : Number of boundary conditions for current block
C IBCDAT         : Data controlling boundary conditions:
C                  IBCDAT(1,ibc) = Boundary Condition Type
C                         2      = Direction (1=i, 2=j, 3=k)
C                         3      = i(1 or idim)
C                         4      = istrt
C                         5      = iend
C                         6      = jstrt
C                         7      = jend
C                         8      = kstrt
C                         9      = kend
C                         10     = profile identifier
C
      include '../header/common.h'
C
      DIMENSION DQ    (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION IBCDAT(10,NBCS)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
C
C 1. Setup memory for looping.
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
      IJKDIM(1) = IDIM
      IJKDIM(2) = JDIM
      IJKDIM(3) = KDIM
C
C 2. Loop over boundary conditions
C
      DO 200 IBCS = 1, NBCS
         IBCTYP   = IBCDAT(1,IBCS)
C
         IF (IBCTYP .EQ. ITAN .OR. IBCTYP .EQ. IWALFN) THEN
C
C Store boundary condition data in mnemonic names
C
            IDIR     = IBCDAT(2,IBCS)
            INROUT   = IBCDAT(3,IBCS)
            ISTRT(1) = IBCDAT(4,IBCS)
            IEND (1) = IBCDAT(5,IBCS)
            ISTRT(2) = IBCDAT(6,IBCS)
            IEND (2) = IBCDAT(7,IBCS)
            ISTRT(3) = IBCDAT(8,IBCS)
            IEND (3) = IBCDAT(9,IBCS)
C
C Modify beginning, ending indices to extend boundary condition to edge/corner
C
            DO 10 L = 1, 3
               IF (L .NE. IDIR) THEN
                  IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
                  IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
               END IF
   10       CONTINUE
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCN   is the point to the inside of the block from IBCPT
C
            IF (INROUT .EQ. 1) THEN
               IBCPT  = 1
               IBCN   = 2
            ELSE
               IBCPT  = IJKDIM(IDIR) + 1
               IBCN   = IJKDIM(IDIR)
            ENDIF
C
            IC1 = ICY (IDIR, 1)
            IC2 = ICY (IDIR, 2)
            IC3 = ICY (IDIR, 3)
C
C 2-a. Replace cell face data with equivalent ghostcell data
C
            DO 130 L = 1, NQ
               DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
                  DO 110 I2 = ISTRT(IC2), IEND(IC2)
C
                     IBC  = IBCPT*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JBC  = IBCPT*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KBC  = IBCPT*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                     IN   = IBCN *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JN   = IBCN *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KN   = IBCN *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
c--                  DQ(JBC,KBC,IBC,L) = DQ(JN,KN,IN,L)
                     DQ(JBC,KBC,IBC,L) = 0.E0
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
         ENDIF
  200 CONTINUE
C
C Finished modifying the boundary conditions for interpolation (prolongation)
C
      RETURN
      END





      SUBROUTINE QBCMOD (IDIM, JDIM, KDIM, Q, NBCS, IBCDAT)
C
C Routine to replace boundary condition data on a computational cell face
C with equivalent data at the ghost cell center to then use in the 
C interpolation (prolongation) operator.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Loop over all boundary conditions
C        If boundary condition is a half cell, then
C        a. Replace cell face data with equivalent ghostcell data
C
C IDIM,JDIM,KDIM : Dimensions of current block
C Q              : Primitive variables at cell centers
C NBCS           : Number of boundary conditions for current block
C IBCDAT         : Data controlling boundary conditions:
C                  IBCDAT(1,ibc) = Boundary Condition Type
C                         2      = Direction (1=i, 2=j, 3=k)
C                         3      = i(1 or idim)
C                         4      = istrt
C                         5      = iend
C                         6      = jstrt
C                         7      = jend
C                         8      = kstrt
C                         9      = kend
C                         10     = profile identifier
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION IBCDAT(10,NBCS)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
C
C 1. Setup memory for looping.
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
      IJKDIM(1) = IDIM
      IJKDIM(2) = JDIM
      IJKDIM(3) = KDIM
C
C 2. Loop over boundary conditions
C
      DO 200 IBCS = 1, NBCS
         IBCTYP   = IBCDAT(1,IBCS)
C
         IF (IBCTYP .EQ. ITAN .OR. IBCTYP .EQ. IWALL .OR. 
     1       IBCTYP .EQ. IWALFN) THEN
C
C Store boundary condition data in mnemonic names
C
            IDIR     = IBCDAT(2,IBCS)
            INROUT   = IBCDAT(3,IBCS)
            ISTRT(1) = IBCDAT(4,IBCS)
            IEND (1) = IBCDAT(5,IBCS)
            ISTRT(2) = IBCDAT(6,IBCS)
            IEND (2) = IBCDAT(7,IBCS)
            ISTRT(3) = IBCDAT(8,IBCS)
            IEND (3) = IBCDAT(9,IBCS)
C
C Modify beginning, ending indices to extend boundary condition to edge/corner
C
            DO 10 L = 1, 3
               IF (L .NE. IDIR) THEN
                  IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
                  IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
               END IF
   10       CONTINUE
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCN   is the point to the inside of the block from IBCPT
C
            IF (INROUT .EQ. 1) THEN
               IBCPT  = 1
               IBCN   = 2
            ELSE
               IBCPT  = IJKDIM(IDIR) + 1
               IBCN   = IJKDIM(IDIR)
            ENDIF
C
            IC1 = ICY (IDIR, 1)
            IC2 = ICY (IDIR, 2)
            IC3 = ICY (IDIR, 3)
C
C 2-a. Replace cell face data with equivalent ghostcell data
C
            DO 130 L = 1, NQ
               DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
                  DO 110 I2 = ISTRT(IC2), IEND(IC2)
C
                     IBC  = IBCPT*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JBC  = IBCPT*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KBC  = IBCPT*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                     IN   = IBCN *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JN   = IBCN *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KN   = IBCN *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                     DELTAQ = Q(JBC,KBC,IBC,L) - Q(JN ,KN ,IN ,L)
                     Q(JBC,KBC,IBC,L) = Q(JBC,KBC,IBC,L) + DELTAQ
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
         ENDIF
  200 CONTINUE
C
C Finished modifying the boundary conditions for interpolation (prolongation)
C
      RETURN
      END





      SUBROUTINE QBCFIX (IDIM, JDIM, KDIM, Q, NBCS, IBCDAT)
C
C Routine to undo QBCMOD modifications moving ghost cells back to faces
C after interpolation (prolongation) operator.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Loop over all boundary conditions
C        If boundary condition is a half cell, then
C        a. Replace ghostcell data with equivalent cell face data 
C
C IDIM,JDIM,KDIM : Dimensions of current block
C Q              : Primitive variables at cell centers
C NBCS           : Number of boundary conditions for current block
C IBCDAT         : Data controlling boundary conditions:
C                  IBCDAT(1,ibc) = Boundary Condition Type
C                         2      = Direction (1=i, 2=j, 3=k)
C                         3      = i(1 or idim)
C                         4      = istrt
C                         5      = iend
C                         6      = jstrt
C                         7      = jend
C                         8      = kstrt
C                         9      = kend
C                         10     = profile identifier
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION IBCDAT(10,NBCS)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
C
C 1. Setup memory for looping.
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
      IJKDIM(1) = IDIM
      IJKDIM(2) = JDIM
      IJKDIM(3) = KDIM
C
C 2. Loop over boundary conditions
C
      DO 200 IBCS = 1, NBCS
         IBCTYP   = IBCDAT(1,IBCS)
C
         IF (IBCTYP .EQ. ITAN .OR. IBCTYP .EQ. IWALL .OR. 
     1       IBCTYP .EQ. IWALFN) THEN
C
C Store boundary condition data in mnemonic names
C
            IDIR     = IBCDAT(2,IBCS)
            INROUT   = IBCDAT(3,IBCS)
            ISTRT(1) = IBCDAT(4,IBCS)
            IEND (1) = IBCDAT(5,IBCS)
            ISTRT(2) = IBCDAT(6,IBCS)
            IEND (2) = IBCDAT(7,IBCS)
            ISTRT(3) = IBCDAT(8,IBCS)
            IEND (3) = IBCDAT(9,IBCS)
C
C Modify beginning, ending indices to extend boundary condition to edge/corner
C
            DO 10 L = 1, 3
               IF (L .NE. IDIR) THEN
                  IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
                  IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
               END IF
   10       CONTINUE
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCN   is the point to the inside of the block from IBCPT
C
            IF (INROUT .EQ. 1) THEN
               IBCPT  = 1
               IBCN   = 2
            ELSE
               IBCPT  = IJKDIM(IDIR) + 1
               IBCN   = IJKDIM(IDIR)
            ENDIF
C
            IC1 = ICY (IDIR, 1)
            IC2 = ICY (IDIR, 2)
            IC3 = ICY (IDIR, 3)
C
C 2-a. Replace cell face data with equivalent ghostcell data
C
            DO 130 L = 1, NQ
               DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
                  DO 110 I2 = ISTRT(IC2), IEND(IC2)
C
                     IBC  = IBCPT*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JBC  = IBCPT*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KBC  = IBCPT*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                     IN   = IBCN *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                     JN   = IBCN *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                     KN   = IBCN *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                     DELTAQ = Q(JBC,KBC,IBC,L) - Q(JN ,KN ,IN ,L)
                     Q(JBC,KBC,IBC,L) = Q(JBC,KBC,IBC,L) - 0.5E0*DELTAQ
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
         ENDIF
  200 CONTINUE
C
C Finished modifying the boundary conditions for interpolation (prolongation)
C
      RETURN
      END





      SUBROUTINE SMOOTH (IDIM, JDIM, KDIM, PHI, MXSECT, A, B, C, RHS)
C
C Routine to perform the implicit smoothing on the variable phi.
C
C This routine is organized as follows:
C     0. Set up constants
C     1. Implicit operator in I direction
C     2. Implicit operator in J direction
C     3. Implicit operator in K direction
C
C IDIM,JDIM,KDIM : Dimensions of current block
C PHI            : Variable to be smoothed
C MXSECT         : Maximum size of any (ixj, jxk or ixk)
C                  (used in dimensioning of temporary space for solver)
C A, B, C        : Storage for coefficient matrices for tridiagonal system
C RHS            : Storage for right hand side of implicit system
C
      include '../header/common.h'
C
      DIMENSION PHI   (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION A     (MXSECT),
     1          B     (MXSECT),
     2          C     (MXSECT),
     3          RHS   (MXSECT)
C
C 0.  Set up constants
C
      CSMTHI = CSMTH
      CSMTHJ = CSMTH
      CSMTHK = CSMTH
C
C 1.  Implicit operator in I direction
C
      CALL SMTHI  (IDIM, JDIM, KDIM, PHI, CSMTHI, A, B, C, RHS)
C
C 2.  Implicit operator in J direction
C
      CALL SMTHJ  (IDIM, JDIM, KDIM, PHI, CSMTHI, A, B, C, RHS)
C
C 3.  Implicit operator in K direction
C
      IF (THREED) THEN
         CALL SMTHK  (IDIM, JDIM, KDIM, PHI, CSMTHI, A, B, C, RHS)
      ENDIF
C
C Finished
C
      RETURN
      END





      SUBROUTINE SMTHI  (IDIM, JDIM, KDIM, PHI, CSMTHI, A, B, C, RHS)
C
C Routine to perform the I-direction implicit smoothing on the variable phi.
C
C This routine is organized as follows:
C     0. Set up constants
C     1. Constant coefficient implicit smoothing
C        a. Set up implicit operator
C           i.  Interior stencil for implicit operator
C           ii. Boundary conditions for implicit operator
C        b. Solve system
C        c. Update smoothed solution
C
C IDIM,JDIM,KDIM : Dimensions of current block
C PHI            : Variable to be smoothed
C CSMTHI         : Coefficient for I direction smoothing
C A, B, C        : Storage for coefficient matrices for tridiagonal system
C RHS            : Storage for right hand side of implicit system
C
      include '../header/common.h'
C
      DIMENSION PHI   (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION A     (2:JDIM,2:IDIM),
     1          B     (2:JDIM,2:IDIM),
     2          C     (2:JDIM,2:IDIM),
     3          RHS   (2:JDIM,2:IDIM)
C
C 0.   Set up constants
C
      NSIM = JDIM - 1
      NPTS = IDIM - 1
C
C 1.   Constant coefficient implicit smoothing
C
      DO L = 1, NQ
         DO K = 2, KDIM
C
C 1-a. Set up implicit operator
C 1-a-i.  Interior stencil for implicit operator
C
            DO J = 2, JDIM
               DO I = 2, IDIM
                  A  (J,I) = - CSMTHI
                  B  (J,I) = 1.E0 + 2.E0 * CSMTHI
                  C  (J,I) = - CSMTHI
                  RHS(J,I) = PHI(J,K,I,L)
               END DO
C
C 1-a-ii. Boundary conditions for implicit operator
C
            END DO
C
C 1-b. Solve system
C
            CALL VSTRI (NPTS, NSIM, A, B, C, RHS)
C
C 1-c. Update smoothed solution
C
            DO J = 2, JDIM
               DO I = 2, IDIM
                  PHI(J,K,I,L) = RHS(J,I)
               END DO
            END DO
C
C      Next K
C
         END DO
      END DO
C
C Finished
C
      RETURN
      END





      SUBROUTINE SMTHJ  (IDIM, JDIM, KDIM, PHI, CSMTHJ, A, B, C, RHS)
C
C Routine to perform the J-direction implicit smoothing on the variable phi.
C
C This routine is organized as follows:
C     0. Set up constants
C     1. Constant coefficient implicit smoothing
C        a. Set up implicit operator
C           i.  Interior stencil for implicit operator
C           ii. Boundary conditions for implicit operator
C        b. Solve system
C        c. Update smoothed solution
C
C IDIM,JDIM,KDIM : Dimensions of current block
C PHI            : Variable to be smoothed
C CSMTHJ         : Coefficient for J direction smoothing
C A, B, C        : Storage for coefficient matrices for tridiagonal system
C RHS            : Storage for right hand side of implicit system
C
      include '../header/common.h'
C
      DIMENSION PHI   (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION A     (2:IDIM,2:JDIM),
     1          B     (2:IDIM,2:JDIM),
     2          C     (2:IDIM,2:JDIM),
     3          RHS   (2:IDIM,2:JDIM)
C
C 0.   Set up constants
C
      NSIM = IDIM - 1
      NPTS = JDIM - 1
C
C 1.   Constant coefficient implicit smoothing
C
      DO L = 1, NQ
         DO K = 2, KDIM
C
C 1-a. Set up implicit operator
C 1-a-i.  Interior stencil for implicit operator
C
            DO I = 2, IDIM
               DO J = 2, JDIM
                  A  (I,J) = - CSMTHJ
                  B  (I,J) = 1.E0 + 2.E0 * CSMTHJ
                  C  (I,J) = - CSMTHJ
                  RHS(I,J) = PHI(J,K,I,L)
               END DO
C
C 1-a-ii. Boundary conditions for implicit operator
C
            END DO
C
C 1-b. Solve system
C
            CALL VSTRI (NPTS, NSIM, A, B, C, RHS)
C
C 1-c. Update smoothed solution
C
            DO I = 2, IDIM
               DO J = 2, JDIM
                  PHI(J,K,I,L) = RHS(I,J)
               END DO
            END DO
C
C      Next K
C
         END DO
      END DO
C
C Finished
C
      RETURN
      END





      SUBROUTINE SMTHK  (IDIM, JDIM, KDIM, PHI, CSMTHK, A, B, C, RHS)
C
C Routine to perform the K-direction implicit smoothing on the variable phi.
C
C This routine is organized as follows:
C     0. Set up constants
C     1. Constant coefficient implicit smoothing
C        a. Set up implicit operator
C           i.  Interior stencil for implicit operator
C           ii. Boundary conditions for implicit operator
C        b. Solve system
C        c. Update smoothed solution
C
C IDIM,JDIM,KDIM : Dimensions of current block
C PHI            : Variable to be smoothed
C CSMTHK         : Coefficient for K direction smoothing
C A, B, C        : Storage for coefficient matrices for tridiagonal system
C RHS            : Storage for right hand side of implicit system
C
      include '../header/common.h'
C
      DIMENSION PHI   (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION A     (2:JDIM,2:KDIM),
     1          B     (2:JDIM,2:KDIM),
     2          C     (2:JDIM,2:KDIM),
     3          RHS   (2:JDIM,2:KDIM)
C
C 0.   Set up constants
C
      NSIM = JDIM - 1
      NPTS = KDIM - 1
C
C 1.   Constant coefficient implicit smoothing
C
      DO L = 1, NQ
         DO I = 2, IDIM
C
C 1-a. Set up implicit operator
C 1-a-i.  Interior stencil for implicit operator
C
            DO J = 2, JDIM
               DO K = 2, KDIM
                  A  (J,K) = - CSMTHK
                  B  (J,K) = 1.E0 + 2.E0 * CSMTHK
                  C  (J,K) = - CSMTHK
                  RHS(J,K) = PHI(J,K,I,L)
               END DO
C
C 1-a-ii. Boundary conditions for implicit operator
C
            END DO
C
C 1-b. Solve system
C
            CALL VSTRI (NPTS, NSIM, A, B, C, RHS)
C
C 1-c. Update smoothed solution
C
            DO J = 2, JDIM
               DO K = 2, KDIM
                  PHI(J,K,I,L) = RHS(J,K)
               END DO
            END DO
C
C      Next I
C
         END DO
      END DO
C
C Finished
C
      RETURN
      END






      SUBROUTINE DQLIMT (NPTS, Q, DQ)
C
C Routine to perform limiting of the coarse grid correction.
C
C This routine is organized as follows:
C     1. Limit mean flow coarse grid correction
C     2. Limit turbulent equations coarse grid correction
C
C NPTS           : Number of points to limit
C Q              : Solution vector
C DQ             : Coarse grid corrections to be limited
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          DQ    (NPTS,NQ)
C
C 1.  Limit mean flow coarse grid correction
C     If coarse grid update for DENSITY exceeds CMXCHG*DENSITY, then
C     set the coarse grid correction for all mean flow variables to zero.
C
      DO I = 1, NPTS
         IF (ABS(DQ(I,1)) .GT. CMXCHG*Q(I,1)) THEN
            DQ(I,1) = 0.E0
            DQ(I,2) = 0.E0
            DQ(I,3) = 0.E0
            DQ(I,4) = 0.E0
            DQ(I,5) = 0.E0
         ENDIF
      END DO
C
C 2.  Limit turbulent equations coarse grid correction
C
c-N.B. Temporary for testing: set DQ for turbulence to zero
      DO L = 6, NQ
         DO I = 1, NPTS
            DQ(I,L) = 0.E0
         END DO
      END DO
C
C Finished
C
      RETURN
      END
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE xDQLIMT (idim, jdim, kdim, Q, DQ)
C
C Routine to perform limiting of the coarse grid correction.
C
C This routine is organized as follows:
C     1. Limit mean flow coarse grid correction
C     2. Limit turbulent equations coarse grid correction
C
C NPTS           : Number of points to limit
C Q              : Solution vector
C DQ             : Coarse grid corrections to be limited
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          DQ    (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
C 1.  Limit mean flow coarse grid correction
C     If coarse grid update for DENSITY exceeds CMXCHG*DENSITY, then
C     set the coarse grid correction for all mean flow variables to zero.
C
      do k = 2, kdim
      DO I = 2, idim
      do j = 2, jdim
         IF (ABS(DQ(J,K,I,1)) .GT. CMXCHG*Q(J,K,I,1)) THEN
            DQ(J,K,I,1) = 0.E0
            DQ(J,K,I,2) = 0.E0
            DQ(J,K,I,3) = 0.E0
            DQ(J,K,I,4) = 0.E0
            DQ(J,K,I,5) = 0.E0
         ENDIF
      END DO
      END DO
      END DO
C
C 2.  Limit turbulent equations coarse grid correction
C
c-N.B. Temporary for testing: set DQ for turbulence to zero
      DO L = 6, NQ
         do k = 2, kdim
         DO I = 2, idim
         do j = 2, jdim
            DQ(J,K,I,L) = 0.E0
         END DO
         END DO
         END DO
      END DO
C
C Finished
C
      RETURN
      END





      SUBROUTINE CMGFF  (NPTS, RESCRS, RES)
C
C Routine to calculate the coarse grid forcing function as 
C       RESCRS = sum RES(fine) - RES(coarse)
C RESCRS contains sum RES(fine) on entry.
C
C NPTS   : Number of points in array
C RESCRS : MG forcing function
C RES    : Coarse grid residual
C
      include '../header/common.h'
C
      INTEGER   NPTS
C
      DIMENSION RESCRS(NPTS),
     1          RES   (NPTS)
C
C Local variables
C
      INTEGER   I
C
C Calculate RESCRS =  sum RES(fine) - RES(coarse)
C
      DO I = 1, NPTS
         RESCRS(I) = RESCRS(I) - RES(I)
      END DO
C
      RETURN
C
C Finished
C
      END





      SUBROUTINE CRESFF (NPTS, RESCRS, RES)
C
C Routine to calculate the residual plus forcing function for the 
C coarse grid update
C       RES = [sum RES(fine) - RES(coarse)] + RES(coarse)
C RESCRS contains [sum RES(fine) - RES(coarse)] on entry.
C
C NPTS   : Number of points in array
C RESCRS : MG forcing function
C RES    : Coarse grid residual
C
      include '../header/common.h'
C
      INTEGER   NPTS
C
      DIMENSION RESCRS(NPTS),
     1          RES   (NPTS)
C
C Local variables
C
      INTEGER   I
C
C Calculate RESCRS =  sum RES(fine) - RES(coarse)
C
      DO I = 1, NPTS
         RES(I) = RES(I) + RESCRS(I)
      END DO
C
      RETURN
C
C Finished
C
      END

      SUBROUTINE xCMGFF  (idim, jdim, kdim, RESCRS, RES)
C
C Routine to calculate the coarse grid forcing function as 
C       RESCRS = sum RES(fine) - RES(coarse)
C RESCRS contains sum RES(fine) on entry.
C
C NPTS   : Number of points in array
C RESCRS : MG forcing function
C RES    : Coarse grid residual
C
      include '../header/common.h'
C
      INTEGER   NPTS
C
      DIMENSION RESCRS(jdim+1,kdim+1,idim+1,nf),
     1          RES   (jdim+1,kdim+1,idim+1,nf)
C
C Local variables
C
      INTEGER   I, j, k, l
C
C Calculate RESCRS =  sum RES(fine) - RES(coarse)
C
      do l = 1, nf
      do k = 2, kdim
      do i = 2, idim
      do j = 2, jdim
         RESCRS(j,k,I,l) = RESCRS(j,k,i,l) - RES(j,k,i,l)
      END DO
      END DO
      END DO
      END DO
C
      RETURN
C
C Finished
C
      END





      SUBROUTINE xCRESFF (idim, jdim, kdim, RESCRS, RES)
C
C Routine to calculate the residual plus forcing function for the 
C coarse grid update
C       RES = [sum RES(fine) - RES(coarse)] + RES(coarse)
C RESCRS contains [sum RES(fine) - RES(coarse)] on entry.
C
C NPTS   : Number of points in array
C RESCRS : MG forcing function
C RES    : Coarse grid residual
C
      include '../header/common.h'
C
      INTEGER   NPTS
C
      DIMENSION RESCRS(jdim+1,kdim+1,idim+1,nf),
     1          RES   (jdim+1,kdim+1,idim+1,nf)
C
C Local variables
C
      INTEGER   I, j, k, l
C
C Calculate RESCRS =  sum RES(fine) - RES(coarse)
C
      do l = 1, nf
      do k = 2, kdim
      do i = 2, idim
      do j = 2, jdim
         RES(j,k,i,l) = RES(j,k,i,l) + RESCRS(j,k,i,l)
      END DO
      END DO
      END DO
      END DO
C
      RETURN
C
C Finished
C
      END
