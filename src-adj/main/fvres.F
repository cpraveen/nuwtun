      SUBROUTINE FVRESID (NCFS, NBLKS, IDIM, JDIM, KDIM, NXSECT, Q0, QN,
     1                    NQNPTS,
     4                    PROPS, NPPTS, S, VOL, R, DTJ, RES, 
     3                    QBC, FHAT, DQP, DQM, QP, QM, TAU, F,
     4                    FM, STEMP, WORKX, IBCDAT, NBCS, 
     5                    IORDER, LIMIT, RKAP, IFLUX, VISCOS, 
     6                    IMODEL, INITRS, IFDELQ, MRCHNG, 
     7                    VGNRN,SIGMA, NPRLIM, NPRNEG, IFL2DQ,
     8                    R2DQDT, CL, CD, ITURB, NXPPTS,
     9                    IOFF, IOFFS, IOFFQ, IOFFP, IOFFQN, IOFFQC,
     A                    IOFQBC, ICUTS, NCUTS, NXPROF, NPROFL,
     B                    EPS, QD, RESD, NSIZEQ, ITER, NITFO)
C
C     Commons
C
      include '../../src-flo/header/maxdim.h'
      include '../../src-flo/header/common.h'
      include '../../src-flo/header/histry.h'
C
C     PERMANENT STORAGE FOR ALL ZONES/BLOCKS
C
C     Storage for flow variables, turbulent eddy viscosity,
C     metrics, volumes, jacobians and grid.
C
      REAL Q0(*), QN(*), PROPS(*), S(*), VOL(*), DTJ(*), R(*), QBC(*)
C
C     Dimensions of domains
C
      DIMENSION IDIM  (MXBLKS,MXLVLS), 
     1          JDIM  (MXBLKS,MXLVLS),
     2          KDIM  (MXBLKS,MXLVLS)
C
C     Offsets of domains used to calculate location in total storage
C     of Q,S,... of a given domain.
C     IOFFQN offset for previous time level storage for QN
C     IOFFQC offset for coarse grid storage for QCRS and RESCRS.
C     IOFQBC offset for multiple profile boundary condition segments
C
      DIMENSION IOFF  (MXBLKS,MXLVLS), 
     1          IOFFS (MXBLKS,MXLVLS),
     2          IOFFQ (MXBLKS,MXLVLS),
     3          IOFFP (MXBLKS,MXLVLS),
     4          IOFFQN(MXBLKS,MXLVLS),
     5          IOFFQC(MXBLKS,MXLVLS),
     6          IOFQBC(MXPROF,MXLVLS)
C
C     Boundary condition data.
C     IBCPRF is re-ordered BC data to read in profile BC segments.
C
      DIMENSION IBCDAT(10,MXBCS,MXBLKS,MXLVLS),  NBCS(MXBLKS),
     1          IBCPRF(10,MXPROF,MXLVLS)
C
C     Cut (domain to domain communication) data.
C
      DIMENSION ICUTS (21,MXCUTS,MXLVLS)
      INTEGER   NCUTS
C
C     Data controlling the accuracy and scheme used.
C
      DIMENSION IORDER(3,MXBLKS,MXLVLS,MXCFS), 
     1          LIMIT (3,MXBLKS), 
     2          RKAP  (3,MXBLKS), 
     3          IFLUX (MXBLKS), 
     4          ITORDR(3)
      LOGICAL   DIAG  (3,MXBLKS)
C
C     Data controlling the viscous/inviscid calculation
C
      LOGICAL   VISCOS(3,MXBLKS)
C
C     Data controlling calculation of delq
C
      DIMENSION IFDELQ(MAXQ)
C
C     Data controlling the turbulence model
C
      DIMENSION IMODEL(MAXMDL),
     1          IMODLX(MAXMDL)
C
C     Data controlling the specification of transition location
C
      DIMENSION ITRSEG(7,MXTRSG)
C
C     Storage for forces
C     FORCE   -> Forces integrated in x, y, z directions
C     RCNTR   -> Location to calculate point vortex and moments about
C
      DIMENSION FORCE(3),
     1          RCNTR(3)
C
C     Data for diagnostic output
C
      REAL      RESMAX
      INTEGER   IFRSMX,
     1          IRESMX(0:3)
      INTEGER   NFIXQ (MAXQ),
     1          NRELIZ(MAXQ),
     2          NPRLIM,
     3          NPRNEG
C
C     STORAGE REQUIRED FOR A ZONE/BLOCK
C
C     Storage for residuals
C          For (RK, AF)     : store residual over entire block
C          For (LGS, MARCH) : store residuals over a plane of block
C
      REAL RES(*)
C
C     Storage for fluxes for fourth order
C
      REAL FHAT(*)
C
C     Storage for MUSCL and flux evaluation
C
      REAL DQP(*), DQM(*), QP(*), QM(*), TAU(*), F(*), FM(*),
     1                                STEMP(*), WORKX(*)
C     Switches
C          RESTRT -> If true then do a restart
C          VGNRN  -> If true then use Vigneron technique when marching
C          CONTRB -> If true then converting turbulence models on restart
C          MRCHNG -> If true then use marching fluxes in residual
C          INITRS -> If true then initialize residuals to zero (used in MG)
C          IFL2DQ -> If true then calculate L2 norm of dQ/dt term (unsteady)
C          YNCALC -> If true then recalculate Y_normal on restart
C          IFMGFF -> If true then calculate the MG forcing function from the
C                    coarse grid residual in the MG cycle
C
      LOGICAL RESTRT, VGNRN, CONTRB, MRCHNG, INITRS, IFL2DQ, YNCALC, 
     1        IFMGFF

      REAL    EPS, QD(*), RESD(*)
C
C     Local variable
C
      REAL    Q(NSIZEQ) 
C
C     Copy Q0 into local variable Q
C
      Q(1:NSIZEQ) = Q0(1:NSIZEQ)
C
C     Set switch to test for storage errors
C
      IERROR = 0
C
C     Initialize force coefficients
C
      CL    = 0.E0
      CD    = 0.E0
C
C 12. TIME DEPENDENT CALCULATION
C
      MRCHNG = .FALSE.
C
C     Loop over coarse to fine sequences
C
      ITRUN  = 0
      DO 900 ICFS = 1, NCFS
        LVLONE = NCFS - ICFS + 1
        ILVL   = LVLONE
C
C     Counter for limiting of turbulence production
C
        NPRLIM = 0
        NPRNEG = 0
C
C     IF EPS>0 then perturb Q to Q + EPS*QD, QD is random
C
        IF(EPS.GT.0.0)THEN
          DO IBLK = 1, NBLKS
            IQ   = IOFFQ (IBLK,ILVL) * NQ          + 1
            CALL INTQD(EPS, IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                 KDIM(IBLK,ILVL), Q(IQ), QD(IQ))
          ENDDO   
        ENDIF
C
C     Do block boundary conditions (CUTS)
C
        DO 530 ICUT = 1, NCUTS
          IBLK1 = ICUTS( 1,ICUT,ILVL)
          IBLK2 = ICUTS(10,ICUT,ILVL)
          IQ1   = IOFFQ(IBLK1,ILVL) * NQ    + 1
          IQ2   = IOFFQ(IBLK2,ILVL) * NQ    + 1
          CALL CUT (IDIM(IBLK1,ILVL), JDIM(IBLK1,ILVL),
     2              KDIM(IBLK1,ILVL), Q(IQ1),
     3              IDIM(IBLK2,ILVL), JDIM(IBLK2,ILVL),
     4              KDIM(IBLK2,ILVL), Q(IQ2), 
     5              ICUTS(1,ICUT,ILVL), NQ)
          IF (ITURB .GE. ITLMNR) THEN
            IP1   = IOFFP(IBLK1,ILVL) * NP    + 1
            IP2   = IOFFP(IBLK2,ILVL) * NP    + 1
            CALL CUT (IDIM(IBLK1,ILVL), JDIM(IBLK1,ILVL),
     2                KDIM(IBLK1,ILVL), PROPS(IP1),
     3                IDIM(IBLK2,ILVL), JDIM(IBLK2,ILVL),
     4                KDIM(IBLK2,ILVL), PROPS(IP2), 
     5                ICUTS(1,ICUT,ILVL), NP)
          ENDIF
  530   CONTINUE
C
C     Loop over all blocks.
C
        DO 500 IBLK = 1, NBLKS
          IQ   = IOFFQ (IBLK,ILVL) * NQ          + 1
          IQN  = IOFFQN(IBLK,ILVL) * NQ * NTMLVL + 1
          IP   = IOFFP (IBLK,ILVL) * NP          + 1
          IS   = IOFFS (IBLK,ILVL) * 4 * 3       + 1
          IV   = IOFFS (IBLK,ILVL)               + 1
          IR   = IOFFS (IBLK,ILVL) * 3           + 1
          ID   = IOFF  (IBLK,ILVL)               + 1
          IABC = IOFF  (IBLK,ILVL)*NRANK*NRANK*3*3+1
C
C     Storage for QN, calculate L2 norm of dQ/dt only if first stage of
C     multistage scheme and on the finest mesh level
C
c-        IF (ITIMED .EQ. IUNSTD) THEN
c-          NQNPTS = (IDIM(IBLK,ILVL)+3) *
c-   1                                    (JDIM(IBLK,ILVL)+3) *
c-   2                                    (KDIM(IBLK,ILVL)+3)
c-          IF (ILVL .EQ. LVLONE .AND.
c-   1                               ISTAGE .EQ. 1) THEN
c-            IFL2DQ = .TRUE.
c-          ELSE
c-            IFL2DQ = .FALSE.
c-          ENDIF
c-        ELSE
          NQNPTS = NXQNPT
          IFL2DQ = .FALSE.
c-        ENDIF
C
C     Storage for PROPS
C
          IF (ITURB .GE. ITLMNR) THEN
            NPPTS = (IDIM(IBLK,ILVL)+3) * 
     1              (JDIM(IBLK,ILVL)+3) *
     2              (KDIM(IBLK,ILVL)+3)
          ELSE
            NPPTS = NXPPTS
          ENDIF
C
C     Set up switch from first order scheme to higher order
C
          ITORDR(1) = IORDER(1,IBLK,ILVL,ICFS)
          ITORDR(2) = IORDER(2,IBLK,ILVL,ICFS)
          ITORDR(3) = IORDER(3,IBLK,ILVL,ICFS)
          IF (ITER .LE. NITFO) THEN
            ITORDR(1) = 1
            ITORDR(2) = 1
            ITORDR(3) = 1
          ENDIF
C
C     Set up use of algebraic turbulent eddy viscosity with higher order
C     models for the initial transient
C
c-        IF (ITER .LE. NITALG) THEN
c-          MUTALG = .TRUE.
c-        ELSE
          MUTALG = .FALSE.
c-        ENDIF
C
C     If running a viscous case, calculate the PROPS array
C     Use QP, QM, F and DQP as temporary storage.
C
          IF (ITURB  .GE. ITLMNR ) THEN
            CALL PRPRTY (IDIM(IBLK,ILVL), 
     1                              JDIM(IBLK,ILVL), KDIM(IBLK,ILVL), 
     2                              1, IDIM(IBLK,ILVL)+1, NXSECT,
     3                              Q(IQ), PROPS(IP), S(IS), VOL(IV),
     4                              IBCDAT(1,1,IBLK,ILVL), NBCS(IBLK), 
     5                              IMODEL, IFDELQ,
     6                              STEMP, QP, F, FM, QM, DQP,DQM,WORKX)
          ENDIF
C
C     Update Boundary Conditions
C
          CALL BC(IDIM(IBLK,ILVL),JDIM(IBLK,ILVL),
     1            KDIM(IBLK,ILVL), 1, 
     2            IDIM(IBLK,ILVL)+1,Q(IQ),PROPS(IP),
     3            NPPTS, S(IS), VOL(IV), R(IR),
     4            NBCS(IBLK),IBCDAT(1,1,IBLK,ILVL), 
     5            NXPROF,NPROFL,IOFQBC(1,ILVL),QBC,
     6            VISCOS(1,IBLK),IMODEL,QP, NXSECT,
     7            CL, RCNTR)
C
C      Calculate the Residual for the current block and level
C
          INITRS = .TRUE.
          IPSTRT = 2
          IPEND  = IDIM(IBLK,ILVL)
          CALL RESID  (IDIM(IBLK,ILVL), 
     1                 JDIM(IBLK,ILVL),KDIM(IBLK,ILVL), 
     2                 IPSTRT, IPEND, NXSECT,
     3                 Q(IQ), QN(IQN), NQNPTS,
     4                 PROPS(IP), NPPTS, S(IS), 
     5                 VOL(IV), R(IR), DTJ(ID), RES,
     6                 FHAT, DQP, DQM, QP, QM, TAU, F, FM, 
     7                 STEMP, WORKX,
     8                 IBCDAT(1,1,IBLK,ILVL), NBCS(IBLK), 
     9                 ITORDR,
     A                 LIMIT(1,IBLK), RKAP(1,IBLK), 
     B                 IFLUX(IBLK), VISCOS(1,IBLK), IMODEL,
     C                 INITRS, IFDELQ, MRCHNG, VGNRN,SIGMA,
     D                 NPRLIM, NPRNEG, IFL2DQ, R2DQDT)
          CALL CPRESD(EPS, IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                KDIM(IBLK,ILVL), RES, RESD(IQ))
  500   CONTINUE
C
C     Calculate forces: update on each subiteration for 2D farfield vortex bc
C
c-      IF (IFFORC .EQ. IYES) THEN
c-        CALL FRCINT (ITURB,NBLKS,IDIM(1,ILVL),
c-   1                               JDIM(1,ILVL),  KDIM(1,ILVL),
c-   2                               IOFF(1,ILVL),  IOFFS(1,ILVL),
c-   3                               IOFFQ(1,ILVL), IOFFP(1,ILVL),
c-   4                               Q, PROPS, S, R, NXBCS, NBCS,
c-   5                               IBCDAT(1,1,1,ILVL), FORCE, CL, CD)
c-      ENDIF
  900 CONTINUE

      RETURN
      END
C
C     MUST BE MODIFIED FOR VISCOUS CASE
C
      SUBROUTINE INTQD(EPS, IDIM, JDIM, KDIM, Q, QD)
      include '../../src-flo/header/common.h'
      DIMENSION Q  (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QD (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)

      DIMENSION QC(NQ)

      QD(:,:,:,:) = 0.0

      DO J=2,JDIM
         DO K=2,KDIM
            DO I=2,IDIM
C           Convert from primitive to conserved
               QC(1) = Q(J,K,I,1)
               QC(2) = Q(J,K,I,1)*Q(J,K,I,2)
               QC(3) = Q(J,K,I,1)*Q(J,K,I,3)
               QC(4) = Q(J,K,I,1)*Q(J,K,I,4)
               QC(5) = Q(J,K,I,5)/GAMM1 + 
     1                 0.5*Q(J,K,I,1)*( Q(J,K,I,2)**2 +
     2                                  Q(J,K,I,3)**2 + 
     3                                  Q(J,K,I,4)**2 )
C           Perturb QC
               DO L=1,NQ
                  QD(J,K,I,L) = 2.0*rand() - 1.0
                  QC(L)       = QC(L) + EPS*QD(J,K,I,L)
c$$$                  Q(J,K,I,L) = Q(J,K,I,L) + EPS*QD(J,K,I,L)
               ENDDO
C           Convert from conserved to primitive
               Q(J,K,I,1) = QC(1)
               Q(J,K,I,2) = QC(2)/QC(1)
               Q(J,K,I,3) = QC(3)/QC(1)
               Q(J,K,I,4) = QC(4)/QC(1)
               Q(J,K,I,5) = GAMM1*( QC(5) - 0.5*(QC(2)**2 + QC(3)**2 +
     1                               QC(4)**2)/QC(1) )
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
C
C     Computes RESD FROM RES
C
      SUBROUTINE CPRESD(EPS, IDIM, JDIM, KDIM, RES, RESD)
      include '../../src-flo/header/common.h'
      DIMENSION RES(JDIM+1,KDIM+1,IDIM+1,NF),
     1          RESD(0:JDIM+2,0:KDIM+2,0:IDIM+2,NF)

      IF(EPS.EQ.0.0)THEN

        RESD(:,:,:,:) = 0.0
        DO L=1,NF
          DO J=2,JDIM
            DO K=2,KDIM
              DO I=2,IDIM
                RESD(J,K,I,L) = -RES(J,K,I,L)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ELSE

        DO L=1,NF
          DO J=2,JDIM
            DO K=2,KDIM
              DO I=2,IDIM
                RESD(J,K,I,L) = (RESD(J,K,I,L) + RES(J,K,I,L))/EPS
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDIF

      RETURN
      END
