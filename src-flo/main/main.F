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
      PROGRAM ISAAC
C
C     Program: ISAAC 
C              (Integrated Solution Algorithm for Arbitrary Configurations)
C     Author:  Joseph H. Morrison
C
C     Original code developed under contract to the NASA Langley Research Center
C     while the author worked for:
C              Analytical Services and Materials, Inc.
C              107 Research Drive
C              Hampton, Virginia 23666
C
C
C     Main routine for ISAAC
C
C     Commons
C
      include '../header/maxdim.h'
      include '../header/common.h'
      include '../header/histry.h'
C
C     PERMANENT STORAGE FOR ALL ZONES/BLOCKS
C
C     Storage for flow variables, turbulent eddy viscosity,
C     metrics, volumes, jacobians and grid.
C
      REAL,DIMENSION(:),ALLOCATABLE:: Q, QN, PROPS, S, VOL, DTJ, R, QBC
C
C     Dimensions of domains
C
      DIMENSION IDIM  (MXBLKS,MXLVLS), 
     1          JDIM  (MXBLKS,MXLVLS),
     2          KDIM  (MXBLKS,MXLVLS)
C
C     Extra dimensions for input of grids
C
      DIMENSION IDIMIN(MXBLKS), 
     1          JDIMIN(MXBLKS),
     2          KDIMIN(MXBLKS)
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
      CHARACTER*20 CUTNAM(MXCUTS)
C
C     Printout control data.
C     IPRINT, NPRINT control output at end of iteration cycle
C     IPRMON, NPRMON control output during iteration cycle
C
      DIMENSION IPRINT(10,MXPRT,MXBLKS), NPRINT(MXBLKS),
     1          IPRMON(10,MXPRT,MXBLKS), NPRMON(MXBLKS)
C
C     Data controlling PLOT3D output
C
      LOGICAL IFPLT3, BINPLT
C
C     Data controlling grid type
C
      LOGICAL GRDBIN
C
C     Data controlling the time integration 
C
      DIMENSION NITS  (MXCFS), MGLVLS(MXCFS), ITSLVL(MXLVLS,MXCFS)
      DIMENSION RKALPH(MXSTAG)
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
C     Storage for L2 Norms 
C     For Perfect Gas:
C          R2NORM(1)   = L2 Norm of the mass conservation equation
C          R2NORM(2-4) = L2 Norm of the momentum equations
C          R2NORM(5)   = L2 Norm of the energy equation
C          R2NORM(6)   = L2 Norm of all equations
C     R2NORM  -> L2 Norm of transport equations (includes dQ/dt for iterative,
C                                                implicit scheme)
C     R2DQDT  -> L2 Norm of dQ/dt terms in iterative, implicit scheme
C
      DIMENSION R2NORM(MAXF+1),
     1          R2DQDT(MAXF+1),
     2          NBRPTS(MXLVLS)
C
C     Storage for forces
C     FORCE   -> Forces integrated in x, y, z directions
C     RCNTR   -> Location to calculate point vortex and moments about
C
      DIMENSION FORCE(3),
     1          RCNTR(3),
     2          FMOMENT(100)
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
      DIMENSION QMIN  (MAXQ)
      CHARACTER FIXTYP*7
C
C     STORAGE REQUIRED FOR RUNGE-KUTTA
C     Required on all zones/blocks for the size of the finest grid
C     N.B. This storage is not needed unless using the RK routine
C
      REAL,DIMENSION(:),ALLOCATABLE:: W
C
C     STORAGE REQUIRED FOR MG FORCING FUNCTION AND STORED RESTRICTED Q
C     Required only for coarse grid levels
C
      REAL,DIMENSION(:),ALLOCATABLE:: QCRS, RESCRS
C
C     STORAGE REQUIRED FOR A ZONE/BLOCK
C
C     Storage for residuals
C          For (RK, AF)     : store residual over entire block
C          For (LGS, MARCH) : store residuals over a plane of block
C
      REAL,DIMENSION(:),ALLOCATABLE:: RES
C
C     Storage for fluxes for fourth order
C
      REAL,DIMENSION(:),ALLOCATABLE:: FHAT
C
C     Storage for MUSCL and flux evaluation 
C
      REAL,DIMENSION(:),ALLOCATABLE:: DQP, DQM, QP, QM, TAU, F, FM,
     1                                STEMP, WORKX
C
C     Temporary storage for implicit terms
C          For (RK)         : not needed
C          For (AF3F)       : NPL * max(IDIM*JDIM, JDIM*KDIM)
C                             where NPL is the number of secondary
C                             planes to vectorize over
C          For (MARCH)      : ABC: JDIM*KDIM*2 (store LU of the current plane)
C                             RTEMP,AT: JDIM*KDIM
C
      REAL,DIMENSION(:),ALLOCATABLE:: ABC, RTEMP, AT
C
C     Temporary storage for the I contribution to the implicit scheme 
C     for LGS and MARCH.  
C
      REAL,DIMENSION(:),ALLOCATABLE:: RM
C
C     Identity Matrix required in implicit work
C
      REAL,DIMENSION(:),ALLOCATABLE:: AI
C
C     Filenames
C     FRDRST  File to read restart from
C     FWRRST  File to write restart to
C     GRDFIL  Array of files to read grid from
C     FPLT3G  File to write PLOT3D Grid file to
C     FPLT3Q  File to write PLOT3D Q (conserved variables) to
C     FPLT3F  File to write PLOT3D function file to
C     FPLT3FN File to write PLOT3D function name file to
C     FPROFL  Array of files to read Profile boundary data from
C     FERROR  File name output if error
C
      CHARACTER*80 FRDRST, FWRRST, GRDFIL(MXBLKS),
     1             FPLT3G, FPLT3Q, FPLT3F, FPLT3FN, FPROFL(MXPROF),
     2             FERROR
C
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
C

#ifdef CPU_TIME
C
C     Storage to calculate execution time
C
      REAL*4  TARRAY(2), RUNTIME
#endif

C
C     Version Number for Output
C
      CHARACTER*48 VERSN
C
C     Local storage to output number of grid points 
C
      DIMENSION NNODES(MXLVLS,2)
C
C 0.  PRELIMINARIES
C     Output Version Number
C
C     VERSN  = '$Revision: xxx$'
C     WRITE (IOUT,1200) VERSN
C1200 FORMAT (' ','NUWTUN Flow Solver - ',A48,/)
C
C     Set switch to test for storage errors
C
      IERROR = 0
C
C     Set storage for IMODEL to maximum
C
      NMDL   = MAXMDL
C
C     Set storage for WORKX to maximum
C
      NW     = MAXW
C
C     Set up the Characteristic Array to define the characteristics 
C     of the available flux evaluation schemes.
C
      IFCHAR(IFROE,1) = IFDS
C
C     Zero appropriate storage
C
      ITTOT = 0
C
C     Initialize time to zero
C
      TIME  = 0.E0
C
C     Initialize force coefficients
C
      CL    = 0.E0
      CD    = 0.E0
C
C 1.  READ INPUT
C     Read input data controlling the run
C
      NXQ    = MAXQ
      NXCFS  = MXCFS
      NXLVLS = MXLVLS
      NXBLKS = MXBLKS
      NXBCS  = MXBCS
      NXCUTS = MXCUTS
      NXPRT  = MXPRT
      NXSTAG = MXSTAG
      NXPROF = MXPROF
      NXTRSG = MXTRSG
      CALL RDDATA (NXQ, NXCFS, NXLVLS, NXBLKS, NXBCS, NXCUTS, NXPRT,
     1             NCFS, MGLVLS, NBLKS, IDIM, JDIM, KDIM,
     2             NBCS, IBCDAT, NCUTS, ICUTS, CUTNAM,
     3             NPRINT, IPRINT, NPRMON, IPRMON, NITMON,
     4             ITUPJ, ITSLVL, NITS, NSUBIT, NITFO, NITBBC, NITALG,
     5             NITRST, NXSTAG, NSTAGE, RKALPH,
     6             RESTRT, FRDRST, FWRRST, GRDFIL, IGRDTP, GRDBIN,
     7             NXPROF, NPROFL, FPROFL, 
     8             IFPLT3, BINPLT, FPLT3G, FPLT3Q, FPLT3F, FPLT3FN,
     9             IORDER, LIMIT, RKAP, IFLUX, VISCOS, DIAG, IMODEL,
     A             IMSTRT, IMEND, TOLER, VGNRN, SIGMA, IFDELQ,
     B             CONTRB, IMODLX, QMIN, ITRSMX, IFFORC, IFRCPR, 
     C             YNCALC, NXTRSG, NTRSEG, ITRSEG, RCNTR)
C
      ITURB = IMODEL( 1)
      IASM  = IMODEL(10)
C
C 2.  VERIFY INPUT
C
C     Check input data for BC, CUTS, Print statements
C
      CALL CHKDAT (NBLKS, IDIM, JDIM, KDIM, VISCOS,
     1             NXBCS, NBCS, IBCDAT, NCUTS, ICUTS, CUTNAM,
     2             NXPRT, NPRINT, IPRINT, IERRCD)
      IF (IERRCD .GT. INO) THEN
         IERROR = IERROR + 1
      ENDIF
C
C     Check to verify that NQ >= NP.  This is required due to temporary
C     arrays in DELQ routines.  If NP is ever greater than NQ, need to
C     change the temporary storage.
C
      IF (NP .GT. NQ) THEN
         WRITE (IOUT,8000) 
         IERROR = IERROR + 1
      ENDIF
 8000 FORMAT (' ','MAIN  : ERROR-> NP is greater than NQ.',
     1       /' ',15X,'This is not allowed due to temporary storage',
     2       /' ',15X,'allocation in DELQ routines.',
     3       /' ',15X,'The temporary storage for PROPSI and PROPSC ',
     4       /' ',15X,'must be modified.')
C
C     Check the values of constants to make sure that the storage exists
C
      IF (NQ     .GT. MAXQ  ) THEN
         WRITE (IOUT,1505) NQ, MAXQ
         IERROR = IERROR + 1
      ENDIF
C
      IF (NF     .GT. MAXF  ) THEN
         WRITE (IOUT,1510) NF, MAXF
         IERROR = IERROR + 1
      ENDIF
C
      IF (NP     .GT. MAXP  ) THEN
         WRITE (IOUT,1515) NP, MAXP
         IERROR = IERROR + 1
      ENDIF
C
      IF (NRANK  .GT. MXRANK) THEN
         WRITE (IOUT,1520) NRANK, MXRANK
         IERROR = IERROR + 1
      ENDIF
C
      IF (NTMLVL .GT. MXTMLV) THEN
         WRITE (IOUT,1525) NTMLVL, MXTMLV
         IERROR = IERROR + 1
      ENDIF
C
 1505 FORMAT (/' ','MAIN  : Number of variables required  ',I5,
     1             ' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MAXQ   and re-submit.')
 1510 FORMAT (/' ','MAIN  : Number of fluxes required     ',I5,
     1             ' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MAXF   and re-submit.')
 1515 FORMAT (/' ','MAIN  : Number of properties required ',I5,
     1             ' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MAXP   and re-submit.')
 1520 FORMAT (/' ','MAIN  : Size of matrices required     ',I5,
     1             ' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MXRANK and re-submit.')
 1525 FORMAT (/' ','MAIN  : Number of time levels required',I5,
     1             ' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MXTMLV and re-submit.')
C
C     Set up the identity matrix used in implicit terms
C
      ALLOCATE( AI (NRANK*NRANK) )
      CALL IDENT (AI)
C
C 3.  CALCULATE POINTERS
C     a. Calculate dimensions for coarser grids and offset pointers.
C     b. Calculate offset pointers and boundary data.
C     c. Calculate ICUTS  for coarser grids.
C     d. Verify coarse grid BC and CUT data
C     e. Set up IORDER for multigrid
C     f. Output grid summary
C
      NDSIZE = 0
      NDCRS  = 0
      NDSECT = 0
      NDRM   = 0
      NLVLS  = 1
      DO 10 ICFS = 1, NCFS
         ILVL = ICFS + MGLVLS(ICFS) - 1
         IF (ILVL .GT. NLVLS) NLVLS = ILVL
 10   CONTINUE
      IOFFT  = 0
      IOFFST = 0
      IOFFQT = 0
      IOFFPT = 0
      IOFQNT = 0
      IOFQCT = 0
C
C 3-a. Calculate dimensions for coarser grids and test allowable grid levels.
C
      MSTLVL = NLVLS
      DO 20 ILVL = 2, NLVLS
         DO 15 IBLK = 1, NBLKS
            IDIM(IBLK,ILVL) = IDIM(IBLK,ILVL-1) / 2 + 1
            JDIM(IBLK,ILVL) = JDIM(IBLK,ILVL-1) / 2 + 1
            KDIM(IBLK,ILVL) = KDIM(IBLK,ILVL-1) / 2 + 1
C
            IDIMF           = 2 * (IDIM(IBLK,ILVL) - 1) + 1
            JDIMF           = 2 * (JDIM(IBLK,ILVL) - 1) + 1
            IF (THREED) THEN
               KDIMF        = 2 * (KDIM(IBLK,ILVL) - 1) + 1
            ELSE
               KDIMF        = 2
            ENDIF
C
            IF (IDIMF .NE. IDIM(IBLK,ILVL-1) .OR.
     1          JDIMF .NE. JDIM(IBLK,ILVL-1) .OR.
     2          KDIMF .NE. KDIM(IBLK,ILVL-1) ) THEN
               MSTLVL = MIN (MSTLVL, (ILVL-1))
            ENDIF
   15    CONTINUE
   20 CONTINUE
C
      IF (NLVLS  .GT. MSTLVL) THEN
         WRITE (IOUT,1530) NLVLS, MSTLVL
         NLVLS = MSTLVL
      ENDIF
      IF (NCFS   .GT. NLVLS ) THEN
         WRITE (IOUT,1535) NCFS,  NLVLS
         NCFS  = NLVLS 
      ENDIF
 1530 FORMAT (/' ','MAIN  : Number of grid levels requested ',I3,
     1             ' exceeds allowable ',I3,
     2        /' ', 8X,'=> Coarser levels are being deleted and ',
     3                 'run continuing.')
 1535 FORMAT (/' ','MAIN  : Number of coarsenings requested ',I3,
     1             ' exceeds allowable ',I3,
     2        /' ', 8X,'=> Coarsenings    are being deleted and ',
     3                 'run continuing.'//)
C
C      Verify all coarsenings fall within NLVLS number of grid levels
C
      DO 25 ICFS = 1, NCFS
         LVLTST = ICFS + MGLVLS(ICFS) - 1
         IF (LVLTST .GT. NLVLS) THEN
            MGLVLS(ICFS) = NLVLS - ICFS + 1
         ENDIF
   25 CONTINUE
C
C 3-b. Calculate offset pointers and boundary data.
C
      DO 50 ILVL = 1, NLVLS
         NBRPTS(ILVL) = 0
         DO 40 IBLK = 1, NBLKS
            IF (ILVL .NE. 1) THEN
C
C      Calculate storage required on coarser grids for forcing function
C
               NDCRS = NDCRS + (IDIM(IBLK,ILVL)+3)*(JDIM(IBLK,ILVL)+3)*
     1                         (KDIM(IBLK,ILVL)+3)
C
C      Calculate IBCDAT for coarser grids.
C
               DO 35 IBC = 1, NBCS(IBLK)
                  IBCDAT( 1,IBC,IBLK,ILVL)=IBCDAT(1,IBC,IBLK,ILVL-1)
                  IBCDAT( 2,IBC,IBLK,ILVL)=IBCDAT(2,IBC,IBLK,ILVL-1)
                  IBCDAT( 3,IBC,IBLK,ILVL)=IBCDAT(3,IBC,IBLK,ILVL-1)/2+1
                  IBCDAT( 4,IBC,IBLK,ILVL)=IBCDAT(4,IBC,IBLK,ILVL-1)/2+1
                  IBCDAT( 5,IBC,IBLK,ILVL)=IBCDAT(5,IBC,IBLK,ILVL-1)/2+1
                  IBCDAT( 6,IBC,IBLK,ILVL)=IBCDAT(6,IBC,IBLK,ILVL-1)/2+1
                  IBCDAT( 7,IBC,IBLK,ILVL)=IBCDAT(7,IBC,IBLK,ILVL-1)/2+1
                  IBCDAT( 8,IBC,IBLK,ILVL)=IBCDAT(8,IBC,IBLK,ILVL-1)/2+1
                  IBCDAT( 9,IBC,IBLK,ILVL)=IBCDAT(9,IBC,IBLK,ILVL-1)/2+1
                  IBCDAT(10,IBC,IBLK,ILVL)=IBCDAT(10,IBC,IBLK,ILVL-1)
   35          CONTINUE
            ENDIF
            NBRPTS(ILVL) = NBRPTS(ILVL) + (IDIM(IBLK,ILVL)-1)*
     1                          (JDIM(IBLK,ILVL)-1)*(KDIM(IBLK,ILVL)-1)
            IOFF  (IBLK,ILVL) = IOFFT
            IOFFS (IBLK,ILVL) = IOFFST
            IOFFQ (IBLK,ILVL) = IOFFQT
            IOFFP (IBLK,ILVL) = IOFFPT
            IOFFQN(IBLK,ILVL) = IOFQNT
            IOFFQC(IBLK,ILVL) = IOFQCT
            IOFFT  = IOFFT  + IDIM(IBLK,ILVL)*JDIM(IBLK,ILVL)*
     1                        KDIM(IBLK,ILVL)
            IOFFST = IOFFST + (IDIM(IBLK,ILVL)+2)*(JDIM(IBLK,ILVL)+2)*
     1                        (KDIM(IBLK,ILVL)+2)
            IOFFQT = IOFFQT + (IDIM(IBLK,ILVL)+3)*(JDIM(IBLK,ILVL)+3)*
     1                        (KDIM(IBLK,ILVL)+3)
            IF (ILVL .NE. 1) THEN
               IOFQCT = IOFQCT+(IDIM(IBLK,ILVL)+3)*(JDIM(IBLK,ILVL)+3)*
     1                        (KDIM(IBLK,ILVL)+3)
            ENDIF
C
            NDSIZE = MAX(NDSIZE,(IDIM(IBLK,ILVL)+3)*(JDIM(IBLK,ILVL)+3)*
     1                          (KDIM(IBLK,ILVL)+3))
            NDSECT = MAX(NDSECT,(IDIM(IBLK,ILVL)+3)*(JDIM(IBLK,ILVL)+3))
            NDSECT = MAX(NDSECT,(JDIM(IBLK,ILVL)+3)*(KDIM(IBLK,ILVL)+3))
C
            IF (ISOLVR .EQ. ILGS .OR. ISOLVR .EQ. IMARCH) THEN
               NDRM  = MAX(NDRM,(JDIM(IBLK,ILVL)-1)*(KDIM(IBLK,ILVL)-1))
            ELSE
               NDRM  = 0
            ENDIF
C
C     If running inviscid case (ITURB < ITLMNR), then don't
C     increment the pointer to the PROPS array.  This allows to not
C     store this data for the inviscid case to cut down on memory.
C
            INCP   = (IDIM(IBLK,ILVL)+3)*(JDIM(IBLK,ILVL)+3)*
     1               (KDIM(IBLK,ILVL)+3)
            IF (ITURB .LT. ITLMNR) THEN
               INCP = 0
            ENDIF
            IOFFPT = IOFFPT + INCP
C
C     If running steady case, then don't increment pointer to previous
C     time level data to allow minimum memory configuration for steady cases.
C
            INCQN  = (IDIM(IBLK,ILVL)+3)*(JDIM(IBLK,ILVL)+3)*
     1               (KDIM(IBLK,ILVL)+3)
            IF (ITIMED .EQ. ISTDY) THEN
               INCQN = 0
            ENDIF
            IOFQNT = IOFQNT + INCQN
   40    CONTINUE
C
C 3-c. Calculate ICUTS  for coarser grids.
C
         IF (ILVL .NE. 1) THEN
            DO 45 ICUT = 1, NCUTS
               ICUTS( 1,ICUT,ILVL) = ICUTS( 1,ICUT,ILVL-1)
               ICUTS( 2,ICUT,ILVL) = ICUTS( 2,ICUT,ILVL-1)
               ICUTS( 3,ICUT,ILVL) = ICUTS( 3,ICUT,ILVL-1) / 2 + 1
               ICUTS( 4,ICUT,ILVL) = ICUTS( 4,ICUT,ILVL-1) / 2 + 1
               ICUTS( 5,ICUT,ILVL) = ICUTS( 5,ICUT,ILVL-1) / 2 + 1
               ICUTS( 6,ICUT,ILVL) = ICUTS( 6,ICUT,ILVL-1) / 2 + 1
               ICUTS( 7,ICUT,ILVL) = ICUTS( 7,ICUT,ILVL-1) / 2 + 1
               ICUTS( 8,ICUT,ILVL) = ICUTS( 8,ICUT,ILVL-1) / 2 + 1
               ICUTS( 9,ICUT,ILVL) = ICUTS( 9,ICUT,ILVL-1) / 2 + 1
               ICUTS(10,ICUT,ILVL) = ICUTS(10,ICUT,ILVL-1)
               ICUTS(11,ICUT,ILVL) = ICUTS(11,ICUT,ILVL-1)
               ICUTS(12,ICUT,ILVL) = ICUTS(12,ICUT,ILVL-1) / 2 + 1
               ICUTS(13,ICUT,ILVL) = ICUTS(13,ICUT,ILVL-1) / 2 + 1
               ICUTS(14,ICUT,ILVL) = ICUTS(14,ICUT,ILVL-1) / 2 + 1
               ICUTS(15,ICUT,ILVL) = ICUTS(15,ICUT,ILVL-1) / 2 + 1
               ICUTS(16,ICUT,ILVL) = ICUTS(16,ICUT,ILVL-1) / 2 + 1
               ICUTS(17,ICUT,ILVL) = ICUTS(17,ICUT,ILVL-1) / 2 + 1
               ICUTS(18,ICUT,ILVL) = ICUTS(18,ICUT,ILVL-1) / 2 + 1
               ICUTS(19,ICUT,ILVL) = ICUTS(19,ICUT,ILVL-1)
               ICUTS(20,ICUT,ILVL) = ICUTS(20,ICUT,ILVL-1)
               ICUTS(21,ICUT,ILVL) = ICUTS(21,ICUT,ILVL-1)
   45       CONTINUE
         ENDIF
   50 CONTINUE
C
C 3-d. Verify coarse grid BC and CUT data
C
      IF (NLVLS .GT. 1) THEN
         CALL CHKCRS (NLVLS,  NXBLKS, NBLKS,  IDIM, JDIM, KDIM, 
     1                NXBCS,  NBCS,   IBCDAT, 
     2                NXCUTS, NCUTS,  ICUTS,  CUTNAM, IERRCD)
         IF (IERRCD .GT. INO) THEN
            IERROR = IERROR + 1
         ENDIF
      ENDIF
C
C 3-e. Set up IORDER for multigrid
C
      IHIGH = 2
      ILOW  = 1
      IF (FOURTH) THEN
         IHIGH = 4
         ILOW  = 4
      ENDIF
C
      DO IBLK = 1, NBLKS
         DO ICFS = 1, NCFS
            LVLONE = NCFS - ICFS + 1
            DO ILVL = 1, NLVLS
               IF (ILVL .EQ. LVLONE) THEN
                  IORDER(1,IBLK,ILVL,ICFS) = IHIGH
                  IORDER(2,IBLK,ILVL,ICFS) = IHIGH
                  IORDER(3,IBLK,ILVL,ICFS) = IHIGH
               ELSE
                  IORDER(1,IBLK,ILVL,ICFS) = ILOW
                  IORDER(2,IBLK,ILVL,ICFS) = ILOW
                  IORDER(3,IBLK,ILVL,ICFS) = ILOW
               ENDIF
            END DO
         END DO
      END DO
C
C 3-f. Output grid summary
C
      WRITE (IOUT,1540) (ILVL, ILVL = 1, NLVLS)
      DO ILVL = 1, NLVLS
         NNODES(ILVL,2) = 0
      END DO
C
      DO IBLK = 1, NBLKS
         DO ILVL = 1, NLVLS
            IF (THREED) THEN
               K = KDIM(IBLK,ILVL)
            ELSE
               K = 1
            ENDIF
            NNODES(ILVL,1) = IDIM(IBLK,ILVL) * JDIM(IBLK,ILVL) * K
            NNODES(ILVL,2) = NNODES(ILVL,2) + NNODES(ILVL,1)
         END DO
         WRITE (IOUT,1541) IBLK,IDIM(IBLK,1),JDIM(IBLK,1),KDIM(IBLK,1),
     1                     (NNODES(ILVL,1), ILVL = 1, NLVLS)
      END DO
      WRITE (IOUT,1542) (NNODES(ILVL,2), ILVL = 1, NLVLS)
      WRITE (IOUT,1543) 
C
 1540 FORMAT (///' ','Grid Summary',
     1          /' ',27X,'Number of Grid Points by Mesh Level',
     2          /' ',11X,'IDIM',1X,'JDIM',1X,'KDIM',6(5x,I2,2X))
 1541 FORMAT (   ' ','Block ',I3,':',3I5,6(1X,I8))
 1542 FORMAT (  /' ','Total ',19X,6(1X,I8))
 1543 FORMAT (///' ')
C
C 4.  Allocate memory
C (a) MXBCPT should be replaced with NDBCPT
C (b) According to comments it should be PROPS(IOFFPT*NP) but this does
C     not work, segmentation fault in subroutine YNABS
C
      ALLOCATE( Q      (IOFFQT*NQ)        )
      ALLOCATE( QN     (IOFQNT*NQ*NTMLVL) )
      ALLOCATE( PROPS  (IOFFQT*NP)        )
      ALLOCATE( S      (IOFFQT*4*3)       )
      ALLOCATE( VOL    (IOFFQT)           )
      ALLOCATE( DTJ    (IOFFQT)           )
      ALLOCATE( R      (IOFFQT*3)         )
      ALLOCATE( QBC    (MXBCPT*NQ)        )
      ALLOCATE( W      (IOFFQT*NQ)        )
      ALLOCATE( QCRS   (NDCRS*NQ)         )
      ALLOCATE( RESCRS (NDCRS*NQ)         )
      ALLOCATE( RES    (NDSIZE*NF)        )
      ALLOCATE( FHAT   (NDSIZE*NF)        )
      ALLOCATE( DQP    (NDSECT*NQ)        )
      ALLOCATE( DQM    (NDSECT*NQ)        )
      ALLOCATE( QP     (NDSECT*NQ)        )
      ALLOCATE( QM     (NDSECT*NQ)        )
      ALLOCATE( TAU    (NDSECT*6)         )
      ALLOCATE( F      (NDSECT*NF)        )
      ALLOCATE( FM     (NDSECT*NF)        )
      ALLOCATE( STEMP  (NDSECT*4)         )
      ALLOCATE( WORKX  (NDSECT*NW)        )
      ALLOCATE( ABC    (NDSECT*NQ*NQ*3)   )
      ALLOCATE( RTEMP  (NDSECT*NQ)        )
      ALLOCATE( AT     (NDSECT*NQ*NQ*3)   )
      ALLOCATE( RM     (NDRM*NQ*NQ)       )
C
C     Used by subroutine WRREST and others
C
      NXSECT = NDSECT
      NXPPTS = IOFFQT
      NXQNPT = IOFQNT
      NXNODE = IOFFQT
      NXABC  = NDSECT
C
C     NOTE: WE MAY NEED ABC(2*NDSECT*NQ*NQ*3) and NXABC = 2*NDSECT
C           MUST CLARIFY THIS POINT.
C
C 5.  INPUT GRID
C     a. Read in grid for all blocks.  
C     b. Loop over cuts and check the validity of the grid cuts
C     c. Extend each block grid from 1:DIM,... to 0:DIM+1,... and calculate
C        coarser grids as every other grid point from next finer grid.
C     d. Loop over cuts on all levels blocks to set image points of the grid.
C     e. Calculate metrics and volumes, reset volumes to be non-zero for
C        boundary conditions that do not lie on a face
C
C 5-a. Read in grid for all blocks, finest grid level.
C
      CALL RDGRID (NBLKS, IDIM, JDIM, KDIM, IOFFS, GRDFIL, R,
     1             IDIMIN, JDIMIN, KDIMIN, IGRDTP, GRDBIN, IEOF, IERRCD)
      IF (IEOF .NE. INO .OR. IERRCD .NE. INO) THEN
         STOP
      ENDIF
C
C 5-b. Loop over cuts and check the validity of the grid cuts
C
      IERRCD = INO
      ILVL   = 1
      DO 55 ICUT = 1, NCUTS
         IBLK1 = ICUTS( 1,ICUT,ILVL)
         IBLK2 = ICUTS(10,ICUT,ILVL)
         IR1   = IOFFS(IBLK1,ILVL) * 3    + 1
         IR2   = IOFFS(IBLK2,ILVL) * 3    + 1
         CALL CHKCUT (IBLK1, IDIM(IBLK1,ILVL), JDIM(IBLK1,ILVL), 
     1                KDIM(IBLK1,ILVL), R(IR1),
     2                IBLK2, IDIM(IBLK2,ILVL), JDIM(IBLK2,ILVL), 
     3                KDIM(IBLK2,ILVL), R(IR2), ICUTS(1,ICUT,ILVL),
     4                CUTNAM(ICUT), IERRCD)
   55 CONTINUE
      IF (IERRCD .NE. INO) THEN
         STOP
      ENDIF
C
C 5-c. Extend each block grid and calculate coarser meshes.
C
      DO 65 IBLK = 1, NBLKS
         ILVL = 1
         IR   = IOFFS(IBLK,ILVL) * 3     + 1
         CALL GDXTND (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                KDIM(IBLK,ILVL), R(IR))
C
         DO 60 ILVL = 2, NLVLS
            IRF  = IOFFS(IBLK,ILVL-1) * 3   + 1
            IR   = IOFFS(IBLK,ILVL  ) * 3   + 1
            CALL CRSGRD (IDIM(IBLK,ILVL-1), JDIM(IBLK,ILVL-1),
     1                   KDIM(IBLK,ILVL-1), R(IRF),
     2                   IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     3                   KDIM(IBLK,ILVL), R(IR))
            CALL GDXTND (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                   KDIM(IBLK,ILVL), R(IR))
 60      CONTINUE
 65   CONTINUE
C
C 5-d. Loop over cuts setting grid image points
C
      DO 75 ILVL = 1, NLVLS
         DO 70 ICUT = 1, NCUTS
            IBLK1 = ICUTS( 1,ICUT,ILVL)
            IBLK2 = ICUTS(10,ICUT,ILVL)
            IR1   = IOFFS(IBLK1,ILVL) * 3    + 1
            IR2   = IOFFS(IBLK2,ILVL) * 3    + 1
            CALL GRDCUT (IDIM(IBLK1,ILVL), JDIM(IBLK1,ILVL), 
     1                   KDIM(IBLK1,ILVL), R(IR1),
     2                   IDIM(IBLK2,ILVL), JDIM(IBLK2,ILVL), 
     3                   KDIM(IBLK2,ILVL), R(IR2), ICUTS(1,ICUT,ILVL))
   70    CONTINUE
   75 CONTINUE
C
C 5-e. Calculate metrics and volumes and reset volumes for some bcs.
C      Reset interior volumes on coarser meshes to be the collection of the
C      fine grid volumes to ensure that the coarse grid transfer operator
C      generates a smooth and accurate coarse grid representation of the fine
C      grid solution.
C
      IERRCD = INO
      DO 85 IBLK = 1, NBLKS
         DO 80 ILVL = 1, NLVLS
            IR   = IOFFS(IBLK,ILVL) * 3     + 1
            IS   = IOFFS(IBLK,ILVL) * 4 * 3 + 1
            IV   = IOFFS(IBLK,ILVL)         + 1
            CALL METRIC (IBLK, ILVL, IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                   KDIM(IBLK,ILVL), R(IR), S(IS), VOL(IV), IERRCD)
            CALL METBC  (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                   KDIM(IBLK,ILVL), S(IS), VOL(IV), 
     2                   NBCS(IBLK), IBCDAT(1,1,IBLK,ILVL))
            IF (ILVL .GT. 1) THEN
               IVF  = IOFFS(IBLK,ILVL-1)    + 1
               CALL COLV2H (IDIM(IBLK,ILVL-1), JDIM(IBLK,ILVL-1),
     1                      KDIM(IBLK,ILVL-1), VOL(IVF),
     2                      IDIM(IBLK,ILVL),   JDIM(IBLK,ILVL),  
     3                      KDIM(IBLK,ILVL),   VOL(IV),  QP)
            ENDIF
   80    CONTINUE
   85 CONTINUE
      IF (IERRCD .NE. INO) THEN
         STOP
      ENDIF
C
C 6.  INITIALIZATION
C     Initialize each block and each grid level to freestream.
C     Initialize each block and each grid level properties.
C     Initialize each block and each grid level at additional time levels.
C     Initialize each block and each grid level wall functions.
C
      DO 100 ILVL = 1, NLVLS
         DO 90 IBLK = 1, NBLKS
            IQ  = IOFFQ(IBLK,ILVL) * NQ    + 1
            IV  = IOFFS(IBLK,ILVL)         + 1
            CALL INIT   (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                   KDIM(IBLK,ILVL), Q(IQ), ITURB)
C
            IF (IASM   .EQ. IASMGS) THEN
               IP  = IOFFP(IBLK,ILVL) * NP    + 1
               CALL INITPR (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                      KDIM(IBLK,ILVL), PROPS(IP), ITURB)
            ENDIF
C
            IF (ITIMED .EQ. IUNSTD) THEN
               IQN    = IOFFQN(IBLK,ILVL) * NQ * NTMLVL + 1
               CALL INITQN (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                      KDIM(IBLK,ILVL), Q(IQ), QN(IQN), 
     2                      QP, QM, IMODEL)
            ENDIF
C
            IF (IFWALF) THEN
               IP  = IOFFP(IBLK,ILVL) * NP    + 1
               CALL INITTQ (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                      KDIM(IBLK,ILVL), PROPS(IP), 
     2                      NBCS(IBLK), IBCDAT(1,1,IBLK,ILVL))
            ENDIF
  90     CONTINUE
 100  CONTINUE
C
C     Loop over all boundary conditions to check for a PROFILE or 
C     SUBSONIC PROFILE bc, then read in the boundary data for the 
C     PROFILE or SUBSONIC PROFILE bcs on the fine grid
C     and calculate the coarsened profile data for coarse grids.
C     Set all offsets to zero in case there are no profile boundary segments.
C     i.  Re-order profile BC segment data based on profile number
C     ii. Loop over all profile segments in order and calculate offsets,
C         read in profile data and coarsen data
C
C     i.  Re-order profile BC segment data based on profile number
C
      DO 115 IBLK = 1, NBLKS
         DO 110 IBC = 1, NBCS(IBLK)
            IBCTYP = IBCDAT(1,IBC,IBLK,1)
            IF (IBCTYP .EQ. IPROFL .OR. IBCTYP .EQ. IPROSB) THEN
               INPROF = IBCDAT(10,IBC,IBLK,1)
               DO 105 ILVL = 1, NLVLS
                  DO 104 L = 1, 10
                     IBCPRF(L,INPROF,ILVL) = IBCDAT(L,IBC,IBLK,ILVL)
  104             CONTINUE
  105          CONTINUE
            ENDIF
  110    CONTINUE
  115 CONTINUE
C
C     ii. Loop over all profile segments in order and calculate offsets,
C         read in profile data and coarsen data
C
      NDBCPT = 0
      IERRCD = INO
      IERRC2 = INO
      DO 125 ILVL = 1, NLVLS
         IOFQBC(1,ILVL) = 0
         DO 120 INPROF = 1, NPROFL
            ISTRT  = IBCPRF( 4,INPROF,ILVL)
            IEND   = IBCPRF( 5,INPROF,ILVL)
            JSTRT  = IBCPRF( 6,INPROF,ILVL)
            JEND   = IBCPRF( 7,INPROF,ILVL)
            KSTRT  = IBCPRF( 8,INPROF,ILVL)
            KEND   = IBCPRF( 9,INPROF,ILVL)
C
            IOFQBC(INPROF,ILVL) = NDBCPT
            IQBC   = IOFQBC(INPROF,ILVL) * NQ + 1
C
            NPTS   = (IEND-ISTRT+1) * (JEND-JSTRT+1) * 
     1               (KEND-KSTRT+1)
            NDBCPT = NDBCPT + NPTS
C
            IF (NDBCPT .GT. MXBCPT) THEN
               IERRCD = IYES
            ENDIF
C
            IF (IERRCD .EQ. INO) THEN
               IF (ILVL .EQ. 1) THEN
                  WRITE (IOUT,1410) FPROFL(INPROF)
                  FERROR = FPROFL(INPROF)
                  OPEN (IRDBC,FILE=FPROFL(INPROF),STATUS='OLD',
     1                        FORM='FORMATTED',ERR=990)
                  CALL RDPRFL (NPTS, QBC(IQBC),
     2                         IBCPRF(1,INPROF,ILVL), IERRC2)
               ELSE
                  ISTRTF  = IBCPRF( 4,INPROF,ILVL-1)
                  IENDF   = IBCPRF( 5,INPROF,ILVL-1)
                  JSTRTF  = IBCPRF( 6,INPROF,ILVL-1)
                  JENDF   = IBCPRF( 7,INPROF,ILVL-1)
                  KSTRTF  = IBCPRF( 8,INPROF,ILVL-1)
                  KENDF   = IBCPRF( 9,INPROF,ILVL-1)
                  NPTSF   = (IENDF-ISTRTF+1) * (JENDF-JSTRTF+1) * 
     1                      (KENDF-KSTRTF+1)
                  IQBCF   = IOFQBC(INPROF,ILVL-1) * NQ + 1
                  CALL CRSPRF (NPTSF, QBC(IQBCF),
     1                         IBCPRF(1,INPROF,ILVL-1), 
     2                         NPTS,  QBC(IQBC),
     3                         IBCPRF(1,INPROF,ILVL), IERRCD)
               ENDIF
            ENDIF
  120    CONTINUE
  125 CONTINUE
C
      IF (IERRCD .NE. INO) THEN
         WRITE (IOUT,1420) NDBCPT, MXBCPT
      ENDIF
C
      IF (IERRCD .NE. INO .OR. IERRC2 .NE. INO) THEN
         STOP
      ENDIF
C
 1410 FORMAT (//' ','MAIN  : Reading input profile from the file:',
     1         /' ',10X,A80)
 1420 FORMAT (' ','MAIN  : ERROR-> Number of boundary points  ',I10,
     1       /' ','                exceeds allocated (MXBCPT) ',I10,
     2       /' ','        Increase MXBCPT in MAIN and resubmit.',
     3       /' ', 8X,'RUN ABORTING!')
C
C 7.  RESTART
C     If this is a restart solution, then read in restart data.
C     The logic here is a little confusing so pay attention:
C     (The restart/conversion process plays with NQ, NP.  It must reset NQ, NP
C      to be the correct value before leaving this section)
C
C     Begin{Restart/Conversion}
C        Open Files
C        Set NQTO,NFTO, NPTO         <- the number of variables of TO   
C                                       turbulence model
C                                       (this is NQ,NP as set in RDDATA)
C        Set NQFROM,NFFROM,NPFROM    <- the number of variables of FROM 
C                                       turbulence model
C        Set NQ <- NQFROM, NF <- NFFROM, NP <- NPFROM
C        For Each Block
C           If NORMAL Restart ( .not. contrb ) 
C              begin
C                 1. Read restart data into Q, PROPS
C                 2. If unsteady, convert QN to conserved variables
C                    N.B. storage for tke in convert is over entire block!
C              end
C           Else (doing a restart with turbulence model conversion)
C              begin
C                 1. Read restart data of NQFROM variables  into RES, 
C                                         NPFROM properties into W
C                    (RES is dimensioned  to hold 1 block of Q,
C                     W   is large enough to hold 1 block of PROPS)
C                 2. Call conversion routine to convert from one turbulence 
C                    model to another.  Returns the new model in Q, PROPS.
C                    Note: The PROPS will be updated after this section.
C                 3. If unsteady, output error as this is not yet supported
C              end
C           Endif
C        Next Block
C        Set NQ <- NQTO, NF <- NFTO, NP <- NPTO
C        Close Files
C     End{Restart/Conversion}
C
      IF ( RESTRT ) THEN
         WRITE (IOUT,1500)
 1500 FORMAT (' ','MAIN  : Reading Restart File'/)
         FERROR = FRDRST
         OPEN (IRDRST,FILE=FRDRST,STATUS='OLD',FORM='UNFORMATTED',
     1                ERR=990)
         REWIND (IRDRST)
C
         NQTO   = NQ
         NFTO   = NF
         NPTO   = NP
C
         IF ( .NOT. CONTRB ) THEN
            NQFROM = NQ
            NFFROM = NF
            NPFROM = NP
         ELSE
            ITURBX = IMODLX( 1)
            IF (ITURBX .EQ. ITLMNR .OR. ITURBX .EQ. ITBLMX) THEN
               NQFROM =  5
               NFFROM =  5
               NPFROM =  5
            ELSE IF (ITURBX .EQ. ITKE .OR. ITURBX .EQ. ITKW) THEN
               NQFROM =  7
               NFFROM =  7
               NPFROM =  5
            ELSE IF (ITURBX .EQ. ITRS) THEN
               NQFROM = 12
               NFFROM = 12
               NPFROM =  5
            ELSE
               WRITE (IOUT,1600) ITURBX
 1600 FORMAT (' ','MAIN  : ERROR-> Invalid FROM turbulence model = ',I5,
     1       /' ', 8X,'RUN ABORTING!')
               STOP
            ENDIF
         ENDIF
C
C     Now reset NQ, NF, NP to be NQFROM, NFFROM, NPFROM
C
         NQ     = NQFROM
         NF     = NFFROM
         NP     = NPFROM
C
C     Loop over all grid levels that are written in the restart file
C
         ILVL = 1
  130    CONTINUE
         DO 135 IBLK = 1, NBLKS
            IQ  = IOFFQ (IBLK,ILVL) * NQTO          + 1
            IQN = IOFFQN(IBLK,ILVL) * NQTO * NTMLVL + 1
            IP  = IOFFP (IBLK,ILVL) * NPTO          + 1
            IS  = IOFFS (IBLK,ILVL) * 4 * 3         + 1
            IV  = IOFFS (IBLK,ILVL)                 + 1
            IF ( .NOT. CONTRB ) THEN
               CALL RDREST (ITURB , IBLK, ILVL, IDIM(IBLK,ILVL), 
     1                      JDIM(IBLK,ILVL), KDIM(IBLK,ILVL), 
     2                      Q(IQ), QN(IQN), PROPS(IP), IEOF, IERRCD)
C
               IF (IEOF   .EQ. IYES) THEN
                  GO TO 139
               ENDIF
               IF (IERRCD .NE. INO ) THEN
                  WRITE (IOUT,1545)
 1545 FORMAT (' ','MAIN  : ERROR-> Error in restart grid level ',I3,
     1       /' ', 8X,'RUN ABORTING!')
                  STOP
               ENDIF
C
               IF (ITIMED .EQ. IUNSTD) THEN
                  IF (INITTD) THEN
                     IQ    = IOFFQ (IBLK,ILVL) * NQ          + 1
                     IQN   = IOFFQN(IBLK,ILVL) * NQ * NTMLVL + 1
                     CALL INITQN (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                            KDIM(IBLK,ILVL), Q(IQ), QN(IQN), 
     2                            QP, QM, IMODEL)
                  ELSE
                     NNPTS = (IDIM(IBLK,ILVL)+3) * (JDIM(IBLK,ILVL)+3) *
     1                       (KDIM(IBLK,ILVL)+3) * (NTIME-1)
                     CALL CONP2C (ITURB, NNPTS, QN(IQN), RES)
                  ENDIF
               ENDIF
C
            ELSE
               ITURBX = IMODLX( 1)
               CALL RDREST (ITURBX, IBLK, ILVL, IDIM(IBLK,ILVL), 
     1                      JDIM(IBLK,ILVL), KDIM(IBLK,ILVL), 
     2                      RES,   QN(IQN), W, IEOF, IERRCD)
C
               IF (IEOF   .EQ. IYES) THEN
                  GO TO 139
               ENDIF
               IF (IERRCD .NE. INO ) THEN
                  STOP
               ENDIF
               IF (ITIMED .EQ. IUNSTD) THEN
                  WRITE (IOUT,1555)
 1555 FORMAT (' ','MAIN  : ERROR-> Conversion of unsteady calculation ',
     1            'is not supported.',
     2       /' ', 8X,'RUN ABORTING!')
                  STOP
               ENDIF
C
C     Temporary variables passed to the conversion routine
C         QP     for QC (variables in one-dimensional array)
C         FM     for PROPSC (properties in one-dimensional array - this is 
C                            passed in using FM.  Must always be sure that 
C                            NP <= NQ or else will need to fix this)
C         QM     for DQDX (derivatives of Q in x direction at cell centers)
C         DQP    for DQDY (derivatives of Q in y direction at cell centers)
C         DQM    for DQDZ (derivatives of Q in z direction at cell centers)
C
               WRITE (IOUT,1550)
 1550 FORMAT (/' ','MAIN  : Converting Solution to New Turbulence ',
     1             'Model')
               CALL CONVRT (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                      KDIM(IBLK,ILVL), 
     2                      IMODLX, RES,   W,         NQFROM, NPFROM,
     3                      IMODEL, Q(IQ), PROPS(IP), NQTO,   NPTO,
     4                      S(IS), VOL(IV), 
     5                      QP, FM, QM, DQP, DQM, TAU, F, 
     6                      IBCDAT(1,1,IBLK,ILVL), NBCS(IBLK), IFDELQ)
            ENDIF
  135    CONTINUE
C
C     Output message, Increment grid level counter, Return to read next 
C     coarser grid.  
C     Do not read coarser grids if current run does not require it.
C
         WRITE (IOUT,1570) ILVL
 1570 FORMAT (' ','MAIN  :         Read Restart Grid Level ',I3)
         ILVL = ILVL + 1
         IF (ILVL .GT. NLVLS) GO TO 139
         GO TO 130
C
C     Finished reading all grid levels
C
  139    CONTINUE
         WRITE (IOUT,*) 
C
C     Now reset NQ, NF, NP back to NQTO, NFTO, NPTO
C
         NQ     = NQTO
         NF     = NFTO
         NP     = NPTO
         CLOSE (IRDRST)
C
C     Check that ITTOT + NITS(all mesh sequences) is less than MXHIST
C
         NITTOT = ITTOT
         DO 138 ICFS = 1, NCFS
            NITTOT = NITTOT + NITS(ICFS)
  138    CONTINUE
         IF (NITTOT .GT. MXHIST) THEN
            WRITE (IOUT,1580) NITTOT, MXHIST
            STOP
         ENDIF
      ENDIF
 1580 FORMAT (' ','MAIN  : ERROR-> Length of residual history    ',
     1            ' (ITTOT + NITS)      ',I6,
     2       /' ',15X,' exceeds maximum allocated      (MXHIST ',
     3            'in hstry.h) ',I6,
     4       /' ',15X,' Increase MXHIST to be greater than ',
     5                'ITTOT + NITS and recompile.',
     6       /' ', 8X,'RUN ABORTING!')
C
C 8.  CALCULATE YNORML
C     a. Calculate y_normal for initial run or if user forces
C     b. Calculate R(cell centers);  
C        Set y_normal = abs(y_normal) to allow for re-setting of transition
C        locations
C     c. Set y_normal = - (y_normal) to specify transition location
C     d. If RESTART calculate forces if needed for the 2D farfield vortex bc
C
C 8-a. Calculate y_normal - Do not recalculate if restarting, unless
C      user forces from input.
C         Use DTJ as temporary storage for y_normal.
C
      IF ( (ITURB .GE. ITLMNR .AND. (.NOT. RESTRT))
     1     .OR. YNCALC ) THEN
         CALL YNORM  (NXLVLS, NXBLKS, NXBCS, NXNODE,
     1                IOFFS, IOFFP,
     2                NLVLS,  NBLKS,  NBCS,  IBCDAT,
     3                IDIM, JDIM, KDIM, R, PROPS, DTJ)
      END IF
C
C 8-b. Calculate R(Cell Centers); Set y_normal = abs(y_normal)
C      Set R to contain position vector of cell center for the normal dependent
C      source terms of the Reynolds stress turbulent models.  R(cell centers)
C      is also needed in the output subroutine.
C      Use RES as temporary storage of one block, then copy the R(cell centers)
C      back into R.
C
      DO 88 ILVL = 1, NLVLS
         DO 87 IBLK = 1, NBLKS
            IR   = IOFFS(IBLK,ILVL) * 3     + 1
            CALL RCENTR (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                   KDIM(IBLK,ILVL), R(IR), RES, 
     2                   IBCDAT(1,1,IBLK,ILVL), NBCS(IBLK))
            CALL RCCOPY (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                   KDIM(IBLK,ILVL), RES, R(IR))
C
            IP   = IOFFP(IBLK,ILVL) * NP    + 1
            NPTS = (IDIM(IBLK,ILVL)+3) * (JDIM(IBLK,ILVL)+3) *
     1             (KDIM(IBLK,ILVL)+3)
            CALL YNABS  (NPTS, PROPS(IP))
   87    CONTINUE
   88 CONTINUE
C
C 8-c. Set y_normal = - (y_normal) to specify transition location
C      Coarsen and set for all mesh levels
C N.B. This destroys the information in the ITRSEG array.
C
      DO ISEG = 1, NTRSEG
         IBLK = ITRSEG(1,ISEG)
         ILVL = 1
         IP   = IOFFP(IBLK,ILVL) * NP    + 1
         CALL YNTRAN (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL), 
     1                KDIM(IBLK,ILVL), PROPS(IP), ITRSEG(1,ISEG))
         DO ILVL = 2, NLVLS
            DO L = 2, 7
               ITRSEG(L,ISEG) = ITRSEG(L,ISEG) / 2 + 1
            END DO
            IP   = IOFFP(IBLK,ILVL) * NP    + 1
            CALL YNTRAN (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL), 
     1                   KDIM(IBLK,ILVL), PROPS(IP), ITRSEG(1,ISEG))
         END DO
      END DO
C
C 8-d. If RESTART then calculate forces if needed for the 2D farfield vortex bc
C      Calculate forces on the grid level to begin the restart on.
C
      IF (RESTRT) THEN
         ILVL = NCFS
         IF (IFFORC .EQ. IYES) THEN
            CALL FRCINT (ITURB,NBLKS,IDIM(1,ILVL), 
     1                   JDIM(1,ILVL),  KDIM(1,ILVL), 
     2                   IOFF(1,ILVL),  IOFFS(1,ILVL),
     3                   IOFFQ(1,ILVL), IOFFP(1,ILVL),
     4                   Q, PROPS, S, R, NXBCS, NBCS, 
     5                   IBCDAT(1,1,1,ILVL), FORCE, CL, CD, FMOMENT)
         ENDIF
      ENDIF
C
C 9.  INITIALIZE PROPERTIES
C     Initialize the PROPERTIES array.
C     Use QP, QM, F and DQP as temporary arrays.
C
      DO 150 ILVL = 1, NLVLS
         DO 140 IBLK = 1, NBLKS
            IQ = IOFFQ(IBLK,ILVL) * NQ    + 1
            IP = IOFFP(IBLK,ILVL) * NP    + 1
            IS = IOFFS(IBLK,ILVL) * 4 * 3 + 1
            IV = IOFFS(IBLK,ILVL)         + 1
            IPSTRT = 1
            IPEND  = IDIM(IBLK,ILVL) + 1
            IF (ITURB .GE. ITLMNR) THEN
               CALL PRPRTY (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                      KDIM(IBLK,ILVL), IPSTRT, IPEND, NXSECT,
     2                      Q(IQ), PROPS(IP), S(IS), VOL(IV),
     3                      IBCDAT(1,1,IBLK,ILVL), NBCS(IBLK), 
     4                      IMODEL, IFDELQ,
     5                      STEMP, QP, F, FM, QM, DQP, DQM, WORKX)
            ENDIF
  140    CONTINUE
  150 CONTINUE
C
C 10. BOUNDARY CONDITION INITIALIZATION AND VERIFICATION
C     Do boundary conditions for all blocks and all levels.  Use QP for work.
C     a. Initialize boundary data to infinity to check input bc data
C     b. Update boundary conditions
C     c. Update Cuts
C     d. Update properties array (synchronize properties with mean flow)
C     e. Update boundary conditions (ensure turbulent bcs use updated props)
C     f. Check boundary data to verify input bc data
C
      DO 170 ILVL = 1, NLVLS
         DO 160 IBLK = 1, NBLKS
            IQ = IOFFQ(IBLK,ILVL) * NQ    + 1
            IP = IOFFP(IBLK,ILVL) * NP    + 1
            IS = IOFFS(IBLK,ILVL) * 4 * 3 + 1
            IV = IOFFS(IBLK,ILVL)         + 1
            IR = IOFFS(IBLK,ILVL) * 3     + 1
C
C     Storage for PROPS
C
            IF (ITURB .GE. ITLMNR) THEN
               NPPTS = (IDIM(IBLK,ILVL)+3) * (JDIM(IBLK,ILVL)+3) *
     1                 (KDIM(IBLK,ILVL)+3)
            ELSE
               NPPTS = NXPPTS
            ENDIF
C
C 10-a. Initialize boundary data to infinity to check input bc data
C       Use WORKX for QTEST.
C
            CALL INITBC (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                   KDIM(IBLK,ILVL), Q(IQ), WORKX, ITURB)
C
C 10-b. Update boundary conditions
C
            CALL BC (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL), KDIM(IBLK,ILVL),
     1               1, IDIM(IBLK,ILVL)+1, Q(IQ), PROPS(IP), NPPTS,
     2               S(IS), VOL(IV), R(IR),
     3               NBCS(IBLK), IBCDAT(1,1,IBLK,ILVL), 
     4               NXPROF, NPROFL, IOFQBC(1,ILVL), QBC,
     5               VISCOS(1,IBLK), IMODEL, QP, NXSECT, CL, RCNTR)
 160     CONTINUE
 170  CONTINUE
C
C 10-c. Update Cuts
C       Do block boundary conditions (CUTS) for all levels
C
      DO 190 ILVL = 1, NLVLS
         DO 180 ICUT = 1, NCUTS
            IBLK1 = ICUTS( 1,ICUT,ILVL)
            IBLK2 = ICUTS(10,ICUT,ILVL)
            IQ1   = IOFFQ(IBLK1,ILVL) * NQ    + 1
            IQ2   = IOFFQ(IBLK2,ILVL) * NQ    + 1
            CALL CUT (IDIM(IBLK1,ILVL), JDIM(IBLK1,ILVL), 
     1                KDIM(IBLK1,ILVL), Q(IQ1),
     2                IDIM(IBLK2,ILVL), JDIM(IBLK2,ILVL), 
     3                KDIM(IBLK2,ILVL), Q(IQ2), ICUTS(1,ICUT,ILVL), NQ)
            IF (ITURB .GE. ITLMNR) THEN
               IP1   = IOFFP(IBLK1,ILVL) * NP    + 1
               IP2   = IOFFP(IBLK2,ILVL) * NP    + 1
               CALL CUT (IDIM(IBLK1,ILVL), JDIM(IBLK1,ILVL), 
     1                   KDIM(IBLK1,ILVL), PROPS(IP1),
     2                   IDIM(IBLK2,ILVL), JDIM(IBLK2,ILVL), 
     3                   KDIM(IBLK2,ILVL), PROPS(IP2), 
     4                   ICUTS(1,ICUT,ILVL), NP)
            ENDIF
  180    CONTINUE
  190 CONTINUE
C
C 10-d,e. Recalculate properties and boundary conditions as properties
C         are functions of quantities just set in BC and some BCs are
C         functions of the properties.
C
      DO 200 ILVL = 1, NLVLS
         DO 195 IBLK = 1, NBLKS
            IQ = IOFFQ(IBLK,ILVL) * NQ    + 1
            IP = IOFFP(IBLK,ILVL) * NP    + 1
            IS = IOFFS(IBLK,ILVL) * 4 * 3 + 1
            IV = IOFFS(IBLK,ILVL)         + 1
            IR = IOFFS(IBLK,ILVL) * 3     + 1
C
C     Storage for PROPS
C
            IF (ITURB .GE. ITLMNR) THEN
               NPPTS = (IDIM(IBLK,ILVL)+3) * (JDIM(IBLK,ILVL)+3) *
     1                 (KDIM(IBLK,ILVL)+3)
            ELSE
               NPPTS = NXPPTS
            ENDIF
C
C 10-b. Update boundary conditions
C
            CALL BC (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL), KDIM(IBLK,ILVL),
     1               1, IDIM(IBLK,ILVL)+1, Q(IQ), PROPS(IP), NPPTS,
     2               S(IS), VOL(IV), R(IR),
     3               NBCS(IBLK), IBCDAT(1,1,IBLK,ILVL), 
     4               NXPROF, NPROFL, IOFQBC(1,ILVL), QBC,
     5               VISCOS(1,IBLK), IMODEL, QP, NXSECT, CL, RCNTR)
C
            IF (ITURB .GE. ITLMNR) THEN
               IPSTRT = 1
               IPEND  = IDIM(IBLK,ILVL) + 1
C
C 10-d. Update properties array (synchronize properties with mean flow)
C
               CALL PRPRTY (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                      KDIM(IBLK,ILVL), IPSTRT, IPEND, NXSECT,
     2                      Q(IQ), PROPS(IP), S(IS), VOL(IV),
     3                      IBCDAT(1,1,IBLK,ILVL), NBCS(IBLK), 
     4                      IMODEL, IFDELQ,
     5                      STEMP, QP, F, FM, QM, DQP, DQM, WORKX)
C
C 10-e. Update boundary conditions (ensure turbulent bcs use updated props)
C
               CALL BC (IDIM(IBLK,ILVL),JDIM(IBLK,ILVL),KDIM(IBLK,ILVL),
     1                  1, IDIM(IBLK,ILVL)+1, Q(IQ), PROPS(IP), NPPTS,
     2                  S(IS), VOL(IV), R(IR),
     3                  NBCS(IBLK), IBCDAT(1,1,IBLK,ILVL), 
     4                  NXPROF, NPROFL, IOFQBC(1,ILVL), QBC,
     5                  VISCOS(1,IBLK), IMODEL, QP, NXSECT, CL, RCNTR)
            ENDIF
 195     CONTINUE
 200  CONTINUE
C
C 10-f. Check boundary data to verify input bc data
C       Check each block on each level to verify that all boundary conditions
C       have been initialized.
C
      IBCERR = 0
      DO 220 ILVL = 1, NLVLS
         DO 210 IBLK = 1, NBLKS
            IQ = IOFFQ(IBLK,ILVL) * NQ    + 1
            CALL CHCKBC (ILVL, IBLK, IDIM(IBLK,ILVL), JDIM(IBLK,ILVL), 
     1                   KDIM(IBLK,ILVL), Q(IQ), IBCERR)
  210    CONTINUE
  220 CONTINUE
      IF (IBCERR .GT. 0) THEN
         WRITE (IOUT,1400) IBCERR
         STOP
      ENDIF
 1400 FORMAT (' ','MAIN  : ERROR-> ',I8,' Unset Boundary Conditions.',
     1       /' ', 8X,'RUN ABORTING!')
C
C 11. MARCHING CALCULATION
C     If Marching, then call marching procedure
C
      IF (ISOLVR .EQ. IMARCH) THEN
         IBLK  = 1
         ILVL  = 1
         ICFS  = 1
         IQ    = IOFFQ (IBLK,ILVL) * NQ          + 1
         IQN   = IOFFQN(IBLK,ILVL) * NQ * NTMLVL + 1
         IP    = IOFFP (IBLK,ILVL) * NP          + 1
         IS    = IOFFS (IBLK,ILVL) * 4 * 3       + 1
         IV    = IOFFS (IBLK,ILVL)               + 1
         IR    = IOFFS (IBLK,ILVL) * 3           + 1
         ID    = IOFF  (IBLK,ILVL)               + 1
C
C     Temporary(?) fix for the global storage of the LU when marching
C
         NABC = NXABC / 2
C
C     Storage for QN - Marching shouldn't be unsteady!
C
         IF (ITIMED .EQ. IUNSTD) THEN
            NQNPTS = (IDIM(IBLK,ILVL)+3) * (JDIM(IBLK,ILVL)+3) *
     1               (KDIM(IBLK,ILVL)+3)
            IFL2DQ = .TRUE.
         ELSE
            NQNPTS = NXQNPT
            IFL2DQ = .FALSE.
         ENDIF
C
C     Storage for PROPS
C
         IF (ITURB .GE. ITLMNR) THEN
            NPPTS = (IDIM(IBLK,ILVL)+3) * (JDIM(IBLK,ILVL)+3) *
     1              (KDIM(IBLK,ILVL)+3)
         ELSE
            NPPTS = NXPPTS
         ENDIF
C
C     Initialize error code
C
         IERRCD = INO
C
         CALL MARCH (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL), KDIM(IBLK,ILVL), 
     1          NXSECT, NABC, Q(IQ), QN(IQN), NQNPTS, PROPS(IP), NPPTS,
     2          S(IS), VOL(IV), R(IR),
     3          DTJ(ID), RES, FHAT, DQP,DQM,QP,QM,TAU,F,FM,STEMP,WORKX,
     4          IBCDAT(1,1,IBLK,ILVL), NBCS(IBLK), 
     5          NXPROF, NPROFL, IOFQBC(1,ILVL), QBC,
     6          ICUTS(1,1,ILVL), NCUTS, R2NORM, ITRSMX, RESMAX, IRESMX,
     7          IORDER(1,IBLK,1,1), LIMIT(1,IBLK), RKAP(1,IBLK),
     8          IFLUX(IBLK), VISCOS(1,IBLK), IMODEL, IFDELQ,
     9          ABC, RM, RTEMP, AT, AI,
     A          IMSTRT, IMEND, NITS(ICFS), NITFO, NITALG, ITUPJ,
     B          TOLER, VGNRN, SIGMA, ITER, QMIN, NFIXQ, NRELIZ, 
     C          NPRLIM, NPRNEG, IFL2DQ, R2DQDT, CL, RCNTR, IERRCD)
C
         IF (IERRCD .NE. INO) THEN
            WRITE (IOUT,1405) QMIN(1), NFIXQ(1), QMIN(5), NFIXQ(5)
            GO TO 910
         ENDIF
 1405 FORMAT (///' ','MAIN  : FIXQ  : ERROR-> Density  < ',1PE9.2,
     1               ' at',I7,' locations'
     2          /' ','                        Pressure < ',1PE9.2,
     3               ' at',I7,' locations'
     4          /' ','RUN ABORTING! Attempting to write output ',
     5               'and restart files.',
     6        ///' ')
      ELSE
C
C 12. TIME DEPENDENT CALCULATION
C
C     Calculate execution time of run
C
#ifdef CPU_TIME
         RUNTIME = DTIME (TARRAY)
#endif
C
         MRCHNG = .FALSE.
C
C     Initialize error code
C
         IERRCD = INO
C
C     Loop over coarse to fine sequences
C
         ITRUN  = 0
         DO 900 ICFS = 1, NCFS
            LVLONE = NCFS - ICFS + 1
C
C     Loop over number of iterations
C
            DO 800 ITER = 1, NITS(LVLONE)
C
C     Check for presence of STOP file; if it exists, write restart and quit
C
               OPEN (ISTOPFL,FILE='STOP',STATUS='OLD',ERR=701)
               WRITE (IOUT,1490)
               GO TO 905
 1490 FORMAT (//' ','MAIN  : STOP file encountered! Run stopping.'//)
C
C     STOP does not exist; continue run
C
  701          CONTINUE
C
C     If unsteady calculation then 
C       i.   Update TIME
C       ii.  Collect Q to coarser levels at time level n
C       iii. Shift Q, QN for next time step
C
               IF (ITIMED .EQ. IUNSTD) THEN
C
C     i.   Update TIME
C
                  TIME = TIME + TIMEDT
C
C     ii.  Collect Q to coarser levels at time level n
C
                  DO 710 ILVL = LVLONE+1, LVLONE + MGLVLS(LVLONE) - 1
                     DO 705 IBLK = 1, NBLKS
                        IQF = IOFFQ(IBLK,ILVL-1) * NQ    + 1
                        IVF = IOFFS(IBLK,ILVL-1)         + 1
                        IQ  = IOFFQ(IBLK,ILVL)   * NQ    + 1
                        IV  = IOFFS(IBLK,ILVL)           + 1
C
C     Temporary storage: QP for QS, RES for TKE
C
                        CALL COLQ2H (IDIM(IBLK,ILVL-1),
     1                              JDIM(IBLK,ILVL-1),KDIM(IBLK,ILVL-1),
     2                              Q(IQF), VOL(IVF),
     3                              IDIM(IBLK,ILVL),  JDIM(IBLK,ILVL),
     4                              KDIM(IBLK,ILVL),  Q(IQ), VOL(IV), 
     5                              QP, ITURB, RES)
  705                CONTINUE
  710             CONTINUE
C
C     iii. Shift Q, QN for next time step
C
                  DO 720 ILVL = LVLONE, LVLONE + MGLVLS(LVLONE) - 1
                     DO 715 IBLK = 1, NBLKS
                        IQ     = IOFFQ (IBLK,ILVL) * NQ          + 1
                        IQN    = IOFFQN(IBLK,ILVL) * NQ * NTMLVL + 1
C
                        CALL QNSHFT (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                               KDIM(IBLK,ILVL), Q(IQ), QN(IQN), 
     2                               QP, QM, IMODEL)
  715                CONTINUE
  720             CONTINUE
               ENDIF
C
C     End of data shift for next time-accurate time step
C
C     Q is Q^(p) and is initialized to Q^(n) to begin subiterations, 
C     QN(1) = Q^(n) and QN(2) = Q^(n-1)
C
C     Perform sub-iterations for each time step
C
               DO 790 ITSUB = 1, NSUBIT
C
C     Set starting multigrid level to be NCFS for coarsest mesh and
C     starting level to be 1 for finest mesh and loop over
C     the multigrid levels.
C
                  DO 700 ILVL = LVLONE, LVLONE + MGLVLS(LVLONE) - 1
C
C     Restrict from fine to coarse level if current level is not finest level
C     Begin Restriction:
C
                     IF (ILVL .NE. LVLONE) THEN
C
C     Collect Q
C
                        DO 300 IBLK = 1, NBLKS
                           IQF = IOFFQ(IBLK,ILVL-1) * NQ    + 1
                           IVF = IOFFS(IBLK,ILVL-1)         + 1
                           IQ  = IOFFQ(IBLK,ILVL)   * NQ    + 1
                           IP  = IOFFP(IBLK,ILVL)   * NP    + 1
                           IS  = IOFFS(IBLK,ILVL)   * 4 * 3 + 1
                           IV  = IOFFS(IBLK,ILVL)           + 1
                           IR  = IOFFS(IBLK,ILVL)   * 3     + 1
C
                           IF (ITURB .GE. ITLMNR) THEN
                              NPPTS = (IDIM(IBLK,ILVL)+3) * 
     1                                (JDIM(IBLK,ILVL)+3) *
     2                                (KDIM(IBLK,ILVL)+3)
                           ELSE
                              NPPTS = NXPPTS
                           ENDIF
C
C     Temporary storage: QP for QS, RES for TKE
C
                           CALL COLQ2H (IDIM(IBLK,ILVL-1),
     1                              JDIM(IBLK,ILVL-1),KDIM(IBLK,ILVL-1),
     2                              Q(IQF), VOL(IVF),
     3                              IDIM(IBLK,ILVL),  JDIM(IBLK,ILVL),
     4                              KDIM(IBLK,ILVL),  Q(IQ), VOL(IV),
     5                              QP, ITURB, RES)
C
                           IPSTRT = 1
                           IPEND  = IDIM(IBLK,ILVL) + 1
C
                           CALL BC     (IDIM(IBLK,ILVL),JDIM(IBLK,ILVL),
     1                           KDIM(IBLK,ILVL), IPSTRT, IPEND, 
     2                           Q(IQ), PROPS(IP), NPPTS, 
     3                           S(IS), VOL(IV), R(IR),
     4                           NBCS(IBLK), IBCDAT(1,1,IBLK,ILVL), 
     5                           NXPROF,NPROFL,IOFQBC(1,ILVL),QBC,
     6                           VISCOS(1,IBLK), IMODEL, QP, NXSECT, 
     7                           CL, RCNTR)
C
                           IF (ITURB .GE. ITLMNR) THEN
                              CALL PRPRTY (IDIM(IBLK,ILVL), 
     1                              JDIM(IBLK,ILVL), KDIM(IBLK,ILVL), 
     2                              IPSTRT, IPEND, NXSECT,
     3                              Q(IQ), PROPS(IP), S(IS), VOL(IV),
     4                              IBCDAT(1,1,IBLK,ILVL), NBCS(IBLK), 
     5                              IMODEL, IFDELQ,
     6                              STEMP, QP, F, FM, QM, DQP,DQM,WORKX)
                           ENDIF
 300                    CONTINUE
C
C      Update CUTS
C
                        DO 320 ICUT = 1, NCUTS
                           IBLK1 = ICUTS( 1,ICUT,ILVL)
                           IBLK2 = ICUTS(10,ICUT,ILVL)
                           IQ1   = IOFFQ(IBLK1,ILVL) * NQ    + 1
                           IQ2   = IOFFQ(IBLK2,ILVL) * NQ    + 1
                           CALL CUT (IDIM(IBLK1,ILVL), JDIM(IBLK1,ILVL),
     2                               KDIM(IBLK1,ILVL), Q(IQ1),
     3                               IDIM(IBLK2,ILVL), JDIM(IBLK2,ILVL),
     4                               KDIM(IBLK2,ILVL), Q(IQ2), 
     5                               ICUTS(1,ICUT,ILVL), NQ)
                           IF (ITURB .GE. ITLMNR) THEN
                              IP1   = IOFFP(IBLK1,ILVL) * NP    + 1
                              IP2   = IOFFP(IBLK2,ILVL) * NP    + 1
                              CALL CUT 
     1                              (IDIM(IBLK1,ILVL), JDIM(IBLK1,ILVL),
     2                               KDIM(IBLK1,ILVL), PROPS(IP1),
     3                               IDIM(IBLK2,ILVL), JDIM(IBLK2,ILVL),
     4                               KDIM(IBLK2,ILVL), PROPS(IP2), 
     5                               ICUTS(1,ICUT,ILVL), NP)
                           ENDIF
  320                   CONTINUE
C
C     Copy QC and store
C
                        DO 340 IBLK = 1, NBLKS
                           IQ    = IOFFQ (IBLK,ILVL)   * NQ    + 1
                           IQCRS = IOFFQC(IBLK,ILVL)   * NQ    + 1
C
                           NPTS  = (IDIM(IBLK,ILVL)+3) * 
     1                             (JDIM(IBLK,ILVL)+3) *
     2                             (KDIM(IBLK,ILVL)+3) * NQ
C
                           CALL CPVCTR (NPTS, Q(IQ), QCRS(IQCRS))
  340                   CONTINUE
C
C     Calculate coarse and fine grid residuals and coarse grid forcing function
C
                        NPRLIM = 0
                        NPRNEG = 0
                        DO 360 IBLK = 1, NBLKS
                           IQF    = IOFFQ (IBLK,ILVL-1) * NQ         + 1
                           IQNF   = IOFFQN(IBLK,ILVL-1) * NQ *NTMLVL + 1
                           IPF    = IOFFP (IBLK,ILVL-1) * NP         + 1
                           ISF    = IOFFS (IBLK,ILVL-1) * 4 * 3      + 1
                           IVF    = IOFFS (IBLK,ILVL-1)              + 1
                           IDF    = IOFF  (IBLK,ILVL-1)              + 1
                           IRF    = IOFFS (IBLK,ILVL-1) * 3          + 1
C
                           IQC    = IOFFQ (IBLK,ILVL)   * NQ         + 1
                           IQNC   = IOFFQN(IBLK,ILVL)   * NQ *NTMLVL + 1
                           IPC    = IOFFP (IBLK,ILVL)   * NP         + 1
                           ISC    = IOFFS (IBLK,ILVL)   * 4 * 3      + 1
                           IVC    = IOFFS (IBLK,ILVL)                + 1
                           IDC    = IOFF  (IBLK,ILVL)                + 1
                           IRC    = IOFFS (IBLK,ILVL)   * 3          + 1
C
                           IRSCRS = IOFFQC(IBLK,ILVL)   * NF         + 1
C
                           IF (ITIMED .EQ. IUNSTD) THEN
                              NQNPTF = (IDIM(IBLK,ILVL-1)+3) * 
     1                                 (JDIM(IBLK,ILVL-1)+3) *
     2                                 (KDIM(IBLK,ILVL-1)+3)
                              NQNPTC = (IDIM(IBLK,ILVL)  +3) * 
     1                                 (JDIM(IBLK,ILVL)  +3) *
     2                                 (KDIM(IBLK,ILVL)  +3)
                           ELSE
                              NQNPTF = NXQNPT
                              NQNPTC = NXQNPT
                           ENDIF
C
                           IF (ITURB .GE. ITLMNR) THEN
                              NPPTSF = (IDIM(IBLK,ILVL-1)+3) * 
     1                                 (JDIM(IBLK,ILVL-1)+3) *
     2                                 (KDIM(IBLK,ILVL-1)+3)
                              NPPTSC = (IDIM(IBLK,ILVL)  +3) * 
     1                                 (JDIM(IBLK,ILVL)  +3) *
     2                                 (KDIM(IBLK,ILVL)  +3)
                           ELSE
                              NPPTSF = NXPPTS
                              NPPTSC = NXPPTS
                           ENDIF
C
                           IF (ILVL .EQ. (LVLONE + 1) ) THEN
                              INITRS = .TRUE.
                           ELSE
                              INITRS = .FALSE.
                              IRSCRF = IOFFQC(IBLK,ILVL-1)   * NF  + 1
                              NPTRSF = (IDIM(IBLK,ILVL-1)+1) *
     1                                 (JDIM(IBLK,ILVL-1)+1) *
     2                                 (KDIM(IBLK,ILVL-1)+1) * NF
                              CALL CPVCTR (NPTRSF, RESCRS(IRSCRF), RES)
                           ENDIF
                           MRCHNG = .FALSE.
                           IFL2DQ = .FALSE.
C
                           IPSTRT = 2
                           IPEND  = IDIM(IBLK,ILVL-1)
                           CALL RESID  (IDIM(IBLK,ILVL-1), 
     1                              JDIM(IBLK,ILVL-1),KDIM(IBLK,ILVL-1),
     2                              IPSTRT, IPEND, NXSECT,
     3                              Q(IQF), QN(IQNF), NQNPTF,
     4                              PROPS(IPF), NPPTSF, S(ISF), 
     5                              VOL(IVF), R(IRF), DTJ(IDF), RES,
     6                              FHAT, DQP, DQM, QP, QM, TAU, F, FM, 
     7                              STEMP, WORKX,
     8                              IBCDAT(1,1,IBLK,ILVL-1), NBCS(IBLK),
     9                              IORDER(1,IBLK,ILVL-1,ICFS), 
     A                              LIMIT(1,IBLK), RKAP(1,IBLK), 
     B                              IFLUX(IBLK), VISCOS(1,IBLK), IMODEL,
     C                              INITRS, IFDELQ, MRCHNG, VGNRN,SIGMA,
     D                              NPRLIM, NPRNEG, IFL2DQ, R2DQDT)
C
C     Set logical to calculate multigrid forcing function on first resid call
C
                           IFMGFF = .TRUE.
C
                           CALL COLR2H (IDIM(IBLK,ILVL-1),
     1                              JDIM(IBLK,ILVL-1),KDIM(IBLK,ILVL-1),
     2                              RES, 
     3                              IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     4                              KDIM(IBLK,ILVL), RESCRS(IRSCRS), QP)
  360                   CONTINUE
                     ENDIF
C
C     End   Restriction:
C
C     Update time step.
C     Use DQP as temporary storage for PFACTR.
C  
                     IF ( MOD((ITER-1), ITUPDT) .EQ. 0 .AND.
     1                    ITSUB .EQ. 1 ) THEN
                        DTMIN =  999.E0
                        DTMAX = -999.E0
                        DO 400 IBLK = 1, NBLKS
                           IQ    = IOFFQ(IBLK,ILVL) * NQ    + 1
                           IP    = IOFFP(IBLK,ILVL) * NP    + 1
                           ID    = IOFF (IBLK,ILVL)         + 1
                           IV    = IOFFS(IBLK,ILVL)         + 1
                           IS    = IOFFS(IBLK,ILVL) * 4 * 3 + 1
                           CALL DELTAT (IBLK,        IDIM(IBLK,ILVL), 
     1                              JDIM(IBLK,ILVL), KDIM(IBLK,ILVL), 
     2                              2, IDIM(IBLK,ILVL),
     3                              Q(IQ), PROPS(IP), DTJ(ID), VOL(IV),
     4                              S(IS), DQP,
     5                              VISCOS(1,IBLK), ITURB, ITER, CFLCUR,
     6                              DTMIN,IBLKMN,IDTMIN,JDTMIN,KDTMIN,
     7                              DTMAX,IBLKMX,IDTMAX,JDTMAX,KDTMAX)
 400                    CONTINUE
C
C     Output minimum and maximum time step/CFL
C
                        IF (DT .GT. 0.0E0) THEN
                           WRITE (IOUT,1440) CFLCUR,
     1                               DTMIN,IBLKMN,IDTMIN,JDTMIN,KDTMIN,
     2                               DTMAX,IBLKMX,IDTMAX,JDTMAX,KDTMAX
                        ELSE
                           WRITE (IOUT,1450) CFLCUR,
     1                               DTMIN,IBLKMN,IDTMIN,JDTMIN,KDTMIN,
     2                               DTMAX,IBLKMX,IDTMAX,JDTMAX,KDTMAX
                        ENDIF
                     ENDIF
 1440 FORMAT (' ','MAIN  : DT =',1PE10.3,'  CFLMIN =',1PE10.3,
     1            ' At BLOCK:I,J,K=(',I3,':',I4,',',I4,',',I4,')',
     2       /' ',22X,'  CFLMAX =',1PE10.3,
     3            ' At BLOCK:I,J,K=(',I3,':',I4,',',I4,',',I4,')')
 1450 FORMAT (' ','MAIN  : CFL =',1PE10.3,'  DTMIN =',1PE10.3,
     1            ' At BLOCK:I,J,K=(',I3,':',I4,',',I4,',',I4,')',
     2       /' ',23X,'  DTMAX =',1PE10.3,
     3            ' At BLOCK:I,J,K=(',I3,':',I4,',',I4,',',I4,')')
C
C     Perform ITSLVL iterations (time step integrations) per level.
C
                     DO 600 INRIT = 1, ITSLVL(ILVL,LVLONE)
C
C     Reset norms to zero
C
                        DO 420 L = 1, NF + 1
                           R2NORM(L) = 0.E0
                           R2DQDT(L) = 0.E0
 420                    CONTINUE
C
C     Perform NSTAGE stages of multistage scheme
C
                        DO 550 ISTAGE = 1, NSTAGE
C
C     Reset counters to zero
C
                           DO 425 L = 1, NQ
                              NFIXQ (L) = 0
                              NRELIZ(L) = 0
 425                       CONTINUE
C
C     Counter for limiting of turbulence production
C
                           NPRLIM = 0
                           NPRNEG = 0
C
C     Set up control parameter for calculating maximum residual
C
                           RESMAX = 0.E0
                           IF (INRIT  .EQ. 1 .AND. 
     1                         ISTAGE .EQ. 1 .AND.
     2                         MOD((ITER),ITRSMX).EQ.0) THEN
                              IFRSMX = IYES
                           ELSE
                              IFRSMX = INO
                           ENDIF
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
                              IF (ITIMED .EQ. IUNSTD) THEN
                                 NQNPTS = (IDIM(IBLK,ILVL)+3) * 
     1                                    (JDIM(IBLK,ILVL)+3) *
     2                                    (KDIM(IBLK,ILVL)+3)
                                 IF (ILVL   .EQ. LVLONE .AND.
     1                               ISTAGE .EQ. 1) THEN
                                    IFL2DQ = .TRUE.
                                 ELSE
                                    IFL2DQ = .FALSE.
                                 ENDIF
                              ELSE
                                 NQNPTS = NXQNPT
                                 IFL2DQ = .FALSE.
                              ENDIF
C
C     Storage for PROPS
C
                              IF (ITURB .GE. ITLMNR) THEN
                                 NPPTS = (IDIM(IBLK,ILVL)+3) * 
     1                                   (JDIM(IBLK,ILVL)+3) *
     2                                   (KDIM(IBLK,ILVL)+3)
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
                              IF (ITER .LE. NITALG) THEN
                                 MUTALG = .TRUE.
                              ELSE
                                 MUTALG = .FALSE.
                              ENDIF
C
C     If running a viscous case, calculate the PROPS array
C     Use QP, QM, F and DQP as temporary storage.
C
                              IF (ITURB  .GE. ITLMNR .AND.
     1                            ISTAGE .EQ. 1) THEN
                                 CALL PRPRTY (IDIM(IBLK,ILVL), 
     1                              JDIM(IBLK,ILVL), KDIM(IBLK,ILVL), 
     2                              1, IDIM(IBLK,ILVL)+1, NXSECT,
     3                              Q(IQ), PROPS(IP), S(IS), VOL(IV),
     4                              IBCDAT(1,1,IBLK,ILVL), NBCS(IBLK), 
     5                              IMODEL, IFDELQ,
     6                              STEMP, QP, F, FM, QM, DQP,DQM,WORKX)
                              ENDIF
C
C      Calculate the Residual for the current block and level
C
                              INITRS = .TRUE.
                              IPSTRT = 2
                              IPEND  = IDIM(IBLK,ILVL)
                              CALL RESID  (IDIM(IBLK,ILVL), 
     1                              JDIM(IBLK,ILVL),KDIM(IBLK,ILVL), 
     2                              IPSTRT, IPEND, NXSECT,
     3                              Q(IQ), QN(IQN), NQNPTS,
     4                              PROPS(IP), NPPTS, S(IS), 
     5                              VOL(IV), R(IR), DTJ(ID), RES,
     6                              FHAT, DQP, DQM, QP, QM, TAU, F, FM, 
     7                              STEMP, WORKX,
     8                              IBCDAT(1,1,IBLK,ILVL), NBCS(IBLK), 
     9                              ITORDR,
     A                              LIMIT(1,IBLK), RKAP(1,IBLK), 
     B                              IFLUX(IBLK), VISCOS(1,IBLK), IMODEL,
     C                              INITRS, IFDELQ, MRCHNG, VGNRN,SIGMA,
     D                              NPRLIM, NPRNEG, IFL2DQ, R2DQDT)
C
C     Calculate L2 norms only on the finest mesh
C          If first stage, call l2norm to calculate L2 norm of residual,
C                          call l2max  to calculate maximum residual
C
                              IF (ILVL .EQ. LVLONE) THEN
                                 IF (ISTAGE .EQ. 1) THEN
                                    CALL L2NORM (IDIM(IBLK,ILVL), 
     1                                  JDIM(IBLK,ILVL),KDIM(IBLK,ILVL),
     2                                  IPSTRT, IPEND,
     3                                  RES, VOL(IV), R2NORM)
                                    IF (IFRSMX .EQ. IYES) THEN
                                       CALL L2MAX (IBLK,IDIM(IBLK,ILVL),
     1                                  JDIM(IBLK,ILVL),KDIM(IBLK,ILVL),
     2                                  IPSTRT, IPEND,
     3                                  RES, VOL(IV), RESMAX, IRESMX)
                                    ENDIF
                                 ENDIF
                              ENDIF
C
C     For MG: i. For the first residual evaluation after coarsening,
C                calculate the MG forcing function
C            ii. Set the residual to RES + FF before time integration
C
                              IF (ILVL .NE. LVLONE) THEN
                                 IRSCRS = IOFFQC(IBLK,ILVL)   * NF  + 1
                                 NPTS   = (IDIM(IBLK,ILVL)+1) *
     1                                    (JDIM(IBLK,ILVL)+1) *
     2                                    (KDIM(IBLK,ILVL)+1) * NF
                                 IF (IFMGFF) THEN
c%%%%%                              CALL CMGFF (NPTS,RESCRS(IRSCRS),RES)
                                    CALL xCMGFF (IDIM(IBLK,ILVL),
     1                                 JDIM(IBLK,ILVL), KDIM(IBLK,ILVL),
     2                                 RESCRS(IRSCRS),RES)
                                 ENDIF
c%%%%%                           CALL CRESFF (NPTS, RESCRS(IRSCRS), RES)
                                 CALL xCRESFF (IDIM(IBLK,ILVL),
     1                                 JDIM(IBLK,ILVL), KDIM(IBLK,ILVL),
     2                                 RESCRS(IRSCRS),RES)
                              ENDIF
C
                              IF (ISOLVR .EQ. IRKN) THEN
                                 CALL RK (IBLK, IDIM(IBLK,ILVL), 
     1                          JDIM(IBLK,ILVL), KDIM(IBLK,ILVL),NXSECT,
     2                          NSTAGE, ISTAGE, RKALPH(ISTAGE),
     3                          Q(IQ), QN(IQN), NQNPTS,
     4                          PROPS(IP), NPPTS, S(IS), VOL(IV), 
     5                          R(IR),DTJ(ID), RES, FHAT, DQP,DQM,QP,QM,
     6                          TAU, F, FM, STEMP, WORKX,
     7                          IBCDAT(1,1,IBLK,ILVL), NBCS(IBLK),
     8                          ITORDR, LIMIT(1,IBLK), RKAP(1,IBLK),
     9                          IFLUX(IBLK), VISCOS(1,IBLK),
     A                          IMODEL, IFDELQ, W, 
     B                          QMIN, NFIXQ, NRELIZ, NPRLIM, NPRNEG,
     C                          IERRCD)
C
                              ELSE IF (ISOLVR .EQ. IAF3F) THEN
                                 CALL AF3F (IBLK, IDIM(IBLK,ILVL), 
     1                       JDIM(IBLK,ILVL), KDIM(IBLK,ILVL), NXSECT, 
     2                       NXABC, Q(IQ), QN(IQN), NQNPTS,
     3                       PROPS(IP), NPPTS, S(IS), VOL(IV), 
     4                       R(IR), DTJ(ID), RES, FHAT, DQP, DQM, QP,QM,
     5                       TAU, F, FM, STEMP, WORKX,
     6                       IBCDAT(1,1,IBLK,ILVL), NBCS(IBLK),
     7                       ITORDR, LIMIT(1,IBLK), RKAP(1,IBLK),
     8                       IFLUX(IBLK), VISCOS(1,IBLK), DIAG(1,IBLK),
     9                       IMODEL, IFDELQ, ABC, RTEMP, AT, AI,
     A                       QMIN, NFIXQ, NRELIZ, NPRLIM, NPRNEG,IERRCD)
C
                              ENDIF
C
C     Update Boundary Conditions
C
                              CALL BC(IDIM(IBLK,ILVL),JDIM(IBLK,ILVL),
     1                                KDIM(IBLK,ILVL), 1, 
     2                                IDIM(IBLK,ILVL)+1,Q(IQ),PROPS(IP),
     3                                NPPTS, S(IS), VOL(IV), R(IR),
     4                                NBCS(IBLK),IBCDAT(1,1,IBLK,ILVL), 
     5                                NXPROF,NPROFL,IOFQBC(1,ILVL),QBC,
     6                                VISCOS(1,IBLK),IMODEL,QP, NXSECT,
     7                                CL, RCNTR)
  500                      CONTINUE
C
C     Reset the logical controlling the calculation of the MG forcing function
C
                           IFMGFF = .FALSE.
C
C     Check error code and abort
C
                           IF (IERRCD .NE. INO) THEN
                              WRITE (IOUT,1405) QMIN(1), NFIXQ(1), 
     1                                          QMIN(5), NFIXQ(5)
                              GO TO 910
                           ENDIF
C
C     Do block boundary conditions (CUTS) 
C
                           DO 530 ICUT = 1, NCUTS
                              IBLK1 = ICUTS( 1,ICUT,ILVL)
                              IBLK2 = ICUTS(10,ICUT,ILVL)
                              IQ1   = IOFFQ(IBLK1,ILVL) * NQ    + 1
                              IQ2   = IOFFQ(IBLK2,ILVL) * NQ    + 1
                              CALL CUT 
     1                              (IDIM(IBLK1,ILVL), JDIM(IBLK1,ILVL),
     2                               KDIM(IBLK1,ILVL), Q(IQ1),
     3                               IDIM(IBLK2,ILVL), JDIM(IBLK2,ILVL),
     4                               KDIM(IBLK2,ILVL), Q(IQ2), 
     5                               ICUTS(1,ICUT,ILVL), NQ)
                              IF (ITURB .GE. ITLMNR) THEN
                                 IP1   = IOFFP(IBLK1,ILVL) * NP    + 1
                                 IP2   = IOFFP(IBLK2,ILVL) * NP    + 1
                                 CALL CUT 
     1                              (IDIM(IBLK1,ILVL), JDIM(IBLK1,ILVL),
     2                               KDIM(IBLK1,ILVL), PROPS(IP1),
     3                               IDIM(IBLK2,ILVL), JDIM(IBLK2,ILVL),
     4                               KDIM(IBLK2,ILVL), PROPS(IP2), 
     5                               ICUTS(1,ICUT,ILVL), NP)
                              ENDIF
  530                      CONTINUE
C
C     End of stage loop
C
  550                   CONTINUE
C
C     Output monitoring data on finest mesh level only
C
                        IF (ILVL .EQ. LVLONE) THEN
C
C     Calculate iteration for output
C
                           ITCUR     = ITER  + ITRUN
                           IT        = ITCUR + ITTOT
C
C     i.  Output monitoring data on subiterations
C
                           IF (ITSUB .NE. 1) THEN
C
C     Output residual data scaled by number of cell centers on current level
C
                              DO 570 LL = 1, NF + 1
                                 R2NORM(LL) = SQRT(R2NORM(LL)) /
     1                                        NBRPTS(ILVL)
                                 R2DQDT(LL) = SQRT(R2DQDT(LL)) /
     1                                        NBRPTS(ILVL)
  570                         CONTINUE
                              WRITE (IOUT,2200) ITCUR, ITSUB, 
     1                                          (R2NORM(LL),LL=1,NF+1)
                              IF (ITIMED .EQ. IUNSTD) THEN
                                 WRITE (IOUT,2210) ITCUR, ITSUB,
     1                                      (R2DQDT(LL),LL=1,NF+1)
                              ENDIF
 2200 FORMAT (' ','Iter =',I8,' Subiteration =',I4,8X,
     1            ' Res = ',13(1PE14.7,2X))
 2210 FORMAT (' ','Iter =',I8,' Subiteration =',I4,2X,
     1            ' L2(dQ/dt) = ',13(1PE14.7,2X))
C
C     ii. Output monitoring data only on main iterations
C
                           ELSE
C
C     Output residual data scaled by number of cell centers on current level
C
                              DO 580 LL = 1, NF + 1
                                 R2NORM(LL) = SQRT(R2NORM(LL)) /
     1                                        NBRPTS(ILVL)
                                 R2DQDT(LL) = SQRT(R2DQDT(LL)) /
     1                                        NBRPTS(ILVL)
  580                         CONTINUE
                              R2(IT) = R2NORM(NF+1)
                              IF (ITIMED .EQ. ISTDY) THEN
                                 WRITE (IOUT,2300) ITCUR, IT, 
     1                                      (R2NORM(LL),LL=1,NF+1)
                                 WRITE (IRESID,2301) ITCUR, IT, 
     1                                      (R2NORM(LL),LL=1,NF+1)
                              ELSE
                                 WRITE (IOUT,2305) ITCUR, IT, TIME,
     1                                      (R2NORM(LL),LL=1,NF+1)
                                 WRITE (IRESID,2306) ITCUR, IT, 
     1                                      (R2NORM(LL),LL=1,NF+1)
                                 WRITE (IOUT,2210) ITCUR, ITSUB,
     1                                      (R2DQDT(LL),LL=1,NF+1)
                              ENDIF
 2300 FORMAT (' ','ITER =',2I8,' RES = ',13(1PE14.7,2X))
 2301 FORMAT (2I8,13(1PE14.7,2X))
 2305 FORMAT (' ','ITER =',2I8,' TIME =',1PE12.5,
     1                         ' RES = ',13(1PE14.7,2X))
 2306 FORMAT (2I8,1PE12.5,13(1PE14.7,2X))
C
C     Output maximum residual
C
                              IF (IFRSMX .EQ. IYES) THEN
                                 RESMAX = SQRT (RESMAX)
                                 WRITE (IOUT,2350) ITCUR, IT, RESMAX, 
     1                                          (IRESMX(L),L=0,3)
                              ENDIF
 2350 FORMAT (' ','Iter =',2I8,' Max Res = ',1PE12.5,
     1            ' At BLOCK:I,J,K=(',I3,':',I4,',',I4,',',I4,')')
C
C     Output if reset turbulence production
C
                              IF (NPRNEG .GT. 0) THEN
                                 WRITE (IOUT,2400) NPRNEG
                              ENDIF
 2400 FORMAT (' ','MAIN  : SRCRES: Turb. Prod reset posit',
     1            'ive   at ',I6,' locations')
C
                              IF (NPRLIM .GT. 0) THEN
                                 WRITE (IOUT,2401) PRDLIM, NPRLIM
                              ENDIF
 2401 FORMAT (' ','MAIN  : SRCRES: Turb. Prod limited to ',1PE10.2,
     1            '*diss at ',I6,' locations')
C
C     Output number of times positivity violated
C
                              DO 585 L = 1, NQ
                                 IF (IFIXQ .EQ. IFIXMN) THEN
                                    FIXTYP = 'Minimum'
                                 ELSE IF (IFIXQ .EQ. IFIXAV) THEN
                                    FIXTYP = 'Average'
                                 ENDIF
                                 IF (NFIXQ(L) .GT. 0) THEN
                                    WRITE (IOUT,2410) L, QMIN(L), 
     1                                                NFIXQ(L), FIXTYP
                                 ENDIF
  585                         CONTINUE
 2410 FORMAT (' ','MAIN  : FIXQ  : Q(',I2,') < ',1PE9.2,
     1            ' at ',I7,' locations.  Replaced with ',A7,'.')
C
C     Output number of times realizability violated
C
                              IF (ITURB. EQ. ITRS) THEN
                                 IF (NRELIZ( 9) .GT. 0) THEN
                                    WRITE (IOUT,2420) 'XY', NRELIZ( 9)
                                 ENDIF
C
                                 IF (NRELIZ(10) .GT. 0) THEN
                                    WRITE (IOUT,2420) 'XZ', NRELIZ(10)
                                 ENDIF
C
                                 IF (NRELIZ(11) .GT. 0) THEN
                                    WRITE (IOUT,2420) 'YZ', NRELIZ(11)
                                 ENDIF
                              ENDIF
 2420 FORMAT (' ','MAIN  : REALIZ: Tau_',A2,' violates Schwartz ineq. ',
     1            ' at ',I6,' locations.  Reset.')
C
C     Output Monitor data
C
                              IF ( MOD(ITER, NITMON) .EQ. 0 ) THEN
                                 DO 595 IBLK = 1, NBLKS
                                    IQ = IOFFQ(IBLK,ILVL) * NQ    + 1
                                    CALL MONITR
     1                                (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     2                                 KDIM(IBLK,ILVL), Q(IQ),
     3                                 NPRMON(IBLK), IPRMON(1,1,IBLK))
  595                            CONTINUE
                              ENDIF
                           ENDIF
C
C     Calculate forces: update on each subiteration for 2D farfield vortex bc
C
                           IF (IFFORC .EQ. IYES) THEN
                              CALL FRCINT (ITURB,NBLKS,IDIM(1,ILVL), 
     1                               JDIM(1,ILVL),  KDIM(1,ILVL), 
     2                               IOFF(1,ILVL),  IOFFS(1,ILVL),
     3                               IOFFQ(1,ILVL), IOFFP(1,ILVL),
     4                               Q, PROPS, S, R, NXBCS, NBCS, 
     5                               IBCDAT(1,1,1,ILVL), FORCE, CL, CD,
     6                               FMOMENT)
                           ENDIF
C
C     Output only on LAST sub-iteration
C
                           IF (ITSUB .EQ. NSUBIT) THEN
                              IF (IFRCPR .EQ. IYES) THEN
                                 IF (ITIMED .EQ. ISTDY) THEN
                                    WRITE (IOUT,2480) ITCUR, IT, CL, CD
                                    WRITE (ICLCD,2486)ITCUR, IT, CL, CD,
     1                                 (FMOMENT(IMOM),IMOM=1,NMOMENT)
                                 ELSE
                                    WRITE (IOUT,2485) ITCUR, IT, TIME, 
     1                                                CL, CD
                                    WRITE (ICLCD,2486) ITCUR, IT, TIME,
     1                                              CL, CD, 
     2                                   (FMOMENT(IMOM),IMOM=1,NMOMENT)
                                 ENDIF
                              ENDIF
 2480 FORMAT (' ','Iter =',2I8,' CL = ',1PE13.6,' CD = ',1PE13.6)
 2485 FORMAT (' ','Iter =',2I8,' Time =',1PE12.5,
     1                         ' CL = ',1PE13.6,' CD = ',1PE13.6)
 2486 FORMAT (2I8,100E15.6)
                           ENDIF
                        ENDIF
C
C     End of output
C
  600                CONTINUE
  700             CONTINUE
C
C    Begin Prolongation from coarser meshes
C    1. For each block 
C       a. Interpolate coarse grid correction 
C       b. Update boundary conditions
C    2. Update Cuts
C
                  DO 750 ILVL = LVLONE + MGLVLS(LVLONE) - 2, LVLONE, -1
                     LVLC  = ILVL + 1
                     LVLF  = ILVL
                     DO 730 IBLK = 1, NBLKS
                        IQC   = IOFFQ (IBLK,LVLC) * NQ    + 1
                        IQCRS = IOFFQC(IBLK,LVLC) * NQ    + 1
                        IQF   = IOFFQ (IBLK,LVLF) * NQ    + 1
                        IPF   = IOFFP (IBLK,LVLF) * NP    + 1
                        ISF   = IOFFS (IBLK,LVLF) * 4 * 3 + 1
                        IVF   = IOFFS (IBLK,LVLF)         + 1
                        IRF   = IOFFS (IBLK,ILVL) * 3     + 1
C
C     Use RES as temporary storage for DQF (RES should be big enough)
C
                        CALL PROLNG (
     1                         IDIM(IBLK,LVLC), JDIM(IBLK,LVLC),
     2                         KDIM(IBLK,LVLC), Q(IQC), QCRS(IQCRS),
     3                         NBCS(IBLK), IBCDAT(1,1,IBLK,LVLC), 
     4                         IDIM(IBLK,LVLF), JDIM(IBLK,LVLF),
     5                         KDIM(IBLK,LVLF), Q(IQF), RES,
     6                         IBCDAT(1,1,IBLK,LVLF), 
     7                         ITURB, QMIN, NFIXQ, NRELIZ,
     8                         NXSECT, ABC, RTEMP)
C
                        IF (ITURB .GE. ITLMNR) THEN
                           NPPTS = (IDIM(IBLK,LVLF)+3) * 
     1                             (JDIM(IBLK,LVLF)+3) *
     2                             (KDIM(IBLK,LVLF)+3)
                        ELSE
                           NPPTS = NXPPTS
                        ENDIF
                        CALL BC     (
     1                         IDIM(IBLK,LVLF),JDIM(IBLK,LVLF),
     2                         KDIM(IBLK,LVLF), 1, IDIM(IBLK,LVLF)+1, 
     3                         Q(IQF),PROPS(IPF),NPPTS,S(ISF),VOL(IVF), 
     4                         R(IRF),NBCS(IBLK),IBCDAT(1,1,IBLK,LVLF), 
     5                         NXPROF, NPROFL, IOFQBC(1,LVLF), QBC,
     6                         VISCOS(1,IBLK), IMODEL, QP, NXSECT,
     7                         CL, RCNTR)
C
                        IF (ITURB .GE. ITLMNR) THEN
                           CALL PRPRTY (IDIM(IBLK,LVLF), 
     1                           JDIM(IBLK,LVLF), KDIM(IBLK,LVLF), 
     2                           1, IDIM(IBLK,LVLF)+1, NXSECT,
     3                           Q(IQF), PROPS(IPF), S(ISF), VOL(IVF),
     4                           IBCDAT(1,1,IBLK,LVLF), NBCS(IBLK), 
     5                           IMODEL, IFDELQ,
     6                           STEMP, QP, F, FM, QM, DQP, DQM, WORKX)
                        ENDIF
  730                CONTINUE
C
C     Update CUTS
C
                     DO 740 ICUT = 1, NCUTS
                        IBLK1 = ICUTS( 1,ICUT,LVLF)
                        IBLK2 = ICUTS(10,ICUT,LVLF)
                        IQ1   = IOFFQ(IBLK1,LVLF) * NQ    + 1
                        IQ2   = IOFFQ(IBLK2,LVLF) * NQ    + 1
                        CALL CUT (IDIM(IBLK1,LVLF), JDIM(IBLK1,LVLF), 
     1                            KDIM(IBLK1,LVLF), Q(IQ1),
     2                            IDIM(IBLK2,LVLF), JDIM(IBLK2,LVLF), 
     3                            KDIM(IBLK2,LVLF), Q(IQ2), 
     4                            ICUTS(1,ICUT,LVLF), NQ)
                        IF (ITURB .GE. ITLMNR) THEN
                           IP1   = IOFFP(IBLK1,LVLF) * NP    + 1
                           IP2   = IOFFP(IBLK2,LVLF) * NP    + 1
                           CALL CUT (IDIM(IBLK1,LVLF), JDIM(IBLK1,LVLF),
     1                               KDIM(IBLK1,LVLF), PROPS(IP1),
     2                               IDIM(IBLK2,LVLF), JDIM(IBLK2,LVLF),
     3                               KDIM(IBLK2,LVLF), PROPS(IP2), 
     4                               ICUTS(1,ICUT,LVLF), NP)
                        ENDIF
  740                CONTINUE
  750             CONTINUE
C
C    End   Prolongation
C
  790          CONTINUE
C
C    End   Sub-iterations; Q = Q^(n+1)
C
C    Output intermediate restart file - but not on the last iteration
C
               IF (ITER/NITRST*NITRST .EQ. ITER         .AND.
     1             ITER               .NE. NITS(LVLONE)) THEN
                  CALL WRREST (NXBLKS, NXLVLS, NXBCS,  NXSECT,
     2                         IOFFS,  IOFFQ,  IOFFP,  IOFFQN,
     3                         NBLKS,  NLVLS,  IMODEL, IFDELQ, ITER,
     4                         IDIM,   JDIM,   KDIM,
     5                         Q, QN, PROPS, S, VOL, NBCS, IBCDAT,
     6                         RES, QP, QM, DQP, DQM, F, FM,WORKX,STEMP,
     7                         FWRRST, IERRCD)
                  IF (IERRCD .NE. INO) GO TO 999
               ENDIF
  800       CONTINUE
C
C 12-c. If not finest grid, then interpolate to the next finer grid and 
C       update the properties and boundary conditions.
C
            IF (ICFS .NE. NCFS) THEN
               LVLC = NCFS - ICFS + 1
               LVLF = LVLC - 1
               WRITE (IOUT,8100) LVLC, LVLF
 8100 FORMAT (//' ','MAIN  : Interpolating from grid level ',I3,
     1              ' to grid level ',I3//)
C
C 12-c-i.   Interpolate coarse grid to fine grid for each block.
C           Modify cell face bc data to be at ghost cells then call 
C           interpolation routine.
C
               DO 850 IBLK = 1, NBLKS
                  IQC = IOFFQ(IBLK,LVLC) * NQ    + 1
                  IQF = IOFFQ(IBLK,LVLF) * NQ    + 1
                  CALL QBCMOD (IDIM(IBLK,LVLC), JDIM(IBLK,LVLC),
     1                         KDIM(IBLK,LVLC), Q(IQC),
     2                         NBCS(IBLK), IBCDAT(1,1,IBLK,LVLC))
                  CALL INT2HQ (IDIM(IBLK,LVLC), JDIM(IBLK,LVLC),
     1                         KDIM(IBLK,LVLC), Q(IQC),
     2                         IDIM(IBLK,LVLF), JDIM(IBLK,LVLF),
     3                         KDIM(IBLK,LVLF), Q(IQF))
                  CALL QBCFIX (IDIM(IBLK,LVLC), JDIM(IBLK,LVLC),
     1                         KDIM(IBLK,LVLC), Q(IQC),
     2                         NBCS(IBLK), IBCDAT(1,1,IBLK,LVLC))
  850          CONTINUE
C
C 12-c-ii.  Update boundary conditions and properties on each block
C
               DO 870 IBLK = 1, NBLKS
                  IQ = IOFFQ(IBLK,LVLF) * NQ    + 1
                  IP = IOFFP(IBLK,LVLF) * NP    + 1
                  IS = IOFFS(IBLK,LVLF) * 4 * 3 + 1
                  IV = IOFFS(IBLK,LVLF)         + 1
                  IR = IOFFS(IBLK,LVLF) * 3     + 1
C
                  IF (ITURB .GE. ITLMNR) THEN
                     NPPTS = (IDIM(IBLK,LVLF)+3) * (JDIM(IBLK,LVLF)+3) *
     1                       (KDIM(IBLK,LVLF)+3)
                     IPSTRT = 1
                     IPEND  = IDIM(IBLK,LVLF) + 1
C
C           Update properties array 
C
                     CALL PRPRTY (IDIM(IBLK,LVLF), JDIM(IBLK,LVLF), 
     1                            KDIM(IBLK,LVLF), IPSTRT, IPEND,NXSECT,
     2                            Q(IQ), PROPS(IP), S(IS), VOL(IV),
     3                            IBCDAT(1,1,IBLK,LVLF), NBCS(IBLK), 
     4                            IMODEL, IFDELQ,
     5                            STEMP, QP, F, FM, QM, DQP, DQM, WORKX)
                  ELSE
                     NPPTS = NXPPTS
                  ENDIF
C
C           Update boundary conditions 
C
                  CALL BC (IDIM(IBLK,LVLF),JDIM(IBLK,LVLF),
     1                  KDIM(IBLK,LVLF), 1, IDIM(IBLK,LVLF)+1, 
     2                  Q(IQ), PROPS(IP), NPPTS, S(IS), VOL(IV), R(IR),
     3                  NBCS(IBLK), IBCDAT(1,1,IBLK,LVLF), 
     4                  NXPROF, NPROFL, IOFQBC(1,LVLF), QBC,
     5                  VISCOS(1,IBLK), IMODEL, QP, NXSECT, CL, RCNTR)
 870           CONTINUE
C
C 12-c-iii. Do block boundary conditions (CUTS) 
C
               DO 890 ICUT = 1, NCUTS
                  IBLK1 = ICUTS( 1,ICUT,LVLF)
                  IBLK2 = ICUTS(10,ICUT,LVLF)
                  IQ1   = IOFFQ(IBLK1,LVLF) * NQ    + 1
                  IQ2   = IOFFQ(IBLK2,LVLF) * NQ    + 1
                  CALL CUT (IDIM(IBLK1,LVLF), JDIM(IBLK1,LVLF), 
     1                      KDIM(IBLK1,LVLF), Q(IQ1),
     2                      IDIM(IBLK2,LVLF), JDIM(IBLK2,LVLF), 
     3                      KDIM(IBLK2,LVLF), Q(IQ2), 
     4                      ICUTS(1,ICUT,LVLF), NQ)
                  IF (ITURB .GE. ITLMNR) THEN
                     IP1   = IOFFP(IBLK1,LVLF) * NP    + 1
                     IP2   = IOFFP(IBLK2,LVLF) * NP    + 1
                     CALL CUT (IDIM(IBLK1,LVLF), JDIM(IBLK1,LVLF), 
     1                         KDIM(IBLK1,LVLF), PROPS(IP1),
     2                         IDIM(IBLK2,LVLF), JDIM(IBLK2,LVLF), 
     3                         KDIM(IBLK2,LVLF), PROPS(IP2), 
     4                         ICUTS(1,ICUT,LVLF), NP)
                  ENDIF
  890          CONTINUE
            ENDIF
            ITRUN = ITRUN + NITS(LVLONE)
  900    CONTINUE
C
C     Jump here if STOP file is encountered
C
  905    CONTINUE
C
C     Calculate execution time of run
C

#ifdef CPU_TIME
C
         RUNTIME = DTIME (TARRAY)
C
         WRITE (IOUT,8200) RUNTIME
 8200 FORMAT (//' ','MAIN  : Approximate solver execution time: ',
     1              1PE12.4,' seconds ')
C
#endif

C
      ENDIF
C
C     Write global result
C
      CALL RESULT(ITCUR, R2NORM, CL, CD, APGINT, FMOMENT)
C
C     Calculate adverse pressure gradient integral
C
      ILVL = 1
      CALL APGRAD (ITURB,NBLKS,IDIM(1,ILVL), 
     1             JDIM(1,ILVL),  KDIM(1,ILVL), 
     2             IOFF(1,ILVL),  IOFFS(1,ILVL),
     3             IOFFQ(1,ILVL), IOFFP(1,ILVL),
     4             Q, S, R, NXBCS, NBCS, 
     5             IBCDAT(1,1,1,ILVL), APGINT)
C
C 13. OUTPUT RESTART DATA
C     Update properties.
C     If unsteady, convert QN to primitive (N.B. storage for tke is over nnpts)
C     Write restart file.
C
  910 CONTINUE
C
      CALL WRREST (NXBLKS, NXLVLS, NXBCS,  NXSECT,
     2             IOFFS,  IOFFQ,  IOFFP,  IOFFQN,
     3             NBLKS,  NLVLS,  IMODEL, IFDELQ, ITRUN,
     4             IDIM,   JDIM,   KDIM,
     5             Q, QN, PROPS, S, VOL, NBCS, IBCDAT,
     6             RES, QP, QM, DQP, DQM, F, FM, WORKX, STEMP,
     7             FWRRST, IERRCD)
      IF (IERRCD .NE. INO) GO TO 999
C
C 14. OUTPUT 
C     Write output and PLOT3D output
C     (PLOT3D output in multiple grid format)
C     Use RES array as temporary space to interpolate the grid
C     nodal values to cell centers.
C
      ILVL = 1
      CALL WRTP3D (IFPLT3, BINPLT, FPLT3G, FPLT3Q, FPLT3F, FPLT3FN,
     1             FERROR, NBLKS, IDIM, JDIM, KDIM, IOFFS, IOFFQ,
     2             IOFFP, ITURB, NPRINT, IPRINT, NCUTS, ICUTS,
     3             ILVL, Q, PROPS, R)
C
C Output skin friction and heat transfer data
C
      DO 982 IBLK = 1, NBLKS
         IQ = IOFFQ(IBLK,ILVL) * NQ    + 1
         IP = IOFFP(IBLK,ILVL) * NP    + 1
         IR = IOFFS(IBLK,ILVL) * 3     + 1
         IS = IOFFS(IBLK,ILVL) * 4 * 3 + 1
         IV = IOFFS(IBLK,ILVL)         + 1
C
         INCP   = (IDIM(IBLK,ILVL)+3)*(JDIM(IBLK,ILVL)+3)*
     1            (KDIM(IBLK,ILVL)+3)
         IP2 = IP  + INCP
         IP3 = IP2 + INCP
C
         DO 981 IBC = 1, NBCS(IBLK)
            IF ( (IBCDAT(1,IBC,IBLK,ILVL) .EQ. IWALL .OR.
     1            IBCDAT(1,IBC,IBLK,ILVL) .EQ. IWALFN ) ) THEN
               CALL SFHT  (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                     KDIM(IBLK,ILVL), 2, IDIM(IBLK,ILVL),ITURB,
     2                     Q(IQ), PROPS(IP), S(IS), VOL(IV),
     3                     IBCDAT(1,IBC,IBLK,ILVL), R(IR), IBLK)
            ENDIF
  981    CONTINUE
  982 CONTINUE
C
      GO TO 999
C
C 15.   Error Handling Section
C 15-a. Bad file name
C
  990 CONTINUE
      WRITE (IOUT,9990) FERROR
 9990 FORMAT (/' ','MAIN  : ERROR-> Error opening file ',A80,
     1        /' ', 8X,'RUN ABORTING!')
C
C 16.   Close output and finish
C
  999 CONTINUE
      CLOSE (IOUT)
C
C     Finished with ISAAC run.
C
      STOP
      END
