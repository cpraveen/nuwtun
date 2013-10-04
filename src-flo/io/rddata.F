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
      SUBROUTINE RDDATA (NXQ, NXCFS, NXLVLS, NXBLKS, NXBCS,NXCUTS,NXPRT,
     1                   NCFS, MGLVLS, NBLKS,  IDIM, JDIM, KDIM,
     2                   NBCS, IBCDAT, NCUTS, ICUTS, CUTNAM,
     3                   NPRINT, IPRINT, NPRMON, IPRMON, NITMON,
     4                   ITUPJ, ITSLVL,NITS,NSUBIT,NITFO,NITBBC,NITALG,
     5                   NITRST, NXSTAG, NSTAGE, RKALPH,
     6                   RESTRT, FRDRST, FWRRST, GRDFIL, IGRDTP, GRDBIN,
     7                   NXPROF, NPROFL, FPROFL, 
     8                   IFPLT3, BINPLT, FPLT3G, FPLT3Q, FPLT3F,FPLT3FN,
     9                   IORDER, LIMIT, RKAP, IFLUX, VISCOS,DIAG,IMODEL,
     A                   IMSTRT, IMEND, TOLER, VGNRN, SIGMA, IFDELQ,
     B                   CONTRB, IMODLX, QMIN, ITRSMX, IFFORC, IFRCPR, 
     C                   YNCALC, NXTRSG, NTRSEG, ITRSEG, RCNTR)
C
C Routine to set default values and read input data controlling run.
C
C This routine is organized as follows:
C     1. Set defaults
C     2. Read input deck
C     3. Set values based on input
C     4. Error handling section
C        a. End-Of-File Handler
C        b. Error Handler
C
C NXQ    : Maximum number of Q variables
C NXCFS  :
C NXLVLS :
C NXBLKS :
C NXBCS  :
C NXCUTS :
C NXPRT  :
C NCFS   :
C MGLVLS :
C NBLKS  :
C IDIM   :
C JDIM   :
C KDIM   :
C NBCS   :
C IBCDAT :
C NCUTS  :
C ICUTS  :
C NPRINT :
C IPRINT :
C NPRMON :
C IPRMON :
C NITMON :
C ITUPJ  :
C ITSLVL :
C NITS   : Number of iterations to perform this run.
C NSUBIT : Number of sub-iterations per time step.
C NITFO  : Number of iterations first order.  First order in all directions
C          for time dependent calculations.  First order in the streamwise
C          (marching) direction for marching calculations.
C NITBBC : Number of iterations of the blunt body condition - NOT USED.
C NITALG : Number of iterations to update the eddy viscosity using the
C          algebraic model.  Used for the initial iterations of higher
C          order models to assist convergence.  The higher order models
C          continue to be updated - but the eddy viscosity used for the
C          flow equations is calculated from the algebraic model.
C NITRST : Number of iterations to save intermediate restart file.
C NXSTAG : Maximum number of stages
C NSTAGE : Number of stages for temporal integration scheme
C          RK : Multistage scheme   NSTAGE = 3 (default or typical value)
C          AF : Single stage scheme NSTAGE = 1
C RKALPH : Coefficients for multistage scheme
C RESTRT :
C FRDRST :
C FWRRST :
C GRDFIL :
C IGRDTP :
C GRDBIN :
C NXPROF : Maximum number of profile boundary condition segments
C NPROFL : Number of profile boundary condition segments
C FPROFL : Array of filenames for profile boundary condition segments
C IFPLT3 :
C BINPLT :
C FPLTSG : File name for PLOT3D grid     file
C FPLT3Q : File name for PLOT3D solution file
C FPLT3F : File name for PLOT3D function file
C FPLT3FN: File name for PLOT3D function name file
C IORDER :
C LIMIT  :
C RKAP   :
C IFLUX  :
C VISCOS :
C DIAG   : Control of diagonalized scheme
C IMODEL : Turbulence model array
C IMSTRT :
C IMEND  :
C TOLER  :
C VGNRN  :
C SIGMA  :
C IFDELQ : Array controlling calculation of delq
C CONTRB : Logical.  If true then convert from one turbulence model to
C          another when do restart.
C IMODLX : Turbulence model array defining model to convert FROM if 
C          CONTRB is true.
C QMIN   : Control array for minimum value of Q
C ITRSMX : Number of iterations between output of maximum residual
C IFFORC : Control of calculation on forces
C IFRCPR : Control of output of forces
C YNCALC : Logical to force recalculation of Y_normal
C NXTRSG : Maximum number of transition specification segments
C NTRSEG : Actual  number of transition specification segments
C ITRSEG : Data for transition specification segments
C RCNTR  : Location of point vortex for 2D farfield bc. Also location to 
C          calculate moments about.
C
      include '../header/common.h'
      include '../header/histry.h'
C
      CHARACTER*80 FRDRST, FWRRST, GRDFIL(NXBLKS), 
     1             FPLT3G, FPLT3Q, FPLT3F, FPLT3FN, FPROFL(NXPROF)
      CHARACTER*80 FNAME
      CHARACTER*25 VAR, VAR1, BCTYP, SOLVER
C
      LOGICAL RESTRT, IFPLT3, GRDBIN, BINPLT, VGNRN, CONTRB, YNCALC
C
      DIMENSION IDIM   (NXBLKS,NXLVLS), 
     1          JDIM   (NXBLKS,NXLVLS),
     2          KDIM   (NXBLKS,NXLVLS)
      DIMENSION IBCDAT (10,NXBCS,NXBLKS,NXLVLS), 
     1          NBCS   (NXBLKS)
      DIMENSION ICUTS  (21,NXCUTS,NXLVLS)
      INTEGER   NCUTS
      CHARACTER*20 CUTNAM(NXCUTS)
      DIMENSION IPRINT (10,NXPRT,NXBLKS), 
     1          NPRINT (NXBLKS),
     2          IPRMON (10,NXPRT,NXBLKS), 
     3          NPRMON (NXBLKS)
      DIMENSION NITS   (NXCFS), 
     1          MGLVLS (NXCFS), 
     2          ITSLVL (NXLVLS,NXCFS)
      DIMENSION RKALPH (NXSTAG)
      DIMENSION IORDER (3,NXBLKS,NXLVLS,NXCFS), 
     1          LIMIT  (3,NXBLKS), 
     2          RKAP   (3,NXBLKS), 
     3          IFLUX  (NXBLKS)
      DIMENSION IMODEL (NMDL),
     1          IMODLX (NMDL)
      LOGICAL   VISCOS (3,NXBLKS),
     1          DIAG   (3,NXBLKS)
      DIMENSION IFDELQ (NXQ)
      DIMENSION QMIN   (NXQ)
C
C     Transition specification
C
      DIMENSION ITRSEG (7,NXTRSG)
C
      DIMENSION RCNTR  (3)
C
C     Local variables
C
      DIMENSION ISTRT (3), IEND(3), IORD(3)
      DIMENSION IDUM  (3)
      DIMENSION ITORDR(3)
C
C     Temporary storage to read print data into to avoid overwriting memory
C
      DIMENSION IDUMPR(10)
C
C 1.  Set defaults
C     Initialize INPERR to be 0.  This is incremented if there is an error
C     in the input dek so that the entire input deck is scanned for potential
C     errors without running an incorrect case.  INPERR counts the number
C     of errors.
C
      INPERR  = 0
C
C     Freestream Quantities
C
      RHOINF = 1.E0
      AINF   = 1.E0
C
C     Initialize variables to default values
C
      NQ     = 5
      NF     = 5
      NP     = 5
      NRANK  = 5
      NTIME  = 1
      NTMTAU = 0
C
      RESTRT = .FALSE.
      INITTD = .FALSE.
      VGNRN  = .FALSE.
      SIGMA  = 0.95E0
C
      IFPLT3 = .FALSE.
      BINPLT = .TRUE.
C
      FOURTH = .FALSE.
      NITBBC = 0
      THNLYR = .FALSE.
      GLOBAL = .TRUE.
      THREED = .TRUE.
      AXISYM = .FALSE.
      TWSPE  = .FALSE.
      SIUNIT = .FALSE.
      YNCALC = .FALSE.
C
      AXIDTH = 1.E0
C
C     Initialize number of profile boundary condition segments to zero
C
      NPROFL = 0
C
C     Grid type
C
      IGRDTP = IGTSEP
      GRDBIN = .TRUE.
C
C     Flow Variables
C
      GAMMA  = 1.4E0
      FSMACH = 0.8E0
      ALPHA  = 0.E0
      BETA   = 0.E0
      RE     = 1.E06
      PR     = 0.7E0
      TINF   = 450.E0
      TWALL  = 1060.E0
      AREARF = 1.E0
      ALENRF = 1.E0
C
C     Location of point vortex boundary condition and location for moments
C
      RCNTR(1) = 0.25E0
      RCNTR(2) = 0.E0
      RCNTR(3) = 0.E0
C
C     Subsonic inflow/outflow bc conditions.  Initialize to zero.  If input
C     does not reset to non-zero value, then set based on input freestream.
C     Input T_total as Kelvins or degrees R, 
C           P_back, P_total as p/p_inf 
C     and non-dimensionalize correctly at end of routine.
C
      PBAKBC = 0.E0
      PTOTBC = 0.E0
      TTOTBC = 0.E0
C
C     Jet conditions 
C     Velocities are input non-dimensionalized by ainf
C     Pressure is input as p(jet)/pinf 
C        - this is later converted to p(jet)/rhoinf*ainf**2
C     Temperature is input as T(jet)/Tinf
C
      UJET   = 0.E0
      VJET   = 0.E0
      WJET   = 0.E0
      PJET   = 1.E0
      TJET   = 1.E0
C
      TKEJET = 1.E-04
      EPSJET = 1.E-04
      OMEJET = 1.E-04
C
      TAUJET(1) = 2.E0 / 3.E0 * TKEJET
      TAUJET(2) = 2.E0 / 3.E0 * TKEJET
      TAUJET(3) = 2.E0 / 3.E0 * TKEJET
      TAUJET(4) = 0.E0
      TAUJET(5) = 0.E0
      TAUJET(6) = 0.E0
C
C     MUSCL Variables, diagonalized control
C
      DO 3 IBLK = 1, NXBLKS
         RKAP  (1,IBLK) = 0.E0
         RKAP  (2,IBLK) = 0.E0
         RKAP  (3,IBLK) = 0.E0
         LIMIT (1,IBLK) = ILMNMD
         LIMIT (2,IBLK) = ILMNMD
         LIMIT (3,IBLK) = ILMNMD
         IFLUX (IBLK)   = IFROE
         DIAG  (1,IBLK) = .TRUE.
         DIAG  (2,IBLK) = .TRUE.
         DIAG  (3,IBLK) = .TRUE.
    3 CONTINUE
C
C     Entropy Fix
C
      CNTRPY = 0.E0
C
C     Artificial dissipation for central difference scheme
C
      ADKP46 = 1.E0 / 64.E0
C
C     Constants for dQ/dt difference formula - default to first-order
C
      TDTHET = 1.E0
      TDPHI  = 0.E0
      TDPHIP = 0.E0
      ITMTYP = ITTS
C
C     Viscous/Inviscid
C
      DO 6 IBLK = 1, NXBLKS
         VISCOS(1,IBLK) = .FALSE.
         VISCOS(2,IBLK) = .FALSE.
         VISCOS(3,IBLK) = .FALSE.
    6 CONTINUE
C
C     Turn off all calculations of delq
C
      DO 7 L = 1, NXQ
         IFDELQ(L) = 0
    7 CONTINUE
C
C     Turbulence model - set to inviscid flow, no wall functions
C                        turbulence production limit = 100
C
      ITURB     = ITNVSD
      PRT       = 0.9E0
      CONTRB    = .FALSE.
      IFWALF    = .FALSE.
      IDAMP     = 0
      IEPSEQ    = 0
      IEPSLN    = 0
      ISRCW     = INO
      IADRM     = IADRNO
      IASM      = IASMBU
      IEPSC     = IECNO
      IPDIL     = IPDNO
      ICMUST    = 0
      PRDLIM    = 100.E0
C
C     Set default to no transition specification
C
      NTRSEG    = 0
C
C     Set default to limit the production on the epsilon/omega equation also
C
      PRDE      = 1.E0
      PRDEM1    = 1.E0 - PRDE
C
C     Set default to higher order advection on turbulence equations
C
      FOTURB    = .FALSE.
C
C     Set default to not enforce positivity of production
C
      POSPRD    = .FALSE.
C
C     Set default convert turbulence model to default turbulence model
C
      ITURBX    = ITURB
      IDAMPX    = IDAMP
      IASMX     = IASM
C
C     Storage locations for wall shear stress and heat transfer
C
      ITQTAU    = 3
      ITQDOT    = 4
C
C     Run Control Variables
C          NITS   <= MXHIST (maximum size of history array)
C          ITDFNL >= 2 to avoid divide by zero when evaluate (ITDFNL - 1)
C          NITRST  = MXHIST + 1 to turn off intermediate restart file output
C
C     Multigrid control
C
      NCFS      = 1
      DO 8 ICFS = 1, NXCFS
         MGLVLS(ICFS) = 1
         NITS  (ICFS) = MIN (200, MXHIST)
    8 CONTINUE
      CSMTH     = 0.25E0
      CMXCHG    = 0.25E0
C
      NITRST    = MXHIST + 1
      NSUBIT    = 1
      NSTAGE    = 1
      NITFO     = 0
      NITALG    = 0
      ISRCE     = .TRUE.
      ISOLVR    = IAF3F
      IMPSRC    = ISDIAG
      ITIMED    = ISTDY
      NBLKS     = 0
      NCUTS     = 0
      ITUPDT    = 1
      ITUPJ     = 1
      CFL       = 1.E0
      DT        = 0.E0
      CFLFNL    = 0.E0
      TIMEDT    = - 1.E0
      ITDFNL    = MAX (NITS(1), 2)
      ITRSMX    = MXHIST + 1
      IFFORC    = INO
      IFRCPR    = IYES
      NMOMENT   = 0
C
C     Control of FIX routines for positivity preservation
C
      IFIXQ     = IFIXMN
C
      QMIN( 1)  =   1.E-12
      QMIN( 2)  = - 1.E0
      QMIN( 3)  = - 1.E0
      QMIN( 4)  = - 1.E0
      QMIN( 5)  =   1.E-12
C
C     Output Control Variables
C
      NITMON    = 999999
      DO 9 IBLK = 1, NXBLKS
         NPRINT(IBLK) = 0
         NPRMON(IBLK) = 0
    9 CONTINUE
C
      FWRRST    = 'isaac.rst'
C
      DO 20 ICFS = 1, NXCFS
         DO 10 ILVL = 1, NXLVLS
            ITSLVL(ILVL,ICFS) = 1
   10    CONTINUE
   20 CONTINUE
C
C 2.  Read input deck
C     This is the beginning of the loop to read in and process the input
C
      WRITE (IOUT,1001)
C
C     Read card, set variable based on input
C
   50 CONTINUE
      READ (INPUT,*,END=990,ERR=995) VAR, VALUE
      WRITE (IOUT,1000) VAR, VALUE
      IF (VAR .EQ. 'END') THEN
         GO TO 900
C
C     Restart Control
C
      ELSE IF (VAR .EQ. 'RESTART IN') THEN
         RESTRT = .TRUE.
         READ (INPUT,*,END=990,ERR=995) FRDRST
         WRITE (IOUT,1025) FRDRST
      ELSE IF (VAR .EQ. 'RESTART OUT') THEN
         READ (INPUT,*,END=990,ERR=995) FWRRST
         WRITE (IOUT,1005) FWRRST
      ELSE IF (VAR .EQ. 'SAVE RESTART') THEN
         NITRST = INT(VALUE)
C
C     Plot Output
C
      ELSE IF (VAR .EQ. 'PLOT3D') THEN
         IFPLT3 = .TRUE.
         READ (INPUT,*,END=990,ERR=995) VAR
         WRITE (IOUT,1227) VAR
         IF (VAR .EQ. 'BINARY') THEN
            BINPLT = .TRUE.
         ELSE IF (VAR .EQ. 'FORMATTED') THEN
            BINPLT = .FALSE.
         ELSE
            WRITE (IOUT,1010) VAR
            INPERR = INPERR + 1
         ENDIF
         READ (INPUT,*,END=990,ERR=995) FNAME
         FPLT3G  = FNAME
         FPLT3Q  = FNAME
         FPLT3F  = FNAME
         FPLT3FN = FNAME
         CALL CHAPND (FPLT3G,  '.g')
         CALL CHAPND (FPLT3Q,  '.q')
         CALL CHAPND (FPLT3F,  '.qt')
         CALL CHAPND (FPLT3FN, '.qt.nam')
         IF (BINPLT) THEN
            CALL CHAPND (FPLT3G, '.unf')
            CALL CHAPND (FPLT3Q, '.unf')
            CALL CHAPND (FPLT3F, '.unf')
         ELSE
            CALL CHAPND (FPLT3G, '.fmt')
            CALL CHAPND (FPLT3Q, '.fmt')
            CALL CHAPND (FPLT3F, '.fmt')
         ENDIF
         WRITE (IOUT,1015) FPLT3G, FPLT3Q, FPLT3F
C
C     Solution Procedure
C
      ELSE IF (VAR .EQ. 'SOLVER') THEN
         READ (INPUT,*,END=990,ERR=995) SOLVER, VALUE
         IF (SOLVER .EQ. 'RK') THEN
            ISOLVR    = IRKN
            NSTAGE    = INT (VALUE)
            IF (NSTAGE .LE. 1) THEN
C
C     Default coefficients from Rick Gaffney optimized for kappa = 1/3 scheme
C     give maximum CFL = 1.76
C
               NSTAGE    = 3
               RKALPH(1) = 0.27E0
               RKALPH(2) = 0.50E0
               RKALPH(3) = 1.0E0
            ELSE
C
C     Read input coefficients
C
               DO 96 ISTAGE = 1, NSTAGE
                  READ (INPUT,*,END=990,ERR=995) RKALPH(ISTAGE)
   96          CONTINUE
            ENDIF
            IF (NSTAGE .GT. NXSTAG) THEN
               WRITE (IOUT,1020) NSTAGE, NXSTAG
               INPERR = INPERR + 1
            ENDIF
            WRITE (IOUT,1021) NSTAGE, (RKALPH(ISTAGE),ISTAGE=1,NSTAGE)
         ELSE IF (SOLVER .EQ. 'AF3F') THEN
            ISOLVR = IAF3F
            DO 98 IBLK = 1, NXBLKS
               DIAG(1,IBLK) = .FALSE.
               DIAG(2,IBLK) = .FALSE.
               DIAG(3,IBLK) = .FALSE.
   98       CONTINUE
         ELSE IF (SOLVER .EQ. 'DIAGONAL AF3F') THEN
            ISOLVR = IAF3F
            DO 99 IBLK = 1, NXBLKS
               DIAG(1,IBLK) = .TRUE.
               DIAG(2,IBLK) = .TRUE.
               DIAG(3,IBLK) = .TRUE.
   99       CONTINUE
         ELSE IF (SOLVER .EQ. 'LGS') THEN
            ISOLVR = ILGS
         ELSE IF (SOLVER .EQ. 'MARCH') THEN
            ISOLVR = IMARCH
            READ (INPUT,*,END=990,ERR=995) IMSTRT, IMEND
            WRITE (IOUT,1022) IMSTRT, IMEND
         ELSE
            WRITE (IOUT,1023) SOLVER
            INPERR = INPERR + 1
         ENDIF
      ELSE IF (VAR .EQ. 'INITIALIZE TIME DEPENDENT') THEN
         INITTD = .TRUE.
      ELSE IF (VAR .EQ. 'IMPLICIT TIME DEPENDENT') THEN
         ITIMED = IUNSTD
         TDTHET = 1.0E0
         TDPHI  = 0.5E0
         TDPHIP = 0.E0
         NTIME  = 2
         NSUBIT = 2
         DO 53 IIMP = 1, INT(VALUE)
            READ (INPUT,*,END=990,ERR=995) VAR, VALUE
            WRITE (IOUT,1061) VAR, VALUE
            IF (VAR .EQ. 'TIME STEP') THEN
               TIMEDT = VALUE
            ELSE IF (VAR .EQ. 'SUB ITERATIONS') THEN
               NSUBIT = INT(VALUE)
            ELSE IF (VAR .EQ. 'T-TS') THEN
               ITMTYP = ITTS
            ELSE IF (VAR .EQ. 'TAU-TS') THEN
               ITMTYP = ITAUTS
               TDPHIP = 0.E0
               NTMTAU = 0
            ELSE IF (VAR .EQ. 'TAU-TS 2ND ORDER') THEN
               ITMTYP = ITAUTS
               TDPHIP = 0.5E0
               NTMTAU = 1
            ELSE
               WRITE (IOUT,1029) VAR
               INPERR = INPERR + 1
            ENDIF
   53    CONTINUE
      ELSE IF (VAR .EQ. 'CFL') THEN
         CFL = VALUE
      ELSE IF (VAR .EQ. 'DT') THEN
         DT  = VALUE
      ELSE IF (VAR .EQ. 'CFLFNL') THEN
         CFLFNL = VALUE
      ELSE IF (VAR .EQ. 'ITUPDT') THEN
         ITUPDT = INT(VALUE)
      ELSE IF (VAR .EQ. 'ITDFNL') THEN
         ITDFNL = MAX (INT(VALUE), 2)
      ELSE IF (VAR .EQ. 'ITUPJ') THEN
         ITUPJ = INT(VALUE)
      ELSE IF (VAR .EQ. 'NITS') THEN
         NITS(1) = INT(VALUE)
         IF (NITS(1) .GT. MXHIST) THEN
            WRITE (IOUT,1024) NITS(1), MXHIST
            NITS(1) = MXHIST
         ENDIF
      ELSE IF (VAR .EQ. 'ENTROPY FIX') THEN
         CNTRPY = VALUE
         IF (CNTRPY .LT. 0.E0 .OR. CNTRPY .GT. 1.E0) THEN
            WRITE (IOUT,1124) 
            INPERR = INPERR + 1
         ENDIF
      ELSE IF (VAR .EQ. 'PRINT MAX RESIDUAL') THEN
         ITRSMX = INT(VALUE)
      ELSE IF (VAR .EQ. 'RECALCULATE YNORMAL') THEN
         ITEMP  = INT(VALUE)
         IF (ITEMP .EQ. 0) THEN
            YNCALC = .FALSE.
         ELSE
            YNCALC = .TRUE.
         ENDIF
      ELSE IF (VAR .EQ. 'PRINT FORCES') THEN
         ITEMP  = INT(VALUE)
         IF (ITEMP .EQ. 0) THEN
            IFRCPR = INO
         ELSE IF (ITEMP .EQ. 1) THEN
            IFRCPR = IYES
         ELSE
            WRITE (IOUT,1026) ITEMP
            INPERR = INPERR + 1
         ENDIF
      ELSE IF (VAR .EQ. 'REFERENCE AREA') THEN
         AREARF = VALUE
      ELSE IF (VAR .EQ. 'REFERENCE LENGTH') THEN
         ALENRF = VALUE
      ELSE IF (VAR .EQ. 'MOMENT') THEN
         NMOMENT = INT(VALUE)
         IF(NMOMENT.LE.0)THEN
            WRITE (IOUT,'(1X,"NMOMENT IS NEGATIVE")')
            INPERR = INPERR + 1
         ENDIF
         DO IMOM=1,NMOMENT
            READ (INPUT,*,END=990,ERR=995) AXISNAME(IMOM),
     1         AMOMAXIS(IMOM,1,:), AMOMAXIS(IMOM,2,:)
            WRITE(IOUT,'(10X,A15,3E12.4)') AXISNAME(IMOM),
     1         AMOMAXIS(IMOM,1,:)
            WRITE(IOUT,'(25X,3E12.4)') AMOMAXIS(IMOM,2,:)
            ALEN = SQRT((AMOMAXIS(IMOM,2,1)-AMOMAXIS(IMOM,1,1))**2 +
     1                  (AMOMAXIS(IMOM,2,2)-AMOMAXIS(IMOM,1,2))**2 +
     2                  (AMOMAXIS(IMOM,2,3)-AMOMAXIS(IMOM,1,3))**2)
            AMOMAXIS(IMOM,3,:)=(AMOMAXIS(IMOM,2,:)-AMOMAXIS(IMOM,1,:))/
     1                         ALEN
         ENDDO
      ELSE IF (VAR .EQ. 'MULTIGRID LEVELS') THEN
         DO 54 ICFS = 1, NXCFS
            MGLVLS(ICFS) = MAX ((INT(VALUE) + 1 - ICFS), 1)
   54    CONTINUE
      ELSE IF (VAR .EQ. 'IMPLICIT SMOOTHING') THEN
         CSMTH  = VALUE
      ELSE IF (VAR .EQ. 'MESH SEQUENCE') THEN
         NCFS    = INT(VALUE) + 1
         IF (NCFS .GT. NXCFS) THEN
            WRITE (IOUT,1210) NCFS, NXCFS, NCFS-1
            INPERR = INPERR + 1
C     Skip input lines if error
            DO 55 ICFS = NCFS, 2, -1
               READ (INPUT,*,END=990,ERR=995) 
   55       CONTINUE
         ELSE
            DO 60 ICFS = NCFS, 2, -1
               READ (INPUT,*,END=990,ERR=995) VALUE
               NITS  (ICFS) = INT(VALUE)
               WRITE (IOUT,1220) ICFS, NITS(ICFS)
               IF (NITS(ICFS) .GT. MXHIST) THEN
                  WRITE (IOUT,1024) NITS(ICFS), MXHIST
                  NITS(ICFS) = MXHIST
               ENDIF
   60       CONTINUE
         ENDIF
      ELSE IF (VAR .EQ. 'NITFO') THEN
         NITFO  = INT(VALUE)
      ELSE IF (VAR .EQ. 'FIRST ORDER TURBULENCE') THEN
         FOTURB = .TRUE.
      ELSE IF (VAR .EQ. 'ENFORCE POSITIVE PROD.') THEN
         POSPRD = .TRUE.
      ELSE IF (VAR .EQ. 'NITALG') THEN
         NITALG = INT(VALUE)
      ELSE IF (VAR .EQ. 'SOURCE') THEN
         ITEMP = INT(VALUE)
         IF (ITEMP .EQ. 0) THEN
            ISRCE = .FALSE.
         ELSE IF (ITEMP .EQ. 1) THEN
            ISRCE = .TRUE.
         ELSE
            WRITE (IOUT,1027) ITEMP
            INPERR = INPERR + 1
         ENDIF
      ELSE IF (VAR .EQ. 'IMPLICIT SOURCE') THEN
         READ (INPUT,*,END=990,ERR=995) VAR
         WRITE (IOUT,1065) VAR
         IF (VAR .EQ. 'NONE') THEN
            IMPSRC = INO
         ELSE IF (VAR .EQ. 'DIAGONAL') THEN
            IMPSRC = ISDIAG
         ELSE IF (VAR .EQ. 'BLOCK') THEN
            IMPSRC = ISBLOC
         ELSE
            INPERR = INPERR + 1
            WRITE (IOUT,1028) VAR
         ENDIF
C
C     Control of FIX positivity routines
C
      ELSE IF (VAR .EQ. 'FIXQ') THEN
         READ (INPUT,*,END=990,ERR=995) VAR
         WRITE (IOUT,1227) VAR
         IF (VAR .EQ. 'NONE') THEN
            IFIXQ  = INO
         ELSE IF (VAR .EQ. 'MINIMUM') THEN
            IFIXQ  = IFIXMN
         ELSE IF (VAR .EQ. 'AVERAGE') THEN
            IFIXQ  = IFIXAV
         ELSE 
            WRITE (IOUT,1229) VAR
            INPERR = INPERR + 1
         ENDIF
C
C     Input transition specification segments
C
      ELSE IF (VAR .EQ. 'ENFORCE LAMINAR FLOW') THEN
         NTRSEG = INT (VALUE)
         IF (NTRSEG .GT. NXTRSG) THEN
            WRITE (IOUT,1545) NTRSEG, NXTRSG
            INPERR = INPERR + 1
            GO TO 900
         ENDIF
C
C      Input each segment specification; Check that ending index is >=
C      beginning index; Switch from cell node to cell center.
C
         DO ISEG = 1, NTRSEG
            READ (INPUT,*,END=990,ERR=995) IBLOCK, (ISTRT(L),L=1,3), 
     1                                             (IEND (M),M=1,3)
            WRITE (IOUT,1580) IBLOCK, (ISTRT(L),L=1,3),
     1                                (IEND (M),M=1,3)
            DO L = 1, 3
               IF (IEND(L) .LT. ISTRT(L)) THEN
                  WRITE (IOUT,1590) L, IBLOCK, (ISTRT(LL),LL=1,3),
     1                                         (IEND (MM),MM=1,3)
                  INPERR = INPERR + 1
               ENDIF
            END DO
            ITRSEG(1,ISEG) = IBLOCK
            ITRSEG(2,ISEG) = ISTRT(1) + 1
            ITRSEG(3,ISEG) = ISTRT(2) + 1
            ITRSEG(4,ISEG) = ISTRT(3) + 1
            ITRSEG(5,ISEG) = IEND (1)
            ITRSEG(6,ISEG) = IEND (2)
            ITRSEG(7,ISEG) = IEND (3)
         END DO
C
C     Control of Marching Algorithm
C
      ELSE IF (VAR .EQ. 'VIGNERON') THEN
         VGNRN  = .TRUE.
      ELSE IF (VAR .EQ. 'TOLER') THEN
         TOLER = VALUE
      ELSE IF (VAR .EQ. 'SIGMA') THEN
         SIGMA = VALUE
c     ELSE IF (VAR .EQ. 'BLUNT BODY') THEN
c        BLUNT  = .TRUE.
c        NITBBC = INT(VALUE)
C
C     Viscous Terms
C
      ELSE IF (VAR .EQ. 'THIN LAYER') THEN
         THNLYR = .TRUE.
      ELSE IF (VAR .EQ. 'NAVIER STOKES') THEN
         THNLYR = .FALSE.
C
C     Turbulence Models
C
      ELSE IF (VAR .EQ. 'BALDWIN LOMAX') THEN
         ITURB  = ITBLMX
C
C     SIGK and SIGT2 are defined as the inverse of sigma(k) and sigma(epsilon)
C
      ELSE IF (VAR .EQ. 'K EPSILON') THEN
C
C     Set default k-epsilon model to be Zhang, So, Gatski, and Speziale
C
C     SAA : Set k-epsilon turbulence constants from Speziale, Abid 
C          and Anderson AIAA 90-1481 as
C             sigma k       = 1.36
C             sigma k       = 1.00 - CORRECTED TO THIS 
C                             (Abid had recommended 0.75 for supersonic flows)
C             sigma epsilon = 1.36
C             C2 epsilon    = 1.83
C        ICASE 90-9
C             alf1 kepsilon = 1.0 without pressure dilatation
C             alf1 kepsilon = 0.6    with pressure dilatation (Rick Gaffney)
C        ICASE 90-42
C             alf2 kepsilon = 0.4 
C             alf3 kepsilon = 0.2
C        CPDFRO,CPDFT = 0 as these models are not good
C     ZSGS: Zhang, So, Gatski, Speziale in Near-Wall Turbulent Flows (1993)
C     RNG : Yakhot, Orszag, Thangam, Gatski, Speziale in Phys. Fluids A
C
C     Default values
C
         ITURB  = ITKE
         IEPSEQ = IEEZSG
         IDAMP  = IFMZSG
         IADRM  = IADRNO
         IASM   = IASMBU
         IEPSC  = IECNO
         IPDIL  = IPDNO
         ICMUST = 0
         IEPSLN = 7
         IFDL2K = 1
         NQ     = 7
         NP     = 5
         NF     = 7
         NRANK  = 7
         SIGK   = 1.0E0 / 1.00E0
         SIGT2  = 1.0E0 / 1.45E0
         SIGRHO = 1.0E0 / 0.95E0
         CMU    = 0.096E0
         CEPS1  = 1.50E0
         CEPS2  = 1.83E0
         A2KE   = 4.9E0
         ALF1   = 0.0E0
         ALF2   = 0.0E0
         ALF3   = 0.0E0
C--->    CPDFRO = SIGRHO
C--->    CPDFT  = 1.E0 / PRT
         CPDFRO = 0.E0
         CPDFT  = 0.E0
         CPVELK = 0.E0
C
         PRDLIM = 100.E0
C
C     Freestream levels
C     Set TKEINF to (turbulence intensity=u'/U)^2; Non-dimensionalize later
C         RMUTNF to mu_T / mu
C
         TKEINF = 1.0E-04 
         RMUTNF = 1.0E0
C
         QMIN( 1)  =   1.E-8
         QMIN( 2)  = - 1.E0
         QMIN( 3)  = - 1.E0
         QMIN( 4)  = - 1.E0
         QMIN( 5)  =   1.E-8
         QMIN( 6)  =   1.E-8
         QMIN( 7)  =   1.E-8
C
C     Turn on dui/dxj calculation
C
         IFDELQ(2) = 1
         IFDELQ(3) = 1
         IFDELQ(4) = 1
C
C     Read in values of constants
C
         NCON   = INT (VALUE)
         DO 101 ICON = 1, NCON
            READ (INPUT,*,END=990,ERR=995) VAR, VALUE
            WRITE (IOUT,1061) VAR, VALUE
            IF (VAR .EQ. 'SAA') THEN
               IEPSEQ = IEESAA
               IDAMP  = IFMSAA
               IFDL2K = 0
               SIGK   = 1.0E0 / 1.00E0
               SIGT2  = 1.0E0 / 1.36E0
               SIGRHO = 1.0E0 / 0.95E0
               CMU    = 0.09E0
               CEPS1  = 1.44E0
               CEPS2  = 1.83E0
               A2KE   = 4.9E0
               ALF1   = 0.6E0
               ALF2   = 0.4E0
               ALF3   = 0.2E0
               CPDFRO = 0.E0
               CPDFT  = 0.E0
               CPVELK = 0.E0
            ELSE IF (VAR .EQ. 'ZSGS') THEN
               IEPSEQ = IEEZSG
               IDAMP  = IFMZSG
               IFDL2K = 1
               SIGK   = 1.0E0 / 1.00E0
               SIGT2  = 1.0E0 / 1.45E0
               SIGRHO = 1.0E0 / 0.95E0
               CMU    = 0.096E0
               CEPS1  = 1.50E0
               CEPS2  = 1.83E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               CPDFRO = 0.E0
               CPDFT  = 0.E0
               CPVELK = 0.E0
            ELSE IF (VAR .EQ. 'ZSSL') THEN
               IEPSEQ = IEEZSL
               IDAMP  = IFMZSG
               IFDL2K = 1
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.45E0
               SIGRHO = 1.0E0 / 0.50E0
               CMU    = 0.096E0
               CEPS1  = 1.50E0
               CEPS2  = 1.83E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               CPDFRO = 0.E0
               CPDFT  = 0.E0
               CPVELK = 0.E0
            ELSE IF (VAR .EQ. 'HIGH RE') THEN
               IEPSEQ = IEEHR
               IDAMP  = IFMHR
               IFDL2K = 0
               SIGK   = 1.0E0 / 1.00E0
               SIGT2  = 1.0E0 / 1.30E0
               SIGRHO = 1.0E0 / 0.95E0
               CMU    = 0.09E0
               CEPS1  = 1.44E0
               CEPS2  = 1.92E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               CPDFRO = 0.E0
               CPDFT  = 0.E0
               CPVELK = 0.E0
            ELSE IF (VAR .EQ. 'RNG') THEN
               IEPSEQ = IEERNG
               IDAMP  = IFMHR
               IFDL2K = 0
               SIGK   = 1.39E0
               SIGT2  = 1.39E0
               SIGRHO = 1.0E0 / 0.95E0
               CMU    = 0.085E0
               CEPS1  = 1.42E0
               CEPS2  = 1.68E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               CPDFRO = 0.E0
               CPDFT  = 0.E0
               CPVELK = 0.E0
            ELSE IF (VAR .EQ. 'ASM') THEN
               IASM   = IASMGS
               ICMUST = 5
               IEPSEQ = IEEABD
               IDAMP  = IFMHR
               IFDL2K = 0
               SIGK   = 1.0E0 / 1.00E0
               SIGRHO = 0.0E0
               CMU    = 0.081E0
               CEPS1  = 1.44E0
               CEPS2  = 1.83E0
               A2KE   = 5.5E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               CPDFRO = 0.E0
               CPDFT  = 0.E0
               CPVELK = 0.E0
               SIGT2  = 1.0E0 / ( 0.41e0**2/((CEPS2-CEPS1)*SQRT(CMU)) )
            ELSE IF (VAR .EQ. 'ADRM') THEN
               IADRM  = IADRGS
               IEPSEQ = IEEHR
               IDAMP  = IFMHR
               IFDL2K = 0
               SIGK   = 1.0E0 / 1.00E0
               SIGT2  = 1.0E0 / 1.87E0
               SIGRHO = 0.0E0
               CMU    = 0.094E0
               CEPS1  = 1.20E0
               CEPS2  = 1.83E0
               A2KE   = 0.0E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               CPDFRO = 0.E0
               CPDFT  = 0.E0
               CPVELK = 0.E0
            ELSE IF (VAR .EQ. 'SIGMAK') THEN
               SIGK   = 1.0E0 / VALUE
            ELSE IF (VAR .EQ. 'SIGMAE') THEN
               SIGT2  = 1.0E0 / VALUE
            ELSE IF (VAR .EQ. 'SIGMARHO') THEN
               SIGRHO = 1.0E0 / VALUE
               CPDFRO = SIGRHO
            ELSE IF (VAR .EQ. 'CMU') THEN
               CMU    = VALUE
            ELSE IF (VAR .EQ. 'CEPS1') THEN
               CEPS1  = VALUE
            ELSE IF (VAR .EQ. 'CEPS2') THEN
               CEPS2  = VALUE
            ELSE IF (VAR .EQ. 'A2KE') THEN
               A2KE   = VALUE
            ELSE IF (VAR .EQ. 'ALF1') THEN
               ALF1   = VALUE
            ELSE IF (VAR .EQ. 'ALF2') THEN
               ALF2   = VALUE
            ELSE IF (VAR .EQ. 'ALF3') THEN
               ALF3   = VALUE
            ELSE IF (VAR .EQ. 'PRT') THEN
               PRT    = VALUE
               CPDFT  = 1.E0 / PRT
            ELSE IF (VAR .EQ. 'CPDFRO') THEN
               CPDFRO = VALUE
            ELSE IF (VAR .EQ. 'CPDFT') THEN
               CPDFT  = VALUE
            ELSE IF (VAR .EQ. 'CPVELK') THEN
               CPVELK = VALUE
            ELSE IF (VAR .EQ. 'PRDLIM') THEN
               PRDLIM = VALUE
            ELSE IF (VAR .EQ. 'FULL EPSILON PROD') THEN
               PRDE   = 0.E0
               PRDEM1 = 1.E0 - PRDE
            ELSE IF (VAR .EQ. 'TKEINF') THEN
               TKEINF = VALUE
            ELSE IF (VAR .EQ. 'RMUTNF') THEN
               RMUTNF = VALUE
            ELSE IF (VAR .EQ. 'MINIMUM TKE') THEN
               QMIN(6) = VALUE
            ELSE IF (VAR .EQ. 'MINIMUM EPSILON') THEN
               QMIN(7) = VALUE
            ELSE IF (VAR .EQ. 'COMPRESSIBLE DISS') THEN
               READ (INPUT,*,END=990,ERR=995) VAR
               WRITE (IOUT,1065) VAR
               IF (VAR .EQ. 'NONE') THEN
                  IEPSC  = IECNO
               ELSE IF (VAR .EQ. 'SARKAR') THEN
                  IEPSC  = IECSAR
               ELSE IF (VAR .EQ. 'RISTORCELLI') THEN
                  IEPSC  = IECRIS
               ELSE
                  INPERR = INPERR + 1
                  WRITE (IOUT,1068) VAR
               ENDIF
            ELSE IF (VAR .EQ. 'PRESSURE DILATATION') THEN
               READ (INPUT,*,END=990,ERR=995) VAR
               WRITE (IOUT,1065) VAR
               IF (VAR .EQ. 'NONE') THEN
                  IPDIL  = IPDNO
               ELSE IF (VAR .EQ. 'SARKAR') THEN
                  IPDIL  = IPDSAR
               ELSE IF (VAR .EQ. 'RISTORCELLI') THEN
                  IPDIL  = IPDRIS
               ELSE
                  INPERR = INPERR + 1
                  WRITE (IOUT,1069) VAR
               ENDIF
            ELSE
               WRITE (IOUT,1029) VAR
               INPERR = INPERR + 1
            ENDIF
  101    CONTINUE
C
C     These variables are used only for output
C
         IF (IEPSEQ .EQ. IEESAA) THEN
            WRITE (IOUT,1070) 'Speziale-Abid-Anderson        '
         ELSE IF (IEPSEQ .EQ. IEEHR ) THEN
            WRITE (IOUT,1070) 'High Reynolds Number          '
         ELSE IF (IEPSEQ .EQ. IEEZSG) THEN
            WRITE (IOUT,1070) 'Zhang-So-Gatski-Speziale      '
         ELSE IF (IEPSEQ .EQ. IEEZSL) THEN
            WRITE (IOUT,1070) 'Zhang-So-Speziale-Lai         '
         ELSE IF (IEPSEQ .EQ. IEERNG) THEN
            WRITE (IOUT,1070) 'RNG (Yak-Ors-Than-Gat-Spez)   '
         ELSE IF (IEPSEQ .EQ. IEEABD) THEN
            WRITE (IOUT,1070) 'Abid ASM (Abid-Morr-Gat-Spez) '
         ENDIF
         SK     = 1.E0 / SIGK
         ST2    = 1.E0 / SIGT2
         WRITE (IOUT,1075) CEPS1, CEPS2, CMU, SK, ST2, SIGRHO, A2KE,
     1                     ALF1, ALF2, ALF3, CPDFRO, CPDFT, CPVELK,
     2                     TKEINF, RMUTNF
C
C     k-omega constants from Wilcox, AIAA Journal, Vol. 26, No. 11
C
      ELSE IF (VAR .EQ. 'K OMEGA') THEN
C
C     Default values
C
         ITURB  = ITKW
         IEPSEQ = 0
         IDAMP  = IFMHR
         IADRM  = IADRNO
         IASM   = IASMBU
         IEPSC  = IECNO
         IPDIL  = IPDNO
         ICMUST = 0
         IEPSLN = 7
         IFDL2K = 0
         NQ     = 7
         NP     = 5
         NF     = 7
         NRANK  = 7
         SIGK   = 0.5E0
         SIGT2  = 0.5E0
         SIGRHO = 0.0E0
         BSTRKW = 0.09E0
         BKW    = 3.0E0 / 40.0E0
         CMU    = 1.0E0
         GKW    = 5.0E0 / 9.0E0
         ROUGHK = 0.0E0
         CPDFRO = 0.0E0
         CPDFT  = 0.0E0
         CPVELK = 0.0E0
C
         PRDLIM = 100.E0
C
C     Freestream levels
C     Set TKEINF to (turbulence intensity=u'/U)^2; Non-dimensionalize later
C         RMUTNF to mu_T / mu
C
         TKEINF = 1.0E-06 
         RMUTNF = 1.0E-01
C
         QMIN( 1)  =   1.E-8
         QMIN( 2)  = - 1.E0
         QMIN( 3)  = - 1.E0
         QMIN( 4)  = - 1.E0
         QMIN( 5)  =   1.E-8
         QMIN( 6)  =   1.E-12
         QMIN( 7)  =   1.E-8
C
C     Turn on dui/dxj calculation
C
         IFDELQ(2) = 1
         IFDELQ(3) = 1
         IFDELQ(4) = 1
C
C     Read in values of constants
C
         NCON   = INT (VALUE)
         DO 102 ICON = 1, NCON
            READ (INPUT,*,END=990,ERR=995) VAR, VALUE
            WRITE (IOUT,1061) VAR, VALUE
            IF (VAR .EQ. 'ASM') THEN
               IASM   = IASMGS
               ICMUST = 5
               BSTRKW = 1.0E0
               CMU    = 0.088E0
               SIGK   = 1.E0 / 1.4E0
               SIGT2  = 1.E0 / 2.0E0
               BKW    = 0.83E0
               GKW    = BKW - (0.41E0)**2 / (SQRT(CMU) / SIGT2)
            ELSE IF (VAR .EQ. 'SIGMAK') THEN
               SIGK   = VALUE
            ELSE IF (VAR .EQ. 'SIGMAW') THEN
               SIGT2  = VALUE
            ELSE IF (VAR .EQ. 'SIGMARHO') THEN
               SIGRHO = 1.0E0 / VALUE
               CPDFRO = SIGRHO
            ELSE IF (VAR .EQ. 'CMU') THEN
               CMU    = VALUE
            ELSE IF (VAR .EQ. 'BSTRKW') THEN
               BSTRKW = VALUE
            ELSE IF (VAR .EQ. 'BKW') THEN
               BKW    = VALUE
            ELSE IF (VAR .EQ. 'GKW') THEN
               GKW    = VALUE
            ELSE IF (VAR .EQ. 'ROUGHK') THEN
               ROUGHK = VALUE
            ELSE IF (VAR .EQ. 'PRT') THEN
               PRT    = VALUE
               CPDFT  = 1.E0 / PRT
            ELSE IF (VAR .EQ. 'CPDFRO') THEN
               CPDFRO = VALUE
            ELSE IF (VAR .EQ. 'CPDFT') THEN
               CPDFT  = VALUE
            ELSE IF (VAR .EQ. 'CPVELK') THEN
               CPVELK = VALUE
            ELSE IF (VAR .EQ. 'PRDLIM') THEN
               PRDLIM = VALUE
            ELSE IF (VAR .EQ. 'FULL OMEGA PROD') THEN
               PRDE   = 0.E0
               PRDEM1 = 1.E0 - PRDE
            ELSE IF (VAR .EQ. 'TKEINF') THEN
               TKEINF = VALUE
            ELSE IF (VAR .EQ. 'RMUTNF') THEN
               RMUTNF = VALUE
            ELSE IF (VAR .EQ. 'MINIMUM TKE') THEN
               QMIN(6) = VALUE
            ELSE IF (VAR .EQ. 'MINIMUM OMEGA') THEN
               QMIN(7) = VALUE
            ELSE
               WRITE (IOUT,1029) VAR
               INPERR = INPERR + 1
            ENDIF
  102    CONTINUE
C
C     These variables are used only for output
C
         WRITE (IOUT,1080) SIGK, SIGT2, BSTRKW, BKW, CMU, GKW, ROUGHK,
     1                     SIGRHO, CPDFRO, CPDFT, CPVELK,
     2                     TKEINF, RMUTNF
C
C     Reynolds stress model
C
      ELSE IF (VAR .EQ. 'REYNOLDS STRESS') THEN
C
C     Set Reynolds Stress turbulence constants 
C     AIAA 90-1481 as
C          sigma k       =  1.36
C          sigma k       =  0.75 CORRECTED TO THIS
C          sigma epsilon =  1.36
C          C  epsilon    =  0.10
C          C1 epsilon    =  1.50
C          C2 epsilon    =  1.83
C          C3 epsilon    = -2.00
C          C4 epsilon    = -1.50
C     ICASE 90-9
C          alf1 kepsilon = 1.0 without pressure dilatation
C          alf1 kepsilon = 0.6    with pressure dilatation (Rick Gaffney)
C          alf1          = 0   for Lai-So
C     ICASE 90-42
C          alf2          = 0   for Lai-So
C          alf3          = 0   for Lai-So
C          alf2 kepsilon = 0.4 
C          alf3 kepsilon = 0.2
C     CPDFRO,CPDFT = 0 as these models are not good
C     Pressure-Strain Correlation
C       Lai-So
C          cpi1 = 3.0
C          cpi2 = 0.6
C
C     Default values
C
         ITURB  = ITRS
         IADRM  = IADRNO
         IASM   = IASMBU
         IEPSC  = IECNO
         IPDIL  = IPDNO
         ICMUST = 0
         IEPSLN = 12
         IFDL2K = 1
         NQ     = 12
         NP     =  5
         NF     = 12
         NRANK  = 12
         SIGK   = 1.0E0 / 0.75E0
         SIGT2  = 1.0E0 / 1.36E0
         SIGRHO = 1.0E0 / 0.50E0
         CMU    = 0.09E0
         CEPSLN = 0.10E0
         CEPS1  = 1.50E0
         CEPS2  = 1.83E0
         CEPS3  = - 2.0E0
         CEPS4  = - 1.5E0
         CEPS5  = 0.0E0
         CSUBS  = 0.11E0
         A2KE   = 4.9E0
         ALF1   = 0.0E0
         ALF2   = 0.0E0
         ALF3   = 0.0E0
C--->    CPDFRO = SIGRHO
C--->    CPDFT  = 1.E0 / PRT
         CPDFRO = 0.E0
         CPDFT  = 0.E0
         CPVELK = 0.E0
C
         PRDLIM = 100.E0
C
C     Freestream levels
C     Set TKEINF to (turbulence intensity=u'/U)^2; Non-dimensionalize later
C         RMUTNF to mu_T / mu
C
         TKEINF = 1.0E-04 * SQRT(10.E0)
         RMUTNF = 1.0E0 * 10.E0
C
         QMIN( 1)  =   1.E-8
         QMIN( 2)  = - 1.E0
         QMIN( 3)  = - 1.E0
         QMIN( 4)  = - 1.E0
         QMIN( 5)  =   1.E-8
         QMIN( 6)  =   1.E-8
         QMIN( 7)  =   1.E-8
         QMIN( 8)  =   1.E-8
         QMIN( 9)  = - 1.E0
         QMIN(10)  = - 1.E0
         QMIN(11)  = - 1.E0
         QMIN(12)  =   1.E-8
C
C     Pressure-Strain Correlation
C
C     Lai-So form (Zhang's) Model
C
         IDMMDL = IDMSO
         IDTMDL = IDTHL
         IPSMDL = IPLS
         CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
C
C     Turn on dui/dxj calculation
C
         IFDELQ(2) = 1
         IFDELQ(3) = 1
         IFDELQ(4) = 1
C
C     Read in values of constants
C
         NCON   = INT (VALUE)
         DO 103 ICON = 1, NCON
            READ (INPUT,*,END=990,ERR=995) VAR, VALUE
            WRITE (IOUT,1061) VAR, VALUE
            IF (VAR .EQ. 'LRR') THEN
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.36E0
               CMU    = 0.09E0
               CEPSLN = 0.15E0
               CEPS1  = 1.44E0
               CEPS2  = 1.90E0
               CEPS3  = 0.0E0
               CEPS4  = 0.0E0
               CEPS5  = 0.0E0
               CSUBS  = 0.11E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               IDMMDL = IDMSO
               IDTMDL = IDTHL
               IPSMDL = IPLRR1
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
            ELSE IF (VAR .EQ. 'SSG') THEN
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.36E0
               CMU    = 0.09E0
               CEPSLN = 0.1842E0
               CEPS1  = 1.44E0
               CEPS2  = 1.83E0
               CEPS3  = 0.0E0
               CEPS4  = 0.0E0
               CEPS5  = 0.0E0
               CSUBS  = 0.11E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               IDMMDL = IDMSO
               IDTMDL = IDTHL
               IPSMDL = IPSSG
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
            ELSE IF (VAR .EQ. 'SSG-ADRM') THEN
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.36E0
               CMU    = 0.094E0
               CEPSLN = 0.217E0
               CEPS1  = 1.00E0
               CEPS2  = 1.83E0
               CEPS3  = 0.0E0
               CEPS4  = 0.0E0
               CEPS5  = 0.0E0
               CSUBS  = 0.11E0
               A2KE   = 0.0E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               IDMMDL = IDMSO
               IDTMDL = IDTHL
               IPSMDL = IPSSG
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
               IEPSIJ = IEADGS
            ELSE IF (VAR .EQ. 'FLT') THEN
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.36E0
               CMU    = 0.09E0
               CEPSLN = 0.15E0
               CEPS1  = 1.45E0
               CEPS2  = 1.90E0
               CEPS3  = 0.0E0
               CEPS4  = 0.0E0
               CEPS5  = 0.0E0
               CSUBS  = 0.11E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               IDMMDL = IDMSO
               IDTMDL = IDTHL
               IPSMDL = IPFLT
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
            ELSE IF (VAR .EQ. 'ZSGS') THEN
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.36E0
               CMU    = 0.09E0
               CEPSLN = 0.10E0
               CEPS1  = 1.50E0
               CEPS2  = 1.83E0
               CEPS3  = - 2.0E0
               CEPS4  = - 1.5E0
               CEPS5  = 0.0E0
               CSUBS  = 0.11E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               IDMMDL = IDMSO
               IDTMDL = IDTHL
               IPSMDL = IPLS
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
            ELSE IF (VAR .EQ. 'SO 95') THEN
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.36E0
               CMU    = 0.09E0
               CEPSLN = 0.12E0
               CEPS1  = 1.50E0
               CEPS2  = 1.90E0
               CEPS3  = 7.E0/(3.E0*SQRT(15.E0)) * 0.0066E0
               CEPS4  = 2.95E0
               CEPS5  = 0.0E0
               CSUBS  = 0.11E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               IDMMDL = IDMSO
               IDTMDL = IDTHL
               IPSMDL = IPSO95
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
            ELSE IF (VAR .EQ. 'GIBSON LAUNDER') THEN
               SIGK   = 1.0E0 / 0.75E0
               SIGT2  = 1.0E0 / 1.36E0
               CMU    = 0.09E0
               CEPSLN = 0.16E0
               CEPS1  = 1.44E0
               CEPS2  = 1.92E0
               CEPS3  = 0.0E0
               CEPS4  = 0.0E0
               CEPS5  = 0.0E0
               CSUBS  = 0.11E0
               A2KE   = 4.9E0
               ALF1   = 0.0E0
               ALF2   = 0.0E0
               ALF3   = 0.0E0
               IDMMDL = IDMSO
               IDTMDL = IDTHL
               IPSMDL = IPGL
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
            ELSE IF (VAR .EQ. 'SIGMAK') THEN
               SIGK   = 1.0E0 / VALUE
            ELSE IF (VAR .EQ. 'SIGMAE') THEN
               SIGT2  = 1.0E0 / VALUE
            ELSE IF (VAR .EQ. 'SIGMARHO') THEN
               SIGRHO = 1.0E0 / VALUE
               CPDFRO = SIGRHO
            ELSE IF (VAR .EQ. 'CMU') THEN
               CMU    = VALUE
            ELSE IF (VAR .EQ. 'CEPSLN') THEN
               CEPSLN = VALUE
            ELSE IF (VAR .EQ. 'CEPS1') THEN
               CEPS1  = VALUE
            ELSE IF (VAR .EQ. 'CEPS2') THEN
               CEPS2  = VALUE
            ELSE IF (VAR .EQ. 'CEPS3') THEN
               CEPS3  = VALUE
            ELSE IF (VAR .EQ. 'CEPS4') THEN
               CEPS4  = VALUE
            ELSE IF (VAR .EQ. 'CEPS5') THEN
               CEPS5  = VALUE
            ELSE IF (VAR .EQ. 'CSUBS') THEN
               CSUBS  = VALUE
            ELSE IF (VAR .EQ. 'A2KE') THEN
               A2KE   = VALUE
            ELSE IF (VAR .EQ. 'ALF1') THEN
               ALF1   = VALUE
            ELSE IF (VAR .EQ. 'ALF2') THEN
               ALF2   = VALUE
            ELSE IF (VAR .EQ. 'ALF3') THEN
               ALF3   = VALUE
            ELSE IF (VAR .EQ. 'PRT') THEN
               PRT    = VALUE
               CPDFT  = 1.E0 / PRT
            ELSE IF (VAR .EQ. 'CPDFRO') THEN
               CPDFRO = VALUE
            ELSE IF (VAR .EQ. 'CPDFT') THEN
               CPDFT  = VALUE
            ELSE IF (VAR .EQ. 'CPVELK') THEN
               CPVELK = VALUE
            ELSE IF (VAR .EQ. 'PRDLIM') THEN
               PRDLIM = VALUE
            ELSE IF (VAR .EQ. 'FULL EPSILON PROD') THEN
               PRDE   = 0.E0
               PRDEM1 = 1.E0 - PRDE
            ELSE IF (VAR .EQ. 'PSC1') THEN
               PSC1   = VALUE
            ELSE IF (VAR .EQ. 'PSC2') THEN
               PSC2   = VALUE
            ELSE IF (VAR .EQ. 'PSCW') THEN
               PSCW   = VALUE
            ELSE IF (VAR .EQ. 'PSALFA') THEN
               PSALFA = VALUE
            ELSE IF (VAR .EQ. 'PSBETA') THEN
               PSBETA = VALUE
            ELSE IF (VAR .EQ. 'PSGAMA') THEN
               PSGAMA = VALUE
            ELSE IF (VAR .EQ. 'TKEINF') THEN
               TKEINF = VALUE
            ELSE IF (VAR .EQ. 'RMUTNF') THEN
               RMUTNF = VALUE
            ELSE IF (VAR .EQ. 'MINIMUM TKE') THEN
               QMIN( 6) = 2.E0 / 3.E0 * VALUE
               QMIN( 7) = 2.E0 / 3.E0 * VALUE
               QMIN( 8) = 2.E0 / 3.E0 * VALUE
            ELSE IF (VAR .EQ. 'MINIMUM EPSILON') THEN
               QMIN(12) = VALUE
            ELSE IF (VAR .EQ. 'TURBULENT DIFFUSION') THEN
               READ (INPUT,*,END=990,ERR=995) VAR
               WRITE (IOUT,1065) VAR
               IF (VAR .EQ. 'DH') THEN
                  IDTMDL = IDTDH
               ELSE IF (VAR .EQ. 'HL') THEN
                  IDTMDL = IDTHL
               ELSE IF (VAR .EQ. 'MH') THEN
                  IDTMDL = IDTMH
               ELSE
                  INPERR = INPERR + 1
                  WRITE (IOUT,1066) VAR
               ENDIF
            ELSE IF (VAR .EQ. 'PRESSURE STRAIN') THEN
               READ (INPUT,*,END=990,ERR=995) VAR
               WRITE (IOUT,1065) VAR
               IF (VAR .EQ. 'LRR1') THEN
                  IPSMDL = IPLRR1
               ELSE IF (VAR .EQ. 'LRR2') THEN
                  IPSMDL = IPLRR2
               ELSE IF (VAR .EQ. 'SSG') THEN
                  IPSMDL = IPSSG
               ELSE IF (VAR .EQ. 'FLT') THEN
                  IPSMDL = IPFLT
               ELSE IF (VAR .EQ. 'LAI SO') THEN
                  IPSMDL = IPLS
               ELSE IF (VAR .EQ. 'GIBSON LAUNDER') THEN
                  IPSMDL = IPGL
               ELSE
                  INPERR = INPERR + 1
                  WRITE (IOUT,1067) VAR
               ENDIF
               CALL SETPIJ (IPSMDL, IEPSIJ, IEPSEQ, IDAMP, ISRCW) 
            ELSE IF (VAR .EQ. 'COMPRESSIBLE DISS') THEN
               READ (INPUT,*,END=990,ERR=995) VAR
               WRITE (IOUT,1065) VAR
               IF (VAR .EQ. 'NONE') THEN
                  IEPSC  = IECNO
               ELSE IF (VAR .EQ. 'SARKAR') THEN
                  IEPSC  = IECSAR
               ELSE IF (VAR .EQ. 'RISTORCELLI') THEN
                  IEPSC  = IECRIS
               ELSE
                  INPERR = INPERR + 1
                  WRITE (IOUT,1068) VAR
               ENDIF
            ELSE IF (VAR .EQ. 'PRESSURE DILATATION') THEN
               READ (INPUT,*,END=990,ERR=995) VAR
               WRITE (IOUT,1227) VAR
               IF (VAR .EQ. 'NONE') THEN
                  IPDIL  = IPDNO
               ELSE IF (VAR .EQ. 'SARKAR') THEN
                  IPDIL  = IPDSAR
               ELSE IF (VAR .EQ. 'RISTORCELLI') THEN
                  IPDIL  = IPDRIS
               ELSE
                  INPERR = INPERR + 1
                  WRITE (IOUT,1069) VAR
               ENDIF
            ELSE
               WRITE (IOUT,1029) VAR
               INPERR = INPERR + 1
            ENDIF
  103    CONTINUE
C
C     These variables are used only for output
C
         SK     = 1.E0 / SIGK
         ST2    = 1.E0 / SIGT2
         WRITE (IOUT,1090) CEPSLN, CEPS1, CEPS2, CEPS3, CEPS4, CEPS5, 
     1                     CSUBS, CMU, SK, ST2, SIGRHO, A2KE,
     2                     ALF1, ALF2, ALF3, CPDFRO, CPDFT, CPVELK,
     3                     TKEINF, RMUTNF
         IF (IPSMDL .EQ. IPLRR1) THEN
            WRITE (IOUT,1095) 'LRR (Launder-Reece-Rodi)      '
         ELSE IF (IPSMDL .EQ. IPLRR2) THEN
            WRITE (IOUT,1095) 'LRR2 (Launder-Reece-Rodi 2)   '
         ELSE IF (IPSMDL .EQ. IPSSG ) THEN
            WRITE (IOUT,1095) 'SSG (Speziale-Sarkar-Gatski)  '
         ELSE IF (IPSMDL .EQ. IPFLT ) THEN
            WRITE (IOUT,1095) 'FLT (Fu-Launder-Tselepidakis) '
         ELSE IF (IPSMDL .EQ. IPLS  ) THEN
            WRITE (IOUT,1095) 'LS (LAI-SO)                   '
         ELSE IF (IPSMDL .EQ. IPGL  ) THEN
            WRITE (IOUT,1095) 'GL (Gibson-Launder)           '
         ENDIF
         WRITE (IOUT,1097) PSC1, PSC2, PSCW, PSALFA, PSBETA, PSGAMA
C
C     Large-Eddy Simulation (LES)
C
      ELSE IF (VAR .EQ. 'LES') THEN
C
C     Set default LES model to Smagorinsky high Reynolds number
C
         ITURB = ITLES
         ISGS  = ISGSMG
         IDAMP = IFMHR
         NQ    = 5
         NP    = 5
         NF    = 5
         NRANK = 5
         CMU   = 0.2E0
C
C     Turn on dui/dxj calculation
C
         IFDELQ(2) = 1
         IFDELQ(3) = 1
         IFDELQ(4) = 1
C
C     Read additional input lines
C
         NCON   = INT (VALUE)
         DO 203 ICON = 1, NCON
            READ (INPUT,*,END=990,ERR=995) VAR, VALUE
            WRITE (IOUT,1061) VAR, VALUE
            IF (VAR .EQ. 'SMAGORINSKY') THEN
               ISGS  = ISGSMG
               IDAMP = IFMHR
               CMU   = 0.2E0
            ELSE IF (VAR .EQ. 'STRUCTURE FUNCTION MODEL') THEN
               ISGS  = ISGSFM
               IDAMP = IFMHR
               CKSFM = 1.4E0
               CMU   = 0.105E0 * CKSFM**(-3.E0/2.E0)
            ELSE
               WRITE (IOUT,1029) VAR
               INPERR = INPERR + 1
            ENDIF
  203    CONTINUE
C
         IF (ISGS .EQ. ISGSMG) THEN
            WRITE (IOUT,1098) 'Smagorinsky                   '
         ELSE IF (ISGS .EQ. ISGSFM) THEN
            WRITE (IOUT,1098) 'Structure Function Model      '
         ENDIF
         WRITE (IOUT,1099) CMU
C
C     Convert turbulence model on restart
C
      ELSE IF (VAR .EQ. 'CONVERT') THEN
         CONTRB = .TRUE.
         READ (INPUT,*,END=990,ERR=995) VAR, VALUE
         WRITE (IOUT,1227) VAR
         IF (VAR .EQ. 'LAMINAR') THEN
            ITURBX = ITLMNR
            WRITE (IOUT,1110) VAR
         ELSE IF (VAR .EQ. 'BALDWIN LOMAX') THEN
            ITURBX = ITBLMX
            WRITE (IOUT,1110) VAR
         ELSE IF (VAR .EQ. 'K EPSILON') THEN
            WRITE (IOUT,1110) VAR
            ITURBX = ITKE
            IDAMPX = IFMZSG
            IASMX  = IASMBU
            NCON   = INT(VALUE)
            DO 104 ICON = 1, NCON
               READ (INPUT,*,END=990,ERR=995) VAR
               WRITE (IOUT,1065) VAR
               IF (VAR .EQ. 'SAA') THEN
                  IDAMPX = IFMSAA
               ELSE IF (VAR .EQ. 'ZSGS') THEN
                  IDAMPX = IFMZSG
               ELSE IF (VAR .EQ. 'ZSSL') THEN
                  IDAMPX = IFMZSG
               ELSE IF (VAR .EQ. 'HIGH RE') THEN
                  IDAMPX = IFMHR
               ELSE IF (VAR .EQ. 'RNG') THEN
                  IDAMPX = IFMHR
               ELSE IF (VAR .EQ. 'ASM') THEN
                  ICMUST = 5
                  IASMX  = IASMGS
               ELSE
                  WRITE (IOUT,1105) VAR
                  INPERR = INPERR + 1
               END IF
  104       CONTINUE
         ELSE IF (VAR .EQ. 'K OMEGA') THEN
            WRITE (IOUT,1110) VAR
            ITURBX = ITKW
            IDAMPX = IFMHR
            IASMX  = IASMBU
            NCON   = INT(VALUE)
            DO 105 ICON = 1, NCON
               READ (INPUT,*,END=990,ERR=995) VAR
               IF (VAR .EQ. 'ASM') THEN
                  ICMUST = 5
                  IASMX  = IASMGS
               ELSE
                  WRITE (IOUT,1105) VAR
                  INPERR = INPERR + 1
               END IF
  105       CONTINUE
         ELSE IF (VAR .EQ. 'REYNOLDS STRESS') THEN
            WRITE (IOUT,1110) VAR
            ITURBX = ITRS
            IDAMPX = IFMZSG
            IASMX  = IASMBU
         ELSE
            WRITE (IOUT,1100) VAR
            INPERR = INPERR + 1
         ENDIF
C
C        Set up IMODLX to control conversion model 
C        Only require ITURBX, IDAMPX, and IASMX set all else to error code
C
         DO 106 IM = 1, NMDL
            IMODLX(IM) = - 1
  106    CONTINUE
         IMODLX( 1) = ITURBX
         IMODLX( 3) = IDAMPX
         IMODLX(10) = IASMX
C
C     End of Turbulence Models
C
C     Metric/English units
C
      ELSE IF (VAR .EQ. 'SIUNIT') THEN
         SIUNIT = .TRUE.
      ELSE IF (VAR .EQ. 'ENGLISH') THEN
         SIUNIT = .FALSE.
C
C     Two/Three-Dimensional
C
      ELSE IF (VAR .EQ. 'TWOD') THEN
         THREED = .FALSE.
         AXISYM = .FALSE.
      ELSE IF (VAR .EQ. 'THREED') THEN
         THREED = .TRUE.
         AXISYM = .FALSE.
      ELSE IF (VAR .EQ. 'AXISYMMETRIC') THEN
         THREED = .FALSE.
         AXISYM = .TRUE.
      ELSE IF (VAR .EQ. 'GRID') THEN
         READ (INPUT,*,END=990,ERR=995) GRDFIL(1)
         WRITE (IOUT,1040) GRDFIL(1)
         IGRDTP = IGTP3D
      ELSE IF (VAR .EQ. 'GRID BINARY') THEN
         GRDBIN = .TRUE.
      ELSE IF (VAR .EQ. 'GRID FORMATTED') THEN
         GRDBIN = .FALSE.
C
C     Flow Quantities
C
      ELSE IF (VAR .EQ. 'GAMMA') THEN
         GAMMA  = VALUE
      ELSE IF (VAR .EQ. 'MACH') THEN
         FSMACH = VALUE
      ELSE IF (VAR .EQ. 'ALPHA') THEN
         ALPHA = VALUE
      ELSE IF (VAR .EQ. 'BETA') THEN
         BETA = VALUE
      ELSE IF (VAR .EQ. 'RE') THEN
         RE = VALUE
      ELSE IF (VAR .EQ. 'PR') THEN
         PR = VALUE
      ELSE IF (VAR .EQ. 'PRT') THEN
         PRT = VALUE
      ELSE IF (VAR .EQ. 'TINF') THEN
         TINF = VALUE
      ELSE IF (VAR .EQ. 'TWALL') THEN
         TWSPE = .TRUE.
         TWALL = VALUE
C
C     MUSCL Control
C
      ELSE IF (VAR .EQ. 'KAPPA') THEN
         DO 108 IBLK = 1, NXBLKS
            RKAP(1,IBLK) = VALUE
            RKAP(2,IBLK) = VALUE
            RKAP(3,IBLK) = VALUE
  108    CONTINUE
      ELSE IF (VAR .EQ. 'LIMIT') THEN
         LIM  = INT(VALUE)
         LIMX = INT(VALUE)
         IF (LIMX .EQ. ILNONE) THEN
            LIM = ILNONE
         ELSE IF (LIMX .EQ. ILSMTH) THEN
            LIM = ILSMTH
         ELSE IF (LIMX .EQ. ILMNMD) THEN
            LIM = ILMNMD
         ELSE IF (LIMX .EQ. ILVNKT) THEN
            LIM = ILVNKT
         ELSE
            WRITE (IOUT,1058) VAR, VALUE
            INPERR = INPERR + 1
         ENDIF
         DO 110 IBLK = 1, NXBLKS
            LIMIT(1,IBLK) = LIM
            LIMIT(2,IBLK) = LIM
            LIMIT(3,IBLK) = LIM
  110    CONTINUE
C
C     Fourth order
C
      ELSE IF (VAR .EQ. 'CENTRAL DIFFERENCE') THEN
  114    READ (INPUT,*,END=990,ERR=995) VAR, VALUE
         WRITE (IOUT,2010) VAR, VALUE
         IF (VAR .EQ. 'FOURTH ORDER') THEN
            FOURTH = .TRUE.
         ELSE IF (VAR .EQ. 'ARTIFICIAL DISSIPATION') THEN
            ADKP46 = VALUE
         ELSE IF (VAR .EQ. 'END CENTRAL DIFFERENCE') THEN
            GO TO 118
         ELSE
            WRITE (IOUT,1059) VAR, VALUE
            INPERR = INPERR + 1
         ENDIF
C
C     Loop for next control of central difference
C
         GO TO 114
C
C     Finished with central difference control
C
  118    CONTINUE
C
C     Subsonic Inflow/Outflow BC conditions
C
      ELSE IF (VAR .EQ. 'BACK PRESSURE') THEN
         PBAKBC = VALUE
      ELSE IF (VAR .EQ. 'INFLOW PTOTAL') THEN
         PTOTBC = VALUE
      ELSE IF (VAR .EQ. 'INFLOW TTOTAL') THEN
         TTOTBC = VALUE
C
C     Jet Conditions
C
      ELSE IF (VAR .EQ. 'JET CONDITIONS') THEN
         JET  = INT(VALUE)
C
C     Read in Jet conditions
C
  120    READ (INPUT,*,END=990,ERR=995) VAR, VALUE
         WRITE (IOUT,2010) VAR, VALUE
         IF (VAR .EQ. 'U') THEN
            UJET = VALUE
         ELSE IF (VAR .EQ. 'V') THEN
            VJET = VALUE
         ELSE IF (VAR .EQ. 'W') THEN
            WJET = VALUE
         ELSE IF (VAR .EQ. 'P') THEN
            PJET = VALUE
         ELSE IF (VAR .EQ. 'T') THEN
            TJET = VALUE
         ELSE IF (VAR .EQ. 'K') THEN
            TKEJET = VALUE
         ELSE IF (VAR .EQ. 'EPSILON') THEN
            EPSJET = VALUE
         ELSE IF (VAR .EQ. 'OMEGA') THEN
            OMEJET = VALUE
         ELSE IF (VAR .EQ. 'TAUXX') THEN
            TAUJET(1) = VALUE
         ELSE IF (VAR .EQ. 'TAUYY') THEN
            TAUJET(2) = VALUE
         ELSE IF (VAR .EQ. 'TAUZZ') THEN
            TAUJET(3) = VALUE
         ELSE IF (VAR .EQ. 'TAUXY') THEN
            TAUJET(4) = VALUE
         ELSE IF (VAR .EQ. 'TAUXZ') THEN
            TAUJET(5) = VALUE
         ELSE IF (VAR .EQ. 'TAUYZ') THEN
            TAUJET(6) = VALUE
         ELSE IF (VAR .EQ. 'END JET') THEN
            GO TO 129
         ELSE
            WRITE (IOUT,1060) VAR, VALUE
            INPERR = INPERR + 1
         ENDIF
C
C     Loop for next Jet condition
C
         GO TO 120
C
C     Finished inputing Jet conditions
C
  129    CONTINUE
C
C     Printout control
C
      ELSE IF (VAR .EQ. 'PRINT') THEN
         NPRT = INT(VALUE)
         DO 300 IPRT = 1, NPRT
            READ (INPUT,*,END=990,ERR=995) IBLOCK, (IDUMPR(LL),LL=1,10)
            WRITE (IOUT,2005)              IBLOCK, (IDUMPR(LL),LL=1,10)
            NPRINT(IBLOCK) = NPRINT(IBLOCK) + 1
            IF (NPRINT(IBLOCK) .GT. NXPRT) THEN
               WRITE (IOUT,1540) NPRINT(IBLOCK), IBLOCK, NXPRT
               GO TO 900
            ENDIF
            DO 299 LL = 1, 10
               IPRINT(LL,NPRINT(IBLOCK),IBLOCK) = IDUMPR(LL)
 299        CONTINUE
 300     CONTINUE
      ELSE IF (VAR .EQ. 'MONITOR') THEN
         NPRT = INT(VALUE)
         READ (INPUT,*,END=990,ERR=995) NITMON
         DO 310 IPRT = 1, NPRT
            READ (INPUT,*,END=990,ERR=995) IBLOCK,
     1                (IPRMON(LL,NPRMON(IBLOCK)+1,IBLOCK),LL=1,10)
            WRITE (IOUT,2005) IBLOCK, 
     1                (IPRMON(LL,NPRMON(IBLOCK)+1,IBLOCK),LL=1,10)
            NPRMON(IBLOCK) = NPRMON(IBLOCK) + 1
 310     CONTINUE
C
C     CUT boundary condition data
C
      ELSE IF (VAR .EQ. 'CUT') THEN
         NCUTS = INT(VALUE)
         IF (NCUTS*2 .GT. NXCUTS) THEN
            WRITE (IOUT,1530) NCUTS*2, NXCUTS
            GO TO 900
         ENDIF
         DO 329 ICUT = 1, NCUTS
            IC1 = 2 * ICUT - 1
            IC2 = 2 * ICUT
c           READ (INPUT,*,END=990,ERR=995) VAR, IBLK, (ISTRT(L),L=1,3), 
c    1                                        (IEND(M),M=1,3)
            read (input,*,END=990,ERR=995) var, iblk,istrt(1),istrt(2),
     1           istrt(3), iend(1), iend(2), iend(3)
            WRITE (IOUT,1600) VAR, IBLK, (ISTRT(L),L=1,3), 
     1                                   (IEND(M),M=1,3)
C
C     Determine CUT face. This is defined by constant index on that face.
C     Note that this will fail for the pathological case where have two
C     or more indices equal for a cut. In case an error occurs, set the
C     face direction to 1 to finish checking the CUT data.
C
            IDIR = 0
            DO 320 L = 1, 3
               IF (ISTRT(L) .EQ. IEND(L)) IDIR = L
  320       CONTINUE
            IF (IDIR .LT. 1 .OR. IDIR .GT. 3) THEN
               WRITE (IOUT,1630) 
               INPERR = INPERR + 1
               IDIR   = 1
            ENDIF
C
C     Modify start or end appropriately to account for cut input as grid
C     points but cut to be cell-centers.  Do this only for the indices that
C     the cut loops over, not the idir direction.
C
            DO 321 L = 1, 3
               IF (L .NE. IDIR) THEN
                  IF (IEND(L) .GT. ISTRT(L)) THEN
                     ISTRT(L) = ISTRT(L) + 1
                  ELSE
                     IEND(L)  = IEND(L)  + 1
                  ENDIF
               ENDIF
  321       CONTINUE
C
            CUTNAM(IC1) = VAR
            CUTNAM(IC2) = VAR
C
            ICUTS( 1,IC1,1) = IBLK
            ICUTS( 2,IC1,1) = IDIR
            ICUTS( 3,IC1,1) = ISTRT(IDIR)
            ICUTS( 4,IC1,1) = ISTRT(1)
            ICUTS( 5,IC1,1) = IEND (1)
            ICUTS( 6,IC1,1) = ISTRT(2)
            ICUTS( 7,IC1,1) = IEND (2)
            ICUTS( 8,IC1,1) = ISTRT(3)
            ICUTS( 9,IC1,1) = IEND (3)
C
            ICUTS(10,IC2,1) = IBLK
            ICUTS(11,IC2,1) = IDIR
            ICUTS(12,IC2,1) = ISTRT(IDIR)
            ICUTS(13,IC2,1) = ISTRT(1)
            ICUTS(14,IC2,1) = IEND (1)
            ICUTS(15,IC2,1) = ISTRT(2)
            ICUTS(16,IC2,1) = IEND (2)
            ICUTS(17,IC2,1) = ISTRT(3)
            ICUTS(18,IC2,1) = IEND (3)
C
            READ (INPUT,*,END=990,ERR=995) VAR1, IBLK, (ISTRT(L),L=1,3), 
     1                             (IEND(M),M=1,3), (IORD(N),N=1,3)
c           read (input,*,END=990,ERR=995) var1, iblk,istrt(1),istrt(2),
c    1           istrt(3), iend(1), iend(2), iend(3)
            WRITE (IOUT,1605) VAR1, IBLK, (ISTRT(L),L=1,3), 
     1                             (IEND(M),M=1,3), (IORD(N),N=1,3)
            IF (VAR .NE. VAR1) THEN
               INPERR = INPERR + 1
               WRITE (IOUT,1610) 
            ENDIF
C
            DO 324 L = 1, 3
               IF (IORD(L) .LT. 1 .OR. IORD(L) .GT. 3) THEN
                  INPERR = INPERR + 1
                  WRITE (IOUT,1620) L, IORD(L)
               ENDIF
  324       CONTINUE
C
C     Determine CUT face. This is defined by constant index on that face.
C     Note that this will fail for the pathological case where have two
C     or more indices equal for a cut. In case an error occurs, set the
C     face direction to 1 to finish checking the CUT data.
C
            IDIR = 0
            DO 330 L = 1, 3
               IF (ISTRT(L) .EQ. IEND(L)) IDIR = L
  330       CONTINUE
            IF (IDIR .LT. 1 .OR. IDIR .GT. 3) THEN
               WRITE (IOUT,1630) 
               INPERR = INPERR + 1
               IDIR   = 1
            ENDIF
C
C     Modify start or end appropriately to account for cut input as grid
C     points but cut to be cell-centers.  Do this only for the indices that
C     the cut loops over, not the idir direction.
C
            DO 331 L = 1, 3
               IF (L .NE. IDIR) THEN
                  IF (IEND(L) .GT. ISTRT(L)) THEN
                     ISTRT(L) = ISTRT(L) + 1
                  ELSE
                     IEND(L)  = IEND(L)  + 1
                  ENDIF
               ENDIF
  331       CONTINUE
C
            ICUTS(10,IC1,1) = IBLK
            ICUTS(11,IC1,1) = IDIR
            ICUTS(12,IC1,1) = ISTRT(IDIR)
            ICUTS(13,IC1,1) = ISTRT(1)
            ICUTS(14,IC1,1) = IEND (1)
            ICUTS(15,IC1,1) = ISTRT(2)
            ICUTS(16,IC1,1) = IEND (2)
            ICUTS(17,IC1,1) = ISTRT(3)
            ICUTS(18,IC1,1) = IEND (3)
            ICUTS(19,IC1,1) = IORD (1)
            ICUTS(20,IC1,1) = IORD (2)
            ICUTS(21,IC1,1) = IORD (3)
C
            ICUTS( 1,IC2,1) = IBLK
            ICUTS( 2,IC2,1) = IDIR
            ICUTS( 3,IC2,1) = ISTRT(IDIR)
            ICUTS( 4,IC2,1) = ISTRT(1)
            ICUTS( 5,IC2,1) = IEND (1)
            ICUTS( 6,IC2,1) = ISTRT(2)
            ICUTS( 7,IC2,1) = IEND (2)
            ICUTS( 8,IC2,1) = ISTRT(3)
            ICUTS( 9,IC2,1) = IEND (3)
            ICUTS(19,IC2,1) = ICY(IORD(1),1)
            ICUTS(20,IC2,1) = ICY(IORD(2),1)
            ICUTS(21,IC2,1) = ICY(IORD(3),1)
  329    CONTINUE
         NCUTS = NCUTS * 2
C
C     Data for each grid block/zone 
C
      ELSE IF (VAR .EQ. 'BLOCK') THEN
         IBLOCK = INT(VALUE)
         NBLKS  = NBLKS + 1
         IF (NBLKS .GT. NXBLKS) THEN
            WRITE (IOUT,1525) NBLKS, NXBLKS
            GO TO 900
         ENDIF
  500    CONTINUE
         READ (INPUT,*,END=990,ERR=995) VAR, VALUE
         WRITE (IOUT,2010) VAR, VALUE
         IF (VAR .EQ. 'DIMENSIONS') THEN
            READ (INPUT,*,END=990,ERR=995) IDIM(IBLOCK,1), 
     1                             JDIM(IBLOCK,1), KDIM(IBLOCK,1)
            WRITE (IOUT,1035) IDIM(IBLOCK,1), 
     1                             JDIM(IBLOCK,1), KDIM(IBLOCK,1)
         ELSE IF (VAR .EQ. 'KAPPA') THEN
            READ (INPUT,*,END=990,ERR=995) (RKAP(LL,IBLOCK),LL=1,3)
            WRITE (IOUT,2000) ' Kappa ', (RKAP(LL,IBLOCK),LL=1,3)
C
C N.B. this input needs modified to use the limiter mnemonics
C
         ELSE IF (VAR .EQ. 'LIMIT') THEN
            READ (INPUT,*,END=990,ERR=995) (LIMIT(LL,IBLOCK),LL=1,3)
            WRITE (IOUT,2002) ' Limit ', (LIMIT(LL,IBLOCK),LL=1,3)
C-       ELSE IF (VAR .EQ. 'ORDER') THEN
C-          READ (INPUT,*,END=990,ERR=995) (ITORDR(LL),LL=1,3)
C-          WRITE (IOUT,2002) ' Order ', (ITORDR(LL),LL=1,3)
C-          DO 509 ILVL = 1, NXLVLS
C-             IORDER(1,IBLOCK,ILVL,ILVL) = ITORDR(1)
C-             IORDER(2,IBLOCK,ILVL,ILVL) = ITORDR(2)
C-             IORDER(3,IBLOCK,ILVL,ILVL) = ITORDR(3)
C-509       CONTINUE
         ELSE IF (VAR .EQ. 'FLUX') THEN
            READ (INPUT,*,END=990,ERR=995) VAR
            WRITE (IOUT,2010) VAR
            IF (VAR .EQ. 'ROE') THEN
               IFLUX(IBLOCK) = IFROE
            ELSE
               WRITE (IOUT,1030)
               INPERR = INPERR + 1
            ENDIF
         ELSE IF (VAR .EQ. 'VISCOUS') THEN
            IV = INT(VALUE)
            READ (INPUT,*,END=990,ERR=995) (IDUM(LL),LL=1,IV)
            WRITE (IOUT,1035) (IDUM(LL),LL=1,IV)
            DO 90 LL = 1, IV
               VISCOS(IDUM(LL),IBLOCK) = .TRUE.
 90         CONTINUE
            ITURB = MAX (ITURB, ITLMNR)
         ELSE IF (VAR .EQ. 'DIAGONAL') THEN
            IV = INT(VALUE)
            READ (INPUT,*,END=990,ERR=995) (IDUM(LL),LL=1,IV)
            WRITE (IOUT,1035) (IDUM(LL),LL=1,IV)
            DO 91 LL = 1, IV
               DIAG  (IDUM(LL),IBLOCK) = .TRUE.
 91         CONTINUE
         ELSE IF (VAR .EQ. 'INVISCID') THEN
            DO 95 LL = 1, 3
               VISCOS(LL,IBLOCK) = .FALSE.
 95         CONTINUE
         ELSE IF (VAR .EQ. 'GRID') THEN
            READ (INPUT,*,END=990,ERR=995) GRDFIL(IBLOCK)
            WRITE (IOUT,1040) GRDFIL(IBLOCK)
            IGRDTP = IGTSEP
C
C     Boundary Condition Data
C
         ELSE IF (VAR .EQ. 'BC') THEN
            NBCS(IBLOCK) = INT (VALUE)
            IF (NBCS(IBLOCK) .GT. NXBCS) THEN
               INPERR = INPERR + 1
               WRITE (IOUT,1535) NBCS(IBLOCK), IBLOCK, NXBCS
               GO TO 900
            ENDIF
            DO 200 I = 1, NBCS(IBLOCK)
c              READ (INPUT,*,END=990,ERR=995) VAR, (ISTRT(L),L=1,3), 
c    1                                (IEND(M),M=1,3), BCTYP
               read (input,*,END=990,ERR=995) var, istrt(1), istrt(2),
     1              istrt(3), iend(1), iend(2), iend(3), bctyp
               WRITE (IOUT,1045) VAR, (ISTRT(L),L=1,3), 
     1                                (IEND(M),M=1,3), BCTYP
C
C     Determine BC face. This is defined by constant index on that face.
C     Note that this will fail for the pathological case where have two
C     or more indices equal for a bc. In case an error occurs, set the
C     face direction to 1 to finish checking the BC data.
C
               IDIR = 0
               DO 150 L = 1, 3
                  IF (ISTRT(L) .EQ. IEND(L)) IDIR = L
 150           CONTINUE
               IF (IDIR .LT. 1 .OR. IDIR .GT. 3) THEN
                  WRITE (IOUT,1049) 
                  INPERR = INPERR + 1
                  IDIR   = 1
               ENDIF
C
C     Set profile boundary segment identifier to zero for all other cases
C
               INPROF = 0
C
               IF (BCTYP .EQ. 'WALL') THEN
                  IBCDAT(1,I,IBLOCK,1) = IWALL
               ELSE IF (BCTYP .EQ. 'WALL FUNCTIONS') THEN
                  IBCDAT(1,I,IBLOCK,1) = IWALFN
                  IFWALF = .TRUE.
               ELSE IF (BCTYP .EQ. 'TANGENCY') THEN
                  IBCDAT(1,I,IBLOCK,1) = ITAN
               ELSE IF (BCTYP .EQ. 'FARFIELD') THEN
                  IBCDAT(1,I,IBLOCK,1) = IFAR
               ELSE IF (BCTYP .EQ. '2D FARFIELD VORTEX') THEN
                  IFFORC               = IYES
                  IBCDAT(1,I,IBLOCK,1) = IFAR2D
               ELSE IF (BCTYP .EQ. 'SYMMETRY XY') THEN
                  IBCDAT(1,I,IBLOCK,1) = ISYMXY
               ELSE IF (BCTYP .EQ. 'SYMMETRY XZ') THEN
                  IBCDAT(1,I,IBLOCK,1) = ISYMXZ
               ELSE IF (BCTYP .EQ. 'SYMMETRY YZ') THEN
                  IBCDAT(1,I,IBLOCK,1) = ISYMYZ
               ELSE IF (BCTYP .EQ. 'EXTRAPOLATE') THEN
                  IBCDAT(1,I,IBLOCK,1) = IEXT
               ELSE IF (BCTYP .EQ. 'FIX') THEN
                  IBCDAT(1,I,IBLOCK,1) = IFIX
               ELSE IF (BCTYP .EQ. 'JET') THEN
                  IBCDAT(1,I,IBLOCK,1) = IJET
               ELSE IF (BCTYP .EQ. 'PERIODIC') THEN
                  IBCDAT(1,I,IBLOCK,1) = IPERD
               ELSE IF (BCTYP .EQ. 'HOMOGENEOUS SHEAR') THEN
                  IBCDAT(1,I,IBLOCK,1) = IHSHR
               ELSE IF (BCTYP .EQ. 'SUBSONIC INFLOW') THEN
                  IBCDAT(1,I,IBLOCK,1) = IINFLO
               ELSE IF (BCTYP .EQ. 'SUBSONIC OUTFLOW') THEN
                  IBCDAT(1,I,IBLOCK,1) = IOUTFL
               ELSE IF (BCTYP .EQ. 'AXISYMMETRIC') THEN
                  IBCDAT(1,I,IBLOCK,1) = IAXISM
               ELSE IF (BCTYP .EQ. 'AXISYMMETRIC CENTERLINE') THEN
                  IBCDAT(1,I,IBLOCK,1) = IAXICL
               ELSE IF (BCTYP .EQ. 'PROFILE CELL CENTER') THEN
                  IBCDAT(1,I,IBLOCK,1)  = IPROFL
                  NPROFL                = NPROFL + 1
                  INPROF                = NPROFL
                  READ (INPUT,*,END=990,ERR=995) FPROFL(NPROFL)
                  WRITE (IOUT,1047) FPROFL(NPROFL)
               ELSE IF (BCTYP .EQ. 'SUBSONIC PROFILE CC') THEN
                  IBCDAT(1,I,IBLOCK,1)  = IPROSB
                  NPROFL                = NPROFL + 1
                  INPROF                = NPROFL
                  READ (INPUT,*,END=990,ERR=995) FPROFL(NPROFL)
                  WRITE (IOUT,1048) FPROFL(NPROFL)
               ELSE
                  WRITE (IOUT,1050) BCTYP
                  INPERR = INPERR + 1
               ENDIF
C
C     Modify starting point appropriately to account for bc input as grid
C     points but bc to be cell-centers.  Do this only for the indices that
C     the bc loops over, not the idir direction.  BCs must be input with
C     indices in increasing order.
C     The starts are the max of the input starting point and 2 for the 
C     looping directions on a given boundary so that only loop over the 
C     interior boundary points to fit with the dimensioning of Q(+/-) etc 
C     as used in the MUSCL and FVSTL routines.
C
               DO 151 L = 1, 3
                  IF (L .NE. IDIR) THEN
                     ISTRT(L) = ISTRT(L) + 1
                  ENDIF
  151          CONTINUE
C
               IBCDAT( 2,I,IBLOCK,1) = IDIR
               IBCDAT( 3,I,IBLOCK,1) = ISTRT(IDIR)
               IBCDAT( 4,I,IBLOCK,1) = ISTRT(1)
               IBCDAT( 5,I,IBLOCK,1) = IEND (1)
               IBCDAT( 6,I,IBLOCK,1) = ISTRT(2)
               IBCDAT( 7,I,IBLOCK,1) = IEND (2)
               IBCDAT( 8,I,IBLOCK,1) = ISTRT(3)
               IBCDAT( 9,I,IBLOCK,1) = IEND (3)
               IBCDAT(10,I,IBLOCK,1) = INPROF
C
  200       CONTINUE
C
C     If read an END BLOCK, then have reached the end of the data for
C     the grid block/zone and are back in the main loop
C
         ELSE IF (VAR .EQ. 'END BLOCK') THEN
            GO TO 50
         ELSE
            WRITE (IOUT,1055) VAR, VALUE
            INPERR = INPERR + 1
         ENDIF
C
C     Go back and get the next line of input from the block/zone data
C
         GO TO 500
C
C     End of data for grid block/zone
C
      ELSE
C
C     Error in main loop
C
         WRITE (IOUT,1060) VAR, VALUE
         INPERR = INPERR + 1
      ENDIF
C
C     Go back and read the next line
C
      GO TO 50
C
C     End of input file
C
  900 CONTINUE
      WRITE (IOUT,1002)
C
C 3.  Set values based on input
C
C     Iterative, implicit items
C     Total additional levels of time storage
C
      NTMLVL = NTIME + NTMTAU
C
C     Set calculation of forces: if output of forces is requested or the
C                                2D farfield vortex bc is used then calculate
C
      IF (IFRCPR .EQ. IYES) IFFORC = IYES
C
C     If unsteady, then verify using AF3F or DIAGONAL AF3F scheme.
C     Set TIMEDT to DT if not previously set; if DT not set then output error.
C     t-ts scheme must have DT = TIMEDT.
C
      IF (ITIMED .EQ. IUNSTD) THEN
         IF ( .NOT. (ISOLVR .EQ. IAF3F) ) THEN
            WRITE (IOUT,2015)
            INPERR = INPERR + 1
         ENDIF
C
         IF (TIMEDT .LE. 0.E0) THEN
            IF (DT .GT. 0.E0) THEN
               TIMEDT = DT
               WRITE (IOUT,2020) DT
            ELSE
               WRITE (IOUT,2030)
               INPERR = INPERR + 1
            ENDIF
         ENDIF
C
         IF (ITMTYP .EQ. ITTS) THEN
            IF (DT .NE. TIMEDT) THEN
               WRITE (IOUT,2040)
               INPERR = INPERR + 1
            ENDIF
         ENDIF
      ENDIF
C
C     Check for errors in the input deck
C
      IF (INPERR .GT. 0) THEN
         WRITE (IOUT,2050) INPERR
         STOP
      ENDIF
C
C     Set-up unit normal to use for collapsed cell faces
C     SHAT = (1,1,1) / sqrt(3) for 3D
C            (1,1,0) / sqrt(2) for 3D
C
      IF (THREED) THEN
         DENOM = SQRT (3.E0)
         SHATX = 1.E0 / DENOM
         SHATY = 1.E0 / DENOM
         SHATZ = 1.E0 / DENOM
      ELSE
         DENOM = SQRT (2.E0)
         SHATX = 1.E0 / DENOM
         SHATY = 1.E0 / DENOM
         SHATZ = 0.E0 
      ENDIF
C
C     Set-up IMODEL
C
      IMODEL( 1) = ITURB
      IMODEL( 2) = IEPSEQ
      IMODEL( 3) = IDAMP
      IMODEL( 4) = IADRM
      IMODEL( 5) = IFDL2K
      IF (ITURB .EQ. ITRS) THEN
         IMODEL( 6) = IPSMDL
         IMODEL( 7) = IEPSIJ
         IMODEL( 8) = IDMMDL
         IMODEL( 9) = IDTMDL
      ELSE IF (ITURB .EQ. ITLES) THEN
         IMODEL( 6) = ISGS
      ENDIF
      IMODEL(10) = IASM
      IMODEL(11) = IEPSC
      IMODEL(12) = IPDIL
      IMODEL(13) = ISRCW
C
C     Set values that depend on the entire range of input
C
      IF (ITURB .EQ. ITKE) THEN
         TKEINF = TKEINF * FSMACH * FSMACH
         EPSINF = CMU * RHOINF * TKEINF * TKEINF / RMUTNF * RE / FSMACH
      ELSE IF (ITURB .EQ. ITKW) THEN
         TKEINF = TKEINF * FSMACH * FSMACH
         OMEINF = CMU * RHOINF * TKEINF / RMUTNF * RE / FSMACH
      ELSE IF (ITURB .EQ. ITRS) THEN
         TKEINF = TKEINF * FSMACH * FSMACH
         EPSINF = CMU * RHOINF * TKEINF * TKEINF / RMUTNF * RE / FSMACH
C
C Isotropic turbulence
C
         TAUINF(1) = 2.E0 / 3.E0 * TKEINF
         TAUINF(2) = 2.E0 / 3.E0 * TKEINF
         TAUINF(3) = 2.E0 / 3.E0 * TKEINF
         TAUINF(4) = 0.E0
         TAUINF(5) = 0.E0
         TAUINF(6) = 0.E0
      ENDIF
C
C     Set PBAKBC, PTOTBC, TTOTBC if not input.
C     Convert PBAKBC, PTOTBC, TTOTBC from dimensional to non-dimensional.
C
      ARG    = 1.E0 + 0.5E0 * (GAMMA-1.E0) * FSMACH * FSMACH
      IF (PBAKBC .LE. 0.E0) PBAKBC = 1.E0
      IF (PTOTBC .LE. 0.E0) PTOTBC = ARG ** (GAMMA/(GAMMA-1.E0))
      IF (TTOTBC .LE. 0.E0) TTOTBC = TINF * ARG
C
      PBAKBC = PBAKBC / GAMMA
      PTOTBC = PTOTBC / GAMMA
      TTOTBC = TTOTBC / TINF
C
C     Convert PJET from P(jet) / Pinf to P(jet) / (rhoinf*ainf**2) 
C     Set the value of the Jet density
C
      PJET   = PJET / GAMMA
      RHOJET = GAMMA * PJET / TJET
C
C     Set TDPDTJ, TDPDEL for t-ts or tau-ts time stepping.
C
      IF (ITMTYP .EQ. ITTS) THEN
         TDPDTJ = TDPHI  + 1.E0
         TDPDEL = 0.E0
      ELSE
         TDPDTJ = TDPHIP + 1.E0
         TDPDEL = TDPHI  + 1.E0
      ENDIF
C
C     Sort the BCs into ascending order based on BC type
C
      WRITE (IOUT,2100)
      DO 980 IBLOCK = 1, NBLKS
         DO 970 I = 1, NBCS(IBLOCK)
            IBCTYP = IBCDAT(1,I,IBLOCK,1) 
            DO 960 J = I+1, NBCS(IBLOCK)
               IBCTST = IBCDAT(1,J,IBLOCK,1) 
               IF (IBCTST .LT. IBCTYP) THEN
                  IBCTYP = IBCTST
                  DO 950 L = 1, 10
                     IBCTMP               = IBCDAT(L,I,IBLOCK,1)
                     IBCDAT(L,I,IBLOCK,1) = IBCDAT(L,J,IBLOCK,1)
                     IBCDAT(L,J,IBLOCK,1) = IBCTMP
  950             CONTINUE
               ENDIF
  960       CONTINUE
  970    CONTINUE
  980 CONTINUE
C
C     Finished
C
      RETURN
C
C 4.   Error handling section
C 4-a. End-Of-File Handler
C
  990 CONTINUE
      WRITE (IOUT,9990)
      STOP
C
C 4-b. Error Handler
C
  995 CONTINUE
      WRITE (IOUT,9995)
      STOP
C
 1001 FORMAT (' ','RDDATA: Begin Input Echo')
 1002 FORMAT (' ','RDDATA: End Input Echo')
 1000 FORMAT (5X,A25,5X,E14.7)
 1005 FORMAT (' ',4X,'Restart being written to file:',
     1       /' ',9X,A70)
 1010 FORMAT (' ','Illegal PLOT3D format option = ',A25,
     1       /' ','Please specify either BINARY or FORMATTED ')
 1015 FORMAT (' ',4X,'Output PLOT3D Grid    : ',A50
     1       /' ',4X,'              Solution: ',A50,
     2       /' ',4X,'              Function: ',A50)
 1020 FORMAT (' ','Specified number of stages, NSTAGE ',I2,
     1        ' ','greater than allowed, MXSTAG ',I2,
     2       /' ','Reset MXSTAG in main and recompile.')
 1021 FORMAT (' ',4X,'Runge-Kutta: Number of Stages = ',I2,
     1       /' ',4X,'             Coefficients     = ',6F10.6)
 1022 FORMAT (5X,'IMSTRT = ',I5,' IMEND = ',I5)
 1023 FORMAT (' ','Illegal solution procedure specified ',A25)
 1024 FORMAT (' ','Specified number of iterations, NITS ',I10,
     1        ' ','greater than allowed, MXHIST ',I10,
     2       /' ','Reset to MXHIST and continue.')
 1124 FORMAT (' ','RDDATA: Error-> Invalid coefficient for entropy fix', 
     1       /' ','                CNTRPY must be 0 <= CNTRPY <= 1.')
 1025 FORMAT (' ',4X,'Restarting from file:',
     1       /' ',9X,A70)
 1026 FORMAT (' ','Illegal specification of force calculation = ',I10,
     1       /' ','     PRINT FORCES = 0 DO NOT calculate forces ',
     2       /' ','     PRINT FORCES = 1 Calculate forces ')
 1027 FORMAT (' ','Illegal specification of source term control = ',I10,
     1       /' ','     SOURCE = 0 DO NOT calculate source terms ',
     2       /' ','     SOURCE = 1 Calculate source terms ')
 1028 FORMAT (10X,'Invalid Implicit Source Term treatment ',A25)
 1210 FORMAT (' ','Required mesh levels ',I3,' exceed maximum ',I3,
     1       /' ','     Increase MXCFS in MAIN or decrease mesh ',
     2            'sequencing.'
     3       /' ','     Skipping next ',I3,' input lines.')
 1220 FORMAT (10X,'Mesh Sequence Level:',I4,' Number of Iterations:',I6)
 1227 FORMAT (10X,A25)
 1229 FORMAT (' ','Illegal specification of FIXQ  pos. control = ',I10,
     1       /' ','     FIXQ   = NONE    -> DO NOT call FIXQ ',
     2       /' ','     FIXQ   = MINIMUM -> Call FIXQ with Minimum ',
     3       /' ','     FIXQ   = AVERAGE -> Call FIXQ with Average ')
 1029 FORMAT (' ','Illegal turbulent constant specified ',A25)
 1030 FORMAT (///' ','***** ILLEGAL FLUX EVALUATION ROUTINE ****')
 1035 FORMAT (20X,3I5)
 1040 FORMAT (10X,'Grid File ',5X,A80)
 1045 FORMAT (10X,A15,6(1X,I4),2X,A15)
 1047 FORMAT (15X,'Profile File: ',5X,A80)
 1048 FORMAT (15X,'Subsonic Profile File: ',5X,A80)
 1049 FORMAT (' ','RDDATA: ERROR-> BC Face not specified. One pair of ',
     1            'BC indices must match',
     2       /' ','                to define the BC Face.')
 1050 FORMAT (' ','INVALID BOUNDARY CONDITION TYPE ',5X,A25)
 1055 FORMAT (10X,'Invalid Option in Block/Zone input ',A25,5X,E14.7)
 1058 FORMAT (' ','Invalid Limiter Specified ',5X,A25,5X,E14.7)
 1059 FORMAT (' ','Invalid Central Difference Option ',5X,A25,5X,E14.7)
 1060 FORMAT (' ','INVALID COMMAND',5X,A25,5X,E14.7)
 1061 FORMAT (10X,A25,5X,E14.7)
 1065 FORMAT (15X,A25)
 1066 FORMAT (10X,'Invalid Turbulent-DIffusion model ',A25)
 1067 FORMAT (10X,'Invalid Pressure-Strain model ',A25)
 1068 FORMAT (10X,'Invalid Compressible Dissipation model ',A25)
 1069 FORMAT (10X,'Invalid Pressure-Dilatation  model ',A25)
 1070 FORMAT (' ',5X,' K-Epsilon Model: ',A30)
 1075 FORMAT (' ',5X,' K-Epsilon Modelling Constants ',
     1       /' ',5X,'      CEPS1        = ',F12.8,
     2       /' ',5X,'      CEPS2        = ',F12.8,
     3       /' ',5X,'      CMU          = ',F12.8,
     4       /' ',5X,'      SIGMAK       = ',F12.8,
     5       /' ',5X,'      SIGMAE       = ',F12.8,
     6       /' ',5X,'      SIGMARHO[-1] = ',F12.8,
     7       /' ',5X,'      A2KE         = ',F12.8,
     8       /' ',5X,'      ALF1         = ',F12.8,
     9       /' ',5X,'      ALF2         = ',F12.8,
     A       /' ',5X,'      ALF3         = ',F12.8,
     B       /' ',5X,'      CPDFRO       = ',F12.8,
     C       /' ',5X,'      CPDFT        = ',F12.8,
     D       /' ',5X,'      CPVELK       = ',F12.8,
     E       /' ',5X,'      TKEINF       = ',F12.8,' * Mach^2',
     F       /' ',5X,'      RMUTNF       = ',F12.8)
 1080 FORMAT (' ',5X,' K-Omega Modelling Constants ',
     1       /' ',5X,'      SIGMAK       = ',F12.8,
     2       /' ',5X,'      SIGMAW       = ',F12.8,
     3       /' ',5X,'      BSTRKW       = ',F12.8,
     4       /' ',5X,'      BKW          = ',F12.8,
     5       /' ',5X,'      CMU          = ',F12.8,
     6       /' ',5X,'      GKW          = ',F12.8,
     7       /' ',5X,'      ROUGHK       = ',F12.8,
     8       /' ',5X,'      SIGMARHO[-1] = ',F12.8,
     9       /' ',5X,'      CPDFRO       = ',F12.8,
     A       /' ',5X,'      CPDFT        = ',F12.8,
     B       /' ',5X,'      CPVELK       = ',F12.8,
     C       /' ',5X,'      TKEINF       = ',F12.8,' * Mach^2',
     D       /' ',5X,'      RMUTNF       = ',F12.8)
 1090 FORMAT (' ',5X,' Reynolds Stress Modelling Constants ',
     1       /' ',5X,'      CEPSLN       = ',F12.8,
     2       /' ',5X,'      CEPS1        = ',F12.8,
     3       /' ',5X,'      CEPS2        = ',F12.8,
     4       /' ',5X,'      CEPS3        = ',F12.8,
     5       /' ',5X,'      CEPS4        = ',F12.8,
     6       /' ',5X,'      CEPS5        = ',F12.8,
     7       /' ',5X,'      CSUBS        = ',F12.8,
     8       /' ',5X,'      CMU          = ',F12.8,
     9       /' ',5X,'      SIGMAK       = ',F12.8,
     A       /' ',5X,'      SIGMAE       = ',F12.8,
     B       /' ',5X,'      SIGMARHO[-1] = ',F12.8,
     C       /' ',5X,'      A2KE         = ',F12.8,
     D       /' ',5X,'      ALF1         = ',F12.8,
     E       /' ',5X,'      ALF2         = ',F12.8,
     F       /' ',5X,'      ALF3         = ',F12.8,
     G       /' ',5X,'      CPDFRO       = ',F12.8,
     H       /' ',5X,'      CPDFT        = ',F12.8,
     I       /' ',5X,'      CPVELK       = ',F12.8,
     J       /' ',5X,'      TKEINF       = ',F12.8,' * Mach^2',
     K       /' ',5X,'      RMUTNF       = ',F12.8)
 1095 FORMAT (' ',10X,' Pressure-Strain Correlation: ',A25)
 1097 FORMAT (' ',10X,'      PSC1         = ',F12.8,
     1       /' ',10X,'      PSC2         = ',F12.8,
     2       /' ',10X,'      PSCW         = ',F12.8,
     3       /' ',10X,'      PSALFA       = ',F12.8,
     4       /' ',10X,'      PSBETA       = ',F12.8,
     5       /' ',10X,'      PSGAMA       = ',F12.8)
 1098 FORMAT (' ',5X,' LES Subgrid-Scale Stress Model: ',A30)
 1099 FORMAT (' ',5X,' Large-Eddy Simulation Model Constants ',
     1       /' ',5X,'      C_s          = ',F12.8)
C
 1100 FORMAT (' ',5X,'Invalid conversion turbulence model specified',
     1               ' = ',A25)
 1105 FORMAT (' ',10X,'Invalid model specified for conversion ',A25)
 1110 FORMAT (' ',5X,'Convert on Restart From ',A25)
C
 1525 FORMAT (/' ','RDDATA: Number of blocks required ',I5,
     1             ' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MXBLKS in MAIN and re-submit.',
     3        /' ', 8X,'   Unable to continue reading input.',
     4        /' ', 8X,'   RUN ABORTING!')
 1530 FORMAT (/' ','RDDATA: Number of cuts   required ',I5,
     1             ' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MXCUTS in MAIN and re-submit.',
     3        /' ', 8X,'   (Storage requires 2*number of cuts)',
     4        /' ', 8X,'   Unable to continue reading input.',
     5        /' ', 8X,'   RUN ABORTING!')
 1535 FORMAT (/' ','RDDATA: Number of BCs    required ',I5,
     1             ' for block ',I3,' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MXBCS  in MAIN and re-submit.',
     3        /' ', 8X,'   Unable to continue reading input.',
     4        /' ', 8X,'   RUN ABORTING!')
 1540 FORMAT (/' ','RDDATA: Number of print  required ',I5,
     1             ' for block ',I3,' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MXPRT  in MAIN and re-submit.',
     3        /' ', 8X,'   Unable to continue reading input.',
     4        /' ', 8X,'   RUN ABORTING!')
 1545 FORMAT (/' ','RDDATA: Number of transition segments required ',I5,
     1             ' exceeds allocated ',I5,
     2        /' ', 8X,'=> Increase MXTRSG in MAIN and re-submit.',
     3        /' ', 8X,'   Unable to continue reading input.',
     4        /' ', 8X,'   RUN ABORTING!')
C
 1580 FORMAT (10X,'Block = ',I5,5X,'Start = ',3I4,5X,'End = ',3I4)
 1590 FORMAT (' ','RDDATA: ERROR-> Starting index is greater than ',
     1            'ending index for transition at:'
     2       /' ','IDIR = ',I2,'Block = ',I5,5X,'Start = ',3I4,5X,
     3            'End = ',3I4)
C
 1600 FORMAT (10X,A15,I5,6(1X,I4))
 1605 FORMAT (10X,A15,I5,9(1X,I4))
 1610 FORMAT (' ','CUT names do not match')
 1620 FORMAT (' ','IORD(',I1,') =',I5,' does not equal 1, 2, or 3.')
 1630 FORMAT (' ','RDDATA: ERROR-> CUT Face not specified. One pair ',
     1            'of CUT indices must match',
     2       /' ','                to define the CUT Face.')
 2000 FORMAT (10X,A25,3E16.7)
 2002 FORMAT (10X,A25,3I5)
 2005 FORMAT (5X,'   Block = ',I4,'  IDIR = ',I4,'  Start = ',3I4,
     1           '  End = ',3I4,'  Increment = ',3I4)
 2010 FORMAT (10X,A25,5X,E14.7)
 2015 FORMAT (/' ','RDDATA: ERROR-> Iterative, implicit time ',
     1             'dependent currently only supports the',
     2        /' ','                AF3F or DIAGONAL AF3F solvers.',
     3        /' ','                Please switch to one of these and ',
     4             'resubmit.')
 2020 FORMAT (/' ','RDDATA: WARNING-> Iterative, implicit time ',
     1             'dependent time step set to DT ',1PE14.7)
 2030 FORMAT (/' ','RDDATA: ERROR-> Iterative, implicit time dependent',
     1             ' requires time step.',
     2        /' ',8X,'Input time step with TIME STEP or DT.')
 2040 FORMAT (/' ','RDDATA: ERROR-> t-ts Iterative, implicit time ',
     1             'dependent requires DT = TIME STEP',
     2        /' ',8X,'Check input TIME STEP and DT.')
 2050 FORMAT (/' ','RDDATA: ERROR-> There were ',I5,' errors in the ',
     1             'input deck.',
     2        /' ', 8X,'RUN ABORTING!')
 2100 FORMAT (//' ','RDDATA: Sorting boundary conditions into ',
     1              'ascending order based on type',//)
 9990 FORMAT (' ','RDDATA: ERROR-> Unexpected End-Of-File encountered.',
     1       /' ', 8X,'RUN ABORTING!')
 9995 FORMAT (' ','RDDATA: ERROR-> Syntax error in input line.'
     1       /' ', 8X,'RUN ABORTING!')
C
      END
