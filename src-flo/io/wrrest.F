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
      SUBROUTINE WRREST (NXBLKS, NXLVLS, NXBCS,  NXSECT,
     2                   IOFFS,  IOFFQ,  IOFFP,  IOFFQN,
     3                   NBLKS,  NLVLS,  IMODEL, IFDELQ, ITRUN,
     4                   IDIM,   JDIM,   KDIM, 
     5                   Q, QN, PROPS, S, VOL, NBCS, IBCDAT,
     6                   RES, QP, QM, DQP, DQM, F, FM, WORKX, STEMP,
     7                   FWRRST, IERRCD)
C
C Routine to write in the restart file for a block.
C
C N.B.  QN is converted to PRIMITIVE variables in restart file. 
C       QN is converted back to conserved variables before exiting this
C          routine to continue calculation.
C
C This routine is organized as follows:
C     0. Setup logicals to control additional output to restart file
C     1. Open restart file
C     2. Write (new) Residual History Data
C     3. Write Restart File Header Information 
C        a. Number of variables
C        b. Number of properties
C        c. Number of additional time levels and time for unsteady
C     4. For each level and block Write Flow Variables and Properties
C     5. Close restart file
C     6. Error handling section
C
C NXBLKS         : Maximum number of blocks dimensioned in main
C NXLVLS         : Maximum number of grid levels dimensioned in main
C                  dimensioned in main
C NXBCS          : Maximum number of boundary conditions dimensioned in main
C NXSECT         : Maximum size of a cross section as dimensioned in main
C NBLKS          : Number of blocks
C NLVLS          : Number of grid levels
C IMODEL         : Data controlling turbulence model
C IFDELQ         : Array controlling calculation of delq
C ITRUN          : Number of iterations to write residual history
C IDIM,JDIM,KDIM : Dimensions of each block
C Q              : Primitive variables stored at cell centers
C QN             : Primitive variables at time levels n, n-1, etc.
C                  N.B. Only NTIME-1 time levels are written/read on restart
C PROPS          : Properties at cell centers
C S              : Metrics
C VOL            : Cell volumes
C NBCS           : Number of boundary conditions for each block
C IBCDAT         : Data controlling boundary conditions
C RES,QP,QM,DQP,DQM,F,FM,WORKX,STEMP
C                : Used for temporary space
C FWRRST         : Restart filename
C IERRCD         : Error code
C
      include '../header/common.h'
      include '../header/histry.h'
C
      DIMENSION IDIM  (NXBLKS,NXLVLS),
     1          JDIM  (NXBLKS,NXLVLS),
     2          KDIM  (NXBLKS,NXLVLS)
C
      DIMENSION IOFFS (NXBLKS,NXLVLS),
     1          IOFFQ (NXBLKS,NXLVLS),
     2          IOFFP (NXBLKS,NXLVLS),
     3          IOFFQN(NXBLKS,NXLVLS)
C
      DIMENSION IMODEL(NMDL)
C
      DIMENSION Q     (*),
     1          QN    (*),
     2          PROPS (*),
     3          S     (*),
     4          VOL   (*)
C
      DIMENSION IBCDAT(10,NXBCS,NXBLKS,NXLVLS),  
     1          NBCS  (NXBLKS)
C
      DIMENSION RES   (*),
     1          QP    (*),
     2          QM    (*),
     3          DQP   (*),
     4          DQM   (*),
     5          F     (*),
     6          FM    (*),
     7          WORKX (*),
     8          STEMP (*)
C
      DIMENSION IFDELQ(NQ)
C
      CHARACTER*80 FWRRST,
     1             FERROR
C
C Control of whether or not to write properties data and additional time levels
C
      LOGICAL   WRPROP, WRUNST
C
C 0.   Setup logicals to control additional output to restart file
C
      ITURB = IMODEL(1)
C
      IF (ITURB  .GE. ITLMNR) THEN
         WRPROP = .TRUE.
      ELSE
         WRPROP = .FALSE.
      ENDIF
C
      IF (ITIMED .EQ. IUNSTD) THEN
         WRUNST = .TRUE.
      ELSE
         WRUNST = .FALSE.
      ENDIF
C
C 1.   Open restart file
C
      IERRCD = INO
      FERROR = FWRRST
      OPEN (IWRRST,FILE=FWRRST,STATUS='UNKNOWN',FORM='UNFORMATTED',
     1             ERR=900)
      REWIND (IWRRST)
C
C 2.   Write (new) Residual History Data
C
      DO 200 IT = 1, ITRUN + ITTOT
         WRITE (IWRRST) IT, R2(IT)
  200 CONTINUE
      ITF = -999
      R2F = -999.99
      WRITE (IWRRST) ITF, R2F
C
C 3.   Write Restart File Header Information
C 3-a. Number of flow variables
C
      WRITE (IWRRST) NQ
C
C 3-b. Number of properties
C
      IF (WRPROP) THEN
         WRITE (IWRRST) NP
      ENDIF
C
C 3-c. Number of additional time levels and time for unsteady
C
      IF (WRUNST) THEN
         WRITE (IWRRST) NTIME-1
         WRITE (IWRRST) TIME
      ENDIF
C
C 4.   For each level and block Write Flow Variables and Properties
C
      WRITE (IOUT,1000)
      DO 420 ILVL = 1, NLVLS
         DO 410 IBLK = 1, NBLKS
            IQ  = IOFFQ (IBLK,ILVL) * NQ          + 1
            IQN = IOFFQN(IBLK,ILVL) * NQ * NTMLVL + 1
            IP  = IOFFP (IBLK,ILVL) * NP          + 1
            IS  = IOFFS (IBLK,ILVL) * 4 * 3       + 1
            IV  = IOFFS (IBLK,ILVL)               + 1
            IF (ITURB .GE. ITLMNR) THEN
               IPSTRT = 1
               IPEND  = IDIM(IBLK,ILVL) + 1
               CALL PRPRTY (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL),
     1                      KDIM(IBLK,ILVL), IPSTRT, IPEND, NXSECT,
     2                      Q(IQ), PROPS(IP), S(IS), VOL(IV),
     3                      IBCDAT(1,1,IBLK,ILVL), NBCS(IBLK), 
     4                      IMODEL, IFDELQ,
     5                      STEMP, QP, F, FM, QM, DQP, DQM, WORKX)
            ENDIF
            IF (ITIMED .EQ. IUNSTD) THEN
               NNPTS = (IDIM(IBLK,ILVL)+3) * (JDIM(IBLK,ILVL)+3) *
     1                 (KDIM(IBLK,ILVL)+3) * (NTIME-1)
               CALL CONC2P (ITURB, NNPTS, QN(IQN), RES)
            ENDIF
C
            WRITE (IOUT,1010) IBLK, ILVL
C
            CALL WRRSTB (ITURB, ITRUN, IBLK, ILVL, IDIM(IBLK,ILVL),
     1                  JDIM(IBLK,ILVL), KDIM(IBLK,ILVL),
     2                  Q(IQ), QN(IQN), PROPS(IP), WRPROP, WRUNST)
C
            IF (ITIMED .EQ. IUNSTD) THEN
               CALL CONP2C (ITURB, NNPTS, QN(IQN), RES)
            ENDIF
  410    CONTINUE
  420 CONTINUE
C
C 5.   Close restart file
C
      CLOSE (IWRRST)
C
C     Finished writing restart file
C
      WRITE (IOUT,1000)
      RETURN
C
C 6.   Error handling section
C
  900 CONTINUE
      IERRCD = IYES
      WRITE (IOUT,1100) FERROR
      RETURN
C
 1000 FORMAT (////' ')
 1010 FORMAT (' ','WRREST: Writing Restart File for Block ',I4,
     1            '    Level ',I4)
 1100 FORMAT (/' ','WRREST: ERROR-> Error opening file ',A80,
     1        /' ', 8X,'RUN ABORTING!')
C
      END
