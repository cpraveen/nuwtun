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
      SUBROUTINE RDREST (ITURB, IBLOCK, ILVL, IDIM, JDIM, KDIM, 
     1                   Q, QN, PROPS, IEOF,  IERRCD)
C
C Routine to read in the restart file for a block.
C
C Note: To read a restart file from ISAAC versions earlier than 3.5.15,
C       set RDPROP = FALSE.
C
C N.B.  QN is converted to PRIMITIVE variables for restart file.
C
C This routine is organized as follows:
C     If (IBLOCK = 1 AND ILVL = 1) 
C     {
C        1. Read Residual History Data 
C           a. Read input residual history data
C           b. Copy residual history data to output restart file
C           c. Echo residual history to listing
C        2. Read Restart File Header Information 
C           a. Number of variables
C           b. Number of properties
C           c. Number of additional time levels and time for unsteady
C     }
C     3. Read Flow Variables and Properties
C        a. Read Flow Variables 
C        b. Read Properties
C        c. Read flow variables at additional time levels for unsteady flow
C     4. Error handling section
C        a. End-Of-File Handler
C        b. Error Handler
C
C ITURB          : Turbulence model
C IBLOCK         : Block number
C ILVL           : Grid level
C IDIM,JDIM,KDIM : Dimensions of current block
C Q              : Primitive variables stored at cell centers
C QN             : Primitive variables at time levels n, n-1, etc.
C                  N.B. Only NTIME-1 time levels are written/read on restart
C PROPS          : Properties at cell centers
C IEOF           : End-of-file switch
C                  IEOF   = INO  if End-of-file NOT reached
C                  IEOF   = IYES if End-of-file WAS reached
C IERRCD         : Error Code
C                  IERRCD = INO  if NO error reading file
C                  IERRCD = IYES if    error reading file
C
      include '../header/common.h'
      include '../header/histry.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QN    (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ,NTIME-1),
     2          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP)
C
C Control of whether or not to read properties data and additional time levels
C
      LOGICAL   RDPROP, RDUNST
C
C Initialize error switches
C
      IEOF   = INO
      IERRCD = INO
C
      IF (ITURB  .GE. ITLMNR) THEN
         RDPROP = .TRUE.
      ELSE
         RDPROP = .FALSE.
      ENDIF
C
      IF (ITIMED .EQ. IUNSTD .AND. .NOT.(INITTD)) THEN
         RDUNST = .TRUE.
      ELSE
         RDUNST = .FALSE.
      ENDIF
C
C 1.  Read residual history and file header only for the first block
C
      IF (IBLOCK .EQ. 1 .AND. ILVL .EQ. 1) THEN
C
C     Initial residual output increment
C
         INC   = 1
C
C     Initialize variables in case restart with no residual history data
C
         ITOLD = 0
         R2OLD = 0.E0
         R2ONE = 1.E0
C
C 1-a. Read Residual History Data
C 1-b. Copy Residual History Data to output restart file
C
  100    CONTINUE
         READ (IRDRST,END=900,ERR=910)  IT, R2T
         IF (IT .LE. 0) GO TO 190
         ITOLD  = IT
         R2OLD  = R2T
         IF (IT .LT. MXHIST) THEN
            R2(IT) = R2T
         ELSE
            IERRCD = IYES
         ENDIF
         IF (IT/INC*INC .EQ. IT) THEN
            WRITE (IOUT,1000) IT, R2T
         ENDIF
         IF (IT .GT. 10)  INC = 10
         IF (IT .GT. 100) INC = 100
         GO TO 100
C
C 1-c. Echo the last residual to listing
C
  190    CONTINUE
         ITTOT = ITOLD 
         WRITE (IOUT,1000) ITOLD, R2OLD 
C
         IF (IERRCD .NE. INO) THEN
            WRITE (IOUT,1010) ITTOT, MXHIST
            RETURN
         ENDIF
C
         R2ONE = R2(1)
C
C 2.   Read Restart File Header Information
C 2-a. NQ
C
         READ (IRDRST,END=900,ERR=910) NQIN
C
         IF (NQIN .NE. NQ) THEN
            WRITE (IOUT,1100) NQIN, NQ
            STOP
         ENDIF
C
C 2-b. NP (not if inviscid)
C
         IF (RDPROP) THEN
            READ (IRDRST,END=900,ERR=910) NPIN
C
            IF (NPIN .NE. NP) THEN
               WRITE (IOUT,1105) NPIN, NP
               STOP
            ENDIF
         ENDIF
C
C 2-c. Number of additional time levels and time for unsteady
C
         IF (RDUNST) THEN
            READ (IRDRST,END=900,ERR=910) NTIMIN
            READ (IRDRST,END=900,ERR=910) TIME
C
            IF (NTIMIN .NE. (NTIME-1)) THEN
               WRITE (IOUT,1110) NTIMIN, NTIME-1
               STOP
            ENDIF
         ENDIF
C
C     End of history and file header section
C
      ENDIF
C
C 3.   Read Flow Variables and Properties
C 3-a. Flow Variables (Q)
C
      READ (IRDRST,END=900,ERR=910) Q
C
C 3-b. Properites (PROPS)
C
      IF (RDPROP) THEN
         READ (IRDRST,END=900,ERR=910) PROPS
      ENDIF
C
C 3-c. Read flow variables at additional time levels for unsteady flow
C
      IF (RDUNST) THEN
         READ (IRDRST,END=900,ERR=910) QN
      ENDIF
C
C     Succesful completion reading restart - return
C
      RETURN
C
C 4.   Error handling section
C 4-a. End-Of-File Handler
C
  900 CONTINUE
      IEOF   = IYES
      RETURN
C
C 4-b. Error Handler
C
  910 CONTINUE
      IERRCD = IYES
      WRITE (IOUT,1900)
      RETURN
C
 1000 FORMAT (' ','   ITER = ',I5,'   R2 = ',E15.8)
 1010 FORMAT (' ','RDREST: ERROR-> Length of residual history in ',
     1            'restart file (ITTOT) ',I6,
     2       /' ',15X,' exceeds maximum allocated      (MXHIST ',
     3            'in hstry.h) ',I6,
     4       /' ',15X,' Increase MXHIST to be greater than ',
     5                'ITTOT + NITS and recompile.',
     6       /' ', 8X,'RUN ABORTING!')
 1100 FORMAT (' ','RDREST: ERROR-> Number of variables  (NQ) in ',
     1            'restart file',I5,
     2       /' ',15X,' does not match number specified for run  ',I5,
     3       /' ', 8X,'Check input for incorrect model and resubmit.',
     4       /' ', 8X,'RUN ABORTING!')
 1105 FORMAT (' ','RDREST: ERROR-> Number of properties (NP) in ',
     1            'restart file',I5,
     2       /' ',15X,' does not match number specified for run  ',I5,
     3       /' ', 8X,'Check input for incorrect model and resubmit.',
     4       /' ', 8X,'RUN ABORTING!')
 1110 FORMAT (' ','RDREST: ERROR-> Number of additional time levels ',
     1            '(NTIME-1) in restart file',I5,
     2       /' ',15X,' does not match number specified for run  ',I5,
     3       /' ', 8X,'This may be caused by attempting to read a ',
     4                'steady restart file', 
     5       /' ', 8X,'without specifying the INITIALIZE TIME ',
     6                'DEPENDENT input option to ',
     7       /' ', 8X,'convert to an unsteady solution.',
     7      //' ', 8X,'Check input for incorrect scheme and resubmit.',
     8       /' ', 8X,'RUN ABORTING!')
 1900 FORMAT (' ','RDREST: ERROR-> Unclassifiable error reading ',
     1            'restart file.',
     2       /' ', 8X,'Check specified restart file and retry.',
     3       /' ', 8X,'RUN ABORTING!')
      END
