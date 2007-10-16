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
      SUBROUTINE RDGRID (NBLKS, IDIM, JDIM, KDIM, IOFFS, GRDFIL, R, 
     1                   IDIMIN, JDIMIN, KDIMIN, IGRDTP, GRDBIN, 
     2                   IEOF, IERRCD)
C
C Routine to read in the complete grid data.  This routine will either read
C the grid for each block from a separate file (IGRDTP = IGTSEP) or
C will read the entire grid from a single PLOT3D file in single or multiblock
C format as appropriate (IGRDTP = IGTP3D).
C
C This routine is organized as follows:
C     1. Read Grid
C        a. Read Separate File (IGRDTP = IGTSEP) [formatted or binary]
C           For each block {
C           i.   Open file
C           ii.  Read header information
C           iii. Read block grid 
C           iv.  Close file
C           }
C        b. Read PLOT3D File   (IGRDTP = IGTP3D) [formatted or binary]
C           i.   Open file
C           ii.  Read header information for all blocks
C           iii. Verify dimensions
C           iv.  Loop over all blocks and read grids
C           v.   Close file
C        c. Error: Invalid type
C     2. Error handling section
C        a. End-Of-File Handler
C        b. Error Handler
C
C NBLKS          : Number of blocks
C IDIM,JDIM,KDIM : Dimensions of all blocks on finest level
C IOFFS          : Offset pointer to each block
C GRDFIL         : Grid file name(s)
C R              : Position vector of grid vertices (x,y,z)
C IDIMIN,JDIMIN,
C KDIMIN         : Dimensions of all blocks read in the grid
C IGRDTP         : Grid type (separate files or Plot3D multiblock)
C GRDBIN         : Logical for binary (.TRUE.) or formatted (.FALSE.) file
C IEOF           : End-of-file switch
C                  IEOF   = INO  if End-of-file NOT reached
C                  IEOF   = IYES if End-of-file WAS reached
C IERRCD         : Error Code
C                  IERRCD = INO  if NO error reading file
C                  IERRCD = IYES if    error reading file
C
      include '../header/common.h'
C
      DIMENSION IDIM  (NBLKS),
     1          JDIM  (NBLKS),
     2          KDIM  (NBLKS),
     3          IOFFS (NBLKS)
C
      DIMENSION IDIMIN(NBLKS),
     1          JDIMIN(NBLKS),
     2          KDIMIN(NBLKS)
C
      CHARACTER*80 GRDFIL(NBLKS),
     1             FERROR
C
      DIMENSION R     (*)
C
      LOGICAL GRDBIN
C
      CHARACTER*80 FILTYP
C
C Initialize error handling
C
      IEOF   = INO
      IERRCD = INO
C
C Initialize file type
C
      IF (GRDBIN) THEN
         FILTYP = 'UNFORMATTED'
      ELSE
         FILTYP = 'FORMATTED'
      ENDIF
C
C 1.   Read Grid
C 1-a. Read Separate Files (IGRDTP = IGTSEP) [formatted or binary]
C
      IF (IGRDTP .EQ. IGTSEP) THEN
         DO 190 IBLK = 1, NBLKS
            FERROR = GRDFIL(IBLK)
            OPEN (IGRID,FILE=GRDFIL(IBLK),STATUS='OLD',FORM=FILTYP,
     1                  ERR=920)
            REWIND (IGRID)
            WRITE (IOUT,1000) IBLK, GRDFIL(IBLK)
            IR   = IOFFS(IBLK) * 3     + 1
C
C      Binary  files
C
            IF (GRDBIN) THEN
               IF (THREED) THEN
                  READ (IGRID,END=900,ERR=910) IDIM1, JDIM1, KDIM1
                  IF (IDIM1 .NE. IDIM(IBLK) .OR. 
     1                JDIM1 .NE. JDIM(IBLK) .OR.
     2                KDIM1 .NE. KDIM(IBLK)) THEN
                     WRITE (IOUT,1100) IBLK, 
     1                                 IDIM(IBLK),JDIM(IBLK),KDIM(IBLK),
     2                                 IDIM1,     JDIM1,     KDIM1
                     IERRCD = IYES
                     RETURN
                  ENDIF
                  CALL RDP3D  (GRDBIN, IDIM(IBLK), JDIM(IBLK), 
     1                         KDIM(IBLK), R(IR), IEOF, IERRCD)
               ELSE
                  READ (IGRID,END=900,ERR=910) IDIM1, JDIM1
                  IF (IDIM1 .NE. IDIM(IBLK) .OR. 
     1                JDIM1 .NE. JDIM(IBLK)) THEN
                     WRITE (IOUT,1105) IBLK, IDIM(IBLK), JDIM(IBLK), 
     1                                       IDIM1,      JDIM1
                     IERRCD = IYES
                     RETURN
                  ENDIF
                  CALL RDP3D  (GRDBIN, IDIM(IBLK), JDIM(IBLK), 
     1                         KDIM(IBLK), R(IR), IEOF, IERRCD)
               ENDIF
C
C      Formatted files
C
            ELSE
               IF (THREED) THEN
                  READ (IGRID,*,END=900,ERR=910) IDIM1, JDIM1, KDIM1
                  IF (IDIM1 .NE. IDIM(IBLK) .OR. 
     1                JDIM1 .NE. JDIM(IBLK) .OR.
     2                KDIM1 .NE. KDIM(IBLK)) THEN
                     WRITE (IOUT,1100) IBLK, 
     1                                 IDIM(IBLK),JDIM(IBLK),KDIM(IBLK),
     2                                 IDIM1,     JDIM1,     KDIM1
                     IERRCD = IYES
                     RETURN
                  ENDIF
                  CALL RDP3D  (GRDBIN, IDIM(IBLK), JDIM(IBLK), 
     1                         KDIM(IBLK), R(IR), IEOF, IERRCD)
               ELSE
                  READ (IGRID,*,END=900,ERR=910) IDIM1, JDIM1
                  IF (IDIM1 .NE. IDIM(IBLK) .OR. 
     1                JDIM1 .NE. JDIM(IBLK)) THEN
                     WRITE (IOUT,1105) IBLK, IDIM(IBLK), JDIM(IBLK), 
     1                                       IDIM1,      JDIM1
                     IERRCD = IYES
                     RETURN
                  ENDIF
                  CALL RDP3D  (GRDBIN, IDIM(IBLK), JDIM(IBLK), 
     1                         KDIM(IBLK), R(IR), IEOF, IERRCD)
               ENDIF
            ENDIF
            CLOSE (IGRID)
            IF (IERRCD .NE. INO .OR. IEOF .NE. INO) RETURN
  190    CONTINUE
C
C 1-b. Read PLOT3D File   (IGRDTP = IGTP3D) [formatted or binary]
C 1-b-i.   Open file
C
      ELSE IF (IGRDTP .EQ. IGTP3D) THEN
         FERROR = GRDFIL(1)
         OPEN (IGRID,FILE=GRDFIL(1),STATUS='OLD',FORM=FILTYP,
     1               ERR=920)
         REWIND (IGRID)
C
C 1-b-ii.  Read grid header file
C
C      Binary
C
         IF (GRDBIN) THEN
C
C      Read multiblock data if more than 1 block
C
            IF (NBLKS .GT. 1) THEN
               READ (IGRID,END=900,ERR=910) NBLKIN
               IF (NBLKIN .NE. NBLKS) THEN
                  IERRCD = IYES
                  WRITE (IOUT,1110) NBLKS, NBLKIN
                  RETURN
               ENDIF
            ENDIF
C
C      3D or 2D
C
            IF (THREED) THEN
               READ (IGRID,END=900,ERR=910) (IDIMIN(L), JDIMIN(L), 
     1                                       KDIMIN(L), L=1,NBLKS)
            ELSE
               READ (IGRID,END=900,ERR=910) (IDIMIN(L), JDIMIN(L), 
     1                                       L=1,NBLKS)
            ENDIF
C
C      Formatted
C
         ELSE
C
C      Read multiblock data if more than 1 block
C
            IF (NBLKS .GT. 1) THEN
               READ (IGRID,*,END=900,ERR=910) NBLKIN
               IF (NBLKIN .NE. NBLKS) THEN
                  IERRCD = IYES
                  WRITE (IOUT,1110) NBLKS, NBLKIN
                  RETURN
               ENDIF
            ENDIF
C
C      3D or 2D
C
            IF (THREED) THEN
               READ (IGRID,*,END=900,ERR=910) (IDIMIN(L), JDIMIN(L), 
     1                                         KDIMIN(L), L=1,NBLKS)
            ELSE
               READ (IGRID,*,END=900,ERR=910) (IDIMIN(L), JDIMIN(L), 
     1                                         L=1,NBLKS)
            ENDIF
         ENDIF
C
C 1-b-iii. Verify dimensions
C
         IF (THREED) THEN
            DO 210 IBLK = 1, NBLKS
               IF (IDIM(IBLK) .NE. IDIMIN(IBLK) .OR.
     1             JDIM(IBLK) .NE. JDIMIN(IBLK) .OR.
     2             KDIM(IBLK) .NE. KDIMIN(IBLK) ) THEN
                  IERRCD = IYES
                  WRITE (IOUT,1100) IBLK, 
     1                           IDIM  (IBLK),JDIM  (IBLK),KDIM  (IBLK),
     2                           IDIMIN(IBLK),JDIMIN(IBLK),KDIMIN(IBLK)
               ENDIF
  210       CONTINUE
         ELSE
            DO 220 IBLK = 1, NBLKS
               IF (IDIM(IBLK) .NE. IDIMIN(IBLK) .OR.
     1             JDIM(IBLK) .NE. JDIMIN(IBLK) ) THEN
                  IERRCD = IYES
                  WRITE (IOUT,1105) IBLK, 
     1                              IDIM  (IBLK),JDIM  (IBLK),
     2                              IDIMIN(IBLK),JDIMIN(IBLK)
               ENDIF
  220       CONTINUE
         ENDIF
         IF (IERRCD .NE. INO .OR. IEOF .NE. INO) RETURN
C
C 1-b-iv.  Loop over all blocks reading grid
C
         DO 290 IBLK = 1, NBLKS
            WRITE (IOUT,1000) IBLK, GRDFIL(1)
            IR   = IOFFS(IBLK) * 3     + 1
            CALL RDP3D  (GRDBIN, IDIM(IBLK), JDIM(IBLK), 
     1                   KDIM(IBLK), R(IR), IEOF, IERRCD)
            IF (IERRCD .NE. INO .OR. IEOF .NE. INO) RETURN
  290    CONTINUE
C
C 1-b-v.   Close file
C
         CLOSE (IGRID)
C
C 1-c. Error
C
      ELSE
         IERRCD = IYES
         WRITE (IOUT,1120) IGRDTP
         RETURN
      ENDIF
C
C     Succesful completion reading grid - return
C
      RETURN
C
C 2.   Error handling section
C 2-a. End-Of-File Handler
C
  900 CONTINUE
      IEOF   = IYES
      WRITE (IOUT,1200) 
      RETURN
C
C 2-b. Error Handler
C
  910 CONTINUE
      IERRCD = IYES
      WRITE (IOUT,1210) 
      RETURN
C
C 2-c. Open Error Handler
C
  920 CONTINUE
      IERRCD = IYES
      WRITE (IOUT,1220) FERROR
      RETURN
C
 1000 FORMAT (' ','RDGRID: Reading grid for block ',I3,' from file ',
     1            A80)
 1100 FORMAT (' ','RDGRID: ERROR-> Reading input 3D grid for Block ',I3,
     1       /' ','                     IDIM, JDIM, KDIM = ',3I6,
     2       /' ','                Grid IDIM, JDIM, KDIM = ',3I6,
     3       /' ', 8X,'RUN ABORTING!')
 1105 FORMAT (' ','RDGRID: ERROR-> Reading input 2D grid for Block ',I3,
     1       /' ','                     IDIM, JDIM = ',2I6,
     2       /' ','                Grid IDIM, JDIM = ',2I6,
     3       /' ', 8X,'RUN ABORTING!')
 1110 FORMAT (' ','RDGRID: ERROR-> Reading input grid number of blocks',
     1       /' ','                     NBLKS = ',I6,
     2       /' ','                Grid NBLKS = ',I6,
     3       /' ', 8X,'RUN ABORTING!')
 1120 FORMAT (' ','RDGRID: ERROR-> Invalid grid type specified ',
     1            'IGRDTP = ',I6,
     2       /' ', 8X,'RUN ABORTING!')
C
 1200 FORMAT (' ','RDGRID: ERROR-> Unexpected End-Of-File Encountered ',
     1       /' ', 8X,'RUN ABORTING!')
 1210 FORMAT (' ','RDGRID: ERROR-> Unknown Error Encountered ',
     1       /' ', 8X,'RUN ABORTING!')
 1220 FORMAT (' ','RDGRID: ERROR-> Error opening file ',A80,
     1       /' ', 8X,'RUN ABORTING!')
C
      END
