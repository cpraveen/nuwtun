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

      BLOCK DATA COMMNS
C
C     Block Data routine to initialize the common data.
C
      include '../header/common.h'
C
C     Data statements for no and yes
C
      DATA INO    /0/,
     1     IYES   /1/
C
C     Data statements for constants
C
      DATA RMAX   /9.99E+20/,
     1     RSMALL /1.00E-12/,
     2     RSMASQ /1.00E-24/
C
C     Data statements to set up the units for input and output.
C
      DATA INPUT   /5/,
     1     IOUT    /6/, 
     2     IRDRST  /8/, 
     3     IWRRST  /9/, 
     4     IGRID  /10/,
     5     IPLT3G /11/,
     6     IPLT3Q /12/,
     7     IPLT3F /13/,
     8     IPLT3FN/14/,
     9     IRDBC  /15/,
     A     ISTOPFL/16/,
     B     IRESID /17/,
     C     ICLCD  /18/,
     C     IRESULT/19/
C
C     Data statements to set up the file types supported for grids
C
      DATA IGTSEP  /1/,
     1     IGTP3D  /2/
C
C     Data statements to set up the available boundary conditions.
C     Boundary conditions are sorted based on increasing value of this 
C     variable, therefore walls are done near last, etc.
C
      DATA ITAN   /22/,
     1     IWALL  /23/,
     2     IWALFN /21/,
     3     IFAR    /4/,
     4     IFAR2D  /5/,
     5     IEXT    /6/,
     6     IFIX    /7/,
     7     ISING   /8/,
     8     IJET    /9/,
     9     IPERD  /11/,
     A     IPROFL  /1/,
     B     IPROSB  /2/,
     C     ISYMYZ /12/,
     D     ISYMXZ /13/,
     E     ISYMXY /14/,
     F     IINFLO /15/,
     G     IOUTFL /16/,
     H     IAXISM /17/,
     I     IAXICL /18/,
     J     IHSHR  /10/
C
C     Data statements to set up the global constants defining the
C     available flux evaluation schemes.  The schemes should be
C     numbered from 1 to the maximum number to work correctly with
C     the IFCHAR array.
C
      DATA IFROE  /1/
C
C     Data statements to set the available types of flux evaluation
C     schemes.
C
      DATA IFDS   /1/,
     1     IFVS   /2/,
     2     IFCD   /0/
C
C     Data statements to set up available limiters
C
      DATA ILNONE /0/, 
     1     ILSMTH /1/, 
     2     ILMNMD /2/, 
     3     ILVNKT /3/
C
C     Data statements to set up the global constants defining the 
C     available turbulence models.
C
      DATA ITNVSD     /0/,
     1     ITLMNR    /10/,
     2     ITBLMX   /500/,
     3     ITKE    /2100/,
     4     ITKW    /2200/,
     5     ITRS   /10000/,
     6     ITLES    /100/
C
C     Data statements for Fmu models
C
      DATA IFMHR     /0/,
     1     IFMSAA    /1/,
     2     IFMZSG    /2/,
     3     IFMVD     /3/
C
C     Data statements for epsilon equation models
C
      DATA IEEHR     /0/,
     1     IEESAA    /1/,
     2     IEEZSG    /2/,
     3     IEEZSL    /3/,
     4     IEERNG    /4/,
     5     IEES95    /5/,
     6     IEEABD    /6/
C
C     Data statements for compressible dissipation models
C
      DATA IECNO     /0/,
     1     IECSAR    /1/,
     2     IECRIS    /2/
C
C     Data statements for Pressure-Dilatation models
C
      DATA IPDNO     /0/,
     1     IPDSAR    /1/,
     2     IPDRIS    /2/
C
C     Data statements for Anisotropic Dissipation Rate Models
C
      DATA IADRNO    /0/,
     1     IADRGS    /1/
C
C     Data statements for Algebraic stress models
C
      DATA IASMBU    /0/,
     1     IASMGS    /1/
C
C     Data statements to set up pressure-strain correlation models.
C
      DATA IPLRR1    /1/,
     1     IPLRR2    /2/,
     2     IPSSG     /3/,
     3     IPGL      /4/,
     4     IPLS      /5/,
     5     IPFLT     /6/,
     6     IPSO95    /7/
C
C     Data statements for Epsilon_ij models
C
      DATA IEISO     /0/,
     1     IELAIS    /1/,
     2     IEFLT     /2/,
     3     IEADGS    /3/,
     4     IESO95    /4/
C
C     Data statements to set molecular diffusion models
C
      DATA IDMSO     /1/,
     1     IDMSPZ    /2/
C
C     Data statements to set turbulent diffusion models
C
      DATA IDTDH     /0/,
     1     IDTHL     /1/,
     2     IDTMH     /2/
C
C     Data statements to set SGS models
C
      DATA ISGSMG    /1/,
     1     ISGSFM    /2/
C
C     Data statements to set the FIX routines for positivity preservation
C     IFIXMN, IFIXAV cannot equal INO
C
      DATA IFIXMN    /1/,
     1     IFIXAV    /2/
C
C     Data statements to set the available types of solution procedure.
C
      DATA IRKN    /1/,
     1     IAF3F   /2/,
     2     ILGS    /3/,
     3     IMARCH /10/
C
      DATA ISDIAG  /1/,
     1     ISBLOC  /2/
C
      DATA ISTDY   /1/,
     1     IUNSTD  /2/
C
      DATA ITTS    /1/,
     1     ITAUTS  /2/
C
C     Data statements to set the Kronecker Delta array and Cyclic array
C
      DATA IKD / 1, 0, 0,
     1           0, 1, 0,
     2           0, 0, 1 /
      DATA ICY / 1, 2, 3,
     1           2, 3, 1,
     2           3, 1, 2 /
C
      END






      SUBROUTINE CHKDAT (NBLKS,  IDIM, JDIM, KDIM, VISCOS,
     1                   NXBCS, NBCS, IBCDAT, NCUTS, ICUTS, CUTNAM,
     2                   NXPRT, NPRINT, IPRINT, IERRCD)
C
C Routine to check input data for boundary conditions, cuts, and print output.
C This routine looks at the finest grid level only.  All other data will be
C derived from the fine grid level.  The approach is to output an error, then
C correct the data (the best that can be done), and then pass this back to
C the main routine so that it can continue checking additional data.  This
C routine can only be called after NBLKS, NCUTS, and NPRINT have been tested
C and corrected to fit in memory.  If NBLKS, NCUTS, or NPRINT exceeded their
C dimensions, the output from this routine could be corrupt as RDDATA may have
C overwritten key data.
C
C This routine is organized as follows:
C     0. Verify KDIM for non-3D cases
C     1. Verify Boundary Condition (BC) data
C     2. Verify Cut data
C        a. Dimensions
C        b. Cut lengths
C     3. Verify Print data
C
C NBLKS          : Number of blocks
C IDIM,JDIM,KDIM : Dimensions of blocks
C VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow in block
C NXBCS          : Maximum number of bcs    per block as dimensioned in main
C NBCS           : Number of boundary conditions in each block
C IBCDAT         : Boundary condition data
C NCUTS          : Number of cuts
C ICUTS          : Cut data
C NXPRT          : Maximum number of prints per block as dimensioned in main
C NPRINT         : Number of print statements in each block
C IPRINT         : Print data
C IERRCD         : Error code
C
      include '../header/common.h'
C
      DIMENSION IDIM   (NBLKS), 
     1          JDIM   (NBLKS),
     2          KDIM   (NBLKS)
C
      LOGICAL   VISCOS (3,NBLKS)
C
      DIMENSION IBCDAT (10,NXBCS,NBLKS), 
     1          NBCS   (NBLKS)
C
      DIMENSION ICUTS  (21,NCUTS)
      CHARACTER*20 CUTNAM(NCUTS)
C
      DIMENSION IPRINT (10,NXPRT,NBLKS), 
     1          NPRINT (NBLKS)
C
C Local storage used to test indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),
     2          IEND  (3),
     3          IORD  (3),
     4          ISTRT2(3),
     5          IEND2 (3),
     6          INC   (3),
     7          ICOUNT(3)
C
C Output data
C
      CHARACTER*20 BCNAME
      CHARACTER*1  DIR (3)
      DATA         DIR /'I', 'J', 'K'/
C
C Initialize error code
C
      IERRCD = INO
C
      WRITE (IOUT,1000)
C
C 0.   Verify KDIM for non-3D cases
C
      IF (.NOT. THREED) THEN
         DO 10 IBLK = 1, NBLKS
            IF (KDIM(IBLK) .NE. 2) THEN
               IERRCD = IYES
               WRITE (IOUT,1010) IBLK, IDIM(IBLK), JDIM(IBLK), 
     1                           KDIM(IBLK)
            ENDIF
   10    CONTINUE
      ENDIF
C
C 1.   Verify Boundary Condition (BC) data
C
      DO 190 IBLK = 1, NBLKS
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
         IJKDIM(1) = IDIM (IBLK)
         IJKDIM(2) = JDIM (IBLK)
         IJKDIM(3) = KDIM (IBLK)
C
         DO 170 IBC = 1, NBCS(IBLK)
C
C Store boundary condition data in mnemonic names
C
            IBCTYP   = IBCDAT(1,IBC,IBLK)
            IDIR     = IBCDAT(2,IBC,IBLK)
            INROUT   = IBCDAT(3,IBC,IBLK)
            ISTRT(1) = IBCDAT(4,IBC,IBLK)
            IEND (1) = IBCDAT(5,IBC,IBLK)
            ISTRT(2) = IBCDAT(6,IBC,IBLK)
            IEND (2) = IBCDAT(7,IBC,IBLK)
            ISTRT(3) = IBCDAT(8,IBC,IBLK)
            IEND (3) = IBCDAT(9,IBC,IBLK)
C
C Modify starting point to switch back to input form for testing
C
            DO 140 L = 1, 3
               IF (L .NE. IDIR) THEN
                  ISTRT(L) = ISTRT(L) - 1
               ENDIF
  140       CONTINUE
C
            DO 150 L = 1, 3
               IF (ISTRT(L) .LT. 1 .OR. IEND(L) .GT. IJKDIM(L)) THEN
                  IERRCD = IYES
                  CALL GETNAM (IBCTYP, BCNAME)
                  WRITE (IOUT,1110) IBLK, DIR(L), IJKDIM(L),
     1                     (ISTRT(L1),L1=1,3), (IEND(L2),L2=1,3), BCNAME
               ENDIF
  150       CONTINUE
C
C Verify that viscous terms are enabled in direction normal to all wall or
C wall function boundary conditions.
C
            IF (IBCTYP .EQ. IWALL .OR. IBCTYP .EQ. IWALFN) THEN
               IF (.NOT. VISCOS(IDIR,IBLK)) THEN
                  IERRCD = IYES
                  CALL GETNAM (IBCTYP, BCNAME)
                  WRITE (IOUT,1120) IBLK, DIR(IDIR), 
     1                     (ISTRT(L1),L1=1,3), (IEND(L2),L2=1,3), BCNAME
               ENDIF
            ENDIF
C
C Verify that the TWOD  option is     set for 2D farfield point vortex bc 
C Verify that the MARCH option is not set for 2D farfield point vortex bc
C
            IF (IBCTYP .EQ. IFAR2D) THEN
               IF (THREED .OR. AXISYM) THEN
                  IERRCD = IYES
                  CALL GETNAM (IBCTYP, BCNAME)
                  WRITE (IOUT,1130) IBLK, DIR(IDIR), 
     1                     (ISTRT(L1),L1=1,3), (IEND(L2),L2=1,3), BCNAME
               ENDIF
               IF (ISOLVR .EQ. IMARCH) THEN
                  IERRCD = IYES
                  CALL GETNAM (IBCTYP, BCNAME)
                  WRITE (IOUT,1140) IBLK, DIR(IDIR), 
     1                     (ISTRT(L1),L1=1,3), (IEND(L2),L2=1,3), BCNAME
               ENDIF
            ENDIF
C
  170    CONTINUE
  190 CONTINUE
C
C 2.   Verify Cut data 
C 2-a. Dimensions
C
      DO 240 ICUT = 1, NCUTS
         IBLK     = ICUTS ( 1,ICUT)
         IDIR     = ICUTS ( 2,ICUT)
         INROUT   = ICUTS ( 3,ICUT)
         ISTRT(1) = ICUTS ( 4,ICUT)
         IEND (1) = ICUTS ( 5,ICUT)
         ISTRT(2) = ICUTS ( 6,ICUT)
         IEND (2) = ICUTS ( 7,ICUT)
         ISTRT(3) = ICUTS ( 8,ICUT)
         IEND (3) = ICUTS ( 9,ICUT)
C
C Modify starting/ending point to switch back to input form for testing
C
         DO 220 L = 1, 3
            IF (L .NE. IDIR) THEN
               IF (IEND(L) .GT. ISTRT(L)) THEN
                  ISTRT(L) = ISTRT(L) - 1
               ELSE
                  IEND (L) = IEND (L) - 1
               ENDIF
            ENDIF
  220    CONTINUE
C
C Check the block specified - if invalid, bypass next test
C
         IF (IBLK .GT. NBLKS) THEN
            IERRCD = IYES
            WRITE (IOUT,1200) IBLK, NBLKS, CUTNAM(ICUT),
     1                  (ISTRT(L1),L1=1,3), (IEND(L2),L2=1,3)
            GO TO 239
         ENDIF
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
         IJKDIM(1) = IDIM (IBLK)
         IJKDIM(2) = JDIM (IBLK)
         IJKDIM(3) = KDIM (IBLK)
C
C Test indices
C
         DO 230 L = 1, 3
            IF (ISTRT(L) .LT. 1 .OR. IEND(L) .GT. IJKDIM(L)) THEN
               IERRCD = IYES
               WRITE (IOUT,1210) IBLK, DIR(L), IJKDIM(L), CUTNAM(ICUT),
     1                  (ISTRT(L1),L1=1,3), (IEND(L2),L2=1,3)
            ENDIF
  230    CONTINUE
  239    CONTINUE
  240 CONTINUE
C
C 2-b. Cut lengths
C
      DO 290 ICUT = 1, NCUTS, 2
         IBLK      = ICUTS ( 1,ICUT)
         IDIR      = ICUTS ( 2,ICUT)
         INROUT    = ICUTS ( 3,ICUT)
         ISTRT (1) = ICUTS ( 4,ICUT)
         IEND  (1) = ICUTS ( 5,ICUT)
         ISTRT (2) = ICUTS ( 6,ICUT)
         IEND  (2) = ICUTS ( 7,ICUT)
         ISTRT (3) = ICUTS ( 8,ICUT)
         IEND  (3) = ICUTS ( 9,ICUT)
C
         IBLK2     = ICUTS (10,ICUT)
         IDIR2     = ICUTS (11,ICUT)
         INROUT2   = ICUTS (12,ICUT)
         ISTRT2(1) = ICUTS (13,ICUT)
         IEND2 (1) = ICUTS (14,ICUT)
         ISTRT2(2) = ICUTS (15,ICUT)
         IEND2 (2) = ICUTS (16,ICUT)
         ISTRT2(3) = ICUTS (17,ICUT)
         IEND2 (3) = ICUTS (18,ICUT)
C
         IORD  (1) = ICUTS (19,ICUT)
         IORD  (2) = ICUTS (20,ICUT)
         IORD  (3) = ICUTS (21,ICUT)
C
C Verify order array contains one each of 1, 2, 3
C
         DO 245 L = 1, 3
            ICOUNT(L) = 0
  245    CONTINUE
         DO 250 L = 1, 3
            ICOUNT(IORD(L)) = ICOUNT(IORD(L)) + 1
  250    CONTINUE
C
C If IORD is invalid, bypass next test
C
         IF (ICOUNT(1) .NE. 1 .OR. ICOUNT(2) .NE. 1 .OR. 
     1       ICOUNT(3) .NE. 1) THEN
            IERRCD = IYES
            WRITE (IOUT,1220) (IORD  (L0),L0=1,3),
     1               CUTNAM(ICUT), IBLK,
     2               (ISTRT (L1),L1=1,3), (IEND (L2),L2=1,3),
     3               CUTNAM(ICUT), IBLK2,
     4               (ISTRT2(L3),L3=1,3), (IEND2(L4),L4=1,3)
            GO TO 289
         ENDIF
C
C Modify starting/ending point to switch back to input form for testing
C
         DO 270 L = 1, 3
            IF (L .NE. IDIR) THEN
               IF (IEND(L) .GT. ISTRT(L)) THEN
                  ISTRT(L) = ISTRT(L) - 1
               ELSE
                  IEND (L) = IEND (L) - 1
               ENDIF
            ENDIF
            IF (L .NE. IDIR2) THEN
               IF (IEND2(L) .GT. ISTRT2(L)) THEN
                  ISTRT2(L) = ISTRT2(L) - 1
               ELSE
                  IEND2 (L) = IEND2 (L) - 1
               ENDIF
            ENDIF
  270    CONTINUE
C
C Check the block specified - if invalid, bypass test
C
         IF (IBLK .GT. NBLKS .OR. IBLK2 .GT. NBLKS) THEN
            GO TO 289
         ENDIF
C
C Check length
C
         DO 280 L = 1, 3
            LEN1 = ABS (IEND(L)        - ISTRT(L)       )
            LEN2 = ABS (IEND2(IORD(L)) - ISTRT2(IORD(L)))
            IF (LEN1 .NE. LEN2) THEN
               IERRCD = IYES
               WRITE (IOUT,1230) IBLK, DIR(L), IBLK2, DIR(IORD(L)), 
     1                  CUTNAM(ICUT), 
     2                  (ISTRT (L1),L1=1,3), (IEND (L2),L2=1,3),
     3                  CUTNAM(ICUT), 
     4                  (ISTRT2(L3),L3=1,3), (IEND2(L4),L4=1,3),
     5                  (IORD  (L5),L5=1,3)
            ENDIF
  280    CONTINUE
C
  289    CONTINUE
  290 CONTINUE
C
C 3.   Verify Print data
C
      DO 390 IBLK = 1, NBLKS
C
C Store IDIM, JDIM, KDIM in IJKDIM 
C Increment IJKDIM since it is allowed to print out the boundary data
C
         IJKDIM(1) = IDIM (IBLK) + 1
         IJKDIM(2) = JDIM (IBLK) + 1
         IJKDIM(3) = KDIM (IBLK) + 1
C
         DO 380 IPR = 1, NPRINT(IBLK)
C
C Store print data in mnemonic names
C
            IDIR     = IPRINT( 1,IPR,IBLK)
            ISTRT(1) = IPRINT( 2,IPR,IBLK)
            ISTRT(2) = IPRINT( 3,IPR,IBLK)
            ISTRT(3) = IPRINT( 4,IPR,IBLK)
            IEND (1) = IPRINT( 5,IPR,IBLK)
            IEND (2) = IPRINT( 6,IPR,IBLK)
            IEND (3) = IPRINT( 7,IPR,IBLK)
            INC  (1) = IPRINT( 8,IPR,IBLK)
            INC  (2) = IPRINT( 9,IPR,IBLK)
            INC  (3) = IPRINT(10,IPR,IBLK)
C
            DO 350 L = 1, 3
               IF (ISTRT(L) .LT. 1 .OR. IEND(L) .GT. IJKDIM(L)) THEN
                  IERRCD = IYES
                  WRITE (IOUT,1300) IBLK, DIR(L), IJKDIM(L), IDIR,
     1                     (ISTRT(L1),L1=1,3), (IEND(L2),L2=1,3),
     2                     (INC  (L3),L3=1,3)
               ENDIF
  350       CONTINUE
  380    CONTINUE
  390 CONTINUE
C
C Finished
C
      WRITE (IOUT,1900)
      RETURN
C
 1000 FORMAT (' ','CHKDAT: Begin    checking BC, CUT, PRINT data')
 1010 FORMAT (' ','CHKDAT: BC    ERROR-> Two-Dimensional option (TWOD)',
     1            ' requires KDIM = 2',
     2       /' ',13X,' Block =',I3,' Dimensions ',3I5)
 1110 FORMAT (' ','CHKDAT: BC    ERROR-> Block =',I3,' Index ',A1,
     1            ' outside allowable range [1:',I4,']',
     2       /' ',13X,' BC Data ',3I4,3X,3I4,3X,A20)
 1120 FORMAT (' ','CHKDAT: BC    ERROR-> Block =',I3,' Index ',A1,
     1            ' Wall normal viscous terms not enabled',
     2       /' ',13X,' BC Data ',3I4,3X,3I4,3X,A20)
 1130 FORMAT (' ','CHKDAT: BC    ERROR-> Block =',I3,' Index ',A1,
     1            ' 2D FARFIELD VORTEX BC must be used',
     2       /' ',22X,'with the TWOD option',
     3       /' ',13X,' BC Data ',3I4,3X,3I4,3X,A20)
 1140 FORMAT (' ','CHKDAT: BC    ERROR-> Block =',I3,' Index ',A1,
     1            ' 2D FARFIELD VORTEX BC cannot be used',
     2       /' ',22X,'with the MARCH option',
     3       /' ',13X,' BC Data ',3I4,3X,3I4,3X,A20)
 1200 FORMAT (' ','CHKDAT: CUT   ERROR-> Invalid Block =',I3,
     1            ' Number of blocks specified is ',I4,
     2       /' ',13X,' CUT Data ',A20,1X,3I4,3X,3I4)
 1210 FORMAT (' ','CHKDAT: CUT   ERROR-> Block =',I3,' Index ',A1,
     1            ' outside allowable range [1:',I4,']',
     2       /' ',13X,' CUT Data ',A20,1X,3I4,3X,3I4)
 1220 FORMAT (' ','CHKDAT: CUT   ERROR-> IORD = ',3I2,
     1            ' has duplicate directions specified',
     2       /' ',13X,' CUT Data ',A20,1X,I3,1X,3I4,3X,3I4,
     3       /' ',13X,'          ',A20,1X,I3,1X,3I4,3X,3I4)
 1230 FORMAT (' ','CHKDAT: CUT   ERROR-> Block',I3,' Index ',A1,
     1            ' different length from Block',I3,' Index ',A1,
     2       /' ',13X,' CUT Data ',A20,1X,3I4,3X,3I4,
     3       /' ',13X,'          ',A20,1X,3I4,3X,3I4,2X,3I2)
 1300 FORMAT (' ','CHKDAT: PRINT ERROR-> Block =',I3,' Index ',A1,
     1            ' outside allowable range [1:',I4,']',
     2       /' ',13X,' Print Data ',I2,3X,3I4,3X,3I4,3X,3I4)
 1900 FORMAT (' ','CHKDAT: Finished checking BC, CUT, PRINT data')
C
      END





      SUBROUTINE CHKCRS (NLVLS,  NXBLKS, NBLKS,  IDIM, JDIM, KDIM, 
     1                   NXBCS,  NBCS,   IBCDAT, 
     2                   NXCUTS, NCUTS,  ICUTS,  CUTNAM, IERRCD)
C
C Routine to check coarsened boundary conditions and cuts for consistency.
C An error is returned in IERRCD.
C
C This routine is organized as follows:
C     1. Verify Boundary Condition (BC) data
C     2. Verify Cut data
C
C NLVLS          : Number of grid levels
C NXBLKS         : Maximum number of blocks as dimensioned in main
C NBLKS          : Number of blocks
C IDIM,JDIM,KDIM : Dimensions of blocks
C NXBCS          : Maximum number of bcs    per block as dimensioned in main
C NBCS           : Number of boundary conditions in each block
C IBCDAT         : Boundary condition data
C NXCUTS         : Maximum number of cuts as dimensioned in main
C NCUTS          : Number of cuts
C ICUTS          : Cut data
C IERRCD         : Error code
C
      include '../header/common.h'
C
      DIMENSION IDIM   (NXBLKS,NLVLS), 
     1          JDIM   (NXBLKS,NLVLS),
     2          KDIM   (NXBLKS,NLVLS)
C
      DIMENSION IBCDAT (10,NXBCS,NXBLKS,NLVLS), 
     1          NBCS   (NXBLKS)
C
      DIMENSION ICUTS  (21,NXCUTS,NLVLS)
      CHARACTER*20 CUTNAM(NXCUTS)
C
C Local storage used to test indices
C
      DIMENSION IJKDIM(3),
     1          ISTRTC(3),
     2          IENDC (3),
     3          ISTRTF(3),
     4          IENDF (3),
     5          IORD  (3),
     6          ISTR2C(3),
     7          IEND2C(3),
     8          ISTR2F(3),
     9          IEND2F(3) 
C
C Output data
C
      CHARACTER*20 BCNAME
      CHARACTER*1  DIR (3)
      DATA         DIR /'I', 'J', 'K'/
C
C Initialize error code
C
      IERRCD = INO
C
      IF (THREED) THEN
         NDIM = 3
      ELSE
         NDIM = 2
      ENDIF
C
      WRITE (IOUT,1000)
C
C 1.   Verify Boundary Condition (BC) data
C
      DO 190 IBLK = 1, NBLKS
C
         DO 170 IBC = 1, NBCS(IBLK)
C
            DO 160 ILVL = 2, NLVLS
               LVLF = ILVL - 1
               LVLC = ILVL
C
C Store boundary condition data in mnemonic names
C
               ISTRTF(1) = IBCDAT(4,IBC,IBLK,LVLF)
               IENDF (1) = IBCDAT(5,IBC,IBLK,LVLF)
               ISTRTF(2) = IBCDAT(6,IBC,IBLK,LVLF)
               IENDF (2) = IBCDAT(7,IBC,IBLK,LVLF)
               ISTRTF(3) = IBCDAT(8,IBC,IBLK,LVLF)
               IENDF (3) = IBCDAT(9,IBC,IBLK,LVLF)
C
               IBCTYP    = IBCDAT(1,IBC,IBLK,LVLC)
               IDIR      = IBCDAT(2,IBC,IBLK,LVLC)
               INROUT    = IBCDAT(3,IBC,IBLK,LVLC)
               ISTRTC(1) = IBCDAT(4,IBC,IBLK,LVLC)
               IENDC (1) = IBCDAT(5,IBC,IBLK,LVLC)
               ISTRTC(2) = IBCDAT(6,IBC,IBLK,LVLC)
               IENDC (2) = IBCDAT(7,IBC,IBLK,LVLC)
               ISTRTC(3) = IBCDAT(8,IBC,IBLK,LVLC)
               IENDC (3) = IBCDAT(9,IBC,IBLK,LVLC)
C
C Modify starting point to switch back to input form for testing
C
               DO 140 L = 1, 3
                  IF (L .NE. IDIR) THEN
                     ISTRTF(L) = ISTRTF(L) - 1
                     ISTRTC(L) = ISTRTC(L) - 1
                  ENDIF
  140          CONTINUE
C
               DO 150 L = 1, NDIM
                  ISTEST = 2 * (ISTRTC(L) - 1) + 1
                  IETEST = 2 * (IENDC (L) - 1) + 1
                  IF (ISTRTF(L) .NE. ISTEST  .OR.
     1                IENDF (L) .NE. IETEST) THEN
                     IERRCD = IYES
                     CALL GETNAM (IBCTYP, BCNAME)
                     WRITE (IOUT,1100) IBLK, LVLC, LVLF, DIR(L),
     1                  (ISTRTC(L1),L1=1,3), (IENDC(L2),L2=1,3), BCNAME,
     2                  (ISTRTF(L3),L3=1,3), (IENDF(L4),L4=1,3), BCNAME
                  ENDIF
  150          CONTINUE
C
  160       CONTINUE
  170    CONTINUE
  190 CONTINUE
C
C 2.   Verify Cut data 
C
      DO 250 ICUT = 1, NCUTS
         DO 240 ILVL = 2, NLVLS
            LVLF = ILVL - 1
            LVLC = ILVL
C
            ISTRTF(1) = ICUTS ( 4,ICUT,LVLF)
            IENDF (1) = ICUTS ( 5,ICUT,LVLF)
            ISTRTF(2) = ICUTS ( 6,ICUT,LVLF)
            IENDF (2) = ICUTS ( 7,ICUT,LVLF)
            ISTRTF(3) = ICUTS ( 8,ICUT,LVLF)
            IENDF (3) = ICUTS ( 9,ICUT,LVLF)
C
            ISTR2F(1) = ICUTS (13,ICUT,LVLF)
            IEND2F(1) = ICUTS (14,ICUT,LVLF)
            ISTR2F(2) = ICUTS (15,ICUT,LVLF)
            IEND2F(2) = ICUTS (16,ICUT,LVLF)
            ISTR2F(3) = ICUTS (17,ICUT,LVLF)
            IEND2F(3) = ICUTS (18,ICUT,LVLF)
C
            IBLK      = ICUTS ( 1,ICUT,LVLC)
            IDIR      = ICUTS ( 2,ICUT,LVLC)
            INROUT    = ICUTS ( 3,ICUT,LVLC)
            ISTRTC(1) = ICUTS ( 4,ICUT,LVLC)
            IENDC (1) = ICUTS ( 5,ICUT,LVLC)
            ISTRTC(2) = ICUTS ( 6,ICUT,LVLC)
            IENDC (2) = ICUTS ( 7,ICUT,LVLC)
            ISTRTC(3) = ICUTS ( 8,ICUT,LVLC)
            IENDC (3) = ICUTS ( 9,ICUT,LVLC)
C
            IBLK2     = ICUTS (10,ICUT,LVLC)
            IDIR2     = ICUTS (11,ICUT,LVLC)
            INROUT2   = ICUTS (12,ICUT,LVLC)
            ISTR2C(1) = ICUTS (13,ICUT,LVLC)
            IEND2C(1) = ICUTS (14,ICUT,LVLC)
            ISTR2C(2) = ICUTS (15,ICUT,LVLC)
            IEND2C(2) = ICUTS (16,ICUT,LVLC)
            ISTR2C(3) = ICUTS (17,ICUT,LVLC)
            IEND2C(3) = ICUTS (18,ICUT,LVLC)
C
            IORD  (1) = ICUTS (19,ICUT,LVLC)
            IORD  (2) = ICUTS (20,ICUT,LVLC)
            IORD  (3) = ICUTS (21,ICUT,LVLC)
C
C Modify starting/ending point to switch back to input form for testing
C
            DO 220 L = 1, 3
               IF (L .NE. IDIR) THEN
                  IF (IENDC(L) .GT. ISTRTC(L)) THEN
                     ISTRTF(L) = ISTRTF(L) - 1
                     ISTRTC(L) = ISTRTC(L) - 1
                  ELSE
                     IENDF (L) = IENDF (L) - 1
                     IENDC (L) = IENDC (L) - 1
                  ENDIF
               ENDIF
               IF (L .NE. IDIR2) THEN
                  IF (IEND2C(L) .GT. ISTR2C(L)) THEN
                     ISTR2F(L) = ISTR2F(L) - 1
                     ISTR2C(L) = ISTR2C(L) - 1
                  ELSE
                     IEND2F(L) = IEND2F(L) - 1
                     IEND2C(L) = IEND2C(L) - 1
                  ENDIF
               ENDIF
  220       CONTINUE
C
C Test indices
C
            DO 230 L = 1, NDIM
               ISTEST = 2 * (ISTRTC(L) - 1) + 1
               IETEST = 2 * (IENDC (L) - 1) + 1
               IS2TST = 2 * (ISTR2C(L) - 1) + 1
               IE2TST = 2 * (IEND2C(L) - 1) + 1
               IF (ISTRTF(L) .NE. ISTEST .OR.
     1             IENDF (L) .NE. IETEST .OR.
     2             ISTR2F(L) .NE. IS2TST .OR.
     3             IEND2F(L) .NE. IE2TST ) THEN
                  IERRCD = IYES
                  WRITE (IOUT,1200) ICUT, LVLC, LVLF, DIR(L), 
     1                  CUTNAM(ICUT), IBLK,
     2                  (ISTRTC(L1),L1=1,3), (IENDC (L2),L2=1,3),
     3                  CUTNAM(ICUT), IBLK2,
     4                  (ISTR2C(L3),L3=1,3), (IEND2C(L4),L4=1,3),
     5                  (IORD  (L5),L5=1,3),
     6                  CUTNAM(ICUT), IBLK,
     7                  (ISTRTF(L1),L1=1,3), (IENDF (L2),L2=1,3),
     8                  CUTNAM(ICUT), IBLK2,
     9                  (ISTR2F(L3),L3=1,3), (IEND2F(L4),L4=1,3),
     A                  (IORD  (L5),L5=1,3)
               ENDIF
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE
C
C Finished
C
      WRITE (IOUT,1900)
      RETURN
C
 1000 FORMAT (' ','CHKCRS: Begin    checking coarsened BC and CUT data')
 1100 FORMAT (' ','CHKCRS: BC    ERROR-> Block =',I3,' Coarse grid ',
     1            'boundary points (Level = ',I2,')',
     2       /' ',22X,'do not correspond to fine grid (Level = ',I2,')',
     3       /' ',22X,'in the ',A1,' direction:',
     4       /' ',13X,' Coarse BC  Data ',3I4,3X,3I4,3X,A20,
     5       /' ',13X,' Fine   BC  Data ',3I4,3X,3I4,3X,A20)
 1200 FORMAT (' ','CHKCRS: CUT   ERROR-> Cut   =',I3,' Coarse grid ',
     1            'cut points (Level = ',I2,')',
     2       /' ',22X,'do not correspond to fine grid (Level = ',I2,')',
     3       /' ',22X,'in the ',A1,' direction:',
     4       /' ',13X,' Coarse CUT Data ',A20,I2,3X,3I4,3X,3I4,
     5       /' ',13X,'                 ',A20,I2,3X,3I4,3X,3I4,3X,3I2,
     6       /' ',13X,' Fine   CUT Data ',A20,I2,3X,3I4,3X,3I4,
     7       /' ',13X,'                 ',A20,I2,3X,3I4,3X,3I4,3X,3I2)
 1900 FORMAT (' ','CHKCRS: Finished checking coarsened BC and CUT data')
C
      END




      SUBROUTINE CONC2P (ITURB, NPTS, Q, TKE)
C
C Routine to convert Conserved variables to Primitive variables.
C
C This routine is organized as follows:
C     1. Divide conserved variables by density
C     2. Calculate turbulent kinetic energy
C     3. Convert total energy to pressure 
C
C ITURB          : Turbulence model
C NPTS           : Number of points to convert
C Q              : Flow variables (Input - conserved; Output - primitive)
C TKE            : Turbulent kinetic energy
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          TKE   (NPTS)
C
      INTEGER   I, L
C
C 1.  Multiply primitive variables by density
C
      DO 110 L = 2, NQ
         IF (L .NE. 5) THEN
            DO 100 I = 1, NPTS
               Q(I,L) = Q(I,L) / Q(I,1)
  100       CONTINUE
         ENDIF
  110 CONTINUE
C
C 2.  Calculate turbulent kinetic energy
C
      CALL GETTKE (ITURB, NPTS, Q, TKE)
C
C 3.  Convert from (rho*Etotal) to pressure 
C
      DO 300 I = 1, NPTS
         RHO    = Q(I,1)
         U      = Q(I,2)
         V      = Q(I,3)
         W      = Q(I,4)
         RHOE   = Q(I,5)
         Q(I,5) = GAMM1 * ( RHOE 
     1          - 0.5E0 * RHO * (U*U + V*V + W*W) 
     2          - RHO * TKE(I) )
  300 CONTINUE
C
C     Finished converting conserved variables to primitive variables
C
      RETURN
      END




      SUBROUTINE CONP2C (ITURB, NPTS, Q, TKE)
C
C Routine to convert Primitive variables to Conserved variables.
C
C This routine is organized as follows:
C     1. Calculate turbulent kinetic energy
C     2. Multiply primitive variables by density
C     3. Convert pressure to total energy
C
C ITURB          : Turbulence model
C NPTS           : Number of points to convert
C Q              : Flow variables (Input - conserved; Output - primitive)
C TKE            : Turbulent kinetic energy
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          TKE   (NPTS)
C
      INTEGER   I, L
C
C 1.  Calculate turbulent kinetic energy
C
      CALL GETTKE (ITURB, NPTS, Q, TKE)
C
C 2.  Multiply primitive variables by density
C
      DO 210 L = 2, NQ
         IF (L .NE. 5) THEN
            DO 200 I = 1, NPTS
               Q(I,L) = Q(I,L) * Q(I,1)
  200       CONTINUE
         ENDIF
  210 CONTINUE
C
C 3.  Convert from pressure to (rho*Etotal)
C
      DO 300 I = 1, NPTS
         RHO    = Q(I,1)
         RHOU   = Q(I,2)
         RHOV   = Q(I,3)
         RHOW   = Q(I,4)
         P      = Q(I,5)
         Q(I,5) = P / GAMM1
     1          + 0.5E0 * (RHOU*RHOU + RHOV*RHOV + RHOW*RHOW) / RHO
     2          + RHO * TKE(I)
  300 CONTINUE
C
C     Finished converting primitive variables to conserved variables
C
      RETURN
      END




      SUBROUTINE CPVCTR (NPTS, VIN, VOUT)
C
C Routine to copy a vector from one array to another.
C
C NPTS  : Number of points in array
C VIN   : Input vector to copy
C VOUT  : Output for copy
C
      include '../header/common.h'
C
      INTEGER   NPTS
C
      DIMENSION VIN   (NPTS),
     1          VOUT  (NPTS)
C
C Local variables
C
      INTEGER   I
C
C Copy
C
      DO 100 I = 1, NPTS
         VOUT(I) = VIN(I)
  100 CONTINUE
C
      RETURN
C
C Finished
C
      END






      SUBROUTINE GETNAM (IBCTYP, BCNAME)
C
C Function to give boundary condition name given the type
C
      include '../header/common.h'
C
      CHARACTER*20 BCNAME
C
      IF (IBCTYP .EQ. IWALL  ) THEN
         BCNAME = 'WALL'
      ELSE IF (IBCTYP .EQ. IWALFN) THEN
         BCNAME = 'WALL FUNCTIONS'
      ELSE IF (IBCTYP .EQ. ITAN  ) THEN
         BCNAME = 'TANGENCY'
      ELSE IF (IBCTYP .EQ. IFAR  ) THEN
         BCNAME = 'FARFIELD'
      ELSE IF (IBCTYP .EQ. IFAR2D) THEN
         BCNAME = '2D FARFIELD VORTEX'
      ELSE IF (IBCTYP .EQ. ISYMXY) THEN
         BCNAME = 'SYMMETRY XY'
      ELSE IF (IBCTYP .EQ. ISYMXZ) THEN
         BCNAME = 'SYMMETRY XZ'
      ELSE IF (IBCTYP .EQ. ISYMYZ) THEN
         BCNAME = 'SYMMETRY YZ'
      ELSE IF (IBCTYP .EQ. IEXT  ) THEN
         BCNAME = 'EXTRAPOLATE'
      ELSE IF (IBCTYP .EQ. IFIX  ) THEN
         BCNAME = 'FIX'
      ELSE IF (IBCTYP .EQ. IJET  ) THEN
         BCNAME = 'JET'
      ELSE IF (IBCTYP .EQ. IPERD ) THEN
         BCNAME = 'PERIODIC'
      ELSE IF (IBCTYP .EQ. IHSHR ) THEN
         BCNAME = 'HOMOGENEOUS SHEAR'
      ELSE IF (IBCTYP .EQ. IPROFL) THEN
         BCNAME = 'PROFILE CELL CENTER'
      ELSE IF (IBCTYP .EQ. IPROSB) THEN
         BCNAME = 'SUBSONIC PROFILE CC'
      ELSE
         BCNAME = '** UNKNOWN **'
      ENDIF
C
      RETURN
      END






      SUBROUTINE CHCKBC (ILVL, IBLK, IDIM, JDIM, KDIM, Q, NERRBC)
C
C     Routine to asure that all faces of each computational zone are
C     set with a boundary condition.
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ) 
C
      LOGICAL   ERRFAC
C
C     Set constant for maximum errors to print out
C
      MAXERR = 25
      ERRFAC = .FALSE.
C
C     I = 0 Face Boundary Conditions
C
      NERFAC = 0
      DO 140 I = 0, 1
         DO 130 K = 2, KDIM
            DO 120 J = 2, JDIM
               DO 110 L = 1, NQ
                  IF (Q(J,K,I,L) .GE. RMAX) THEN
                     ERRFAC = .TRUE.
                  ENDIF
  110          CONTINUE
               IF (ERRFAC) THEN
                  IF (NERFAC .LE. MAXERR) THEN
                     WRITE (IOUT,1000) ILVL, IBLK, I, J, K
                  ENDIF
                  ERRFAC = .FALSE.
                  NERRBC = NERRBC + 1
                  NERFAC = NERFAC + 1
               ENDIF
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
C
C     Output total number of errors on current face
C
      IF (NERFAC .GT. 0) THEN
         WRITE (IOUT,1010) NERFAC, ILVL, IBLK, 'I = 0 Face    '
      ENDIF
C
C     I = IDIM Face Boundary Conditions
C
      NERFAC = 0
      DO 180 I = IDIM+1, IDIM+2
         DO 170 K = 2, KDIM
            DO 160 J = 2, JDIM
               DO 150 L = 1, NQ
                  IF (Q(J,K,I,L) .GE. RMAX) THEN
                     ERRFAC = .TRUE.
                  ENDIF
  150          CONTINUE
               IF (ERRFAC) THEN
                  IF (NERFAC .LE. MAXERR) THEN
                     WRITE (IOUT,1000) ILVL, IBLK, I, J, K
                  ENDIF
                  ERRFAC = .FALSE.
                  NERRBC = NERRBC + 1
                  NERFAC = NERFAC + 1
               ENDIF
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
C
C     Output total number of errors on current face
C
      IF (NERFAC .GT. 0) THEN
         WRITE (IOUT,1010) NERFAC, ILVL, IBLK, 'I = IDIM Face '
      ENDIF
C
C     J = 0 Face Boundary Conditions
C
      NERFAC = 0
      DO 240 J = 0, 1
         DO 230 I = 2, IDIM
            DO 220 K = 2, KDIM
               DO 210 L = 1, NQ
                  IF (Q(J,K,I,L) .GE. RMAX) THEN
                     ERRFAC = .TRUE.
                  ENDIF
  210          CONTINUE
               IF (ERRFAC) THEN
                  IF (NERFAC .LE. MAXERR) THEN
                     WRITE (IOUT,1000) ILVL, IBLK, I, J, K
                  ENDIF
                  ERRFAC = .FALSE.
                  NERRBC = NERRBC + 1
                  NERFAC = NERFAC + 1
               ENDIF
  220       CONTINUE
  230    CONTINUE
  240 CONTINUE
C
C     Output total number of errors on current face
C
      IF (NERFAC .GT. 0) THEN
         WRITE (IOUT,1010) NERFAC, ILVL, IBLK, 'J = 0 Face    '
      ENDIF
C
C     J = JDIM Face Boundary Conditions
C
      NERFAC = 0
      DO 280 J = JDIM+1, JDIM+2
         DO 270 I = 2, IDIM
            DO 260 K = 2, KDIM
               DO 250 L = 1, NQ
                  IF (Q(J,K,I,L) .GE. RMAX) THEN
                     ERRFAC = .TRUE.
                  ENDIF
  250          CONTINUE
               IF (ERRFAC) THEN
                  IF (NERFAC .LE. MAXERR) THEN
                     WRITE (IOUT,1000) ILVL, IBLK, I, J, K
                  ENDIF
                  ERRFAC = .FALSE.
                  NERRBC = NERRBC + 1
                  NERFAC = NERFAC + 1
               ENDIF
  260       CONTINUE
  270    CONTINUE
  280 CONTINUE
C
C     Output total number of errors on current face
C
      IF (NERFAC .GT. 0) THEN
         WRITE (IOUT,1010) NERFAC, ILVL, IBLK, 'J = JDIM Face '
      ENDIF
C
C     Check K = 0 and K = KDIM boundary conditions for 3D or axisymmetric
C
      IF (THREED .OR. AXISYM) THEN
C
C     K = 0 Face Boundary Conditions
C
         NERFAC = 0
         DO 340 K = 0, 1
            DO 330 I = 2, IDIM
               DO 320 J = 2, JDIM
                  DO 310 L = 1, NQ
                     IF (Q(J,K,I,L) .GE. RMAX) THEN
                        ERRFAC = .TRUE.
                     ENDIF
  310             CONTINUE
                  IF (ERRFAC) THEN
                     IF (NERFAC .LE. MAXERR) THEN
                        WRITE (IOUT,1000) ILVL, IBLK, I, J, K
                     ENDIF
                     ERRFAC = .FALSE.
                     NERRBC = NERRBC + 1
                     NERFAC = NERFAC + 1
                  ENDIF
  320          CONTINUE
  330       CONTINUE
  340    CONTINUE
C
C     Output total number of errors on current face
C
         IF (NERFAC .GT. 0) THEN
            WRITE (IOUT,1010) NERFAC, ILVL, IBLK, 'K = 0 Face    '
         ENDIF
C
C     K = KDIM Face Boundary Conditions
C
         NERFAC = 0
         DO 380 K = KDIM+1, KDIM+2
            DO 370 I = 2, IDIM
               DO 360 J = 2, JDIM
                  DO 350 L = 1, NQ
                     IF (Q(J,K,I,L) .GE. RMAX) THEN
                        ERRFAC = .TRUE.
                     ENDIF
  350             CONTINUE
                  IF (ERRFAC) THEN
                     IF (NERFAC .LE. MAXERR) THEN
                        WRITE (IOUT,1000) ILVL, IBLK, I, J, K
                     ENDIF
                     ERRFAC = .FALSE.
                     NERRBC = NERRBC + 1
                     NERFAC = NERFAC + 1
                  ENDIF
  360          CONTINUE
  370       CONTINUE
  380    CONTINUE
C
C     Output total number of errors on current face
C
         IF (NERFAC .GT. 0) THEN
            WRITE (IOUT,1010) NERFAC, ILVL, IBLK, 'K = KDIM Face '
         ENDIF
      ENDIF
C
      RETURN
 1000 FORMAT (' ','CHCKBC: ERROR-> Boundary data not set at ',
     1            '(lvl,block,i,j,k)',I3,':',I3,':',3I4)
 1010 FORMAT (' ','CHCKBC: ERROR-> ',I8,' BC Errors for level =',I3,
     1            ', block =',I3,', ',A//)
      END





      SUBROUTINE CHKCUT (IBLK1, IDIM1, JDIM1, KDIM1, R1,
     1                   IBLK2, IDIM2, JDIM2, KDIM2, R2, ICUTS, 
     2                   CUTNAM, IERRCD)
C
C Routine to check cut boundary conditions.  Grid points in R1 are verified 
C against grid points in R2 for each cut boundary condition.
C
C IBLK1             : Block 1
C IDIM1,JDIM1,KDIM1 : Dimensions of current block 1
C R1                : Position vector of grid vertices (x,y,z) for block 1
C IBLK2             : Block 2
C IDIM2,JDIM2,KDIM2 : Dimensions of current block 2
C R2                : Position vector of grid vertices (x,y,z) for block 2
C ICUTS             : Data controlling transfer from block 2 to block 1
C CUTNAM            : Name of cut
C IERRCD            : Error code
C
      include '../header/common.h'
      DIMENSION R1   (0:JDIM1+1,0:KDIM1+1,0:IDIM1+1,3),
     1          R2   (0:JDIM2+1,0:KDIM2+1,0:IDIM2+1,3)
      DIMENSION ICUTS(21)
      CHARACTER*20 CUTNAM
C
      DIMENSION IJKDM1(3), IJKDM2(3)
      DIMENSION ISTR1 (3), IEND1 (3), 
     1          ISTR2 (3), IEND2 (3), 
     2          ISGN1 (3), ISGN2 (3),
     3          IORD  (3), LEN   (3)
C
C Control
C
      MAXPRT = 25
      GTOLER = 2.5E-05
C
C Store IDIM, JDIM, KDIM in IJKDIM for looping
C
      IJKDM1(1) = IDIM1
      IJKDM1(2) = JDIM1
      IJKDM1(3) = KDIM1
      IJKDM2(1) = IDIM2
      IJKDM2(2) = JDIM2
      IJKDM2(3) = KDIM2
C
C Store boundary condition data in mnemonic names
C
      IDIR1    = ICUTS( 2)
      INOUT1   = ICUTS( 3)
      ISTR1(1) = ICUTS( 4)
      IEND1(1) = ICUTS( 5)
      ISTR1(2) = ICUTS( 6)
      IEND1(2) = ICUTS( 7)
      ISTR1(3) = ICUTS( 8)
      IEND1(3) = ICUTS( 9)
C
      IBLK2    = ICUTS(10)
      IDIR2    = ICUTS(11)
      INOUT2   = ICUTS(12)
      ISTR2(1) = ICUTS(13)
      IEND2(1) = ICUTS(14)
      ISTR2(2) = ICUTS(15)
      IEND2(2) = ICUTS(16)
      ISTR2(3) = ICUTS(17)
      IEND2(3) = ICUTS(18)
C
      IORD(1)  = ICUTS(19)
      IORD(2)  = ICUTS(20)
      IORD(3)  = ICUTS(21)
C
C Translate ICUTS data from cell-centered storage to grid point storage.
C This requires decrementing the minimum ('starting') index. Do this only 
C for the indices that the cut loops over, not the idir direction.
C
      DO 10 L = 1, 3
C
         IF (L .NE. IDIR1) THEN
            IF (IEND1(L) .GE. ISTR1(L)) THEN
               ISTR1(L) = ISTR1(L) - 1
            ELSE
               IEND1(L) = IEND1(L) - 1
            ENDIF
         ENDIF
C
         IF (L .NE. IDIR2) THEN
            IF (IEND2(L) .GE. ISTR2(L)) THEN
               ISTR2(L) = ISTR2(L) - 1
            ELSE
               IEND2(L) = IEND2(L) - 1
            ENDIF
         ENDIF
   10 CONTINUE
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INOUT1 = 1) or the outer boundary (INOUT1 > 1)
C      IBCPT1 = grid point on grid face of block 1
C      IBCPT2 = corresponding grid point on grid face of block 2
C
      IF (INOUT1 .EQ. 1) THEN
         IBCPT1 = 1
      ELSE
         IBCPT1 = IJKDM1(IDIR1)
      ENDIF
C
      IF (INOUT2 .EQ. 1) THEN
         IBCPT2 = 1
      ELSE
         IBCPT2 = IJKDM2(IDIR2)
      ENDIF
C
C Find the length of the two outer loops and loop over these using
C offsets and delta function to set the two cut data points to the
C two interior data points of block 2.
C
      DO 100 L = 1, 3
C
         LEN(L) = ABS ( IEND1(L) - ISTR1(L) )
C
C     Increment/Decrement 
C
         IF ( IEND1(L) .GT. ISTR1(L) ) THEN
            ISGN1(L) =   1
         ELSE
            ISGN1(L) = - 1
         ENDIF
C
C     Increment/Decrement 
C
         IF ( IEND2(L) .GT. ISTR2(L) ) THEN
            ISGN2(L) =   1
         ELSE
            ISGN2(L) = - 1
         ENDIF
  100 CONTINUE
C
C     II1,JJ1,KK1 grid point on boundary of block 1
C     II2,JJ2,KK2 corresponding grid point on boundary of block 2
C
      IC1 = ICY (IDIR1, 1)
      IC2 = ICY (IDIR1, 2)
      IC3 = ICY (IDIR1, 3)
C
      JC1 = IORD (IC1)
      JC2 = IORD (IC2)
      JC3 = IORD (IC3)
C
C Initialize counter
C
      ICOUNT = 0
C
      DO 120 I3 = 0, LEN(IC3)
         DO 110 I2 = 0, LEN(IC2)
            II1 = IBCPT1                       * IKD (IC1, 1) +
     1            (ISTR1(IC2) + ISGN1(IC2)*I2) * IKD (IC2, 1) +
     2            (ISTR1(IC3) + ISGN1(IC3)*I3) * IKD (IC3, 1) 
            JJ1 = IBCPT1                       * IKD (IC1, 2) +
     1            (ISTR1(IC2) + ISGN1(IC2)*I2) * IKD (IC2, 2) +
     2            (ISTR1(IC3) + ISGN1(IC3)*I3) * IKD (IC3, 2) 
            KK1 = IBCPT1                       * IKD (IC1, 3) +
     1            (ISTR1(IC2) + ISGN1(IC2)*I2) * IKD (IC2, 3) +
     2            (ISTR1(IC3) + ISGN1(IC3)*I3) * IKD (IC3, 3) 
C
            II2 = IBCPT2                       * IKD (JC1, 1) +
     1            (ISTR2(JC2) + ISGN2(JC2)*I2) * IKD (JC2, 1) +
     2            (ISTR2(JC3) + ISGN2(JC3)*I3) * IKD (JC3, 1) 
            JJ2 = IBCPT2                       * IKD (JC1, 2) +
     1            (ISTR2(JC2) + ISGN2(JC2)*I2) * IKD (JC2, 2) +
     2            (ISTR2(JC3) + ISGN2(JC3)*I3) * IKD (JC3, 2) 
            KK2 = IBCPT2                       * IKD (JC1, 3) +
     1            (ISTR2(JC2) + ISGN2(JC2)*I2) * IKD (JC2, 3) +
     2            (ISTR2(JC3) + ISGN2(JC3)*I3) * IKD (JC3, 3) 
C
            DX = R1(JJ1,KK1,II1,1) - R2(JJ2,KK2,II2,1)
            DY = R1(JJ1,KK1,II1,2) - R2(JJ2,KK2,II2,2)
            DZ = R1(JJ1,KK1,II1,3) - R2(JJ2,KK2,II2,3)
C
            DIST = SQRT (DX*DX + DY*DY + DZ*DZ)
C
C If the grid points don't match output an error.  Output the header on
C the first point that doesn't match and then output only the data up to
C a maximum number before a summary.
C
            IF (DIST .GT. GTOLER) THEN
               ICOUNT = ICOUNT + 1
               IF (ICOUNT .EQ. 1) THEN
                  WRITE (IOUT,1000) CUTNAM, IBLK1, IBLK2
               ENDIF
               IF (ICOUNT .LE. MAXPRT) THEN
                  WRITE (IOUT,1010) 
     1              II1, JJ1, KK1, R1(JJ1,KK1,II1,1), R1(JJ1,KK1,II1,2), 
     2                             R1(JJ1,KK1,II1,3),
     3              II2, JJ2, KK2, R2(JJ2,KK2,II2,1), R2(JJ2,KK2,II2,2), 
     4                             R2(JJ2,KK2,II2,3)
               ENDIF
            ENDIF
 110     CONTINUE
 120  CONTINUE
C
C Output summary
C
      IF (ICOUNT .GT. 0) THEN
         IERRCD = IYES
         WRITE (IOUT,1020) CUTNAM, IBLK1, IBLK2, ICOUNT
      ENDIF
C
C Finished Checking Cut 
C
      RETURN
 1000 FORMAT (' ','CHKCUT: ERROR-> Grid Points do not match on cut ',
     1            A20,
     2       /' ',20X,'Block = ',I3,28X,'Block = ',I3,
     3       /' ','  I   J   K     X        Y        Z    ',
     4            '  I   J   K     X        Y        Z    ')
 1010 FORMAT (' ',2(3I4,3(1PE9.2)))
 1020 FORMAT (' ','CHKCUT: ERROR-> Cut ',A20,' from Block ',I3,
     1            ' to Block ',I3,' had ',I10,' mismatched points.'/)
      END





      SUBROUTINE CONVRT (IDIM, JDIM, KDIM, 
     1                   IMODLX, QFROM, PROPSF, NQFROM, NPFROM,
     2                   IMODEL, QTO,   PROPST, NQTO,   NPTO,   S, VOL, 
     3                   QC, PROPSC, DQDX, DQDY, DQDZ, TAU, SRC, 
     4                   IBCDAT, NBCS, IFDELQ)
C
C Routine to convert from one equation set to another.  Perform the conversion 
C for the interior of the domain.  The call to the boundary condition routine 
C will then set the boundaries accordingly.
C
C This routine converts the following data:
C      Part:   FROM:                    TO:
C        0.    Mean Flow Equations
C        1.a.  Laminar               -> K-Epsilon
C        1.b.  Laminar               -> K-Omega
C        1.c.  Laminar               -> Reynolds Stress
C        3.    K-Epsilon             -> K-Omega
C        4.    K-Omega               -> K-Epsilon
C        5.    K-Epsilon or K-Omega  -> Reynolds Stress
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IMODLX         : Turbulence model to be converted FROM
C QFROM          : Primitive variables at cell centers converted FROM
C PROPSF         : Properties at cell centers converted FROM
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity 
C                  PROPS(3) = Y+
C                  PROPS(4) = Ynorm
C                  PROPS(5) = C_mu
C NQFROM         : Number of variables  for the model converted FROM
C NPFROM         : Number of properties for the model converted FROM
C IMODEL         : Turbulence model to be converted TO
C QTO            : Primitive variables at cell centers converted TO
C PROPST         : Properties at cell centers converted TO
C NQTO           : Number of variables  for the model converted TO
C NPTO           : Number of properties for the model converted TO
C VOL            : Cell volume
C S              : Metrics
C VOL            : Cell volumes
C RES            : Residual
C QC             : Storage for Q in one-dimensional array
C PROPSC         : Storage for PROPS array in one-dimensional array
C DQDX,DQDY,DQDZ : Storage for derivatives of Q at cell centers
C TAU            : Storage for Reynolds stress terms
C SRC            : Source terms for cell (J,K)
C IBCDAT         : Data controlling boundary conditions
C NBCS           : Number of boundary conditions for current block
C IFDELQ         : Array controlling calculation of delq
C
      include '../header/common.h'
C
      DIMENSION IMODLX(NMDL),
     1          IMODEL(NMDL)
C
      DIMENSION QFROM (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQFROM), 
     1          PROPSF(0:JDIM+2,0:KDIM+2,0:IDIM+2,NPFROM),
     2          QTO   (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQTO),
     3          PROPST(0:JDIM+2,0:KDIM+2,0:IDIM+2,NPTO)
C
      DIMENSION S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     1          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION QC    (2:JDIM,2:KDIM,NQ), 
     1          PROPSC(2:JDIM,2:KDIM,NP),
     2          DQDX  (2:JDIM,2:KDIM,NQ),
     3          DQDY  (2:JDIM,2:KDIM,NQ), 
     4          DQDZ  (2:JDIM,2:KDIM,NQ),
     5          TAU   (2:JDIM,2:KDIM,6),
     6          SRC   (2:JDIM,2:KDIM,NF)
C
      DIMENSION IBCDAT(10,NBCS),
     1          IFDELQ(NQ)
C
      ITURBX = IMODLX( 1)
      ITURB  = IMODEL( 1)
C
C     Number of points in a cross-section to calculate Reynolds stresses at
C
      NPTS = (JDIM - 1) * (KDIM - 1)
C
C 0. Convert the Mean Flow Variables and Properties
C    The mean flow variables are identical for TO and FROM for all conversions
C    Properties 1-5 are the same for all models.
C
      DO 30 I = 0, IDIM+2
         DO 20 K = 0, KDIM+2
            DO 10 J = 0, JDIM+2
               QTO   (J,K,I,1) = QFROM (J,K,I,1)
               QTO   (J,K,I,2) = QFROM (J,K,I,2)
               QTO   (J,K,I,3) = QFROM (J,K,I,3)
               QTO   (J,K,I,4) = QFROM (J,K,I,4)
               QTO   (J,K,I,5) = QFROM (J,K,I,5)
               PROPST(J,K,I,1) = PROPSF(J,K,I,1)
               PROPST(J,K,I,2) = PROPSF(J,K,I,2)
               PROPST(J,K,I,3) = PROPSF(J,K,I,3)
               PROPST(J,K,I,4) = PROPSF(J,K,I,4)
               PROPST(J,K,I,5) = PROPSF(J,K,I,5)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
C
C 1. Initialize additional turbulence variables to freestream values.
C
      IF (ITURBX .EQ. ITLMNR) THEN
C
C 1.a. K-Epsilon
C
         IF (ITURB .EQ. ITKE) THEN
            DO 130 I = 2, IDIM
               DO 120 K = 2, KDIM
                  DO 110 J = 2, JDIM
                     QTO(J,K,I,6) = TKEINF
                     QTO(J,K,I,7) = EPSINF
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
C
C 1.b. K-Omega
C
         ELSE IF (ITURB .EQ. ITKW) THEN
            DO 160 I = 2, IDIM
               DO 150 K = 2, KDIM
                  DO 140 J = 2, JDIM
                     QTO(J,K,I,6) = TKEINF
                     QTO(J,K,I,7) = OMEINF
  140             CONTINUE
  150          CONTINUE
  160       CONTINUE
C
C 1.c. Reynolds Stress
C
         ELSE IF (ITURB .EQ. ITRS) THEN
            TXXINF = 2.E0 / 3.E0 * TKEINF
            TYYINF = 2.E0 / 3.E0 * TKEINF
            TZZINF = 2.E0 / 3.E0 * TKEINF
            DO 190 I = 2, IDIM
               DO 180 K = 2, KDIM
                  DO 170 J = 2, JDIM
                     QTO(J,K,I, 6) = TXXINF
                     QTO(J,K,I, 7) = TYYINF
                     QTO(J,K,I, 8) = TZZINF
                     QTO(J,K,I, 9) = 0.E0
                     QTO(J,K,I,10) = 0.E0
                     QTO(J,K,I,11) = 0.E0
                     QTO(J,K,I,12) = EPSINF
  170             CONTINUE
  180          CONTINUE
  190       CONTINUE
         ENDIF
C
C Do conversion for additional turbulence variables.
C
C 3. Convert K-Epsilon to K-Omega
C
      ELSE IF (ITURBX .EQ. ITKE .AND. ITURB .EQ. ITKW) THEN
         DO 330 I = 2, IDIM
            DO 320 K = 2, KDIM
               DO 310 J = 2, JDIM
                  QTO(J,K,I,6) = QFROM(J,K,I,6)
                  QTO(J,K,I,7) = QFROM(J,K,I,7)/(CMU*QFROM(J,K,I,6))
  310          CONTINUE
  320       CONTINUE
  330    CONTINUE
C
C 4. Convert K-Omega to K-Epsilon
C
      ELSE IF (ITURBX .EQ. ITKW .AND. ITURB .EQ. ITKE) THEN
         DO 430 I = 2, IDIM
            DO 420 K = 2, KDIM
               DO 410 J = 2, JDIM
                  QTO(J,K,I,6) = QFROM(J,K,I,6)
                  QTO(J,K,I,7) = CMU*QFROM(J,K,I,6)*QFROM(J,K,I,7)
  410          CONTINUE
  420       CONTINUE
  430    CONTINUE
C
C 5.  Convert K-Epsilon or K-Omega to Reynolds Stress
C     For Each Cross Section I
C        a. Calculate derivatives at cell centers
C        b. Calculate (negative of) Reynolds Stresses from 2-equation model
C        c. Store Reynolds stresses into Q
C        d. Calculate epsilon from epsilon or omega
C     Next I
C
      ELSE IF ( (ITURBX.EQ.ITKE .OR. ITURBX.EQ.ITKW) 
     1           .AND. ITURB.EQ.ITRS) THEN
C
C     Loop over I performing the conversion at each I station
C
         DO 590 I = 2, IDIM
C
C 5-a. Calculate derivatives of Q at cell centers.
C
C     Temporary variables used internally in DELQC:
C        TAU(1-3) for QJ  
C        SRC(1-3) for QK  
C        TAU(5->) for QPI 
C        SRC(4-6) for QMI 
C     QJ is dimensioned as JDIM,2:KDIM therefore need more room for it 
C     than 1 plane of TAU. Similarly for QK.  QPI and QMI require only 
C     one plane as they are dimensioned 2:JDIM,2:KDIM.
C     Use SRC for some of the temporary space to fix a memory allocation
C     problem.
C
C
C     Calculate derivatives at interior cell centers only.
C
            NBOUND = 0
            CALL DELQC (IDIM, JDIM, KDIM, I, IFDELQ, NBOUND,
     1                  QFROM, PROPSF, S, VOL, IBCDAT, NBCS,
     2                  TAU(2,2,1), SRC(2,2,1), TAU(2,2,5), SRC(2,2,4),
     3                  QC, PROPSC, DQDX, DQDY, DQDZ)
C
C 5-b.  Calculate the Reynold's stresses at the cell centers.  
C                                      ~
C     The array TAU contains {- rho (u" u")}.
C                                     i  j
C     For laminar flow:            TAU_ij => zero.
C         algebraic models:        TAU_ij => modeled with Boussinesq approx.
C         two-equation models:     TAU_ij => modeled with Boussinesq approx.
C         Reynold's stress models: TAU_ij => computed from PDE.
C     Pass in RMUT and DUDX,DUDY,DUDZ for Boussinesq approximation.
C     Pass in QC for Q data at same location as DQDX,... derivative data.
C     The Reynold's stresses are stored as:
C         TAU(1) = TAUXX, TAU(2) = TAUYY, TAU(3) = TAUZZ,
C         TAU(4) = TAUXY, TAU(5) = TAUXZ, TAU(6) = TAUYZ
C
            CALL TAUIJ (IMODLX, NPTS, QC, PROPSC, DQDX, DQDY, DQDZ, TAU)
C
C 5-c.  Calculate the TO variables
C
            DO 520 K = 2, KDIM
               DO 510 J = 2, JDIM
                  QTO(J,K,I, 6) = - TAU(J,K,1) / QTO(J,K,I,1)
                  QTO(J,K,I, 7) = - TAU(J,K,2) / QTO(J,K,I,1)
                  QTO(J,K,I, 8) = - TAU(J,K,3) / QTO(J,K,I,1)
                  QTO(J,K,I, 9) = - TAU(J,K,4) / QTO(J,K,I,1)
                  QTO(J,K,I,10) = - TAU(J,K,5) / QTO(J,K,I,1)
                  QTO(J,K,I,11) = - TAU(J,K,6) / QTO(J,K,I,1)
  510          CONTINUE
  520       CONTINUE
C
C 5-d.  Calculate Epsilon
C
C ITURBX = K-Epsilon
C
            IF (ITURBX .EQ. ITKE) THEN
               DO 540 K = 2, KDIM
                  DO 530 J = 2, JDIM
                     QTO(J,K,I,12) = QFROM(J,K,I,7)
  530             CONTINUE
  540          CONTINUE
C
C ITURBX = K-Omega
C        epsilon = C_mu * K * omega
C
            ELSE IF (ITURBX .EQ. ITKW) THEN
               DO 560 K = 2, KDIM
                  DO 550 J = 2, JDIM
                     QTO(J,K,I,12) = CMU*QFROM(J,K,I,6)*QFROM(J,K,I,7)
  550             CONTINUE
  560          CONTINUE
C
C  Error
C
            ELSE
               WRITE (IOUT,1005) ITURBX
            ENDIF
C
  590    CONTINUE
C
C Error - Unsupported Conversion
C
      ELSE
         WRITE (IOUT,1000) ITURBX, ITURB
         STOP
      ENDIF
C
C Finished converting restart file
C
      RETURN
 1000 FORMAT (' ','CONVRT: ERROR-> Unsupported turbulence model ',
     1            'conversion specified.',
     2       /' ',15X,'Turbulence model FROM = ',I5,
     3        ' ',15X,'Turbulence model TO   = ',I5,
     4       /' ', 8X,'RUN ABORTING!')
 1005 FORMAT (' ','CONVRT: ERROR-> Unsupported two equation model to ',
     1            'convert to Reynolds Stress: ITURBX = ',I5)
      END






      SUBROUTINE DTFACT (IDIM, JDIM, KDIM, I, Q, PFACTR)
C
C     Routine to calculate time step factor to reduce time step at a shock.
C     Based on discussions and code of Jeff White.
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     I              : I Plane (cell center) to perform calculations at 
C     Q              : Primitive variables at cell centers
C     PFACTR         : Factor to reduce CFL by based on pressure gradient
C                      PFACTR = 1 : Maintain time step
C                      PFACTR = 0 : Set time step to minimum value
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ), 
     1          PFACTR(JDIM,KDIM) 
C
C     Calculate factor based on pressure switch to reduce time step at a shock
C
      IF ( .true. ) THEN
         DO 120 K = 2, KDIM
            DO 110 J = 2, JDIM
               PDX  = 0.5E0 * (Q(J,K,I+1,5) - Q(J,K,I-1,5)) / Q(J,K,I,5)
               PDY  = 0.5E0 * (Q(J+1,K,I,5) - Q(J-1,K,I,5)) / Q(J,K,I,5)
               PDZ  = 0.5E0 * (Q(J,K+1,I,5) - Q(J,K-1,I,5)) / Q(J,K,I,5)
C
C     If marching scheme set pdx to zero and only look at shocks in crossflow
C
               IF (ISOLVR .EQ. IMARCH) THEN
                  PDX = 0.E0
               ENDIF
C
               PDEL = SQRT ( PDX*PDX + PDY*PDY + PDZ*PDZ )
               PFACTR(J,K) = (1.E0 - TANH(PDEL))**2
  110       CONTINUE
  120    CONTINUE
C
C     Set up to maintain time step
C
      ELSE
         DO 220 K = 2, KDIM
            DO 210 J = 2, JDIM
               PFACTR(J,K) = 1.E0
  210       CONTINUE
  220    CONTINUE
      ENDIF
C
C     Finished
C
      RETURN
      END




      SUBROUTINE DELTAT (IBLOCK, IDIM, JDIM, KDIM, IPSTRT, IPEND,
     1                   Q, PROPS, DTJ, VOL, S, PFACTR, 
     2                   VISCOS, ITURB, ITER, DTCFLC,
     3                   DTMIN, IBLKMN, IDTMIN, JDTMIN, KDTMIN,
     4                   DTMAX, IBLKMX, IDTMAX, JDTMAX, KDTMAX)
C
C Routine to calculate time step divided by the cell volume and the
C minimum/maximum value and location of DT/CFL
C IF DTINPT > 0   DT = DTINPT             - CONSTANT TIME STEP
C    DTINPT < 0   DT = DT(CFL=ABS(DT))   - CONSTANT CFL
C
C Use average of metrics in calculation of the cfl number
C
C This routine is organized as follows:
C     1. Set up control of time step increment and viscous effects
C     2. Calculate time step and min/max location
C        a. Constant time step (Time Accurate)
C        b. Constant CFL
C     3. Scale time step by Jacobian (inverse of volume) DTJ = dt / vol
C
C IBLOCK         : Current block
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at 
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C DTJ            : Time step divided by cell volume DTJ = DT * J = DT / VOL.
C VOL            : Cell volume
C S              : Metrics
C PFACTR         : Factor to reduce time step by at a shock
C VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C ITURB          : Turbulence model
C ITER           : Current iteration used to vary time step
C DTCFLC         : Value of DT/CFL for current iteration 
C DTMIN          : Minimum time step or CFL
C IBLKMN,IDTMIN, : Location of minimum time step or CFL (Block,I,J,K)
C  JDTMIN,KDTMIN
C DTMAX          : Maximum time step or CFL
C IBLKMX,IDTMAX, : Location of maximum time step or CFL (Block,I,J,K)
C  JDTMAX,KDTMAX
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ), 
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP), 
     2          DTJ   (JDIM,KDIM,IDIM), 
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1), 
     4          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3)
C
      DIMENSION PFACTR(JDIM,KDIM)
C
      LOGICAL VISCOS(3)
C
C 1.  Set up control of time step increment and viscous effects
C
C     Set up factor to calculate new time step for ramping
C
      IF (CFLFNL .EQ. 0) THEN
         FACT = 0.E0
      ELSE
         IF (ITER .GT. ITDFNL) THEN
            FACT = 1.E0
         ELSE 
            FACT = FLOAT(ITER-1)/FLOAT(ITDFNL-1)
         ENDIF
      ENDIF
C
C     Set up constants to include viscous effects and three-dimensional 
C     effectsin the CFL
C
      IF ( VISCOS(1) ) THEN
         VISI = 1.E0
      ELSE
         VISI = 0.E0
      ENDIF
C
      IF ( VISCOS(2) ) THEN
         VISJ = 1.E0
      ELSE
         VISJ = 0.E0
      ENDIF
C
      IF ( VISCOS(3) ) THEN
         VISK = 1.E0
      ELSE
         VISK = 0.E0
      ENDIF
C
      IF (THREED .OR. AXISYM) THEN
         R3D = 1.E0
      ELSE
         R3D = 0.E0
      ENDIF
C
C 2.   Calculate time step and min/max location
C
C 2-a. Constant time step (Time Accurate)
C
      IF (DT .GT. 0.0E0) THEN
         DTCFLC  = DT + FACT*(CFLFNL-DT)
         DO 230 I = IPSTRT, IPEND
            DO 220 K = 2, KDIM
               DO 210 J = 2, JDIM
                  DTJ(J,K,I) = DTCFLC  
C
                  RHO  = Q(J,K,I,1)
                  U    = Q(J,K,I,2)
                  V    = Q(J,K,I,3)
                  W    = Q(J,K,I,4)
                  P    = Q(J,K,I,5)
C
                  A    = SQRT (GAMMA*P/RHO)
C
                  RMU  = PROPS(J,K,I,1)
                  RMUT = PROPS(J,K,I,2)
C
C Calculate the average metric quantities at the cell center.
C
                  DXI = 0.5E0 * (S(J,K,I-1,1,1) + S(J,K,I,1,1))
                  DYI = 0.5E0 * (S(J,K,I-1,2,1) + S(J,K,I,2,1))
                  DZI = 0.5E0 * (S(J,K,I-1,3,1) + S(J,K,I,3,1))
                  DSI = 0.5E0 * (S(J,K,I-1,4,1) + S(J,K,I,4,1))
C
                  DXJ = 0.5E0 * (S(J-1,K,I,1,2) + S(J,K,I,1,2))
                  DYJ = 0.5E0 * (S(J-1,K,I,2,2) + S(J,K,I,2,2))
                  DZJ = 0.5E0 * (S(J-1,K,I,3,2) + S(J,K,I,3,2))
                  DSJ = 0.5E0 * (S(J-1,K,I,4,2) + S(J,K,I,4,2))
C
C Remove the three-dimensional influence in CFL if a 2D case
C
                  DXK = 0.5E0 * (S(J,K-1,I,1,3) + S(J,K,I,1,3)) * R3D
                  DYK = 0.5E0 * (S(J,K-1,I,2,3) + S(J,K,I,2,3)) * R3D
                  DZK = 0.5E0 * (S(J,K-1,I,3,3) + S(J,K,I,3,3)) * R3D
                  DSK = 0.5E0 * (S(J,K-1,I,4,3) + S(J,K,I,4,3)) * R3D
C
C Calculate the invscid and viscous eigenvalues in the three coordinates
C
                  RLAMII     = DSI * ( ABS(DXI*U + DYI*V + DZI*W) + A ) 
                  RLAMIJ     = DSJ * ( ABS(DXJ*U + DYJ*V + DZJ*W) + A )
                  RLAMIK     = DSK * ( ABS(DXK*U + DYK*V + DZK*W) + A )
C
                  RMUI       = 0.25E0 * (PROPS(J,K,I-1,1) + 2.E0*RMU  +
     1                                   PROPS(J,K,I+1,1))
                  RMUTI      = 0.25E0 * (PROPS(J,K,I-1,2) + 2.E0*RMUT + 
     1                                   PROPS(J,K,I+1,2))
                  RLAMVI     = 4.E0 * VISI * DSI*DSI * GAMMA*FSMACH/RE*
     1                       (RMUI/PR + RMUTI/PRT) / (RHO * VOL(J,K,I))
C
                  RMUJ       = 0.25E0 * (PROPS(J-1,K,I,1) + 2.E0*RMU  +
     1                                   PROPS(J+1,K,I,1))
                  RMUTJ      = 0.25E0 * (PROPS(J-1,K,I,2) + 2.E0*RMUT + 
     1                                   PROPS(J+1,K,I,2))
                  RLAMVJ     = 4.E0 * VISJ * DSJ*DSJ * GAMMA*FSMACH/RE*
     1                       (RMUJ/PR + RMUTJ/PRT) / (RHO * VOL(J,K,I))
C
                  RMUK       = 0.25E0 * (PROPS(J,K-1,I,1) + 2.E0*RMU  +
     1                                   PROPS(J,K+1,I,1))
                  RMUTK      = 0.25E0 * (PROPS(J,K-1,I,2) + 2.E0*RMUT + 
     1                                   PROPS(J,K+1,I,2))
                  RLAMVK     = 4.E0 * VISK * DSK*DSK * GAMMA*FSMACH/RE*
     1                       (RMUK/PR + RMUTK/PRT) / (RHO * VOL(J,K,I))
C
C Calculate CFL as per Anderson, Thomas and van Leer for the inviscid
C portion (identical to Swanson et al) and include the viscous
C part due to Swanson et al.
C
                  CFLT       = DTCFLC * ( RLAMII + RLAMIJ + RLAMIK +
     1                                    RLAMVI + RLAMVJ + RLAMVK ) /
     2                                  VOL(J,K,I)
C
C Find minimum and maximum CFL
C
                  IF (CFLT .GT. DTMAX) THEN
                     DTMAX  = CFLT
                     IBLKMX = IBLOCK
                     IDTMAX = I
                     JDTMAX = J
                     KDTMAX = K
                  ENDIF
                  IF (CFLT .LT. DTMIN) THEN
                     DTMIN  = CFLT
                     IBLKMN = IBLOCK
                     IDTMIN = I
                     JDTMIN = J
                     KDTMIN = K
                  ENDIF
 210           CONTINUE
 220        CONTINUE
 230     CONTINUE
C
C 2-b. Constant CFL
C
      ELSE
         DTCFLC   = CFL + FACT*(CFLFNL-CFL)
         CFLSHK = MIN (1.E0, DTCFLC)
         DO 260 I = IPSTRT, IPEND
C
C      Calculate factor to reduce time step by
C
            CALL DTFACT (IDIM, JDIM, KDIM, I, Q, PFACTR)
C
            DO 250 K = 2, KDIM
               DO 240 J = 2, JDIM
                  RHO  = Q(J,K,I,1)
                  U    = Q(J,K,I,2)
                  V    = Q(J,K,I,3)
                  W    = Q(J,K,I,4)
                  P    = Q(J,K,I,5)
C
                  A    = SQRT (GAMMA*P/RHO)
C
                  RMU  = PROPS(J,K,I,1)
                  RMUT = PROPS(J,K,I,2)
C
C Calculate the average metric quantities at the cell center.
C
                  DXI = 0.5E0 * (S(J,K,I-1,1,1) + S(J,K,I,1,1))
                  DYI = 0.5E0 * (S(J,K,I-1,2,1) + S(J,K,I,2,1))
                  DZI = 0.5E0 * (S(J,K,I-1,3,1) + S(J,K,I,3,1))
                  DSI = 0.5E0 * (S(J,K,I-1,4,1) + S(J,K,I,4,1))
C
                  DXJ = 0.5E0 * (S(J-1,K,I,1,2) + S(J,K,I,1,2))
                  DYJ = 0.5E0 * (S(J-1,K,I,2,2) + S(J,K,I,2,2))
                  DZJ = 0.5E0 * (S(J-1,K,I,3,2) + S(J,K,I,3,2))
                  DSJ = 0.5E0 * (S(J-1,K,I,4,2) + S(J,K,I,4,2))
C
C Remove the three-dimensional influence in CFL if a 2D case
C
                  DXK = 0.5E0 * (S(J,K-1,I,1,3) + S(J,K,I,1,3)) * R3D
                  DYK = 0.5E0 * (S(J,K-1,I,2,3) + S(J,K,I,2,3)) * R3D
                  DZK = 0.5E0 * (S(J,K-1,I,3,3) + S(J,K,I,3,3)) * R3D
                  DSK = 0.5E0 * (S(J,K-1,I,4,3) + S(J,K,I,4,3)) * R3D
C
C Calculate local CFL
C
                  CFLLCL     = CFLSHK + PFACTR(J,K) * (DTCFLC - CFLSHK)
C
C Calculate the invscid and viscous eigenvalues in the three coordinates
C
                  RLAMII     = DSI * ( ABS(DXI*U + DYI*V + DZI*W) + A ) 
                  RLAMIJ     = DSJ * ( ABS(DXJ*U + DYJ*V + DZJ*W) + A )
                  RLAMIK     = DSK * ( ABS(DXK*U + DYK*V + DZK*W) + A )
C
                  RMUI       = 0.25E0 * (PROPS(J,K,I-1,1) + 2.E0*RMU  +
     1                                   PROPS(J,K,I+1,1))
                  RMUTI      = 0.25E0 * (PROPS(J,K,I-1,2) + 2.E0*RMUT + 
     1                                   PROPS(J,K,I+1,2))
                  RLAMVI     = 4.E0 * VISI * DSI*DSI * GAMMA*FSMACH/RE*
     1                       (RMUI/PR + RMUTI/PRT) / (RHO * VOL(J,K,I))
C
                  RMUJ       = 0.25E0 * (PROPS(J-1,K,I,1) + 2.E0*RMU  +
     1                                   PROPS(J+1,K,I,1))
                  RMUTJ      = 0.25E0 * (PROPS(J-1,K,I,2) + 2.E0*RMUT + 
     1                                   PROPS(J+1,K,I,2))
                  RLAMVJ     = 4.E0 * VISJ * DSJ*DSJ * GAMMA*FSMACH/RE*
     1                       (RMUJ/PR + RMUTJ/PRT) / (RHO * VOL(J,K,I))
C
                  RMUK       = 0.25E0 * (PROPS(J,K-1,I,1) + 2.E0*RMU  +
     1                                   PROPS(J,K+1,I,1))
                  RMUTK      = 0.25E0 * (PROPS(J,K-1,I,2) + 2.E0*RMUT + 
     1                                   PROPS(J,K+1,I,2))
                  RLAMVK     = 4.E0 * VISK * DSK*DSK * GAMMA*FSMACH/RE*
     1                       (RMUK/PR + RMUTK/PRT) / (RHO * VOL(J,K,I))
C
C Calculate CFL as per Anderson, Thomas and van Leer for the inviscid
C portion (identical to Swanson et al) and include the viscous
C part due to Swanson et al.
C
                  DTJ(J,K,I) = CFLLCL * VOL(J,K,I) / 
     1                         ( RLAMII + RLAMIJ + RLAMIK +
     2                           RLAMVI + RLAMVJ + RLAMVK )
 240           CONTINUE
 250        CONTINUE
 260     CONTINUE
C
C Calculate the minimum and maximum DT for constant CFL time-step.  
C
         DO 290 I = IPSTRT, IPEND
            DO 280 K = 2, KDIM
               DO 270 J = 2, JDIM
                  IF (DTJ(J,K,I) .GT. DTMAX) THEN
                     DTMAX  = DTJ(J,K,I)
                     IBLKMX = IBLOCK
                     IDTMAX = I
                     JDTMAX = J
                     KDTMAX = K
                  ENDIF
                  IF (DTJ(J,K,I) .LT. DTMIN) THEN
                     DTMIN  = DTJ(J,K,I)
                     IBLKMN = IBLOCK
                     IDTMIN = I
                     JDTMIN = J
                     KDTMIN = K
                  ENDIF
  270          CONTINUE
  280       CONTINUE
  290    CONTINUE
C
      ENDIF
C
C 3. Scale time step by Jacobian (inverse of volume) DTJ = dt / vol
C
      DO 330 I = IPSTRT, IPEND
         DO 320 K = 2, KDIM
            DO 310 J = 2, JDIM
               DTJ(J,K,I) = DTJ(J,K,I) / VOL(J,K,I)
  310       CONTINUE
  320    CONTINUE
  330 CONTINUE
C
C Finished
C
      RETURN
      END




      SUBROUTINE FIXQ (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, Q, 
     1                 QMIN, NFIXQ, IERRCD)
C
C Routine to ensure positivity of the variables controlled by the QMIN array.  
C If QMIN(L) >= 0, then perform the positivity correction for the L variable.
C
C Two types of positivty preservation are provided:
C     1. Set Q = QMIN
C     2. Set Q = average of neighboring points
C
C This routine is organized as follows
C     1. Check for positivity of density and pressure; return error code
C     2. Loop over all turbulence variables and enforce positivity 
C        a. Calculate Q = QMIN positivity preservation
C        b. Calculate Q = average of neighboring points
C           i.  MARCHing routine (average only in the plane)
C           ii. Time dependent routine
C        c. Error - invalid option
C
C ITURB          : Turbulence model
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at 
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C Q              : Primitive variables at cell centers
C QMIN           : Control array for minimum value of Q
C NFIXQ          : Number of times reset Q
C IERRCD         : Error code
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION QMIN  (NQ)
C
      DIMENSION NFIXQ (NQ)
C
C 1.  Check for positivity of density and pressure; return error code
C
      DO I = IPSTRT, IPEND
         DO K = 2, KDIM
            DO J = 2, JDIM
               IF (Q(J,K,I,1) .LT. QMIN(1)) THEN
                  NFIXQ(1)   = NFIXQ(1) + 1
                  Q(J,K,I,1) = QMIN(1)
                  IERRCD     = IYES
               ENDIF
               IF (Q(J,K,I,5) .LT. QMIN(5)) THEN
                  NFIXQ(5)   = NFIXQ(5) + 1
                  Q(J,K,I,5) = QMIN(5)
                  IERRCD     = IYES
               ENDIF
            END DO
         END DO
      END DO
C
      IF (IERRCD .EQ. IYES) THEN
         RETURN
      END IF
C
C 2.  Loop over all turbulence variables and perform the positivity fix if
C     QMIN(L) >= 0
C
      DO 900 L = 6, NQ
         IF (QMIN(L) .GE. 0.E0) THEN
C
C 2-a. Set Q = QMIN
C
            IF (IFIXQ .EQ. IFIXMN) THEN
               DO 90 I = IPSTRT, IPEND
                  DO 80 K = 2, KDIM
                     DO 70 J = 2, JDIM
                        IF (Q(J,K,I,L) .LT. QMIN(L)) THEN
                           NFIXQ(L)   = NFIXQ(L) + 1
                           Q(J,K,I,L) = QMIN(L)
                        ENDIF
   70                CONTINUE
   80             CONTINUE
   90          CONTINUE
            ELSE IF (IFIXQ .EQ. IFIXAV) THEN
C
C 2-b. Set Q = average of neighboring points
C
C 2-b-i. MARCHing routine: average data only in the plane
C
               IF (ISOLVR .EQ. IMARCH) THEN
                  IF (THREED) THEN
                     DO 130 I = IPSTRT, IPEND
                        DO 120 K = 2, KDIM
CDIR$ IVDEP
                           DO 110 J = 2, JDIM
                              IF (Q(J,K,I,L) .LT. QMIN(L)) THEN
                                 NFIXQ(L)   = NFIXQ(L) + 1
                                 QNEW       = 
     1                              ( MAX(QMIN(L), Q(J+1,K,I,L))
     2                               + MAX(QMIN(L), Q(J-1,K,I,L))
     3                               + MAX(QMIN(L), Q(J,K+1,I,L))
     4                               + MAX(QMIN(L), Q(J,K-1,I,L)) )/4.E0
                                 Q(J,K,I,L) = QNEW
                              ENDIF 
  110                      CONTINUE
  120                   CONTINUE
  130                CONTINUE
                  ELSE
                     K = 2
                     DO 150 I = IPSTRT, IPEND
CDIR$ IVDEP
                        DO 140 J = 2, JDIM
                           IF (Q(J,K,I,L) .LT. QMIN(L)) THEN
                              NFIXQ(L)   = NFIXQ(L) + 1
                              QNEW = ( MAX(QMIN(L), Q(J+1,K,I,L))
     1                               + MAX(QMIN(L), Q(J-1,K,I,L)) )/2.E0
                              Q(J,K,I,L) = QNEW
                           ENDIF 
  140                   CONTINUE
  150                CONTINUE
                  ENDIF
C
C 2-b-ii. Time dependent section - average all neighboring points
C
               ELSE
                  IF (THREED) THEN
                     DO 230 I = IPSTRT, IPEND
                        DO 220 K = 2, KDIM
CDIR$ IVDEP
                           DO 210 J = 2, JDIM
                              IF (Q(J,K,I,L) .LT. QMIN(L)) THEN
                                 NFIXQ(L)   = NFIXQ(L) + 1
                                 QNEW       = 
     1                               ( MAX(QMIN(L), Q(J+1,K,I,L))
     2                               + MAX(QMIN(L), Q(J-1,K,I,L))
     3                               + MAX(QMIN(L), Q(J,K+1,I,L))
     4                               + MAX(QMIN(L), Q(J,K-1,I,L))
     5                               + MAX(QMIN(L), Q(J,K,I+1,L))
     6                               + MAX(QMIN(L), Q(J,K,I-1,L)) )/6.E0
                                 Q(J,K,I,L) = QNEW
                              ENDIF 
  210                      CONTINUE
  220                   CONTINUE
  230                CONTINUE
                  ELSE
                     K = 2
                     DO 250 I = IPSTRT, IPEND
CDIR$ IVDEP
                        DO 240 J = 2, JDIM
                           IF (Q(J,K,I,L) .LT. QMIN(L)) THEN
                              NFIXQ(L)   = NFIXQ(L) + 1
                              QNEW = ( MAX(QMIN(L), Q(J+1,K,I,L))
     1                               + MAX(QMIN(L), Q(J-1,K,I,L))
     2                               + MAX(QMIN(L), Q(J,K,I+1,L))
     3                               + MAX(QMIN(L), Q(J,K,I-1,L)) )/4.E0
                              Q(J,K,I,L) = QNEW
                           ENDIF 
  240                   CONTINUE
  250                CONTINUE
                  ENDIF
               ENDIF
            ELSE
C
C 2-c. Error - Invalid positivity type input
C
               WRITE (IOUT,1000) IFIXQ
               IERRCD = IYES
               RETURN
            ENDIF
C
C      Finished with this variable, do the next
C
         ENDIF
  900 CONTINUE
C
      RETURN
 1000 FORMAT (' ','FIXQ  : ERROR-> Invalid FIXQ option chosen ',I10,
     1       /' ','RUN ABORTING!')
      END




      SUBROUTINE GETTKE (ITURB, NPTS, Q, TKE)
C
C     Routine to calculate turbulent kinetic energy from Q
C
C     ITURB          : Turbulence model
C     NPTS           : Size of array to calculate fluxes on.
C     Q              : Flow variables
C     TKE            : Turbulent kinetic energy
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          TKE   (NPTS)
C
C     Calculate the Turbulence Kinetic Energy (TKE)
C
C     b. Turbulence Kinetic Energy for Two Equation Turbulence models
C
      IF (ITURB .EQ. ITKE .OR. ITURB .EQ. ITKW) THEN
         DO 100 I = 1, NPTS
            TKE(I) = Q(I,6)
  100    CONTINUE
C
C     c. Turbulence Kinetic Energy for Reynolds Stress Turbulence model
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         DO 110 I = 1, NPTS
            TAUXX  = Q(I,6)
            TAUYY  = Q(I,7)
            TAUZZ  = Q(I,8)
            TKE(I) = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
  110    CONTINUE
C
C     a. Turbulence Kinetic Energy = 0 for Inviscid, Laminar and
C          Algebraic Turbulence models
C
      ELSE
         DO 120 I = 1, NPTS
            TKE(I) = 0.E0
  120    CONTINUE
      ENDIF
C
      RETURN
      END






      SUBROUTINE IDENT (AI)
C
C     Routine to initialize the matrix AI to the identity
C     matrix of rank NRANK.
C
C
      include '../header/common.h'
C
      DIMENSION AI    (NRANK,NRANK)
C
C     Set Identity Matrix up
C     first initialize the entire matrix to zero and then reset the
C     diagonals to one
C
      DO 120 L = 1, NRANK
         DO 110 M = 1, NRANK
            AI(L,M) = 0.E0
  110    CONTINUE
  120 CONTINUE
C
      DO 130 L = 1, NRANK
         AI(L,L) = 1.E0
  130 CONTINUE
C
      RETURN
      END






      SUBROUTINE INIT (IDIM, JDIM, KDIM, Q, ITURB)
C
C Routine to initialize the flowfield to non-dimensionalized values.
C
C This routine is organized as follows:
C     1. Initialize all points to freestream values
C        a. Mean Flow
C        b. Turbulence Quantities
C           i.   k-epsilon
C           ii.  k-omega
C           iii. Reynolds stress
C
C IDIM,JDIM,KDIM : Dimensions of current block
C Q              : Primitive variables at cell centers
C ITURB          : Turbulence model
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
C 1.   Initialize all points to freestream data
C 1-a. Mean Flow
C      Initialize gas variables
C
      GAMM1 = GAMMA - 1.E0
      GSQM1 = GAMMA * GAMMA - 1.E0
C
C Calculate freestream quantities based on the following
C storage and non-dimensionalization:
C      Q(1) = rho / rhoinf
C      Q(2) =   u / ainf
C      Q(3) =   v / ainf
C      Q(4) =   w / ainf
C      Q(5) =   p / pref
C where p = (gamm1 - 1) * (e - 0.5 * rho * q**2)
C       pref = rhoinf * ainf**2
C
      PI     = ACOS(-1.E0)
C
C Alpha in radians
C
      ALPRAD = ALPHA*PI/180.E0
C
C If doing a two-dimensional case, then use the I-J plane
C
      IF ( THREED ) THEN
         UINF   = FSMACH * (COS(ALPRAD))
         VINF   = 0.E0
         WINF   = FSMACH * (SIN(ALPRAD))
      ELSE
         UINF   = FSMACH * (COS(ALPRAD))
         VINF   = FSMACH * (SIN(ALPRAD))
         WINF   = 0.E0
      ENDIF
      PREF   = 1.E0 / GAMMA
      EINF   = PREF / (GAMM1 * RHOINF) +
     1         0.5E0 * (UINF*UINF + VINF*VINF + WINF*WINF)
      WRITE(IOUT,100)UINF, VINF, WINF
  100 FORMAT(" Freestream velocity:", 3(e18.10))
C
C Initialize the flow field
C
      DO 230 I = 0, IDIM+2
         DO 220 K = 0, KDIM+2
            DO 210 J = 0, JDIM+2
               Q(J,K,I,1)  = RHOINF
               Q(J,K,I,2)  = UINF
               Q(J,K,I,3)  = VINF
               Q(J,K,I,4)  = WINF
               Q(J,K,I,5)  = PREF
 210        CONTINUE
 220     CONTINUE
 230  CONTINUE
C
C 1-b-i.   Initialization for k-epsilon two equation model
C          Turbulence intensity = sqrt (k) / u_inf -> Charles Speziale
C
      IF (ITURB .EQ. ITKE) THEN
         EINF = EINF + TKEINF
         DO 413 I = 0, IDIM+2
            DO 412 K = 0, KDIM+2
               DO 411 J = 0, JDIM+2
                  Q(J,K,I,6) = TKEINF
                  Q(J,K,I,7) = EPSINF
  411          CONTINUE
  412       CONTINUE
  413    CONTINUE
C
C 1-b-ii.  Initialization for k-omega two equation model
C          Turbulence intensity = sqrt (k) / u_inf -> Charles Speziale
C
      ELSE IF (ITURB .EQ. ITKW) THEN
         EINF = EINF + TKEINF
         DO 423 I = 0, IDIM+2
            DO 422 K = 0, KDIM+2
               DO 421 J = 0, JDIM+2
                  Q(J,K,I,6) = TKEINF
                  Q(J,K,I,7) = OMEINF
  421          CONTINUE
  422       CONTINUE
  423    CONTINUE
C
C 1-b-iii. Initialization for Reynolds Stress Turbulence model
C          Turbulence intensity = sqrt (k) / u_inf -> Charles Speziale
C
      ELSE IF (ITURB .EQ. ITRS) THEN
         EINF = EINF + TKEINF
         DO 453 I = 0, IDIM+2
            DO 452 K = 0, KDIM+2
               DO 451 J = 0, JDIM+2
                  Q(J,K,I, 6) = TAUINF(1)
                  Q(J,K,I, 7) = TAUINF(2)
                  Q(J,K,I, 8) = TAUINF(3)
                  Q(J,K,I, 9) = TAUINF(4)
                  Q(J,K,I,10) = TAUINF(5)
                  Q(J,K,I,11) = TAUINF(6)
                  Q(J,K,I,12) = EPSINF
  451          CONTINUE
  452       CONTINUE
  453    CONTINUE
      ENDIF
C
C Finished with the flowfield initialization.
C
      RETURN
      END




      SUBROUTINE INITBC (IDIM, JDIM, KDIM, Q, QTEST, ITURB)
C
C Routine to initialize the boundary data to infinity to verify input bc data.
C
C This routine is organized as follows:
C     1. Initialize QTEST
C     2. Initialize boundary to QTEST to test bc's set
C        a. I faces
C        b. J faces
C        c. K faces - If three dimensional
C
C IDIM,JDIM,KDIM : Dimensions of current block
C Q              : Primitive variables at cell centers
C QTEST          : Test array to initialize Q on boundaries
C ITURB          : Turbulence model
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QTEST (NQ)
C
C 1.   Initialize QTEST - density, pressure, turbulent quantities to infinity;
C                         velocities to zero
C
      DO 100 L = 1, NQ
         QTEST(L) = RMAX
  100 CONTINUE
C
      QTEST(2) = 0.E0
      QTEST(3) = 0.E0
      QTEST(4) = 0.E0
C
C 2.   Initialize boundary points to QTEST
C 2-a. I-faces
C
      DO 230 L = 1, NQ
         DO 220 K = 2, KDIM
            DO 210 J = 2, JDIM
               Q(J,K,     0,L) = QTEST(L)
               Q(J,K,     1,L) = QTEST(L)
               Q(J,K,IDIM+1,L) = QTEST(L)
               Q(J,K,IDIM+2,L) = QTEST(L)
  210       CONTINUE
  220    CONTINUE
  230 CONTINUE
C
C 2-b. J-faces
C
      DO 260 L = 1, NQ
         DO 250 K = 2, KDIM
            DO 240 I = 2, IDIM
               Q(     0,K,I,L) = QTEST(L)
               Q(     1,K,I,L) = QTEST(L)
               Q(JDIM+1,K,I,L) = QTEST(L)
               Q(JDIM+2,K,I,L) = QTEST(L)
  240       CONTINUE
  250    CONTINUE
  260 CONTINUE
C
C 2-c. K-faces - If Three Dimensional
C
      IF (THREED .OR. AXISYM) THEN
         DO 290 L = 1, NQ
            DO 280 I = 2, IDIM
               DO 270 J = 2, JDIM
                  Q(J,     0,I,L) = QTEST(L)
                  Q(J,     1,I,L) = QTEST(L)
                  Q(J,KDIM+1,I,L) = QTEST(L)
                  Q(J,KDIM+2,I,L) = QTEST(L)
  270          CONTINUE
  280       CONTINUE
  290    CONTINUE
      ENDIF
C
C Finished with the boundary data initialization.
C
      RETURN
      END




      SUBROUTINE INITPR (IDIM, JDIM, KDIM, PROPS, ITURB)
C
C Routine to initialize the Properties array for the flowfield to CMU
C
C IDIM,JDIM,KDIM : Dimensions of current block
C PROPS          : Properties at cell centers
C ITURB          : Turbulence model
C
      include '../header/common.h'
C
      DIMENSION PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP)
C
      DO 130 I = 0, IDIM+2
         DO 120 K = 0, KDIM+2
            DO 110 J = 0, JDIM+2
               PROPS (J,K,I,ICMUST) = CMU
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C Finished initializing CMUSTR
C
      RETURN
      END





      SUBROUTINE INITQN (IDIM, JDIM, KDIM, Q, QN, QCP, TKE, IMODEL)
C
C Routine to initialize QN to Q.
C
C N.B. This routine stores the previous time level solutions as conserved 
C      variables rather than primitive variables.
C
C This routine is organized as follows:
C     1. Copy Q to work array
C     2. Convert Q from primitive to conserved variables
C     3. Initialize QN to Q
C
C IDIM,JDIM,KDIM : Dimensions of current block
C Q              : Primitive variables at iterative time level p
C QN             : Conserved variables at time levels n, n-1, etc.
C QCP            : Work array for conserved Q at level p
C TKE            : Work array for turbulent kinetic energy
C IMODEL         : Turbulence model
C
      include '../header/common.h'
C
      INTEGER   IDIM, JDIM, KDIM
C
      DIMENSION Q   (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QN  (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ,NTMLVL)
C
      DIMENSION QCP (0:JDIM+2,0:KDIM+2,NQ),
     1          TKE (0:JDIM+2,0:KDIM+2)
C
      DIMENSION IMODEL(NMDL)
C
C Local variables
C
      INTEGER   I, J, K, L, NPTS
C
C     Turbulence model
C
      ITURB = IMODEL(1)
C
C     Loop over each cross section and convert Q to conserved then store in QN
C
      DO 900 I = 0, IDIM+2
C
C 1.  Copy Q to work array
C
         DO 130 L = 1, NQ
            DO 120 K = 0, KDIM+2
               DO 110 J = 0, JDIM+2
                  QCP(J,K,L) = Q(J,K,I,L)
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
C
C 2.  Convert Q from primitive to conserved
C
         NPTS = (JDIM+3) * (KDIM+3)
         CALL CONP2C (ITURB, NPTS, QCP, TKE)
C
C 3.  Initialize QN to Q
C
         DO 340 ITIM = 1, NTMLVL
            DO 330 L = 1, NQ
               DO 320 K = 0, KDIM+2
                  DO 310 J = 0, JDIM+2
                     QN(J,K,I,L,ITIM) = QCP(J,K,L)
  310             CONTINUE
  320          CONTINUE
  330       CONTINUE
  340    CONTINUE
C
C     Next cross section
C
  900 CONTINUE
C
C     Finished initializing QN
C
      RETURN
      END





      SUBROUTINE INITTQ (IDIM, JDIM, KDIM, PROPS, NBCS, IBCDAT)
C
C Routine to initialize tau_wall and qdot_wall for wall functions.
C
C This routine is organized as follows:
C     1. Loop over all boundary conditions 
C     2. If wall function BC then set tau_wall and qdot_wall
C
C IDIM,JDIM,KDIM : Dimensions of current block
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
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
      DIMENSION PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP)
C
      DIMENSION IBCDAT(10,NBCS)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
C
C Initial wall values
C
      TAUW = TKEINF * RHOINF * RE / (3.E0 * FSMACH)
      QW   = 0.E0
C
C Setup memory for looping.
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
      IJKDIM(1) = IDIM
      IJKDIM(2) = JDIM
      IJKDIM(3) = KDIM
C
C 1. Loop over boundary conditions
C
      DO 900 IBC = 1, NBCS
C
C Store boundary condition data in mnemonic names
C
         IBCTYP   = IBCDAT(1,IBC)
C
C 2. Set initial values for wall function BCs
C
         IF (IBCTYP .EQ. IWALFN) THEN
            IDIR     = IBCDAT(2,IBC)
            INROUT   = IBCDAT(3,IBC)
            ISTRT(1) = IBCDAT(4,IBC)
            IEND (1) = IBCDAT(5,IBC)
            ISTRT(2) = IBCDAT(6,IBC)
            IEND (2) = IBCDAT(7,IBC)
            ISTRT(3) = IBCDAT(8,IBC)
            IEND (3) = IBCDAT(9,IBC)
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
C      IBCPT2 is the boundary condition location outside the block from IBCPT
C
C
            IF (INROUT .EQ. 1) THEN
               IBCPT2 = 0
            ELSE
               IBCPT2 = IJKDIM(IDIR) + 2
            ENDIF
C
            IC1 = ICY (IDIR, 1)
            IC2 = ICY (IDIR, 2)
            IC3 = ICY (IDIR, 3)
C
            DO 120 I3 = ISTRT(IC3), IEND(IC3)
CDIR$ IVDEP
               DO 110 I2 = ISTRT(IC2), IEND(IC2)
C
                  IBC2 = IBCPT2*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JBC2 = IBCPT2*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KBC2 = IBCPT2*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
C Initialize tau_wall and qdot_wall
C
                  PROPS (JBC2,KBC2,IBC2,ITQTAU) = TAUW
                  PROPS (JBC2,KBC2,IBC2,ITQDOT) = QW
  110          CONTINUE
  120       CONTINUE
         ENDIF
  900 CONTINUE
C
C Finished initializing tau_wall and qdot_wall
C
      RETURN
      END





      SUBROUTINE L2MAX  (IBLK, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
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
      include '../header/common.h'
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
     1          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION IRESMX(0:3)
C
C 1.  Loop to find the maximum residual.
C     Do not initialize the maximum residual to 0.  This is done in the main
C     routine so that the maximum is over the entire flowfield (i.e.,
C     including all zones) or by lines for marching.
C
      DO 130 I = IPSTRT, IPEND
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





      SUBROUTINE L2NORM (IDIM, JDIM, KDIM, IPSTRT, IPEND, 
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
      include '../header/common.h'
C
      DIMENSION RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF),
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
         DO 300 I = IPSTRT, IPEND
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





      SUBROUTINE QNSHFT (IDIM, JDIM, KDIM, Q, QN, QCP, TKE, IMODEL)
C
C Routine to shift Q, QN to prepare for the next time step.
C
C N.B. This routine assumes that the previous time level solutions are
C      passed as conserved variables rather than primitive variables.
C
C This routine is organized as follows:
C     1. Shift QN (QN(NTIME) <- QN(NTIME-1), etc)
C     2. Copy Q to work array
C     3. Convert Q from primitive to conserved variables
C     4. Store Q in QN(1)
C     5. Store Q as Q^(m-1) in QN(NTIME+1) for 2nd order pseudo time 
C        sub-iterations in tau-ts
C
C IDIM,JDIM,KDIM : Dimensions of current block
C Q              : Primitive variables at iterative time level p
C QN             : Conserved variables at time levels n, n-1, etc.
C QCP            : Work array for conserved Q at level p
C TKE            : Work array for turbulent kinetic energy
C IMODEL         : Turbulence model
C
      include '../header/common.h'
C
      INTEGER   IDIM, JDIM, KDIM
C
      DIMENSION Q   (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QN  (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ,NTMLVL)
C
      DIMENSION QCP (0:JDIM+2,0:KDIM+2,NQ),
     1          TKE (0:JDIM+2,0:KDIM+2)
C
      DIMENSION IMODEL(NMDL)
C
C Local variables
C
      INTEGER   I, J, K, L, NPTS
C
C     Turbulence model
C
      ITURB = IMODEL(1)
C
C 1.  Shift QN (QN(NTIME) <- QN(NTIME-1), etc)
C
      DO 150 ITIM = NTIME, 2, -1
         DO 140 L = 1, NQ
            DO 130 I = 0, IDIM+2
               DO 120 K = 0, KDIM+2
                  DO 110 J = 0, JDIM+2
                     QN(J,K,I,L,ITIM) = QN(J,K,I,L,ITIM-1)
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
C
C     Loop over each cross section and convert Q to conserved then store in QN
C
      DO 900 I = 0, IDIM+2
C
C 2.  Copy Q to work array
C
         DO 230 L = 1, NQ
            DO 220 K = 0, KDIM+2
               DO 210 J = 0, JDIM+2
                  QCP(J,K,L) = Q(J,K,I,L)
  210          CONTINUE
  220       CONTINUE
  230    CONTINUE
C
C 3.  Convert Q from primitive to conserved
C
         NPTS = (JDIM+3) * (KDIM+3)
         CALL CONP2C (ITURB, NPTS, QCP, TKE)
C
C 4.  Store Q in QN(1)
C
         DO 430 L = 1, NQ
            DO 420 K = 0, KDIM+2
               DO 410 J = 0, JDIM+2
                  QN(J,K,I,L,1) = QCP(J,K,L)
  410          CONTINUE
  420       CONTINUE
  430    CONTINUE
C
C 5.  Store Q as Q^(m-1) in QN(NTIME+1) for 2nd order pseudo time 
C     sub-iterations in tau-ts
C
         IF (ITMTYP .EQ. ITAUTS .AND. TDPHIP .NE. 0.E0) THEN
            IT = NTIME + 1
            DO 530 L = 1, NQ
               DO 520 K = 0, KDIM+2
                  DO 510 J = 0, JDIM+2
                     QN(J,K,I,L,IT) = QCP(J,K,L)
  510             CONTINUE
  520          CONTINUE
  530       CONTINUE
         ENDIF
C
C     Next cross section
C
  900 CONTINUE
C
C     Finished shifting Q, QN
C
      RETURN
      END






      SUBROUTINE UPDATE (IDIM, JDIM, KDIM, IPSTRT, IPEND, Q, RES, ITURB,
     1                   QC, TKE, QMIN, NFIXQ, NRELIZ, DELTAQ, IERRCD)
C
C Routine to update the Conserved Variables
C
C This routine is organized as follows:
C     1. Convert dQ_cell-average to dQ_cell-centered
C        Begin I loop {
C     2. Convert Q from primitive to conserved variables
C     3. Update conserved variables
C     4. Convert conserved variables back to primitive variables.
C        } End I loop
C     5. Check for positivity of variables
C     6. Enforce realizability of Reynolds stresses for Reynolds stress models
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at 
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C Q              : Primitive variables at cell centers
C RES            : Residual
C ITURB          : Turbulence model
C QC             : Work array to store variable on a plane
C TKE            : Work array to store turbulent kinetic energy on a plane
C QMIN           : Control array for minimum value of Q
C NFIXQ          : Number of times reset Q
C NRELIZ         : Number of times violate realizability
C DELTAQ         : Work area to calculate delta(Q)_cell-centered
C IERRCD         : Error code
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          RES   (JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF)
C
      DIMENSION QC    (2:JDIM,2:KDIM,NQ),
     1          TKE   (2:JDIM,2:KDIM)
C
      DIMENSION QMIN  (NQ)
      DIMENSION NFIXQ (NQ),
     1          NRELIZ(NQ)
C
      DIMENSION DELTAQ(JDIM+1,KDIM+1,IPSTRT-1:IPEND+1,NF)
C
C 1.  Convert dQ_cell-average to dQ_cell-centered
C     Currently this subroutine is empty. I don't know what it is 
C     supposed to do. Defined in fourth.F
C
      CALL DQCC   (IDIM, JDIM, KDIM, IPSTRT, IPEND, RES, DELTAQ)
C
C     Begin I loop {
C
      DO 490 I = IPSTRT, IPEND
C
C 2.   Convert Q from primitive variables to conservative variables
C
         DO 230 L = 1, NQ
            DO 220 K = 2, KDIM
               DO 210 J = 2, JDIM
                  QC(J,K,L) = Q(J,K,I,L)
  210          CONTINUE
  220       CONTINUE
  230    CONTINUE
C
         NPTS = (JDIM - 1) * (KDIM - 1)
         CALL CONP2C (ITURB, NPTS, QC, TKE)
C
C 3.  Update Conservative Variables
C
         DO 330 L = 1, NF
            DO 320 K = 2, KDIM
               DO 310 J = 2, JDIM
                  QC(J,K,L) = QC(J,K,L) + RES(J,K,I,L)
  310          CONTINUE
  320       CONTINUE
  330    CONTINUE
C
C 4.   Return Q to primitive variables
C
         CALL CONC2P (ITURB, NPTS, QC, TKE)
C
         DO 430 L = 1, NQ
            DO 420 K = 2, KDIM
               DO 410 J = 2, JDIM
                  Q(J,K,I,L) = QC(J,K,L)
  410          CONTINUE
  420       CONTINUE
  430    CONTINUE
C
C     } End I loop
C
  490 CONTINUE
C
C 5. Check for positivity of flow variables
C
      IF (IFIXQ .NE. INO) THEN
         CALL FIXQ   (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, Q, 
     1                QMIN,  NFIXQ, IERRCD)
      ENDIF
C
C 6. Enforce realizability of Reynolds stresses for Reynolds stress models
C
      IF (ITURB .EQ. ITRS) THEN
         CALL REALIZ (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, Q, NRELIZ)
      ENDIF
C
C     Finished with conservative variable update
C
      RETURN
      END






      SUBROUTINE YNORM  (NXLVLS, NXBLKS, NXBCS, NXNODE, 
     1                   IOFFS, IOFFP,
     2                   NLVLS,  NBLKS,  NBCS,  IBCDAT,
     3                   IDIM, JDIM, KDIM, R, PROPS, YNORML)
C
C Routine to calculate distance to wall from cell center in the flowfield.
C Initially calculate the distances as the square of the distance to avoid
C unnecessary sqrt calculations. When all of the distances are complete for
C a block/level, then convert from the distance squared to the distance
C (advice from Veer Vatsa).
C
C This routine is organized as follows:
C     1. Loop over all levels doing steps 2-4
C     2. Initialize all values of y_normal to infinity
C     3. For each grid point on a wall or wall function boundary
C        loop through all blocks and set y_normal to min(y_normal, dist)
C     4. Convert y_normal from distance squared to distance and
C        Interpolate y_normal calculated on grid points to cell centers
C     5. Modify y_normal at first 'cell center' for wall functions
C
C NXLVLS         : 
C NXBLKS         : 
C NXBCS          : 
C NXNODE         : 
C IOFFS          : 
C IOFFP          : 
C NLVLS          : Number of grid levels
C NBLKS          : Number of blocks
C NBCS           : Number of boundary conditions per current block
C IBCDAT         : Data controlling boundary conditions
C IDIM,JDIM,KDIM : Dimensions of each block
C R              : Grid
C PROPS(4)       : Normal distance to wall at cell centers
C YNORML         : Temporary storage for distance to wall at grid points
C
      include '../header/common.h'
C
      DIMENSION IOFFS  (NXBLKS,NXLVLS),
     1          IOFFP  (NXBLKS,NXLVLS)
C
      DIMENSION IDIM   (NXBLKS,NXLVLS),
     1          JDIM   (NXBLKS,NXLVLS),
     2          KDIM   (NXBLKS,NXLVLS)
C
      DIMENSION IBCDAT (10,NXBCS,NXBLKS,NXLVLS),
     1          NBCS   (NXBLKS)
C
      DIMENSION PROPS (*),
     1          R     (*),
     2          YNORML(*)
C
      WRITE(IOUT,*)
      WRITE(IOUT,*)'YNORM: Computing distance to wall ...'
C
C 1.  Loop over all grid levels calculating y_normal
C
      DO 900 ILVL = 1, NLVLS
C
C 2.  Initialize y_normal to infinity
C
         DO 200 IBLK = 1, NBLKS
            IY   = IOFFS(IBLK,ILVL)         + 1
            NPTS = (IDIM(IBLK,ILVL)+2) * (JDIM(IBLK,ILVL)+2) *
     1             (KDIM(IBLK,ILVL)+2)
            CALL YNINIT (NPTS, YNORML(IY))
  200    CONTINUE
C
C 3.  For each block, call routines to check each wall or wall function
C     grid point against all points of all blocks
C
         DO 300 IBLK = 1, NBLKS
            IRC  = IOFFS(IBLK,ILVL) * 3     + 1
            CALL YNWALL (NXNODE, IOFFS(1,ILVL), IOFFP(1,ILVL), 
     1                   NBLKS,  IDIM(1,ILVL),  JDIM(1,ILVL),  
     2                   KDIM(1,ILVL),  R,  YNORML,
     3                   IDIM(IBLK,ILVL), JDIM(IBLK,ILVL), 
     4                   KDIM(IBLK,ILVL), R(IRC), 
     5                   NBCS(IBLK), IBCDAT(1,1,IBLK,ILVL))
  300    CONTINUE
C
C 4.  Convert y_normal from distance squared to distance and
C     Interpolate grid point values to cell centers
C
         DO 400 IBLK = 1, NBLKS
            IY   = IOFFS(IBLK,ILVL)         + 1
            IP   = IOFFP(IBLK,ILVL) * NP    + 1
            NPTS = (IDIM(IBLK,ILVL)+2) * (JDIM(IBLK,ILVL)+2) *
     1             (KDIM(IBLK,ILVL)+2)
            CALL YNSQRT (NPTS, YNORML(IY))
            CALL YNCENT (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL), 
     1                   KDIM(IBLK,ILVL), YNORML(IY), PROPS(IP))
  400    CONTINUE
C
C 5.  Modify y_normal for first 'cell center' for wall functions
C
         IF ( IFWALF ) THEN
            DO 510 IBLK = 1, NBLKS
               IY   = IOFFS(IBLK,ILVL)         + 1
               IP   = IOFFP(IBLK,ILVL) * NP    + 1
               DO 500 IBC = 1, NBCS(IBLK)
                  IBCTYP = IBCDAT(1,IBC,IBLK,ILVL)
                  IF (IBCTYP .EQ. IWALFN) THEN
                     CALL YNWALF (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL), 
     1                            KDIM(IBLK,ILVL), PROPS(IP),
     2                            IBCDAT(1,IBC,IBLK,ILVL))
                  ENDIF
  500          CONTINUE
  510       CONTINUE
         ENDIF
C
C     Next grid level
C
  900 CONTINUE
C
C     Finished calculating y_normal
C
      RETURN
      END





      SUBROUTINE YNINIT (NPTS, YNORML)
C
C Routine to initialize y_normal to infinity
C
C NPTS      : Number of points to calculate distance on
C YNORML    : Normal distance to wall at grid points
C
      include '../header/common.h'
C
      DIMENSION YNORML(NPTS)
C
C     Initialize YNORML to <infinity> in case there are no walls
C
      DO 100 I = 1, NPTS
         YNORML(I) = RMAX
  100 CONTINUE
C
C     Finished with initialization
C
      RETURN
      END





      SUBROUTINE YNWALL (NXNODE, 
     1                   IOFFS, IOFFP, NBLKS,  
     2                   IDIM,  JDIM,  KDIM,  R,  YNORML,
     3                   IDIMC, JDIMC, KDIMC, RC, NBCS,  IBCDAT)
C
C This routine loops over all boundary conditions for a given block
C and calculates distances from wall or wall function points to all points
C in all blocks.
C
C This routine is organized as follows:
C     1. Loop over all boundary conditions doing steps 2-3 for WALL and WALFN
C     2. Loop over all grid points on the boundary doing step 3
C     3. Loop over ALL blocks (including current block) calculating
C        y_normal as min (y_normal, distance to wall point)
C
C NXNODE            : 
C IOFFS             : 
C IOFFP             : 
C NBLKS             : Number of blocks
C IDIM,JDIM,KDIM    : Dimensions of each block
C R                 : Grid
C YNORML            : Temporary storage for distance to wall at grid points
C IDIMC,JDIMC,KDIMC :
C RC                :
C NBCS              : Number of boundary conditions in current block
C IBCDAT            : Data controlling boundary conditions
C
      include '../header/common.h'
C
      DIMENSION IOFFS  (NBLKS),
     1          IOFFP  (NBLKS)
C
      DIMENSION IDIM   (NBLKS),
     1          JDIM   (NBLKS),
     2          KDIM   (NBLKS)
C
      DIMENSION IBCDAT (10,NBCS)
C
      DIMENSION R     (*),
     1          YNORML(*)
C
      DIMENSION RC    (0:JDIMC+1,0:KDIMC+1,0:IDIMC+1,3)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),
     2          IEND  (3)
C
      DIMENSION RWALL (3)
C
C     Store IDIMC, JDIMC, KDIMC in IJKDIM for boundary condition locations
C
      IJKDIM(1) = IDIMC
      IJKDIM(2) = JDIMC
      IJKDIM(3) = KDIMC
C
C 1.  Loop over WALL boundary conditions
C
      DO 500 IBC = 1, NBCS
         IBCTYP = IBCDAT(1,IBC)
         IF (IBCTYP .EQ. IWALL .OR. IBCTYP .EQ. IWALFN) THEN
            IDIR     = IBCDAT(2,IBC)
            INROUT   = IBCDAT(3,IBC)
            ISTRT(1) = IBCDAT(4,IBC)
            IEND (1) = IBCDAT(5,IBC)
            ISTRT(2) = IBCDAT(6,IBC)
            IEND (2) = IBCDAT(7,IBC)
            ISTRT(3) = IBCDAT(8,IBC)
            IEND (3) = IBCDAT(9,IBC)
C
C     Decrement beginning index to account for cell-center/grid point storage.
C     Modify begin, end indices to extend boundary condition to edge/corner
C
            DO 10 L = 1, 3
               IF (L .NE. IDIR) THEN
                  ISTRT(L) = ISTRT(L) - 1
                  IF (ISTRT(L) .EQ. 1        ) ISTRT(L) = 0
                  IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
               END IF
   10       CONTINUE
C
C     Set needed variables depending on whether the boundary is
C     the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C         IBCPT  is the grid boundary location
C
            IF (INROUT .EQ. 1) THEN
               IBCPT  = 1
            ELSE
               IBCPT  = IJKDIM(IDIR)
            ENDIF
C
            IC1 = ICY (IDIR, 1)
            IC2 = ICY (IDIR, 2)
            IC3 = ICY (IDIR, 3)
C
C 2.  Loop over all grid points on the current block boundary
C
            DO 220 I3 = ISTRT(IC3), IEND(IC3)
               DO 210 I2 = ISTRT(IC2), IEND(IC2)
C
                  IBCW = IBCPT*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JBCW = IBCPT*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KBCW = IBCPT*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  RWALL(1) = RC(JBCW,KBCW,IBCW,1)
                  RWALL(2) = RC(JBCW,KBCW,IBCW,2)
                  RWALL(3) = RC(JBCW,KBCW,IBCW,3)
C
C 3.  Loop over all blocks calculating distance from EACH boundary point 
C     to ALL grid points
C
                  DO 100 IBLK = 1, NBLKS
                     IR   = IOFFS(IBLK) * 3     + 1
                     IY   = IOFFS(IBLK)         + 1
                     NPTS = (IDIM(IBLK)+2) * (JDIM(IBLK)+2) *
     1                      (KDIM(IBLK)+2)
                     CALL YNDIST (NPTS, RWALL, R(IR), YNORML(IY))
  100             CONTINUE
  210          CONTINUE
  220       CONTINUE
         ENDIF
  500 CONTINUE
C
C     Finished
C
      RETURN
      END





      SUBROUTINE YNDIST (NPTS, RWALL, R, YNORML)
C
C Routine to calculate distance SQUARED from grid points to wall using simple
C distance formula. Using the distance squared speeds up the code by removing
C unnecessary sqrt calls.
C
C NPTS      : Number of points to calculate distance on
C RWALL     : Position of wall point
C R         : Grid
C YNORML    : Normal distance to wall at grid points
C
      include '../header/common.h'
C
      DIMENSION RWALL (3),
     1          R     (NPTS,3),
     2          YNORML(NPTS)
C
C 1.  Loop over all points calculating distance from grid point to wall point
C
      DO 100 I = 1, NPTS
         DX        = R(I,1) - RWALL(1)
         DY        = R(I,2) - RWALL(2)
         DZ        = R(I,3) - RWALL(3)
         YDIST     = (DX*DX + DY*DY + DZ*DZ)
         YNORML(I) = MIN (YNORML(I), YDIST)
  100 CONTINUE
C
C     Finished calculating YNORML
C
      RETURN
      END





      SUBROUTINE YNCENT (IDIM, JDIM, KDIM, YNORML, PROPS)
C
C Routine to interpolate y_normal from grid points to cell centers
C
C IDIM,JDIM,KDIM : Dimensions of current block
C YNORML         : Normal distance to wall at grid points
C PROPS(4)       : Normal distance to wall at cell centers
C
      include '../header/common.h'
C
      DIMENSION YNORML(0:JDIM+1,0:KDIM+1,0:IDIM+1),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP)
C
C     Calculate cell centered normal distance as average of 8 nodes
C
      DO 130 I = 1, IDIM + 1
         DO 120 K = 1, KDIM + 1
            DO 110 J = 1, JDIM + 1
               PROPS(J,K,I,4) = 0.125E0 * 
     1                     ( YNORML(J  ,K  ,I  ) + YNORML(J-1,K  ,I  )
     2                     + YNORML(J  ,K-1,I  ) + YNORML(J  ,K  ,I-1)
     3                     + YNORML(J-1,K-1,I  ) + YNORML(J-1,K  ,I-1)
     4                     + YNORML(J  ,K-1,I-1) + YNORML(J-1,K-1,I-1) )
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C     Finished interpolating to cell-centers
C
      RETURN
      END




      SUBROUTINE YNSQRT (NPTS, YNORML)
C
C Routine to convert y_normal from the distance squared to the distance
C
C NPTS           : Number of points
C YNORML         : Normal distance to wall at grid points
C
      DIMENSION YNORML(NPTS)
C
C     Convert y_normal from distance squared to distance
C
      DO 100 I = 1, NPTS
         YNORML(I) = SQRT (YNORML(I))
  100 CONTINUE
C
C     Finished converting y_normal from distance squared to distance
C
      RETURN
      END





      SUBROUTINE YNWALF (IDIM, JDIM, KDIM, PROPS, IBCDAT)
C
C Routine to modify y_normal at the first 'cell center' for wall functions
C
C IDIM,JDIM,KDIM : Dimensions of current block
C PROPS(4)       : Normal distance to wall at cell centers
C IBCDAT         : Data controlling boundary conditions
C
      include '../header/common.h'
C
      DIMENSION PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     1          IBCDAT(10)
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),
     2          IEND  (3)
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
      IJKDIM(1) = IDIM
      IJKDIM(2) = JDIM
      IJKDIM(3) = KDIM
C
C     Store boundary condition data in mnemonic names
C
      IDIR      = IBCDAT(2)
      INROUT    = IBCDAT(3)
      ISTRT(1)  = IBCDAT(4)
      IEND (1)  = IBCDAT(5)
      ISTRT(2)  = IBCDAT(6)
      IEND (2)  = IBCDAT(7)
      ISTRT(3)  = IBCDAT(8)
      IEND (3)  = IBCDAT(9)
C
C     Decrement beginning index to account for cell-center/grid point storage.
C     Modify begin, end indices to extend boundary condition to edge/corner
C
      DO 10 L = 1, 3
         IF (L .NE. IDIR) THEN
            IF (ISTRT(L) .EQ. 2        ) ISTRT(L) = 1
            IF (IEND (L) .EQ. IJKDIM(L)) IEND (L) = IJKDIM(L) + 1
         END IF
   10 CONTINUE
C
C     Set needed variables depending on whether the boundary is
C     the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C         IBCN   is the 1st point to the inside of the block from IBCPT
C         IBCN2  is the 2nd point to the inside of the block from IBCPT
C
      IF (INROUT .EQ. 1) THEN
         IBCN   = 2
         IBCN2  = 3
      ELSE
         IBCN   = IJKDIM(IDIR) 
         IBCN2  = IJKDIM(IDIR) - 1
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C     Loop over all grid points on the current block boundary
C
            DO 220 I3 = ISTRT(IC3), IEND(IC3)
               DO 210 I2 = ISTRT(IC2), IEND(IC2)
C
                  IN   = IBCN *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JN   = IBCN *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KN   = IBCN *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  IN2  = IBCN2*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JN2  = IBCN2*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KN2  = IBCN2*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  YINTRF = 2.E0 * PROPS(JN,KN,IN,4)
                  DELTAY = PROPS(JN2,KN2,IN2,4) - YINTRF
                  PROPS(JN,KN,IN,4) = YINTRF - DELTAY
  210    CONTINUE
  220 CONTINUE
C
C     Finished modifying first 'cell-center' y_normal
C
      RETURN
      END





      SUBROUTINE YNABS  (NPTS, PROPS)
C
C Routine to set y_normal to abs(y_normal) to allow the modification of
C the transition location on restart.
C
C NPTS      : Number of points to calculate distance on
C YNORML    : Normal distance to wall at cell centers
C
      include '../header/common.h'
C
      DIMENSION PROPS (NPTS,NP)
C
C     Set y_normal = abs(y_normal)
C
      DO 100 I = 1, NPTS
         PROPS(I,4) = ABS (PROPS(I,4))
  100 CONTINUE
C
C     Finished with setting y_normal = abs(y_normal)
C
      RETURN
      END





      SUBROUTINE YNTRAN (IDIM, JDIM, KDIM, PROPS, ITRSEG)
C
C Routine to set y_normal to - (y_normal) to specify the transition location.
C
C IDIM,JDIM,KDIM : Dimensions of blocks
C PROPS          : Properties at cell centers 
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity 
C                  PROPS(3) = Y+
C                  PROPS(4) = Ynorm
C ITRSEG         : Data specifying transition location
C
C This routine is organized as follows:
C     1. Set up transition indices
C     2. Loop over transition segment and set y_normal = - (y_normal)
C
      include '../header/common.h'
C
      DIMENSION PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP)
C
      DIMENSION ITRSEG(7)
C
C 1.  Set up transition indices
C
      ISTRT = ITRSEG(2)
      JSTRT = ITRSEG(3)
      KSTRT = ITRSEG(4)
      IEND  = ITRSEG(5)
      JEND  = ITRSEG(6)
      KEND  = ITRSEG(7)
C
C 2.  Loop over transition segment and set y_normal = - (y_normal)
C
      DO I = ISTRT, IEND
         DO K = KSTRT, KEND
            DO J = JSTRT, JEND
               PROPS(J,K,I,4) = - PROPS(J,K,I,4)
            END DO
         END DO
      END DO
C
C     Finished
C
      RETURN
      END
