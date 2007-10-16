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
      SUBROUTINE WRTP3D(IFPLT3, BINPLT, FPLT3G, FPLT3Q, FPLT3F, FPLT3FN,
     1                  FERROR, NBLKS, IDIM, JDIM, KDIM, IOFFS, IOFFQ,
     2                  IOFFP, ITURB, NPRINT, IPRINT, NCUTS, ICUTS,
     3                  ILVL, Q, PROPS, R)
      include '../header/maxdim.h'
      include '../header/common.h'
      include '../header/histry.h'
C
C     Data controlling PLOT3D output
C
      LOGICAL IFPLT3, BINPLT
C
C     Filenames
C     FPLT3G  File to write PLOT3D Grid file to
C     FPLT3Q  File to write PLOT3D Q (conserved variables) to
C     FPLT3F  File to write PLOT3D function file to
C     FPLT3FN File to write PLOT3D function name file to
C     FERROR  File name output if error
C
      CHARACTER*80 FPLT3G, FPLT3Q, FPLT3F, FPLT3FN, FERROR
C
C     Dimensions of domains
C
      DIMENSION IDIM  (MXBLKS,MXLVLS),
     1          JDIM  (MXBLKS,MXLVLS),
     2          KDIM  (MXBLKS,MXLVLS)
C
C     Offsets of domains used to calculate location in total
C     storage of Q,S,... of a given domain.
C
      DIMENSION IOFFS (MXBLKS,MXLVLS),
     2          IOFFQ (MXBLKS,MXLVLS),
     3          IOFFP (MXBLKS,MXLVLS)
C
C     Printout control data.
C     IPRINT, NPRINT control output at end of iteration cycle
C
      DIMENSION IPRINT(10,MXPRT,MXBLKS), NPRINT(MXBLKS)
C
C     Cut (domain to domain communication) data.
C
      DIMENSION ICUTS (21,MXCUTS,MXLVLS)
      INTEGER   NCUTS
C
C     Which grid level solution to save, should be 1, the finest grid
C
      INTEGER   ILVL
C
C     Arrays of solution, properties and coordinates
C
      DIMENSION  Q     (*), 
     1           PROPS (*),
     2           R     (*)
C
C     Dimensions of domains
C
C 14. OUTPUT 
C     Write output and PLOT3D output
C     (PLOT3D output in multiple grid format)
C     Use RES array as temporary space to interpolate the grid
C     nodal values to cell centers.
C
      WRITE(*,*)'WRTP3D: Writing plot3d output'

      IF ( IFPLT3 ) THEN
C
C Binary PLOT3D output
C
         IF ( BINPLT ) THEN
C Grid
            FERROR = FPLT3G
            OPEN (IPLT3G,FILE=FPLT3G,FORM='UNFORMATTED',ERR=990)
            WRITE (IPLT3G) NBLKS
            IF ( THREED ) THEN
               WRITE (IPLT3G) (IDIM(IBLK,ILVL)+1,JDIM(IBLK,ILVL)+1,
     1                         KDIM(IBLK,ILVL)+1,IBLK=1,NBLKS)
            ELSE
               WRITE (IPLT3G) (IDIM(IBLK,ILVL)+1,JDIM(IBLK,ILVL)+1,
     1                         IBLK=1,NBLKS)
            ENDIF
C Solution
            FERROR = FPLT3Q
            OPEN (IPLT3Q,FILE=FPLT3Q,FORM='UNFORMATTED',ERR=990)
            WRITE (IPLT3Q) NBLKS
            IF ( THREED ) THEN
               WRITE (IPLT3Q) (IDIM(IBLK,ILVL)+1,JDIM(IBLK,ILVL)+1,
     1                         KDIM(IBLK,ILVL)+1,IBLK=1,NBLKS)
            ELSE
               WRITE (IPLT3Q) (IDIM(IBLK,ILVL)+1,JDIM(IBLK,ILVL)+1,
     1                         IBLK=1,NBLKS)
            ENDIF
C Function file
            IF (NQ .GT. 5) THEN
               FERROR = FPLT3F
               NVAR   = (NQ - 5) + NP
               OPEN (IPLT3F,FILE=FPLT3F,FORM='UNFORMATTED',ERR=990)
               WRITE (IPLT3F) NBLKS
               IF ( THREED ) THEN
                  WRITE (IPLT3F) (IDIM(IBLK,ILVL)+1,JDIM(IBLK,ILVL)+1,
     1                            KDIM(IBLK,ILVL)+1,NVAR,IBLK=1,NBLKS)
               ELSE
                  WRITE (IPLT3F) (IDIM(IBLK,ILVL)+1,JDIM(IBLK,ILVL)+1,
     1                            NVAR,IBLK=1,NBLKS)
               ENDIF
            ENDIF
C
C Formatted PLOT3D output
C
          ELSE
C Grid
            FERROR = FPLT3G
            OPEN (IPLT3G,FILE=FPLT3G,FORM='FORMATTED',ERR=990)
            WRITE (IPLT3G,'(I5)')  NBLKS
            IF ( THREED ) THEN
               WRITE (IPLT3G,'(3I5)') (IDIM(IBLK,ILVL)+1,
     1                                 JDIM(IBLK,ILVL)+1,
     2                                 KDIM(IBLK,ILVL)+1,IBLK=1,NBLKS)
            ELSE 
               WRITE (IPLT3G,'(3I5)') (IDIM(IBLK,ILVL)+1,
     1                                 JDIM(IBLK,ILVL)+1,
     2                                 IBLK=1,NBLKS)
            ENDIF
C Solution
            FERROR = FPLT3Q
            OPEN (IPLT3Q,FILE=FPLT3Q,FORM='FORMATTED',ERR=990)
            WRITE (IPLT3Q,'(I5)')  NBLKS
            IF ( THREED ) THEN
               WRITE (IPLT3Q,'(3I5)') (IDIM(IBLK,ILVL)+1,
     1                                 JDIM(IBLK,ILVL)+1, 
     2                                 KDIM(IBLK,ILVL)+1,IBLK=1,NBLKS)
            ELSE
               WRITE (IPLT3Q,'(3I5)') (IDIM(IBLK,ILVL)+1,
     1                                 JDIM(IBLK,ILVL)+1,
     2                                 IBLK=1,NBLKS)
            ENDIF
C Function file
            IF (NQ .GT. 5) THEN
               FERROR = FPLT3F
               NVAR   = (NQ - 5) + NP
               OPEN (IPLT3F,FILE=FPLT3F,FORM='FORMATTED',ERR=990)
               WRITE (IPLT3F,'(I5)') NBLKS
               IF ( THREED ) THEN
                  WRITE (IPLT3F,'(4I5)') (IDIM(IBLK,ILVL)+1,
     1                                    JDIM(IBLK,ILVL)+1,
     2                                    KDIM(IBLK,ILVL)+1,
     3                                    NVAR,IBLK=1,NBLKS)
               ELSE
                  WRITE (IPLT3F,'(4I5)') (IDIM(IBLK,ILVL)+1,
     1                                    JDIM(IBLK,ILVL)+1,
     2                                    NVAR,IBLK=1,NBLKS)
               ENDIF
            ENDIF
          ENDIF
      ENDIF
C
C Plot3d function name file
C
      IF (NQ .GT. 5) THEN
         FERROR = FPLT3FN
         OPEN (IPLT3FN,FILE=FPLT3FN,FORM='FORMATTED',ERR=990)
         CALL PLOT3DFN (ITURB)
         CLOSE (IPLT3FN)
      ENDIF
C
C Finished writing PLOT3D headers 
C
      DO 980 IBLK = 1, NBLKS
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
         CALL OUTPUT (IDIM(IBLK,ILVL),JDIM(IBLK,ILVL),
     1                KDIM(IBLK,ILVL),Q(IQ),PROPS(IP),PROPS(IP2),
     2                PROPS(IP3),ITURB,
     3                R(IR),IBLK,NPRINT(IBLK),IPRINT(1,1,IBLK))
         IF ( IFPLT3 ) THEN
            IF ( BINPLT ) THEN
               WRITE (IPLT3Q) FSMACH, ALPHA, RE, TIME
            ELSE
               WRITE (IPLT3Q,'(4E16.8)') FSMACH, ALPHA, RE, TIME
            ENDIF
            DO 959 ICUT = 1, NCUTS
               IF (ICUTS(1,ICUT,ILVL) .EQ. IBLK) THEN
                  CALL CUT3D (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL), 
     1                        KDIM(IBLK,ILVL), Q(IQ), NQ,
     2                        ICUTS(1,ICUT,ILVL))
                  IF (NP.GT.0) THEN
                     CALL CUT3D (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL), 
     1                           KDIM(IBLK,ILVL), PROPS(IP), NP,
     2                           ICUTS(1,ICUT,ILVL))
                  ENDIF
               ENDIF
  959       CONTINUE
            CALL COREDG (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL), 
     1                   KDIM(IBLK,ILVL), Q(IQ), NQ)
            IF (NP.GT.0) THEN
               CALL COREDG (IDIM(IBLK,ILVL), JDIM(IBLK,ILVL), 
     1                      KDIM(IBLK,ILVL), PROPS(IP), NP)
            ENDIF
            CALL PLOT3D (IDIM(IBLK,ILVL),JDIM(IBLK,ILVL),
     1                   KDIM(IBLK,ILVL),Q(IQ),PROPS(IP),R(IR),BINPLT)
         ENDIF
  980 CONTINUE

C     Success
      RETURN
C
C     Failure: STOP
C
  990 CONTINUE
      WRITE (IOUT,9990) FERROR
 9990 FORMAT (/' ','WRTP3D: ERROR-> Error opening file ',A80,
     1        /' ', 8X,'RUN ABORTING!')
      CLOSE(IOUT)
      STOP
      END
