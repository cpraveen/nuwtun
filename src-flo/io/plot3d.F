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
      SUBROUTINE PLOT3D (IDIM, JDIM, KDIM, Q, PROPS, RC, BINARY)
C
C Routine to write PLOT3D output. The mean flow variables are output in
C the PLOT3D solution (Q) file in 'conserved' form (neglecting TKE in the
C equation of state) and the turbulence variables and properties (if NQ > 5) 
C are output in the function file. The turbulence quantities are output as 
C they are calculated and are not converted to 'conserved' form.
C
C This routine is organized as follows:
C     1. Convert mean flow from primitive to 'conserved' variables 
C        NOTE: do not include TKE in conversion from pressure to total energy
C              to match the equation of state that is assumed in FAST/PLOT3D
C     2. Write PLOT3D output files (FORMATTED or UNFORMATTED)
C        a. Three-dimensional output
C           i.  Binary    { }
C           ii. Formatted { }
C        b. Two  -dimensional output
C           i.  Binary    { }
C           ii. Formatted { }
C               { a. Grid file     (G)
C                 b. Solution file (Q)
C                 c. If (NQ > 5) write function file (F) }
C     3. Convert mean flow from 'conserved' to primitive variables
C
C IDIM,JDIM,KDIM : Dimensions of current block
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C RC             : Calculated position of grid at cell centers
C BINARY         : Logical for PLOT3D output file type
C                  = TRUE  -> binary (unformatted) file
C                    FALSE -> formatted file
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      LOGICAL BINARY
C 
C 1. Convert mean flow from primitive to 'conserved' variables 
C    NOTE: do not include TKE in conversion from pressure to total energy
C          to match the equation of state that is assumed in FAST/PLOT3D
C 
      DO 130 I = 1, IDIM+1 
         DO 120 K = 1, KDIM+1 
            DO 110 J = 1, JDIM+1 
               Q(J,K,I,5) = Q(J,K,I,5) / GAMM1 + 0.5E0*Q(J,K,I,1)*
     1                           ( Q(J,K,I,2)*Q(J,K,I,2)
     2                           + Q(J,K,I,3)*Q(J,K,I,3)
     3                           + Q(J,K,I,4)*Q(J,K,I,4) )
               Q(J,K,I,2) = Q(J,K,I,2) * Q(J,K,I,1)
               Q(J,K,I,3) = Q(J,K,I,3) * Q(J,K,I,1)
               Q(J,K,I,4) = Q(J,K,I,4) * Q(J,K,I,1)
  110       CONTINUE 
  120    CONTINUE 
  130 CONTINUE 
C 
C 2.   Write PLOT3D output files (FORMATTED or UNFORMATTED)
C
C 2-a. Output for three-dimensional case
C
      IF ( THREED ) THEN
         IF ( BINARY ) THEN
C
C 2-a-i.  Binary output
C
            WRITE (IPLT3G) ((((RC(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                     K=1,KDIM+1),L=1,3)
            WRITE (IPLT3Q) ((((Q (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                     K=1,KDIM+1),L=1,5)
            IF (NQ .GT. 5) THEN
               WRITE (IPLT3F) 
     1                     ((((Q    (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                        K=1,KDIM+1),L=6,NQ),
     2                     ((((PROPS(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     3                                        K=1,KDIM+1),L=1,NP)
            ENDIF
         ELSE
C
C 2-a-ii. Formatted output
C
            WRITE (IPLT3G,1000) ((((RC(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                          K=1,KDIM+1),L=1,3)
            WRITE (IPLT3Q,1000) ((((Q (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                          K=1,KDIM+1),L=1,5)
            IF (NQ .GT. 5) THEN
               WRITE (IPLT3F,1000) 
     1                       ((((Q    (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     2                                          K=1,KDIM+1),L=6,NQ),
     3                       ((((PROPS(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     4                                          K=1,KDIM+1),L=1,NP)
            ENDIF
         ENDIF
C
C 2-b. Output for two-dimensional case
C      Output as a three-dimensional case with KDIM=1 for FAST
C
      ELSE
         K = 2
         IF ( BINARY ) THEN
C
C 2-b-i.  Binary output
C
            WRITE (IPLT3G) (((RC(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                    L=1,2)
            WRITE (IPLT3Q) (((Q (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                    L=1,3),
     2                      ((Q (J,K,I,5),I=1,IDIM+1),J=1,JDIM+1)
            IF (NQ .GT. 5) THEN
               WRITE (IPLT3F) 
     1                     (((Q    (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     2                                       L=6,NQ),
     3                     (((PROPS(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     4                                       L=1,NP)
            ENDIF
         ELSE
C
C 2-b-ii. Formatted output
C
            WRITE (IPLT3G,1000) (((RC(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                         L=1,2)
            WRITE (IPLT3Q,1000) (((Q (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     1                                         L=1,3),
     2                           ((Q (J,K,I,5),I=1,IDIM+1),J=1,JDIM+1)
            IF (NQ .GT. 5) THEN
               WRITE (IPLT3F,1000) 
     1                       (((Q    (J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     2                                         L=6,NQ),
     3                       (((PROPS(J,K,I,L),I=1,IDIM+1),J=1,JDIM+1),
     4                                         L=1,NP)
            ENDIF
         ENDIF
      ENDIF
C 
C 3.  Convert mean flow from 'conserved' to primitive variables
C 
      DO 330 I = 1, IDIM+1 
         DO 320 K = 1, KDIM+1 
            DO 310 J = 1, JDIM+1 
               Q(J,K,I,2) = Q(J,K,I,2) / Q(J,K,I,1)
               Q(J,K,I,3) = Q(J,K,I,3) / Q(J,K,I,1)
               Q(J,K,I,4) = Q(J,K,I,4) / Q(J,K,I,1)
               Q(J,K,I,5) = GAMM1*( Q(J,K,I,5) - 0.5E0*Q(J,K,I,1)*
     1                               ( Q(J,K,I,2)*Q(J,K,I,2)
     2                               + Q(J,K,I,3)*Q(J,K,I,3)
     3                               + Q(J,K,I,4)*Q(J,K,I,4) ) )
  310       CONTINUE 
  320    CONTINUE 
  330 CONTINUE 
C
C Finished with PLOT3D output
C
      RETURN
 1000 FORMAT (5E16.8)
      END
