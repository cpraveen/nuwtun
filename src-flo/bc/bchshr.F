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
      SUBROUTINE BCHSHR (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, IBCDAT, QBC, PBC)
C
C Routine to do PERIODIC boundary condition for HOMOGENEOUS SHEAR in the
C J-direction ONLY.
C
C This routine is organized as follows:
C     1. Setup variables for looping
C     2. Set boundary conditions for all variables
C     3. Set boundary conditions for all properties
C
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start BCs at 
C                  (Minimum Value: 1)
C IPEND          : I Plane (cell center) to end BCs at
C                  (Maximum Value: IDIM+1)
C Q              : Primitive variables at cell centers
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
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
C QBC,PBC        : Temporary storage for Q, PROPS for interpolation
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP)
C
      DIMENSION IBCDAT(10)
C
C Temporary storage for Q and PROPS for interpolation
C
      DIMENSION QBC   (3-IDIM:2*IDIM-1,NQ),
     1          PBC   (3-IDIM:2*IDIM-1,NP)
C
C 1. Setup memory for looping.
C
C
C Store boundary condition data in mnemonic names
C
      IBCTYP   = IBCDAT(1)
      IDIR     = IBCDAT(2)
      INROUT   = IBCDAT(3)
C
C This routine is set-up for the J-Direction ONLY!
C
      IF (IDIR .NE. 2) THEN
         WRITE (IOUT,1000) IDIR
         STOP
      END IF
C
C      Calculate the increment in the I coordinate for interpolation
C
      DLENGTH = 1.E0
c--   DUDY    = 0.
      DUDY    = 0.6847688E-01
      DELTAU  = DUDY * DLENGTH
C
      DX      = DLENGTH / FLOAT(IDIM-1)
      DELTA   = DELTAU * TIMEDT
      DELTAI  = DELTA / DX
C
C      Check to see that DELTAI < IDIM; Stop if an error
C
      IF (INT(DELTAI) .GT. IDIM) THEN
         WRITE (IOUT,1005) INT(DELTAI), IDIM
         STOP
      ENDIF
C
C 2.   Boundary conditions for ALL Variables and ALL Properties
C
C      Loop over all K planes
C
      DO K = 1, KDIM+1
C
C 2-a. Set up the arrays to interpolate in
C
         DO JLOOP = 0, 1
C
C      Set needed variables depending on whether the boundary is
C      the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C
            IF (INROUT .EQ. 1) THEN
               DINC   = DELTAI
               INC    = INT(DINC)
               DU     = - DELTAU
               IF (JLOOP .EQ. 0) THEN
                  JBC = 0
                  JPR = JDIM - 1
               ELSE
                  JBC = 1
                  JPR = JDIM
               END IF
            ELSE
               DINC   = - DELTAI
               INC    = INT(DINC)
               DU     = DELTAU
               IF (JLOOP .EQ. 0) THEN
                  JBC = JDIM + 1
                  JPR = 2
               ELSE
                  JBC = JDIM + 2
                  JPR = 3
               END IF
            END IF
C
C      Calculate the initial point in the interpolating stencil as the third
C      grid index to the left of the interpolating point.
C
            IF (DINC .LT. 0) INC = INC - 1
C
C      Calculate the interpolation coefficients. Due to the design of the
C      homogeneous shear case, the interpolation is independent of the
C      location. The interpolation coefficients are then invariant.
C      5th order Lagrange interpolating polynomial. Carnahan, Luther, and
C      Wilkes, p. 27.
C
            DX0 = DINC - FLOAT(INC) + 2.E0
            DX1 = DX0    - 1.E0
            DX2 = DX0    - 2.E0
            DX3 = DX0    - 3.E0
            DX4 = DX0    - 4.E0
            DX5 = DX0    - 5.E0
C
            A0  = -       DX1 * DX2 * DX3 * DX4 * DX5 / 120.E0
            A1  =   DX0       * DX2 * DX3 * DX4 * DX5 /  24.E0
            A2  = - DX0 * DX1       * DX3 * DX4 * DX5 /  12.E0
            A3  =   DX0 * DX1 * DX2       * DX4 * DX5 /  12.E0
            A4  = - DX0 * DX1 * DX2 * DX3       * DX5 /  24.E0
            A5  =   DX0 * DX1 * DX2 * DX3 * DX4       / 120.E0
C
C      Store Q and PROPS in a periodic one-dimensional array to allow
C      for a simple form of the interpolation.
C
            DO I = 2, IDIM
               DO L = 1, NQ
                  QBC(I       ,L) = Q    (JPR,K,I,L)
                  QBC(I+IDIM-1,L) = Q    (JPR,K,I,L)
                  QBC(I-IDIM+1,L) = Q    (JPR,K,I,L)
               END DO
               DO L = 1, NP
                  PBC(I       ,L) = PROPS(JPR,K,I,L)
                  PBC(I+IDIM-1,L) = PROPS(JPR,K,I,L)
                  PBC(I-IDIM+1,L) = PROPS(JPR,K,I,L)
               END DO
            END DO
C
C 2-b. Interpolate
C      5th order Lagrange interpolating polynomial. Carnahan, Luther, and
C      Wilkes, p. 27. The II index is the cell center location of the
C      first point in the interpolating polynomial corresponding to the
C      A0 coefficient.
C
            DO I = 1, IDIM+1
               II = I + INC - 2
C
               DO L = 1, NQ
                  Q    (JBC,K,I,L) = A0 * QBC(II  ,L) + A1 * QBC(II+1,L)
     1                             + A2 * QBC(II+2,L) + A3 * QBC(II+3,L)
     2                             + A4 * QBC(II+4,L) + A5 * QBC(II+5,L)
               END DO
               DO L = 1, NP
                  PROPS(JBC,K,I,L) = A0 * PBC(II  ,L) + A1 * PBC(II+1,L)
     1                             + A2 * PBC(II+2,L) + A3 * PBC(II+3,L)
     2                             + A4 * PBC(II+4,L) + A5 * PBC(II+5,L)
               END DO
            END DO
C
C     Add/Subtract the mean velocity gradient from U
C
            DO I = 1, IDIM+1
               Q(JBC,K,I,2) = Q(JBC,K,I,2) + DU
            END DO
C
         END DO
      END DO
C
C Finished with PERIODIC HOMOGENEOUS SHEAR boundary condition
C
      RETURN
 1000 FORMAT (' ','ERROR-> BCHSHR: This routine is for the J-direction',
     1            ' only.',
     2       /' ',15x,' An invalid direction ',I2,' has been specified',
     3       /' ',15x,' Run Aborting.')
 1005 FORMAT (' ','ERROR-> BCHSHR: The calculated offset for the ',
     1                             'J-direction ',I5,
     2       /' ','                exceeds the allowed (IDIM) ',I5,
     3       /' ',15x,' Run Aborting.')
      END
