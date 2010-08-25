      SUBROUTINE APGRAD(ITURB, NBLKS, IDIM, JDIM, KDIM, 
     1                   IOFF,  IOFFS, IOFFQ, IOFFP,
     2                   Q, S, RC, NXBCS, NBCS, IBCDAT, 
     3                   APGINT)
C
C Routine to loop over all blocks and all surfaces integrating the forces
C on the surface for the current grid level.  
C If INVISCID case then integrate forces on TANGENCY boundaries.
C If VISCOUS  case then integrate forces on WALL or WALLFN boundaries.
C
C Programming note: The arrays (Q, PROPS, S, RC) are passed in as the
C                   complete array since the offset pointers are calculated
C                   from the very beginning of the array.
C
C This routine is organized as follows:
C     0. Initialize forces to zero
C     1. Calculate forces integrated over all blocks
C        For each block {
C           For each surface boundary condition {
C              Call routine to add incremental forces
C           }
C        }
C     2. Non-dimensionalize forces to force coefficients
C
C ITURB          : Turbulence model
C NBLKS          : Number of blocks
C IDIM,JDIM,KDIM : Dimensions of all blocks on finest level
C IOFF           : Offset pointer
C IOFFS          : Offset pointer for S     array
C IOFFQ          : Offset pointer for Q     array
C IOFFP          : Offset pointer for PROPS array
C Q              : Primitive variables at cell centers
C S              : Metrics
C RC             : Position vector to cell center
C NXBCS          : Maximum number of bcs IBCDAT is dimensioned for
C NBCS           : Number of boundary conditions per block
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
C FORCE          : Force coefficients (1 = x, 2 = y, 3 = z)
C CL,CD          : Lift, Drag coefficients
C
      include '../header/common.h'
C
      DIMENSION IDIM  (NBLKS),
     1          JDIM  (NBLKS),
     2          KDIM  (NBLKS),
     3          IOFF  (NBLKS),
     4          IOFFS (NBLKS),
     5          IOFFQ (NBLKS),
     6          IOFFP (NBLKS)
C
      DIMENSION Q     (*),
     2          S     (*),
     3          RC    (*)
C
      DIMENSION NBCS  (NBLKS),
     1          IBCDAT(10,NXBCS,NBLKS)
C
      LOGICAL   VISCOS
C
C Set viscos to false if doing an inviscid case, else true
C
      IF (ITURB .EQ. ITNVSD) THEN
         VISCOS = .FALSE. 
      ELSE
         VISCOS = .TRUE. 
      ENDIF
C
C 0.  Initialize adverse pressure gradient integral to zero
C
      APGINT = 0.0E0
C
C 1.  Calculate forces integrated over all blocks
C
      DO 190 IBLK = 1, NBLKS
         IQ    = IOFFQ(IBLK) * NQ    + 1
         IP    = IOFFP(IBLK) * NP    + 1
         IS    = IOFFS(IBLK) * 4 * 3 + 1
         IR    = IOFFS(IBLK) * 3     + 1
         DO 150 IBC = 1, NBCS(IBLK)
            IBCTYP   = IBCDAT(1,IBC,IBLK)
            IF (.NOT. VISCOS .AND. IBCTYP .EQ. ITAN) THEN
               CALL APGRADBLK (IDIM(IBLK), JDIM(IBLK), KDIM(IBLK), 
     1                         Q(IQ), S(IS), RC(IR), 
     2                         IBCDAT(1,IBC,IBLK), APGINT)
            ENDIF
  150    CONTINUE
  190 CONTINUE
C
C Finished calculating adverse pressure gradient
C
      RETURN
      END





      SUBROUTINE APGRADBLK (IDIM, JDIM, KDIM, Q, S, RC, 
     1                      IBCDAT, APGINT)
C
C Routine to calculate the increment to the forces for the current surface
C element.
C
C This routine is organized as follows:
C     1. Setup memory for looping.
C     2.  Calculate Pressure force increments
C     3.  Calculate Viscous  force increments
C         i.  Wall contribution
C         ii. Wall function contribution
C
C IDIM,JDIM,KDIM : Dimensions of current block
C Q              : Primitive variables at cell centers
C S              : Metrics
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
C RC             : Cell center grid
C VISCOS         : Logical determining whether or not to include viscous forces
C FORCE          : Forces (1 = x, 2 = y, 3 = z)
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     3          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      DIMENSION IBCDAT(10)
C
      LOGICAL   VISCOS
C
C Local storage used to calculate indices
C
      DIMENSION IJKDIM(3),
     1          ISTRT (3),  
     2          IEND  (3)
C
C Two-Thirds
C
      T23    = 2.E0 / 3.E0
C
C Set control for higher order evaluation of du/dn and dT/dn
C     IHIORD = 0 for first  order
C     IHIORD = 1 for second order
C
      IHIORD = 1
C
C 1.  Setup memory for looping.
C
C Store IDIM, JDIM, KDIM in IJKDIM for boundary condition locations
C
      IJKDIM(1) = IDIM
      IJKDIM(2) = JDIM
      IJKDIM(3) = KDIM
C
C Store boundary condition data in mnemonic names
C
      IBCTYP   = IBCDAT(1)
      IDIR     = IBCDAT(2)
      INROUT   = IBCDAT(3)
      ISTRT(1) = IBCDAT(4)
      IEND (1) = IBCDAT(5)
      ISTRT(2) = IBCDAT(6)
      IEND (2) = IBCDAT(7)
      ISTRT(3) = IBCDAT(8)
      IEND (3) = IBCDAT(9)
C
C Set needed variables depending on whether the boundary is
C the inner boundary (INROUT = 1) or the outer boundary (INROUT > 1)
C      IBCPT  is the boundary condition location
C      IBCPT2 is the boundary condition location outside the block from IBCPT
C      IBCN   is the point to the inside of the block from IBCPT
C      IBCM   is the location of the metrics
C
C
      IF (INROUT .EQ. 1) THEN
         IBCPT  = 1
         IBCPT2 = 0
         IBCN   = 2
         IBCN2  = 3
         IBCM   = 1
         SGNM   = 1.E0
      ELSE
         IBCPT  = IJKDIM(IDIR) + 1
         IBCPT2 = IJKDIM(IDIR) + 2
         IBCN   = IJKDIM(IDIR)
         IBCN2  = IJKDIM(IDIR) - 1
         IBCM   = IJKDIM(IDIR)
         SGNM   = - 1.E0
      ENDIF
C
      IC1 = ICY (IDIR, 1)
      IC2 = ICY (IDIR, 2)
      IC3 = ICY (IDIR, 3)
C
C 2.  Calculate Pressure force increments
C
      DO 220 I2 = ISTRT(IC2), IEND(IC2)
         DO 210 I3 = ISTRT(IC3), IEND(IC3)
C
            IBC  = IBCPT *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JBC  = IBCPT *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KBC  = IBCPT *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            IN   = IBCN  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JN   = IBCN  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KN   = IBCN  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            PW   = Q(JBC,KBC,IBC,5) 
            PN   = Q(JN ,KN ,IN ,5)
C
            IF(IDIR.EQ.1)THEN
               PET = (Q(JN+1,KN,IN,5) - Q(JN-1,KN,IN,5))/2.0
               PZE = (Q(JN,KN+1,IN,5) - Q(JN,KN-1,IN,5))/2.0
               PXI = (PN - PW)/0.5
            ELSE IF(IDIR.EQ.2)THEN
               PXI = (Q(JN,KN,IN+1,5) - Q(JN,KN,IN-1,5))/2.0
               PZE = (Q(JN,KN+1,IN,5) - Q(JN,KN-1,IN,5))/2.0
               PET = (PN - PW)/0.5
            ELSE
               PXI = (Q(JN,KN,IN+1,5) - Q(JN,KN,IN-1,5))/2.0
               PET = (Q(JN+1,KN,IN,5) - Q(JN-1,KN,IN,5))/2.0
               PZE = (PN - PW)/0.5
            ENDIF
C
C             PX = PXI*S(JM,KM,IM,1,1)*S(JM,KM,IM,4,1) + 
C      1           PET*S(JM,KM,IM,1,2)*S(JM,KM,IM,4,2) + 
C      2           PZE*S(JM,KM,IM,1,3)*S(JM,KM,IM,4,3)
C             PY = PXI*S(JM,KM,IM,2,1)*S(JM,KM,IM,4,1) + 
C      1           PET*S(JM,KM,IM,2,2)*S(JM,KM,IM,4,2) + 
C      2           PZE*S(JM,KM,IM,2,3)*S(JM,KM,IM,4,3)
C             PZ = PXI*S(JM,KM,IM,3,1)*S(JM,KM,IM,4,1) + 
C      1           PET*S(JM,KM,IM,3,2)*S(JM,KM,IM,4,2) + 
C      2           PZE*S(JM,KM,IM,3,3)*S(JM,KM,IM,4,3)
C            
            PX = PXI*S(JN,KN,IN,1,1)*S(JN,KN,IN,4,1) + 
     1           PET*S(JN,KN,IN,1,2)*S(JN,KN,IN,4,2) + 
     2           PZE*S(JN,KN,IN,1,3)*S(JN,KN,IN,4,3)
            PY = PXI*S(JN,KN,IN,2,1)*S(JN,KN,IN,4,1) + 
     1           PET*S(JN,KN,IN,2,2)*S(JN,KN,IN,4,2) + 
     2           PZE*S(JN,KN,IN,2,3)*S(JN,KN,IN,4,3)
            PZ = PXI*S(JN,KN,IN,3,1)*S(JN,KN,IN,4,1) + 
     1           PET*S(JN,KN,IN,3,2)*S(JN,KN,IN,4,2) + 
     2           PZE*S(JN,KN,IN,3,3)*S(JN,KN,IN,4,3)
            R  = Q(JBC,KBC,IBC,1)
            U  = Q(JBC,KBC,IBC,2)
            V  = Q(JBC,KBC,IBC,3) 
            W  = Q(JBC,KBC,IBC,4) 
            SP = SQRT(U**2 + V**2 + W**2) + 1.0E-8
            P1 = (PX*U + PY*V + PZ*W)/SP
            P2 = MAX(0.0E0, P1)/SP**2/R
            P3 = P2 * S(JM,KM,IM,4,IDIR)
C             PRINT*
C             PRINT*,'INJNKN ',IN,JN,KN
C             PRINT*,'IMJMKM ',IM,JM,KM
C             PRINT*,'gradp1 ',PXI,PET,PZE
C             PRINT*,'gradp2 ',PX,PY,PZ
C             PRINT*,'veloc ',U,V,W
C             PRINT*,'prest ',P1,P2,P3
C             PRINT*
C
            APGINT   = APGINT + P3 
C
  210    CONTINUE
  220 CONTINUE
C
C Finished with force increment
C
      RETURN
      END
