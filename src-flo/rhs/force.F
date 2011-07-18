      SUBROUTINE FRCINT (ITURB, NBLKS, IDIM, JDIM, KDIM, 
     1                   IOFF,  IOFFS, IOFFQ, IOFFP,
     2                   Q, PROPS, S, RC, NXBCS, NBCS, IBCDAT, 
     3                   FORCE, CL, CD, FMOMENT)
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
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
C                  PROPS(4) = YNORML
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
     1          PROPS (*),
     2          S     (*),
     3          RC    (*)
C
      DIMENSION NBCS  (NBLKS),
     1          IBCDAT(10,NXBCS,NBLKS)
C
      DIMENSION FORCE (3)
      DIMENSION FMOMENT (*)
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
C 0.  Initialize forces to zero
C
      FORCE(1) = 0.E0
      FORCE(2) = 0.E0
      FORCE(3) = 0.E0
C
C 0.  Initialize moments to zero
C
      DO IMOM=1,NMOMENT
         FMOMENT(IMOM) = 0.E0
      ENDDO
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
            IF ( (.NOT. VISCOS .AND. IBCTYP .EQ. ITAN)  .OR. 
     1           (VISCOS .AND. 
     2               (IBCTYP .EQ. IWALL .OR. IBCTYP .EQ. IWALFN)) ) THEN
               CALL FRCEL (IDIM(IBLK), JDIM(IBLK), KDIM(IBLK), 
     1                     Q(IQ), PROPS(IP), S(IS), RC(IR), 
     2                     IBCDAT(1,IBC,IBLK), VISCOS, FORCE, FMOMENT)
            ENDIF
  150    CONTINUE
  190 CONTINUE
C
C 2.  Non-dimensionalize forces to force coefficients
C     C_f = F / (1/2 rho_inf U_inf^2 Area_ref) where Area_ref is the 
C                                              reference area of the body
C
      CONST    = 2.E0 / (FSMACH * FSMACH * AREARF)
      FORCE(1) = FORCE(1) * CONST
      FORCE(2) = FORCE(2) * CONST
      FORCE(3) = FORCE(3) * CONST
C
C     Normalize moment
C
      CONST = CONST / ALENRF
      DO IMOM=1,NMOMENT
         FMOMENT(IMOM) = FMOMENT(IMOM) * CONST
      ENDDO
C
C     C_L and C_D: Defined relative to the freestream velocity direction
C
      PI     = ACOS(-1.E0)
      ALPHAR = ALPHA * PI / 180.E0
      IF (THREED) THEN
         CL  = FORCE(3) * COS(ALPHAR) - FORCE(1) * SIN(ALPHAR)
         CD  = FORCE(1) * COS(ALPHAR) + FORCE(3) * SIN(ALPHAR)
      ELSE
         CL  = FORCE(2) * COS(ALPHAR) - FORCE(1) * SIN(ALPHAR)
         CD  = FORCE(1) * COS(ALPHAR) + FORCE(2) * SIN(ALPHAR)
      ENDIF
C
C Finished calculating forces
C
      RETURN
      END





      SUBROUTINE FRCEL (IDIM, JDIM, KDIM, Q, PROPS, S, RC, 
     1                  IBCDAT, VISCOS, FORCE, FMOMENT)
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
C PROPS          : Properties at cell centers
C                  PROPS(1) = molecular viscosity
C                  PROPS(2) = turbulent eddy viscosity
C                  PROPS(3) = Y+
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
     1          PROPS (0:JDIM+2,0:KDIM+2,0:IDIM+2,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     3          RC    (JDIM+1,KDIM+1,IDIM+1,3)
C
      DIMENSION IBCDAT(10)
C
      DIMENSION FORCE (3)
C
      DIMENSION FMOMENT (*)
C
      DIMENSION FMOM(3), RMOM(3)
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
            IM   = IBCM  *IKD(IC1,1) + I2*IKD(IC2,1) + I3*IKD(IC3,1)
            JM   = IBCM  *IKD(IC1,2) + I2*IKD(IC2,2) + I3*IKD(IC3,2)
            KM   = IBCM  *IKD(IC1,3) + I2*IKD(IC2,3) + I3*IKD(IC3,3)
C
            PW   = Q(JBC,KBC,IBC,5) 
C
C Calculate the pressure force 
C Minus sign since pressure acts opposite direction of face normal
C
            FPRESS   = PW * S(JM,KM,IM,4,IDIR) * SGNM
            FX       = FPRESS * S(JM,KM,IM,1,IDIR)
            FY       = FPRESS * S(JM,KM,IM,2,IDIR)
            FZ       = FPRESS * S(JM,KM,IM,3,IDIR)
            FORCE(1) = FORCE(1) - FX
            FORCE(2) = FORCE(2) - FY
            FORCE(3) = FORCE(3) - FZ

            DO IMOM=1,NMOMENT
C Moment arm
               RMOM(1) = RC(JM,KM,IM,1) - AMOMAXIS(IMOM,1,1)
               RMOM(2) = RC(JM,KM,IM,2) - AMOMAXIS(IMOM,1,2)
               RMOM(3) = RC(JM,KM,IM,3) - AMOMAXIS(IMOM,1,3)
C Moment arm x force
C Negative sign required since we need force on wing
C See also computation of FORCE(*) above
               FMOM(1) = -(RMOM(2) * FZ - RMOM(3) * FY)
               FMOM(2) = -(RMOM(3) * FX - RMOM(1) * FZ)
               FMOM(3) = -(RMOM(1) * FY - RMOM(2) * FX)
C Dot product with moment axis
               FMOMENT(IMOM) = FMOMENT(IMOM) + 
     1                         FMOM(1) * AMOMAXIS(IMOM,3,1) +
     2                         FMOM(2) * AMOMAXIS(IMOM,3,2) +
     3                         FMOM(3) * AMOMAXIS(IMOM,3,3)
            ENDDO
  210    CONTINUE
  220 CONTINUE
C
C 3.    Calculate Viscous  force increments
C
      IF (VISCOS) THEN
C
C 3-i.  Wall contribution
C
         IF (IBCTYP .EQ. IWALL) THEN
            DO 320 I2 = ISTRT(IC2), IEND(IC2)
               DO 310 I3 = ISTRT(IC3), IEND(IC3)
C
                  IBC  = IBCPT *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JBC  = IBCPT *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KBC  = IBCPT *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  IBC2 = IBCPT2*IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JBC2 = IBCPT2*IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KBC2 = IBCPT2*IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  IN   = IBCN  *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JN   = IBCN  *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KN   = IBCN  *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  IN2  = IBCN2 *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JN2  = IBCN2 *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KN2  = IBCN2 *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  IM   = IBCM  *IKD(IC1,1)+I2*IKD(IC2,1)+I3*IKD(IC3,1)
                  JM   = IBCM  *IKD(IC1,2)+I2*IKD(IC2,2)+I3*IKD(IC3,2)
                  KM   = IBCM  *IKD(IC1,3)+I2*IKD(IC2,3)+I3*IKD(IC3,3)
C
                  RHOW = Q(JBC,KBC,IBC,1) 
                  PW   = Q(JBC,KBC,IBC,5) 
C
                  UN   = Q(JN ,KN ,IN ,2) 
                  UN2  = Q(JN2,KN2,IN2,2) 
                  VN   = Q(JN ,KN ,IN ,3) 
                  VN2  = Q(JN2,KN2,IN2,3) 
                  WN   = Q(JN ,KN ,IN ,4) 
                  WN2  = Q(JN2,KN2,IN2,4) 
C
C Viscosity (RMUW) at the wall 
C
                  RMUW   = PROPS(JBC,KBC,IBC,1)
C
C Distances
C
                  DN     = SQRT (
     1                      (RC(JN ,KN ,IN ,1) - RC(JBC,KBC,IBC,1))**2
     2                    + (RC(JN ,KN ,IN ,2) - RC(JBC,KBC,IBC,2))**2
     3                    + (RC(JN ,KN ,IN ,3) - RC(JBC,KBC,IBC,3))**2 )
C
                  DN2    = SQRT (
     1                      (RC(JN2,KN2,IN2,1) - RC(JN ,KN ,IN ,1))**2
     2                    + (RC(JN2,KN2,IN2,2) - RC(JN ,KN ,IN ,2))**2
     3                    + (RC(JN2,KN2,IN2,3) - RC(JN ,KN ,IN ,3))**2 )
C
C Calculate first derivative terms using first or second order differencing
C
                  DUDN   = (1-IHIORD) * (UN / DN)
     1                     +  IHIORD  * (-DN*DN*UN2 + (DN+DN2)**2*UN)
     2                                  / (DN*DN2*(DN+DN2))
C
                  DVDN   = (1-IHIORD) * (VN / DN)
     1                     +  IHIORD  * (-DN*DN*VN2 + (DN+DN2)**2*VN)
     2                                  / (DN*DN2*(DN+DN2))
C
                  DWDN   = (1-IHIORD) * (WN / DN)
     1                     +  IHIORD  * (-DN*DN*WN2 + (DN+DN2)**2*WN)
     2                                  / (DN*DN2*(DN+DN2))
C
                  DUDX   = S(JM,KM,IM,1,IDIR) * DUDN * SGNM
                  DUDY   = S(JM,KM,IM,2,IDIR) * DUDN * SGNM
                  DUDZ   = S(JM,KM,IM,3,IDIR) * DUDN * SGNM
C
                  DVDX   = S(JM,KM,IM,1,IDIR) * DVDN * SGNM
                  DVDY   = S(JM,KM,IM,2,IDIR) * DVDN * SGNM
                  DVDZ   = S(JM,KM,IM,3,IDIR) * DVDN * SGNM
C
                  DWDX   = S(JM,KM,IM,1,IDIR) * DWDN * SGNM
                  DWDY   = S(JM,KM,IM,2,IDIR) * DWDN * SGNM
                  DWDZ   = S(JM,KM,IM,3,IDIR) * DWDN * SGNM
C
                  SXX    = T23 * ( 2.E0*DUDX - DVDY - DWDZ )
                  SYY    = T23 * ( 2.E0*DVDY - DUDX - DWDZ )
                  SZZ    = T23 * ( 2.E0*DWDZ - DUDX - DVDY )
                  SXY    = ( DUDY + DVDX )
                  SXZ    = ( DWDX + DUDZ )
                  SYZ    = ( DWDY + DVDZ )
C
                  VISC   = RMUW * FSMACH / RE
                  SIGXX  = VISC * SXX
                  SIGYY  = VISC * SYY
                  SIGZZ  = VISC * SZZ
                  SIGXY  = VISC * SXY
                  SIGXZ  = VISC * SXZ
                  SIGYZ  = VISC * SYZ
C
C Calculate the Shear force
C
                  CNST     = S(JM,KM,IM,4,IDIR) * SGNM
                  FORCE(1) = FORCE(1) + CNST*(SIGXX*S(JM,KM,IM,1,IDIR)
     1              +SIGXY*S(JM,KM,IM,2,IDIR)+SIGXZ*S(JM,KM,IM,3,IDIR))
                  FORCE(2) = FORCE(2) + CNST*(SIGXY*S(JM,KM,IM,1,IDIR)
     1              +SIGYY*S(JM,KM,IM,2,IDIR)+SIGYZ*S(JM,KM,IM,3,IDIR))
                  FORCE(3) = FORCE(3) + CNST*(SIGXZ*S(JM,KM,IM,1,IDIR)
     1              +SIGYZ*S(JM,KM,IM,2,IDIR)+SIGZZ*S(JM,KM,IM,3,IDIR))
  310          CONTINUE
  320       CONTINUE
C
C 3-ii. Wall function contribution
C
         ELSE IF (IBCTYP .EQ. IWALFN) THEN
            write (iout,'(a)') 
     1            ' Wall functions not yet supported in force routine'
         ENDIF
      ENDIF
C
C Finished with force increment
C
      RETURN
      END
