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
      SUBROUTINE BC     (IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, NPPTS, S, VOL, RC,
     2                   NBCS, IBCDAT, NXPROF, NPROFL, IOFQBC, QBC, 
     3                   VISCOS, IMODEL, WORK, MXSECT, CL, RCNTR)
C
C Routine to set boundary conditions. This routine calls the individual
C boundary condition routines to set each type of boundary condition.
C The boundary condition routines extend the BCs from I=2 to I=1 and
C I=IDIM to I=IDIM+1 and similarly for J and K faces to set the corner
C data. Therefore, calling this routine with IPSTRT=2 is identical to
C using IPSTRT=1 and IPEND=IDIM is identical to IPEND=IDIM+1.
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
C NPPTS          : Number of points at which PROPS are stored.  This allows
C                  the storage of PROPS for viscous flows and the avoidance
C                  of the additional storage for inviscid flows.
C S              : Metrics
C VOL            : Cell volumes
C RC             : Position vector at cell centers
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
C NXPROF         : Maximum number of profile boundary condition segments
C NPROFL         : Number of profile boundary condition segments
C IOFQBC         : Offsets for QBC array
C QBC            : Boundary condition array - for PROFILE data
C VISCOS         : Logical for viscous (TRUE) or inviscid (FALSE) flow
C IMODEL         : Turbulence model array
C WORK           : Work array
C MXSECT         : Maximum size of the work array
C CL             : Lift coefficient. Used for 2D farfield vortex bc.
C RCNTR          : Center of airfoil. Used as location of point vortex.
C                  Also used as airfoil location at which to calculate moments.
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          PROPS (NPPTS,NP),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1),
     4          RC    (1:JDIM+1,1:KDIM+1,1:IDIM+1,3)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION IOFQBC(NXPROF)
      DIMENSION QBC   (*) 
C
      DIMENSION IMODEL(NMDL)
C
      DIMENSION WORK  (MXSECT,NQ)
C
      LOGICAL   VISCOS(3)
C
      DIMENSION RCNTR (3)
C
C Need to add Logic to avoid calling the IDIR = 1 (I) direction boundary
C     conditions when in the interior of the domain and using a planar
C     time integration scheme (e.g. LGS, MARCH).  Need to allow for the
C     call for other time integration schemes and when using a planar
C     scheme on the first plane for a conical solution procedure.
C
      ITURB = IMODEL(1)
C
      DO 900 IBC = 1, NBCS
         IBCTYP = IBCDAT(1,IBC)
C
C     Do not call K-direction boundary conditions for 2D flows
C
         IDIR   = IBCDAT(2,IBC)
         IF (.NOT.THREED .AND. .NOT.AXISYM .AND. IDIR .EQ. 3) GO TO 900
C
C     Wall boundary condition
C
         IF (IBCTYP .EQ. IWALL) THEN
            CALL BCWALL (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Wall Functions boundary condition
C
         ELSE IF (IBCTYP .EQ. IWALFN) THEN
            CALL BCWALF (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, VOL, IBCDAT(1,IBC))
C
C     Tangency boundary condition
C
         ELSE IF (IBCTYP .EQ. ITAN) THEN
            CALL BCTAN  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Farfield boundary condition
C
         ELSE IF (IBCTYP .EQ. IFAR) THEN
            CALL BCFAR  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     2D Farfield boundary condition with point vortex correction
C
         ELSE IF (IBCTYP .EQ. IFAR2D) THEN
            CALL BCFR2D (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, RC, IBCDAT(1,IBC), CL, RCNTR)
C
C     Extrapolation boundary condition
C     This is general for any number of equations
C
         ELSE IF (IBCTYP .EQ. IEXT) THEN
            CALL BCEXT  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Fix boundary condition
C
         ELSE IF (IBCTYP .EQ. IFIX) THEN
            CALL BCFIX  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Jet boundary condition
C
         ELSE IF (IBCTYP .EQ. IJET) THEN
            CALL BCJET  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Symmetry boundary condition
C
         ELSE IF (IBCTYP .EQ. ISYMXY .OR. IBCTYP .EQ. ISYMXZ
     1                               .OR. IBCTYP .EQ. ISYMYZ) THEN
            CALL BCSYM  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC), WORK)
C     
C     Periodic boundary condition
C
         ELSE IF (IBCTYP .EQ. IPERD) THEN
            CALL BCPERD (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Profile boundary condition
C
         ELSE IF (IBCTYP .EQ. IPROFL) THEN
            ISTRT  = IBCDAT( 4,IBC)
            IEND   = IBCDAT( 5,IBC)
            JSTRT  = IBCDAT( 6,IBC)
            JEND   = IBCDAT( 7,IBC)
            KSTRT  = IBCDAT( 8,IBC)
            KEND   = IBCDAT( 9,IBC)
            INPROF = IBCDAT(10,IBC)
            IQBC   = IOFQBC(INPROF) * NQ + 1
            NPTS   = (IEND-ISTRT+1) * (JEND-JSTRT+1) * (KEND-KSTRT+1)
            CALL BCPRFL (IDIM, JDIM, KDIM, IPSTRT, IPEND, NPTS,
     1                   Q, QBC(IQBC), S, IBCDAT(1,IBC), IERROR)
C
C     Subsonic Profile boundary condition
C     - call profile bc first and then update the pressure
C
         ELSE IF (IBCTYP .EQ. IPROSB) THEN
            ISTRT  = IBCDAT( 4,IBC)
            IEND   = IBCDAT( 5,IBC)
            JSTRT  = IBCDAT( 6,IBC)
            JEND   = IBCDAT( 7,IBC)
            KSTRT  = IBCDAT( 8,IBC)
            KEND   = IBCDAT( 9,IBC)
            INPROF = IBCDAT(10,IBC)
            IQBC   = IOFQBC(INPROF) * NQ + 1
            NPTS   = (IEND-ISTRT+1) * (JEND-JSTRT+1) * (KEND-KSTRT+1)
            CALL BCPRFL (IDIM, JDIM, KDIM, IPSTRT, IPEND, NPTS,
     1                   Q, QBC(IQBC), S, IBCDAT(1,IBC), IERROR)
            CALL BCPRSB (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, IBCDAT(1,IBC), IERROR)
C
C     Subsonic Inflow  boundary condition
C
         ELSE IF (IBCTYP .EQ. IINFLO) THEN
            CALL BCSIN  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Subsonic Outflow boundary condition
C
         ELSE IF (IBCTYP .EQ. IOUTFL) THEN
            CALL BCSOUT (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Axisymmetric boundary condition
C
         ELSE IF (IBCTYP .EQ. IAXISM) THEN
            CALL BCAXI  (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C
C     Axisymmetric centerline boundary condition
C
         ELSE IF (IBCTYP .EQ. IAXICL) THEN
            CALL BCAXCL (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, S, IBCDAT(1,IBC))
C     
C     Homogeneous Shear Periodic boundary condition
C
C N.B. Need to fix the WORK array here
C
         ELSE IF (IBCTYP .EQ. IHSHR) THEN
            CALL BCHSHR (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, 
     1                   Q, PROPS, IBCDAT(1,IBC), 
     2                   WORK(1,1), WORK(1,2))
C
C     Error
C
         ELSE
            WRITE (IOUT,1000) IBC, IBCTYP
         ENDIF
 900  CONTINUE
C
C     Finished with boundary condition routine
C
      RETURN
 1000 FORMAT (' ','BC    : ERROR-> An illegal boundary condition ',
     1            'type occurs in the ',I3,
     2            ' boundary condition and is = ',I3)
      END
