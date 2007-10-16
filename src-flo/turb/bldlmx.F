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






      SUBROUTINE BLDLMX (JDIM, RHOW, TAUW, RMUW, UMAX, UMIN, Y, RHO, 
     1                   VORT, F, RMUT)
C
C     Routine to calculate the turbulent eddy viscosity using
C     the two-layer algebraic Baldwin-Lomax model of AIAA 78-257.
C
C
C     JDIM     : Dimension
C     RHOW     : Density at the wall
C     TAUW     : Wall shear stress
C     RMUW     : Molecular viscosity at the wall
C     UMAX     : The maximum value of the velocity magnitude in the profile
C     UMIN     : The minimum value of the velocity magnitude in the profile
C     Y        : Distance from the wall
C     RHO      : Density
C     VORT     : Vorticity magnitude
C     F        : Function for the outer layer
C     RMUT     : The turbulent eddy viscosity calculated by this routine
C
C     The above variables are in non-dimensional form as follows:
C        density   non-dimensionalized by freestream density
C        velocity  non-dimensionalized by freestream speed of sound
C        viscosity non-dimensionalized by freestream molecular viscosity
C        distance  non-dimensionalized by characteristic length L
C
      include '../header/common.h'
C
      DIMENSION Y     (JDIM),
     1          RHO   (JDIM),
     2          VORT  (JDIM),
     3          F     (JDIM),
     4          RMUT  (JDIM)
C
C     Baldwin-Lomax constants
C
      RK     = 0.4E0
      APLUS  = 26.E0
      CLAUSR = 0.0168E0
      CCP    = 1.6E0
      CWK    = 0.25E0
      CKLEB  = 0.3E0
C
      YPWALL = SQRT(RE/FSMACH*RHOW*TAUW) / RMUW
C
C     Loop over J defining the inner turbulent viscosity and
C     determining the function F(y).  
C     In the definition of F(y), the exponential term should be
C     set to zero for a wake.  This is not presently done.
C
      DO 100 J = 1, JDIM
         YPLUS   = YPWALL * Y(J)
         ETERM   = Y(J) * (1.E0 - EXP(-YPLUS/APLUS))
         RLNGTH  = RK * ETERM
         F(J)    = ETERM * VORT(J)
         RMUT(J) = RE / FSMACH * RHO(J) * RLNGTH * RLNGTH * VORT(J)
  100 CONTINUE
C
C     Determine the maximum of the function F(y) and the associated y
C
      JSTRT = 1
      JEND  = JDIM
      FMAX  = 0.E0
      DO 200 J = JSTRT, JEND
         IF (F(J) .GT. FMAX) THEN
            FMAX = F(J)
            YMAX = Y(J)
         ENDIF
  200 CONTINUE
C
C     Form the outer layer eddy viscosity values and store them in
C     the F array.
C
      UDIFF = UMAX - UMIN
      FWAKE = MIN( (YMAX*FMAX), (CWK*YMAX*UDIFF*UDIFF/FMAX) )
      DO 300 J = 1, JDIM
         FKLEB = 1.E0 / (1.E0 + 5.5E0*(CKLEB*Y(J)/YMAX)**6)
         F(J) = RE / FSMACH * CLAUSR * CCP * RHO(J) * FWAKE * FKLEB
  300 CONTINUE
C
C     Determine where the switch from the inner to the outer model occurs
C
      JSWITCH = JDIM
      DO 400 J = 1, JDIM
         IF (RMUT(J) .GE. F(J)) THEN
            JSWITCH = J
            GO TO 410
         ENDIF
  400 CONTINUE
  410 CONTINUE
C
C     Store the outer model value of eddy viscosity into the eddy
C     viscosity array for J >= JSWITCH (the crossover point)
C
      DO 500 J = JSWITCH, JDIM
         RMUT(J) = F(J)
  500 CONTINUE
C
C     Finished calculating the algebraic eddy viscosity
C
      RETURN
      END






      SUBROUTINE TURB (IDIM, JDIM, KDIM, MAXDIM, I, ITURB,
     1                 Q, RMUT, S, VOL, IBCDAT, NBCS,
     2                 VORT, TVIS, DIST,
     3                 RHO, VORTJ, RMUTJ, Y, FUN)
C
C     Routine to calculate the eddy viscosity coefficient.
C
C
C     IDIM,JDIM,KDIM : Dimensions of current block
C     MAXDIM         : Maximum dimension of the J-K plane
C     I              : Current I Plane
C     ITURB          : Turbulence model
C     Q              : Primitive variables at cell centers
C     RMUT           : Eddy Viscosity at cell centers
C     S              : Metrics
C     VOL            : Cell volumes
C     IBCDAT         : Data controlling boundary conditions
C     NBCS           : Number of boundary conditions for current block
C     VORT           : Temporary storage of the magnitude of the vorticity
C     TVIS           : Storage for the turbulent eddy viscosity in J-K plane.
C                      This will store the viscosity for up to 4 walls.
C     DIST           : Storage for the normal distance to wall in J-K plane.
C                      This will store the distance for up to 4 walls.
C     RHO            : Density in the profile
C     VORTJ          : Vorticity in the profile
C     RMUTJ          : Turbulent eddy viscosity in the profile
C     Y              : Normal distance from wall in the profile
C     FUN            : Used in BLDLMX to calculate F(y) and store mu(outer)
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          RMUT  (0:JDIM+2,0:KDIM+2,0:IDIM+2),
     2          S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3), 
     3          VOL   (0:JDIM+1,0:KDIM+1,0:IDIM+1)
C
      DIMENSION VORT  (0:JDIM+1,0:KDIM+1), 
     1          TVIS  (0:JDIM+1,0:KDIM+1,4), 
     2          DIST  (0:JDIM+1,0:KDIM+1,4)
C
C     Profile variables
C
      DIMENSION RHO   (MAXDIM+1),
     1          VORTJ (MAXDIM+1),
     2          RMUTJ (MAXDIM+1),
     3          Y     (MAXDIM+1),
     4          FUN   (MAXDIM+1)
C
      DIMENSION IBCDAT(10,NBCS)
C
C     Function for Sutherlands Law
C
      VIS(TT,VNUM,VTINF) = TT * SQRT( TT ) * VNUM / ( TT + VTINF )
C
C     Constants for Sutherlands Law
C
      SREF  = 198.6E0
      IF ( SIUNIT ) SREF = SREF / 1.8E0
      STINF = SREF / TINF
      SNUM  = 1.E0 + STINF
C
C     Initialize J1W, JDIMW, K1W, KDIMW
C     These variables tell if there is a wall on the J = 1, J = JDIM
C     K = 1 or K = KDIM boundaries.  This data tells the turbulent
C     scheme to calculate turbulence eddy viscosities based on all walls
C     present but must calculate the viscosity over the entire plane
C     rather than just of the wall.
C
      J1W   = 0
      JDIMW = 0
      K1W   = 0
      KDIMW = 0
C
C     Calculate the vorticity at the cell centers
C
C     The vorticity calculated at the wall will not be fixed to be exactly
C     correct.  This is not an issue as the wall value of eddy viscosity
C     is identically zero.
C
C     Store the values of U, V, W used to build up the derivative terms
C     in the DIST array as it is not used until after VORT is defined.
C     Store the derivative terms used to build up the vorticity in the
C     TVIS array as it is not used until after VORT is defined.
C
C     Form the I contribution to the magnitude of vorticity
C
      DO 120 K = 1, KDIM + 1
         DO 110 J = 1, JDIM + 1
C
C     Calculate first derivative terms using Divergence Theorem
C**** These terms will have errors in the I + I+1 terms for the
C**** MARCHing algorithm
C
            DUDY    = 0.5E0*(Q(J,K,I+1,2)+Q(J,K,I,2)) * S(J,K,I,2,1)   
     1                                                * S(J,K,I,4,1) 
     2              - 0.5E0*(Q(J,K,I,2)+Q(J,K,I-1,2)) * S(J,K,I-1,2,1)
     3                                                * S(J,K,I-1,4,1)
            DUDZ    = 0.5E0*(Q(J,K,I+1,2)+Q(J,K,I,2)) * S(J,K,I,3,1)   
     1                                                * S(J,K,I,4,1) 
     2              - 0.5E0*(Q(J,K,I,2)+Q(J,K,I-1,2)) * S(J,K,I-1,3,1) 
     3                                                * S(J,K,I-1,4,1)
            DVDX    = 0.5E0*(Q(J,K,I+1,3)+Q(J,K,I,3)) * S(J,K,I,1,1)   
     1                                                * S(J,K,I,4,1) 
     2              - 0.5E0*(Q(J,K,I,3)+Q(J,K,I-1,3)) * S(J,K,I-1,1,1) 
     3                                                * S(J,K,I-1,4,1)
            DVDZ    = 0.5E0*(Q(J,K,I+1,3)+Q(J,K,I,3)) * S(J,K,I,3,1)   
     1                                                * S(J,K,I,4,1) 
     2              - 0.5E0*(Q(J,K,I,3)+Q(J,K,I-1,3)) * S(J,K,I-1,3,1) 
     3                                                * S(J,K,I-1,4,1)
            DWDX    = 0.5E0*(Q(J,K,I+1,4)+Q(J,K,I,4)) * S(J,K,I,1,1)   
     1                                                * S(J,K,I,4,1) 
     2              - 0.5E0*(Q(J,K,I,4)+Q(J,K,I-1,4)) * S(J,K,I-1,1,1) 
     3                                                * S(J,K,I-1,4,1)
            DWDY    = 0.5E0*(Q(J,K,I+1,4)+Q(J,K,I,4)) * S(J,K,I,2,1)   
     1                                                * S(J,K,I,4,1) 
     2              - 0.5E0*(Q(J,K,I,4)+Q(J,K,I-1,4)) * S(J,K,I-1,2,1) 
     3                                                * S(J,K,I-1,4,1)
            TVIS(J,K,1) = DUDY - DVDX
            TVIS(J,K,2) = DVDZ - DWDY
            TVIS(J,K,3) = DWDX - DUDZ
  110    CONTINUE
  120 CONTINUE
C
C     Form the J contribution to the magnitude of vorticity
C
C     Calculate U,V,W at cell interfaces as averages of cell centers
C     Store U in DIST(1), V in DIST(2) and W in DIST(3)
C
      DO 220 K = 2, KDIM
         DO 210 J = 0, JDIM + 1
            DIST(J,K,1) = 0.5E0 * ( Q(J,K,I,2) + Q(J+1,K,I,2) )
            DIST(J,K,2) = 0.5E0 * ( Q(J,K,I,3) + Q(J+1,K,I,3) )
            DIST(J,K,3) = 0.5E0 * ( Q(J,K,I,4) + Q(J+1,K,I,4) )
  210    CONTINUE
  220 CONTINUE
C
C     Loop over boundary conditions and fix the wall data
C
      DO 250 IBC = 1, NBCS
         IF (IBCDAT(1,IBC) .EQ. IWALL .AND. 
     1       IBCDAT(2,IBC) .EQ. 2) THEN
            INROUT = IBCDAT(3,IBC)
            ISTRT  = IBCDAT(4,IBC)
            IEND   = IBCDAT(5,IBC)
            KSTRT  = IBCDAT(8,IBC)
            KEND   = IBCDAT(9,IBC)
            IF (INROUT .EQ. 1) THEN
               JBCPT = 1
               JNTR  = 1
               J1W   = 1
            ELSE
               JBCPT = JDIM + 1
               JNTR  = JDIM
               JDIMW = 1
            ENDIF
            IF (I .GE. ISTRT .AND. I .LE. IEND) THEN
               DO 230 K = KSTRT, KEND
                  DIST(JNTR,K,1) = 0.E0
                  DIST(JNTR,K,2) = 0.E0
                  DIST(JNTR,K,3) = 0.E0
  230          CONTINUE
            ENDIF
         ENDIF
  250 CONTINUE
C
C     Loop over K forming the contribution along the J direction
C     to the vorticity 
C
      DO 270 K = 2, KDIM
         DO 260 J = 1, JDIM
C
C     Calculate first derivative terms using Divergence Theorem
C
            DUDY = DIST(J,K,1)   * S(J,K,I,2,2)   * S(J,K,I,4,2) -
     1             DIST(J-1,K,1) * S(J-1,K,I,2,2) * S(J-1,K,I,4,2)
            DUDZ = DIST(J,K,1)   * S(J,K,I,3,2)   * S(J,K,I,4,2) -
     1             DIST(J-1,K,1) * S(J-1,K,I,3,2) * S(J-1,K,I,4,2)
            DVDX = DIST(J,K,2)   * S(J,K,I,1,2)   * S(J,K,I,4,2) - 
     1             DIST(J-1,K,2) * S(J-1,K,I,1,2) * S(J-1,K,I,4,2)
            DVDZ = DIST(J,K,2)   * S(J,K,I,3,2)   * S(J,K,I,4,2) -
     1             DIST(J-1,K,2) * S(J-1,K,I,3,2) * S(J-1,K,I,4,2)
            DWDX = DIST(J,K,3)   * S(J,K,I,1,2)   * S(J,K,I,4,2) - 
     1             DIST(J-1,K,3) * S(J-1,K,I,1,2) * S(J-1,K,I,4,2)
            DWDY = DIST(J,K,3)   * S(J,K,I,2,2)   * S(J,K,I,4,2) -
     1             DIST(J-1,K,3) * S(J-1,K,I,2,2) * S(J-1,K,I,4,2)
            TVIS(J,K,1) = TVIS(J,K,1) + DUDY - DVDX
            TVIS(J,K,2) = TVIS(J,K,2) + DVDZ - DWDY
            TVIS(J,K,3) = TVIS(J,K,3) + DWDX - DUDZ
  260    CONTINUE
  270 CONTINUE
C
C     Calculate U,V,W at cell interfaces as averages of cell centers
C     Store U in DIST(1), V in DIST(2) and W in DIST(3)
C
      DO 320 J = 2, JDIM
         DO 310 K = 0, KDIM + 1
            DIST(J,K,1) = 0.5E0 * ( Q(J,K,I,2) + Q(J,K+1,I,2) )
            DIST(J,K,2) = 0.5E0 * ( Q(J,K,I,3) + Q(J,K+1,I,3) )
            DIST(J,K,3) = 0.5E0 * ( Q(J,K,I,4) + Q(J,K+1,I,4) )
  310    CONTINUE
  320 CONTINUE
C
C     Loop over boundary conditions and fix the wall data
C
      DO 350 IBC = 1, NBCS
         IF (IBCDAT(1,IBC) .EQ. IWALL .AND. 
     1       IBCDAT(2,IBC) .EQ. 3) THEN
            INROUT = IBCDAT(3,IBC)
            ISTRT  = IBCDAT(4,IBC)
            IEND   = IBCDAT(5,IBC)
            JSTRT  = IBCDAT(6,IBC)
            JEND   = IBCDAT(7,IBC)
            IF (INROUT .EQ. 1) THEN
               KBCPT = 1
               KNTR  = 1
               K1W   = 1
            ELSE
               KBCPT = KDIM + 1
               KNTR  = KDIM
               KDIMW = 1
            ENDIF
            IF (I .GE. ISTRT .AND. I .LE. IEND) THEN
               DO 330 J = JSTRT, JEND
                  DIST(J,KNTR,1)   = 0.E0
                  DIST(J,KNTR,2)   = 0.E0
                  DIST(J,KNTR,3)   = 0.E0
  330          CONTINUE
            ENDIF
         ENDIF
  350 CONTINUE
C
C     Loop over J forming the contribution along the K direction 
C     to the vorticity 
C
C
      DO 370 J = 2, JDIM
         DO 360 K = 1, KDIM
C
C     Calculate first derivative terms using Divergence Theorem
C
            DUDY = DIST(J,K,1)   * S(J,K,I,2,3)   * S(J,K,I,4,3) -
     1             DIST(J,K-1,1) * S(J,K-1,I,2,3) * S(J,K-1,I,4,3)
            DUDZ = DIST(J,K,1)   * S(J,K,I,3,3)   * S(J,K,I,4,3) -
     1             DIST(J,K-1,1) * S(J,K-1,I,3,3) * S(J,K-1,I,4,3)
            DVDX = DIST(J,K,2)   * S(J,K,I,1,3)   * S(J,K,I,4,3) - 
     1             DIST(J,K-1,2) * S(J,K-1,I,1,3) * S(J,K-1,I,4,3)
            DVDZ = DIST(J,K,2)   * S(J,K,I,3,3)   * S(J,K,I,4,3) -
     1             DIST(J,K-1,2) * S(J,K-1,I,3,3) * S(J,K-1,I,4,3)
            DWDX = DIST(J,K,3)   * S(J,K,I,1,3)   * S(J,K,I,4,3) - 
     1             DIST(J,K-1,3) * S(J,K-1,I,1,3) * S(J,K-1,I,4,3)
            DWDY = DIST(J,K,3)   * S(J,K,I,2,3)   * S(J,K,I,4,3) -
     1             DIST(J,K-1,3) * S(J,K-1,I,2,3) * S(J,K-1,I,4,3)
            TVIS(J,K,1) = TVIS(J,K,1) + DUDY - DVDX
            TVIS(J,K,2) = TVIS(J,K,2) + DVDZ - DWDY
            TVIS(J,K,3) = TVIS(J,K,3) + DWDX - DUDZ
  360    CONTINUE
  370 CONTINUE
C
C     Divide by the Volume and form the magnitude of vorticity
C
c***This may cause problems when VOL = 0 (only for j,k = 1 and j or k = dim + 1
      DO 390 K = 1, KDIM + 1
         DO 380 J = 1, JDIM + 1
            VOLI = 1.E0 / (VOL(J,K,I) + RSMASQ)
            VORT(J,K) = SQRT ( TVIS(J,K,1) * TVIS(J,K,1)
     1                       + TVIS(J,K,2) * TVIS(J,K,2)
     2                       + TVIS(J,K,3) * TVIS(J,K,3) ) * VOLI
  380    CONTINUE
  390 CONTINUE
C
C     Calculate the turbulent eddy viscosity normal to each wall
C
      NBRWLS = 0
C
C     Calculate the turbulent eddy viscosity in the J direction
C
C     J = 1 wall
C
      IF (J1W .EQ. 1) THEN
         NBRWLS = NBRWLS + 1
         JBC    = 1
         JN     = 2
         JM     = 1
         JMP1   = 2
C
C     Calculate the directed distance 'normal' to the wall
C
         DO 405 K = 2, KDIM
            DIST(1,K,NBRWLS) = 0.E0
            DIST(2,K,NBRWLS) = 0.E0
  405    CONTINUE
         DO 420 J = 2, JDIM
            DO 410 K = 2, KDIM
               DY = VOL(J,K,I) / (0.5E0*(S(J,K,I,4,2)+S(J-1,K,I,4,2)))
               DIST(J,K,NBRWLS)   = DIST(J,K,NBRWLS) + 0.5E0 * DY
               DIST(J+1,K,NBRWLS) = DIST(J,K,NBRWLS) + 0.5E0 * DY
  410       CONTINUE
  420    CONTINUE
C
C     Loop over each profile and calculate the turbulent eddy viscosity
C     in that profile.
C
         DO 450 K = 2, KDIM 
C
C     Find the maximum and minimum of the velocity magnitude (Umax, Umin)
C     in the profile and get the density and vorticity in an array for
C     the appropriate algebraic eddy viscosity model.
C
            UMAX = - 999999.E0
            UMIN =   999999.E0
            DO 430 J = 1, JDIM + 1
               UMAGN    = SQRT ( Q(J,K,I,2)*Q(J,K,I,2) 
     1                         + Q(J,K,I,3)*Q(J,K,I,3)
     2                         + Q(J,K,I,4)*Q(J,K,I,4) )
               UMAX     = MAX ( UMAX, UMAGN )
               UMIN     = MIN ( UMIN, UMAGN )
               RHO(J)   = Q(J,K,I,1)
               VORTJ(J) = VORT(J,K)
               Y(J)     = DIST(J,K,NBRWLS)
  430       CONTINUE
C
C     Calculate 'wall' values of density, viscosity and shear stress
C
            RHOW  = Q(JBC,K,I,1)
            TW    = GAMMA * Q(JBC,K,I,5) / Q(JBC,K,I,1)
            RMUW  = VIS (TW, SNUM, STINF)
C
            UBARN = S(JM,K,I,1,2) * Q(JN,K,I,2)
     1            + S(JM,K,I,2,2) * Q(JN,K,I,3)
     2            + S(JM,K,I,3,2) * Q(JN,K,I,4)
            VELN  = SQRT ( (Q(JN,K,I,2) - UBARN*S(JM,K,I,1,2))**2
     1                   + (Q(JN,K,I,3) - UBARN*S(JM,K,I,2,2))**2
     2                   + (Q(JN,K,I,4) - UBARN*S(JM,K,I,3,2))**2 )
C
C     Average metrics and cell volume
C
            VOLAUX = 0.5E0 * VOL(JMP1,K,I)
            SPJ1   = 0.5E0 * ( S(JMP1,K,I,1,2)*S(JMP1,K,I,4,2)
     1                       + S(JM,K,I,1,2)*S(JM,K,I,4,2) )
            SPJ2   = 0.5E0 * ( S(JMP1,K,I,2,2)*S(JMP1,K,I,4,2)
     1                       + S(JM,K,I,2,2)*S(JM,K,I,4,2) )
            SPJ3   = 0.5E0 * ( S(JMP1,K,I,3,2)*S(JMP1,K,I,4,2)
     1                       + S(JM,K,I,3,2)*S(JM,K,I,4,2) )
            SMJ1   = S(JM,K,I,1,2)*S(JM,K,I,4,2)
            SMJ2   = S(JM,K,I,2,2)*S(JM,K,I,4,2)
            SMJ3   = S(JM,K,I,3,2)*S(JM,K,I,4,2)
C
C     Calculate first derivative terms using Divergence Theorem
C
            DVELDX = VELN*SPJ1 / VOLAUX
            DVELDY = VELN*SPJ2 / VOLAUX
            DVELDZ = VELN*SPJ3 / VOLAUX
C
C     Tau(wall) 
C
            DVELDN = DVELDX * S(JM,K,I,1,2) + DVELDY * S(JM,K,I,2,2)
     1             + DVELDZ * S(JM,K,I,3,2)
            TAUW   = RMUW * DVELDN 
C
C     Calculate the eddy viscosity using the algebraic method of
C     Baldwin and Lomax (AIAA 78-257)
C
            JDIMP1 = JDIM + 1
            CALL BLDLMX (JDIMP1, RHOW, TAUW, RMUW, UMAX, UMIN, Y, RHO,
     1                   VORTJ, FUN, RMUTJ)
C
C     Store these values of eddy viscosity into arrays used to build 
C     the total eddy viscosity.
C
            DO 440 J = 1, JDIM + 1
               TVIS(J,K,NBRWLS) = RMUTJ(J)
  440       CONTINUE
  450    CONTINUE
      ENDIF
C
C     Build up the total turbulent eddy viscosity from the contribution 
C     from each wall.
C
      DO 920 K = 1, KDIM
         DO 910 J = 1, JDIM + 1
            RMUT(J,K,I) = TVIS(J,K,1)
  910    CONTINUE
  920 CONTINUE
C
C     Finished calculating the turbulent eddy viscosity
C
      RETURN
 1000 FORMAT (' ','RUN ABORTING. Invalid Turbulence Model specified.')
      END
