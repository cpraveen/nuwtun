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




      SUBROUTINE GETBIJ (NPTS, Q, BIJ)
C
C Routine to calculate the anisotropy tensor for Reynolds stresses.
C
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C BIJ            : Anisotropy tensor
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ)
C
C     +---------------+
C     +      BIJ      +
C     +  (i,j)   sub  +
C     +---------------+
C     +   x,x     1   +
C     +   y,y     2   +
C     +   z,z     3   +
C     +   x,y     4   +
C     +   x,z     5   +
C     +   y,z     6   +
C     +---------------+
C
      DIMENSION BIJ   (NPTS,6)
C
C T13
C
      T13 = 1.E0 / 3.E0
C
C Calculate the Anisotropy tensor b_ij
C
      DO 100 I = 1, NPTS
C
         TAUXX = Q(I,6)
         TAUYY = Q(I,7)
         TAUZZ = Q(I,8)
         TAUXY = Q(I,9)
         TAUXZ = Q(I,10)
         TAUYZ = Q(I,11)
C
         TKE   = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
C b_ij, Anisotropy tensor
C
         BIJ(I,1) = 0.5E0 * TAUXX / TKE - T13
         BIJ(I,2) = 0.5E0 * TAUYY / TKE - T13
         BIJ(I,3) = 0.5E0 * TAUZZ / TKE - T13
         BIJ(I,4) = 0.5E0 * TAUXY / TKE
         BIJ(I,5) = 0.5E0 * TAUXZ / TKE
         BIJ(I,6) = 0.5E0 * TAUYZ / TKE
  100 CONTINUE
C
C Finished calculating anisotropy tensor
C
      RETURN
      END






      SUBROUTINE GETSW  (NPTS, Q, DQDX, DQDY, DQDZ, SIJ, WIJ)
C
C Routine to calculate the symmetric and anti-symmetric tensors S_ij, W_ij
C from the mean velocity gradient.
C
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C DQDX,DQDY,DQDZ : Derivatives of Q at cell centers
C SIJ            : Symmetric tensor      (du_i/dx_j + du_j/dx_i)/2
C WIJ            : Anit-Symmetric tensor (du_i/dx_j - du_j/dx_i)/2
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          DQDX  (NPTS,NQ),
     2          DQDY  (NPTS,NQ),
     3          DQDZ  (NPTS,NQ)
C
C     +-------------+---------------+
C     +     SIJ     +       WIJ     +
C     + (i,j)  sub  +  (i,j)   sub  +
C     +-------------+---------------+
C     +  x,x    1   +   y,x     1   +
C     +  y,y    2   +   z,x     2   +
C     +  z,z    3   +   z,y     3   +
C     +  x,y    4   +   x,y     4   +
C     +  x,z    5   +   x,z     5   +
C     +  y,z    6   +   y,z     6   +
C     +-------------+---------------+
C
      DIMENSION SIJ   (NPTS,6),
     1          WIJ   (NPTS,6)
C
C Calculate the symmetric tensor S_ij and anti-symmetric tensor W_ij
C
      DO 100 I = 1, NPTS
C
C S_ij, Symmetric tensor
C
         SIJ(I,1) = DQDX(I,2)
         SIJ(I,2) = DQDY(I,3)
         SIJ(I,3) = DQDZ(I,4)
         SIJ(I,4) = 0.5E0 * (DQDY(I,2) + DQDX(I,3))
         SIJ(I,5) = 0.5E0 * (DQDZ(I,2) + DQDX(I,4))
         SIJ(I,6) = 0.5E0 * (DQDZ(I,3) + DQDY(I,4))
C
C W_ij, Anti-Symmetric tensor
C
         WIJ(I,1) = 0.5E0 * (DQDX(I,3) - DQDY(I,2))
         WIJ(I,2) = 0.5E0 * (DQDX(I,4) - DQDZ(I,2))
         WIJ(I,3) = 0.5E0 * (DQDY(I,4) - DQDZ(I,3))
         WIJ(I,4) = 0.5E0 * (DQDY(I,2) - DQDX(I,3))
         WIJ(I,5) = 0.5E0 * (DQDZ(I,2) - DQDX(I,4))
         WIJ(I,6) = 0.5E0 * (DQDZ(I,3) - DQDY(I,4))
  100 CONTINUE
C
C     Finished calculating tensors
C
      RETURN
      END






      SUBROUTINE SRCRS (IMODEL, NPTS, Q, PROPS, DQDX, DQDY, DQDZ, 
     1                  PIIJ, DEL2K, TKE, ETILDE, XSIW, FDAMP,
     2                  EPSC, PDIL, SIJ, WIJ, SRC, NPRLIM, NPRNEG)
C
C Routine to calculate the source terms for the Reynolds stress
C turbulence model.
C
C This routine is organized into the following sections:
C      1. Calculate damping functions for epsilon equation
C      2. Calculate b_ij
C      3. Calculate compressible dissipation
C      4. Calculate the pressure-strain rate correlation from b, S, W
C      5. Calculate the source terms for the tau and epsilon equations
C      6. Add the epsilon_ij terms 
C
C IMODEL         : Turbulence model array
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C PROPS          : Properties stored at cell centers
C                  PROPS(1) = RMU   molecular viscosity
C                  PROPS(2) = RMUT  turbulent eddy viscosity
C                  PROPS(3) = YPLUS Y+
C DQDX,DQDY,DQDZ : Derivatives of Q at cell centers
C PIIJ           : Pressure Strain Correlation
C                  This array contains TAU=-rho tau_ij on entrance to routine
C DEL2K          : Laplacian of turbulent kinetic energy
C TKE            :
C ETILDE         :
C XSIW           :
C FDAMP          :
C EPSC           :
C PDIL           :
C SIJ            :
C WIJ            :
C SRC            : Source terms for the k-omega model at the cell centers
C                      of the I-plane
C NPRLIM         : Number of limitations to the production term
C NPRNEG         : Number of times invoke positivity preservation for production
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          PROPS (NPTS,NP)
C
      DIMENSION DQDX  (NPTS,NQ),
     1          DQDY  (NPTS,NQ),
     2          DQDZ  (NPTS,NQ),
     3          PIIJ  (NPTS,6),
     4          DEL2K (NPTS),
     5          TKE   (NPTS),
     6          ETILDE(NPTS),
     7          XSIW  (NPTS),
     8          FDAMP (NPTS,2),
     9          EPSC  (NPTS),
     A          PDIL  (NPTS),
     B          SIJ   (NPTS,6),
     C          WIJ   (NPTS,6),
     D          SRC   (NPTS,NF)
C
      DIMENSION IMODEL(NMDL)
C
C Two-Thirds
C
      T23 = 2.E0 / 3.E0
C
      ITURB  = IMODEL(1)
      IPSMDL = IMODEL(6)
      IEPSIJ = IMODEL(7)
C
C 1. Calculate the damping functions for the epsilon equation
C    PIIJ contains TAU = - rho tau_ij at this point
C
      CALL DAMPEE (IMODEL, NPTS, Q, PROPS, DEL2K, TKE, 
     1             PIIJ, DQDX, DQDY, DQDZ, ETILDE, XSIW, FDAMP)
C
C 2. Call routine to calculate b_ij
C    Temporary space:
C       b_ij: PIIJ
C       NOTE: use b_ij: SRC(1,1) in the future
C    this assumes that SRC is at least 12 long and has not been used yet
C
      CALL GETBIJ (NPTS, Q, PIIJ)
C
C 3. Calculate compressible dissipation
C
      CALL EPSCMP (IMODEL, NPTS, Q, PROPS, TKE, SIJ, WIJ, EPSC)
C
C 4. Calculate pressure-dilatation
c***********
c***This won't work as I currently don't have TAU
c***********
C
      CALL PRDIL  (IMODEL, NPTS, Q, PROPS, DQDX, DQDY, DQDZ, TAU,
     1             TKE, SIJ, WIJ, EPSC, PDIL)
C
C 5. Call routine to calculate Pressure-Strain Correlation
C    b_ij is passed in in PIIJ.  The pressure-strain is returned in PIIJ.
C
      CALL PSTRAN (IPSMDL, NPTS, Q, PROPS, SIJ, WIJ, PIIJ)
C
C 6. Form the source terms for the Reynolds stress model 
C
      DO 100 I = 1, NPTS
C
         RHO    = Q(I,1)
         P      = Q(I,5)
         TAUXX  = Q(I,6)
         TAUYY  = Q(I,7)
         TAUZZ  = Q(I,8)
         TAUXY  = Q(I,9)
         TAUXZ  = Q(I,10)
         TAUYZ  = Q(I,11)
         EPSLN  = Q(I,12)
C
         RMU    = PROPS(I,1)
         RMUT   = PROPS(I,2)
         YPLUS  = PROPS(I,3)
         YNORML = ABS (PROPS(I,4))
C
C Transition specification
C
         TRANS = 1.E0
         IF (PROPS(I,4) .LT. 0.E0) TRANS = 0.E0
C
C DUKDXK = DU_k/DX_k / 3 = (DU/DX + DV/DY + DW/DZ) / 3
C
         DUKDXK = (DQDX(I,2) + DQDY(I,3) + DQDZ(I,4)) / 3.E0
C
C Production terms P_ij:
C    P_ij = - rho Tau_ik du_j/dx_k - rho Tau_jk du_i/dx_k
C
         PXX   = - 2.E0 * RHO * ( TAUXX * DQDX(I,2) + 
     1                            TAUXY * DQDY(I,2) + 
     2                            TAUXZ * DQDZ(I,2) ) * TRANS
         PYY   = - 2.E0 * RHO * ( TAUXY * DQDX(I,3) + 
     1                            TAUYY * DQDY(I,3) + 
     2                            TAUYZ * DQDZ(I,3) ) * TRANS
         PZZ   = - 2.E0 * RHO * ( TAUXZ * DQDX(I,4) + 
     1                            TAUYZ * DQDY(I,4) + 
     2                            TAUZZ * DQDZ(I,4) ) * TRANS
         PXY   = - RHO * ( TAUXX * DQDX(I,3) + 
     1                     TAUXY * DQDY(I,3) + 
     2                     TAUXZ * DQDZ(I,3) +
     3                     TAUXY * DQDX(I,2) +
     4                     TAUYY * DQDY(I,2) +
     5                     TAUYZ * DQDZ(I,2) ) * TRANS
         PXZ   = - RHO * ( TAUXX * DQDX(I,4) + 
     1                     TAUXY * DQDY(I,4) + 
     2                     TAUXZ * DQDZ(I,4) +
     3                     TAUXZ * DQDX(I,2) +
     4                     TAUYZ * DQDY(I,2) +
     5                     TAUZZ * DQDZ(I,2) ) * TRANS
         PYZ   = - RHO * ( TAUXY * DQDX(I,4) + 
     1                     TAUYY * DQDY(I,4) + 
     2                     TAUYZ * DQDZ(I,4) +
     3                     TAUXZ * DQDX(I,3) +
     4                     TAUYZ * DQDY(I,3) +
     5                     TAUZZ * DQDZ(I,3) ) * TRANS
C
         PROD  = 0.5E0 * (PXX + PYY + PZZ)
C
C Positivity preservation for production
C
         IF (POSPRD) THEN
            IF (PROD .LT. 0.E0) THEN
               NPRNEG = NPRNEG + 1
               PROD   = ABS (PROD)
            ENDIF
         ENDIF
C
         PRODUL = PROD
C
C Limit production term for robustness:
C
         TSTDIS = PRDLIM * RHO * EPSLN
         IF (PROD .GT. TSTDIS) THEN
            NPRLIM = NPRLIM + 1
            PROD   = TSTDIS
         ENDIF
C
C Production term for epsilon:
C->  PRODE = - CEPS1*epsilon/k*rho*tau_ij*(Du_i/Dx_j-Du_k/Dx_k*delta_ij/3)
C->          - 4/3*rho*epsilon*Du_k/Dx_k
C->  from: Speziale and Sarkar, ICASE 91-9
C
         PRODE  = CEPS1 * FDAMP(I,1) * EPSLN / TKE(I) 
     1            * (PRDE * PROD + PRDEM1 * PRODUL)
C
c->      PRODE  = - CEPS1 * FDAMP(I,1) * RHO * EPSLN / TKE(I) *
c->  1               ( TAUXX * (DQDX(I,2) - DUKDXK)
c->  2               + TAUYY * (DQDY(I,3) - DUKDXK)
c->  3               + TAUZZ * (DQDZ(I,4) - DUKDXK)
c->  4               + TAUXY * (DQDY(I,2) + DQDX(I,3)) 
c->  5               + TAUXZ * (DQDZ(I,2) + DQDX(I,4))
c->  6               + TAUYZ * (DQDZ(I,3) + DQDY(I,4)) )
C
C Additional term from Speziale - factor of 3 is required because of 
C definition DUKDXK
C
         DNUDT  = (1.E0 - 0.7E0 * GAMM1) * DUKDXK * 3.E0
C
C Dissipation
C
         DISSE  = CEPS2 * FDAMP(I,2) * RHO * EPSLN * ETILDE(I) / TKE(I)
C
C Calculate the Pressure-Velocity correlation term:
C    PVEL  = CPVELK * mu_T * sigma_rho / rho^2 * (DrhoDx_i DPDx_i)
C
         PVEL  = CPVELK * RMUT * SIGRHO / (RHO * RHO) * FSMACH / RE *
     1           ( DQDX(I,1) * DQDX(I,5) + DQDY(I,1) * DQDY(I,5) 
     2           + DQDZ(I,1) * DQDZ(I,5) ) 
         PVELXX = 0.e0
         PVELYY = 0.e0
         PVELZZ = 0.e0
         PVELXY = 0.e0
         PVELXZ = 0.e0
         PVELYZ = 0.e0
C
C Calculate Source Term
C
         SRC(I,1)  = 0.E0
         SRC(I,2)  = 0.E0
         SRC(I,3)  = 0.E0
         SRC(I,4)  = 0.E0
         SRC(I,5)  = 0.E0
         SRC(I,6)  = PXX + PIIJ(I,1) - PVELXX + PDIL(I)
         SRC(I,7)  = PYY + PIIJ(I,2) - PVELYY + PDIL(I)
         SRC(I,8)  = PZZ + PIIJ(I,3) - PVELZZ + PDIL(I)
         SRC(I,9)  = PXY + PIIJ(I,4) - PVELXY 
         SRC(I,10) = PXZ + PIIJ(I,5) - PVELXZ 
         SRC(I,11) = PYZ + PIIJ(I,6) - PVELYZ 
         SRC(I,12) = PRODE - DISSE + XSIW(I)
c->  1               + (DNUDT - 4.E0 * DUKDXK) * RHO * EPSLN
  100 CONTINUE
C
C 6. Add source terms due to dissipation
C
      CALL SRCEIJ (IEPSIJ, NPTS, Q, PROPS, DQDX, DQDY, DQDZ, 
     1             EPSC, DEL2K, SRC)
C
C Finished with Reynolds Stress source terms
C
      RETURN
      END






      SUBROUTINE SRCRSW (IMODEL, Q, PROPS, SN, DQDX, DQDY, DQDZ, 
     1                   TKE, EPSC, SRC, SRCW,
     2                   IDIM, JDIM, KDIM, I, NBCS, IBCDAT, 
     3                   S, RC, YNI, YNITOT)
C
C Routine to calculate the wall normal dependent source terms for the 
C Reynolds stress turbulence model.
C
C Assumptions:
C    I.   The contributions due to walls are restricted to walls on the 
C         current block
C    II.  The total contribution due to all walls is calculated as weighted
C         sum of each wall as:
C             a_ij_walls = Sum (a_ij_l / yn_l)   / Sum (1 / yn_l)
C                           l=1,walls               l=1,walls
C    III. The distance to any given wall is sufficiently approximated as the
C         distance from the cell center to the wall point = dist (R - R_wall)
C
C This routine is organized into the following sections:
C      1. Initialize total inverse distance to zero
C      2. Loop over all boundary conditions
C         a. If WALL boundary then
C            i.   Initialize SN, YNI to zero on current plane
C            ii.  Calculate SN, inverse distance to wall along normal coordinate
C                 1. I direction )
C                 2. J direction } For each wall point, calculate SN, YNI
C                 3. K direction )
C                 4. Error
C            iii. Sum inverse distance to total inverse distance
C            iv.  Add weighted pressure-strain to additional source terms
C            v.   Add weighted epsilon_ij to additional source terms
C      3. Weight additional source terms by total inverse distance and add 
C         to residual
C
C IMODEL         : Turbulence model array
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C PROPS          : Properties stored at cell centers
C                  PROPS(1) = RMU   molecular viscosity
C                  PROPS(2) = RMUT  turbulent eddy viscosity
C                  PROPS(3) = YPLUS Y+
C SN             : Wall normals
C DQDX,DQDY,DQDZ : Derivatives of Q at cell centers
C TKE            :
C EPSC           :
C SRC            : Source terms at the cell centers of the I-plane
C SRCW           : Wall normal dependent Reynolds stress source terms
C IDIM,JDIM,KDIM : Dimensions of current block
C I              : I plane to calculate additional source terms at
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
C S              : Metrics
C RC             : Grid position vector to cell center
C YNI            : Storage for inverse normal distance to wall
C YNITOT         : Storage for total inverse normal distance to wall for 
C                  all walls
C
      include '../header/common.h'
C
      DIMENSION Q     (2:JDIM,2:KDIM,NQ),
     1          PROPS (2:JDIM,2:KDIM,NP),
     2          SN    (2:JDIM,2:KDIM,3)
C
      DIMENSION DQDX  (2:JDIM,2:KDIM,NQ),
     1          DQDY  (2:JDIM,2:KDIM,NQ),
     2          DQDZ  (2:JDIM,2:KDIM,NQ),
     3          TKE   (2:JDIM,2:KDIM),
     4          EPSC  (2:JDIM,2:KDIM),
     5          SRC   (2:JDIM,2:KDIM,NF),
     6          SRCW  (2:JDIM,2:KDIM,6)
C
      DIMENSION IBCDAT(10,NBCS)
C
      DIMENSION S     (0:JDIM+1,0:KDIM+1,0:IDIM+1,4,3),
     1          RC    (  JDIM+1,  KDIM+1,  IDIM+1,3)
C
      DIMENSION YNI   (2:JDIM,2:KDIM),
     1          YNITOT(2:JDIM,2:KDIM)
C
      DIMENSION IMODEL(NMDL)
C
      LOGICAL   DOPLAN,
     1          DOWALL
C
C Two-Thirds
C
      T23 = 2.E0 / 3.E0
C
      ITURB  = IMODEL(1)
      IPSMDL = IMODEL(6)
      IEPSIJ = IMODEL(7)
C
C     Number of points in a cross-section to calculate source terms at
C
      NPTS   = (JDIM - 1) * (KDIM - 1)
C
C 1.  Initialize total inverse distance and wall source terms to zero
C
      DO 120 K = 2, KDIM
         DO 110 J = 2, JDIM
            YNITOT(J,K)   = 0.E0
            SRCW  (J,K,1) = 0.E0
            SRCW  (J,K,2) = 0.E0
            SRCW  (J,K,3) = 0.E0
            SRCW  (J,K,4) = 0.E0
            SRCW  (J,K,5) = 0.E0
            SRCW  (J,K,6) = 0.E0
  110    CONTINUE
  120 CONTINUE
C
C     Initialize DOPLAN to false so that the source terms are updated
C     only if walls affect the current plane
C
      DOPLAN = .FALSE.
C
C 2.  Loop over all boundary conditions
C
      DO 290 IBC = 1, NBCS
C
C     Initialize DOWALL to false so that the additional source terms
C     are calculated only if the current wall affects this plane
C
         DOWALL   = .FALSE.
         IBCTYP   = IBCDAT(1,IBC)
C
C 2.a. If WALL boundary then
C Store boundary condition data in mnemonic names
C
         IF (IBCTYP .EQ. IWALL) THEN
            IDIR     = IBCDAT(2,IBC)
            INROUT   = IBCDAT(3,IBC)
            ISTRT    = IBCDAT(4,IBC) 
            IEND     = IBCDAT(5,IBC)
            JSTRT    = IBCDAT(6,IBC)
            JEND     = IBCDAT(7,IBC)
            KSTRT    = IBCDAT(8,IBC)
            KEND     = IBCDAT(9,IBC)
C
C 2.a.i.    Initialize SN, YNI to zero
C
            DO 210 K = 2, KDIM
               DO 205 J = 2, JDIM
                  SN (J,K,1) = 0.E0
                  SN (J,K,2) = 0.E0
                  SN (J,K,3) = 0.E0
                  YNI(J,K)   = 0.E0
  205          CONTINUE
  210       CONTINUE
C
C 2.a.ii.   Calculate SN, inverse distance to wall along normal coordinate
C
C 2.a.ii-1. I direction
C
            IF (IDIR .EQ. 1) THEN
c**            write (iout,'(a)') ' Stopping in SRCRSw'
c**            stop
               DOWALL = .TRUE.
C
               IF (INROUT .EQ. 1) THEN
                  IBCW = 1
               ELSE
                  IBCW = IDIM + 1
               ENDIF
                  DO 230 K = KSTRT, KEND
                     DO 220 J = JSTRT, JEND
                        SN (J,K,1) = S (J,K,IBCW,1,IDIR)
                        SN (J,K,2) = S (J,K,IBCW,2,IDIR)
                        SN (J,K,3) = S (J,K,IBCW,3,IDIR)
                        DX         = RC(J,K,I,1) - RC(J,K,IBCW,1)
                        DY         = RC(J,K,I,2) - RC(J,K,IBCW,2)
                        DZ         = RC(J,K,I,3) - RC(J,K,IBCW,3)
                        YNI(J,K)   = SQRT (DX*DX + DY*DY + DZ*DZ)
  220                CONTINUE
  230             CONTINUE
C
C 2.a.ii-2. J direction
C
            ELSE IF (IDIR .EQ. 2) THEN
               DOWALL = .TRUE.
               IF (I .GE. ISTRT .AND. I .LE. IEND) THEN
                  IF (INROUT .EQ. 1) THEN
                     JBCW = 1
                  ELSE
                     JBCW = JDIM + 1
                  ENDIF
                  DO 250 K = KSTRT, KEND
                     DO 240 J = 2, JDIM
                        SN (J,K,1) = S(JBCW,K,I,1,IDIR)
                        SN (J,K,2) = S(JBCW,K,I,2,IDIR)
                        SN (J,K,3) = S(JBCW,K,I,3,IDIR)
                        DX         = RC(J,K,I,1) - RC(JBCW,K,I,1)
                        DY         = RC(J,K,I,2) - RC(JBCW,K,I,2)
                        DZ         = RC(J,K,I,3) - RC(JBCW,K,I,3)
                        YNI(J,K)   = SQRT (DX*DX + DY*DY + DZ*DZ)
  240                CONTINUE
  250             CONTINUE
               ENDIF
C
C 2.a.ii-3. K direction
C
            ELSE IF (IDIR .EQ. 3) THEN
c**            write (iout,'(a)') ' Stopping in SRCRSw'
c**            stop
               DOWALL = .TRUE.
               IF (I .GE. ISTRT .AND. I .LE. IEND) THEN
                  IF (INROUT .EQ. 1) THEN
                     KBCW = 1
                  ELSE
                     KBCW = KDIM + 1
                  ENDIF
                  DO 270 K = 2, KDIM
                     DO 260 J = JSTRT, JEND
                        SN (J,K,1) = S(J,KBCW,I,1,IDIR)
                        SN (J,K,2) = S(J,KBCW,I,2,IDIR)
                        SN (J,K,3) = S(J,KBCW,I,3,IDIR)
                        DX         = RC(J,K,I,1) - RC(J,KBCW,I,1)
                        DY         = RC(J,K,I,2) - RC(J,KBCW,I,2)
                        DZ         = RC(J,K,I,3) - RC(J,KBCW,I,3)
                        YNI(J,K)   = SQRT (DX*DX + DY*DY + DZ*DZ)
  260                CONTINUE
  270             CONTINUE
               ENDIF
C
C 2.a.ii-4. Error
C
            ELSE
               WRITE (IOUT,1000) IDIR
               STOP
            ENDIF
c---------
c     do k = 2, kdim
c        do j = 2, jdim
c           write (iout,'(a,3i3,4e14.6)') ' i,j,k,s,y ',
c    1             i, j, k, sn(j,k,1), sn(j,k,2), sn(j,k,3), yni(j,k)
c        end do
c     end do
c---------
C
C           Only call the additional source terms if a wall actually affects
C           the current plane
C
            IF (DOWALL) THEN
C
C           Set DOPLAN so that the source terms are updated
C
               DOPLAN = .TRUE.
C
C 2.a.iii.  Sum inverse distance to total inverse distance
C
               DO 285 K = 2, KDIM
                  DO 280 J = 2, JDIM
                     YNITOT(J,K) = YNITOT(J,K) + YNI(J,K)
  280             CONTINUE
  285          CONTINUE
C
C 2.a.iv.  Add weighted pressure-strain to additional source terms
C
               CALL PSTRANW (IPSMDL, NPTS, Q, PROPS, SN, YNI, SRCW)
C
C 2.a.v.   Add weighted epsilon_ij to additional source terms
C
               CALL SRCEIJW (IEPSIJ, NPTS, Q, PROPS, SN, EPSC, YNI,SRCW)
            ENDIF
         ENDIF
  290 CONTINUE
C
C 3. Weight additional source terms by total inverse distance and add 
C    to residual
C
      IF (DOPLAN) THEN
         DO 320 K = 2, KDIM
            DO 310 J = 2, JDIM
               DENOMI    = 1.E0 / (YNITOT(J,K) + RSMALL)
               SRC(J,K, 6) = SRC(J,K, 6) + SRCW(J,K,1) * DENOMI
               SRC(J,K, 7) = SRC(J,K, 7) + SRCW(J,K,2) * DENOMI
               SRC(J,K, 8) = SRC(J,K, 8) + SRCW(J,K,3) * DENOMI
               SRC(J,K, 9) = SRC(J,K, 9) + SRCW(J,K,4) * DENOMI
               SRC(J,K,10) = SRC(J,K,10) + SRCW(J,K,5) * DENOMI
               SRC(J,K,11) = SRC(J,K,11) + SRCW(J,K,6) * DENOMI
  310       CONTINUE
  320    CONTINUE
      ENDIF
C
C Finished with wall normal dependent Reynolds Stress source terms
C
      RETURN
 1000 FORMAT (' ','ERROR - RUN ABORTING',
     1       /' ','Incorrect direction specified in routine SRCRSW = ',
     2            I10)
      END






      SUBROUTINE SRCEIJ (IEPSIJ, NPTS, Q, PROPS, DQDX, DQDY, DQDZ, 
     1                   EPSC, DEL2K, SRC)
C
C Routine to calculate the source terms for the dissipation and epsilon
C equation in the Reynolds stress turbulence model.
C
C This routine is organized into the following sections:
C      1. Isotropic distribution (For isotropic models and the isotropic
C         portion of Lai-So model with the near wall portions in SRCEIJW)
C      2. Fu-Launder-Tselepidakis Model
C      3. Gatski-Speziale Algebraic Anisotropic Dissipation Rate Model
C      4. So et al. 1995 model (Additional near wall in the diffusion terms)
C
C IEPSIJ         : Model for epsilon_ij
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C PROPS          : Properties stored at cell centers
C                  PROPS(1) = RMU   molecular viscosity
C                  PROPS(2) = RMUT  turbulent eddy viscosity
C                  PROPS(3) = YPLUS Y+
C DQDX,DQDY,DQDZ : Derivatives of Q at cell centers
C EPSC           : Compressible dissipation
C DEL2K          : Laplacian of turbulent kinetic energy
C SRC            : Source terms at the cell centers
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          PROPS (NPTS,NP)
C
      DIMENSION DQDX  (NPTS,NQ),
     1          DQDY  (NPTS,NQ),
     2          DQDZ  (NPTS,NQ),
     3          EPSC  (NPTS),
     4          DEL2K (NPTS),
     5          SRC   (NPTS,NF)
C
C One-Third
C Two-Thirds
C
      T13 = 1.E0 / 3.E0
      T23 = 2.E0 / 3.E0
C
C Add the source terms for the dissipation and the epsilon equation to SRC
C
C 1. Isotropic (High Reynolds Number Model)
C    Isotropic portion of Lai-So model (SN (wall normal) dependent portion
C    of Lai-So model is in SRCEIJW.
C    NASA Contractor Report 4369 modifications from Thomas Sommer
C
      IF (IEPSIJ .EQ. IEISO .OR. IEPSIJ .EQ. IELAIS) THEN
         DO 100 I = 1, NPTS
C
            RHO    = Q(I,1)
            EPSLN  = Q(I,12)
C
C Dissipation = Solenoidal dissipation + compressible dissipation
C
            DISS   = T23 * RHO * (EPSLN + EPSC(I))
C
C Calculate Source Term
C
            SRC(I, 6) = SRC(I, 6) - DISS
            SRC(I, 7) = SRC(I, 7) - DISS
            SRC(I, 8) = SRC(I, 8) - DISS
  100    CONTINUE
C
C 2. Fu-Launder-Tselepidakis Model
C
      ELSE IF (IEPSIJ .EQ. IEFLT) THEN
         DO 200 I = 1, NPTS
C
            RHO    = Q(I,1)
            P      = Q(I,5)
            TAUXX  = Q(I,6)
            TAUYY  = Q(I,7)
            TAUZZ  = Q(I,8)
            TAUXY  = Q(I,9)
            TAUXZ  = Q(I,10)
            TAUYZ  = Q(I,11)
            EPSLN  = Q(I,12)
C
            TKE    = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
            RMU    = PROPS(I,1)
            RMUT   = PROPS(I,2)
C
C FLT variables
C
            AXX    = TAUXX / TKE - T23
            AYY    = TAUYY / TKE - T23
            AZZ    = TAUZZ / TKE - T23
            AXY    = TAUXY / TKE
            AXZ    = TAUXZ / TKE
            AYZ    = TAUYZ / TKE
C
C A_2 = a_ik a_ki
C
            A2     =       AXX*AXX + AYY*AYY + AZZ*AZZ
     1             + 2.E0*(AXY*AXY + AXZ*AXZ + AYZ*AYZ)
C
C A_3 = a_ik a_kl a_li
C
            A3     =      AXX*AXX*AXX + AYY*AYY*AYY + AZZ*AZZ*AZZ
     1             + 3.E0*AXY*AXY*(AXX + AYY)
     2             + 3.E0*AXZ*AXZ*(AXX + AZZ)
     3             + 3.E0*AYZ*AYZ*(AYY + AZZ)
     4             + 6.E0*AXY*AXZ*AYZ
C
C A   =  1 - 9(A_2 - A_3)/8
C     => 0 in 2D limit
C     => 1 in isotropic turbulence
C
            ASQRT  = SQRT (1.E0 - 1.125E0*(A2 - A3))
C
C Dissipation
C
            DISS   = RHO  * (EPSLN + EPSC(I))
            DISSXX = DISS * (TAUXX*(1.E0-ASQRT)/TKE + T23*ASQRT)
            DISSYY = DISS * (TAUYY*(1.E0-ASQRT)/TKE + T23*ASQRT)
            DISSZZ = DISS * (TAUZZ*(1.E0-ASQRT)/TKE + T23*ASQRT)
            DISSXY = DISS * (TAUXY*(1.E0-ASQRT)/TKE)
            DISSXZ = DISS * (TAUXZ*(1.E0-ASQRT)/TKE)
            DISSYZ = DISS * (TAUYZ*(1.E0-ASQRT)/TKE)
C
C Calculate Source Term
C
            SRC(I, 6) = SRC(I, 6) - DISSXX
            SRC(I, 7) = SRC(I, 7) - DISSYY
            SRC(I, 8) = SRC(I, 8) - DISSZZ
            SRC(I, 9) = SRC(I, 9) - DISSXY
            SRC(I,10) = SRC(I,10) - DISSXZ
            SRC(I,11) = SRC(I,11) - DISSYZ
  200    CONTINUE
C
C 3. Gatski-Speziale Algebraic Anisotropic Dissipation Rate Model
C
      ELSE IF (IEPSIJ .EQ. IEADGS) THEN
         CE5    = 5.8E0
         AL3    = 0.6E0
         CALFA  = 0.75E0 * (14.E0*AL3/11.E0 - 16.E0/33.E0)
         DO 300 I = 1, NPTS
C
            RHO       = Q(I, 1)
            P         = Q(I, 5)
            TAUXX     = Q(I, 6)
            TAUYY     = Q(I, 7)
            TAUZZ     = Q(I, 8)
            TAUXY     = Q(I, 9)
            TAUXZ     = Q(I,10)
            TAUYZ     = Q(I,11)
            EPSLN     = Q(I,12)
C
            TKE       = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
C Production term:
C    PROD = TAUXX*DUDX + TAUYY*DVDY + TAUZZ*DWDZ +
C           TAUXY*(DUDY+DVDX) + TAUXZ*(DUDZ+DWDX) + TAUYZ*(DVDZ+DWDY)
C    note: TAU_ij is defined as -rho tau_ij
C
            PROD      = - RHO * (  TAUXX * DQDX(I,2)
     1                           + TAUYY * DQDY(I,3)
     2                           + TAUZZ * DQDZ(I,4)
     3                           + TAUXY * (DQDY(I,2) + DQDX(I,3))
     4                           + TAUXZ * (DQDZ(I,2) + DQDX(I,4))
     5                           + TAUYZ * (DQDZ(I,3) + DQDY(I,4)) )
C
            IF (POSPRD) PROD = ABS (PROD)
C
C Limit production term for robustness:
C
            TSTDIS    = PRDLIM * RHO * EPSLN
            PROD      = MIN (PROD, TSTDIS)
C
C     Constants based on SSG pressure-strain
C
            POVERE    = PROD / (RHO*EPSLN)
C
            TIMSCL    = TKE / EPSLN
C
C     S_ij 
C
            SXX       = DQDX(I,2)
            SYY       = DQDY(I,3)
            SZZ       = DQDZ(I,4)
            SXY       = 0.5E0 * (DQDY(I,2) + DQDX(I,3))
            SXZ       = 0.5E0 * (DQDZ(I,2) + DQDX(I,4))
            SYZ       = 0.5E0 * (DQDZ(I,3) + DQDY(I,4))
C
C     S_kk = (S_xx + S_yy + S_zz) / 3
C
            SKK       = (SXX + SYY + SZZ) / 3.E0
C
C     S_kl * S_lk / 3
C
            SKLSLK    =  ( SXX*SXX + SYY*SYY + SZZ*SZZ
     1             + 2.E0*(SXY*SXY + SXZ*SXZ + SYZ*SYZ) ) / 3.E0
C
C     W_ij 
C
            WYX       = 0.5E0 * (DQDX(I,3) - DQDY(I,2))
            WZX       = 0.5E0 * (DQDX(I,4) - DQDZ(I,2))
            WZY       = 0.5E0 * (DQDY(I,4) - DQDZ(I,3))
            WXY       = 0.5E0 * (DQDY(I,2) - DQDX(I,3))
            WXZ       = 0.5E0 * (DQDZ(I,2) - DQDX(I,4))
            WYZ       = 0.5E0 * (DQDZ(I,3) - DQDY(I,4))
C
            ETASQ     =      TIMSCL * TIMSCL *
     1                      ( SXX*SXX + SYY*SYY + SZZ*SZZ
     2                + 2.E0*(SXY*SXY + SXZ*SXZ + SYZ*SYZ) )
            ZETSQ     =      TIMSCL * TIMSCL *
     1                  2.E0*(WXY*WXY + WXZ*WXZ + WYZ*WYZ) 
C
            DENOM     = CE5 + POVERE - 1.E0
C
            B1        = ( 7.E0 * AL3 + 1.E0) / (11.E0 * DENOM)
            B2        = (15.E0 * AL3 - 1.E0) / (11.E0 * DENOM)
C
            CMUEPS    = 1.E0 / ( 15.E0 * DENOM * 
     1                     (1.E0 + 2.E0*B1*B1*ZETSQ - T23*B2*B2*ETASQ) )
C
            B1KESQ    = B1 * TIMSCL * TIMSCL
            B2KESQ    = B2 * TIMSCL * TIMSCL
C
C Scaled Dissipation Anisotropy Tensor
C
            DXX    = - 2.E0 * CMUEPS * ( (SXX-SKK) * TIMSCL
     1                  +      B1KESQ*(        SXY*WYX+SXZ*WZX) 
     2                  + 2.E0*B2KESQ*(SXX*SXX+SXY*SXY+SXZ*SXZ-SKLSLK) )
            DYY    = - 2.E0 * CMUEPS * ( (SYY-SKK) * TIMSCL
     1                  +      B1KESQ*(SXY*WXY        +SYZ*WZY) 
     2                  + 2.E0*B2KESQ*(SXY*SXY+SYY*SYY+SYZ*SYZ-SKLSLK) )
            DZZ    = - 2.E0 * CMUEPS * ( (SZZ-SKK) * TIMSCL
     1                  +      B1KESQ*(SXZ*WXZ+SYZ*WYZ        ) 
     2                  + 2.E0*B2KESQ*(SXZ*SXZ+SYZ*SYZ+SZZ*SZZ-SKLSLK) )
            DXY    = - 2.E0 * CMUEPS * (  SXY * TIMSCL
     1                  +      B1KESQ*(SXX*WXY        +SXZ*WZY
     2                                        +SYY*WYX+SXZ*WZX) 
     3                  + 2.E0*B2KESQ*(SXX*SXY+SXY*SYY+SXZ*SYZ) )
            DXZ    = - 2.E0 * CMUEPS * (  SXZ * TIMSCL
     1                  +      B1KESQ*(SXX*WXZ+SXY*WYZ          
     2                                        +SYZ*WYX+SZZ*WZX) 
     3                  + 2.E0*B2KESQ*(SXX*SXZ+SXY*SYZ+SXZ*SZZ) )
            DYZ    = - 2.E0 * CMUEPS * (  SYZ * TIMSCL
     1                  +      B1KESQ*(SXY*WXZ+SYY*WYZ          
     2                                +SXZ*WXY        +SZZ*WZY) 
     3                  + 2.E0*B2KESQ*(SXY*SXZ+SYY*SYZ+SYZ*SZZ) )
C
C Eps_ij
C
            DISSXX    = 2.E0 * RHO * EPSLN * (DXX + T13)
            DISSYY    = 2.E0 * RHO * EPSLN * (DYY + T13)
            DISSZZ    = 2.E0 * RHO * EPSLN * (DZZ + T13)
            DISSXY    = 2.E0 * RHO * EPSLN * (DXY)
            DISSXZ    = 2.E0 * RHO * EPSLN * (DXZ)
            DISSYZ    = 2.E0 * RHO * EPSLN * (DYZ)
C
C Source term in dissipation transport equation
C           R = - 2 (1 + alpha) rho epsilon d_ij S_ij
C
            REPSLN    = - 2.E0 * (1.E0 + CALFA) * RHO * EPSLN *
     1                            ( DXX*SXX + DYY*SYY + DZZ*SZZ
     2                      + 2.E0*(DXY*SXY + DXZ*SXZ + DYZ*SYZ) )
C
C Calculate Source Term
C
            SRC(I, 6) = SRC(I, 6) - DISSXX
            SRC(I, 7) = SRC(I, 7) - DISSYY
            SRC(I, 8) = SRC(I, 8) - DISSZZ
            SRC(I, 9) = SRC(I, 9) - DISSXY
            SRC(I,10) = SRC(I,10) - DISSXZ
            SRC(I,11) = SRC(I,11) - DISSYZ
            SRC(I,12) = SRC(I,12) + REPSLN
  300    CONTINUE
C
C 4.  So et al. 1995 model (Additional near wall in the diffusion terms)
C
      ELSE IF (IEPSIJ .EQ. IESO95) THEN
         DO 400 I = 1, NPTS
C
            RHO       = Q(I, 1)
            TAUXX     = Q(I, 6)
            TAUYY     = Q(I, 7)
            TAUZZ     = Q(I, 8)
            TAUXY     = Q(I, 9)
            TAUXZ     = Q(I,10)
            TAUYZ     = Q(I,11)
            EPSLN     = Q(I,12)
C
            TKE       = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
            RMU       = PROPS(I,1)
            YPLUS     = PROPS(I,3)
C
C Dissipation = Solenoidal dissipation + compressible dissipation
C
            DISS      = T23 * RHO * (EPSLN + EPSC(I))
C
C Epsilon_ij = Isotropic + Near Wall Correction
C
            RET       = RHO * TKE * TKE / (EPSLN * RMU) * RE / FSMACH
            FW2       = EXP( -(RET /  40.E0)**2 )
C
            TERM      = FW2 * RHO * (EPSLN + EPSC(I)) / (TKE + RSMASQ)
            EPSXX     = DISS + TERM * (TAUXX - T23*TKE)
            EPSYY     = DISS + TERM * (TAUYY - T23*TKE)
            EPSZZ     = DISS + TERM * (TAUZZ - T23*TKE)
            EPSXY     =        TERM *  TAUXY
            EPSXZ     =        TERM *  TAUXZ
            EPSYZ     =        TERM *  TAUYZ
C
C Calculate Source Term
C
            SRC(I, 6) = SRC(I, 6) - EPSXX
            SRC(I, 7) = SRC(I, 7) - EPSYY
            SRC(I, 8) = SRC(I, 8) - EPSZZ
            SRC(I, 9) = SRC(I, 9) - EPSXY
            SRC(I,10) = SRC(I,10) - EPSXZ
            SRC(I,11) = SRC(I,11) - EPSYZ
  400    CONTINUE
C
C Error
C
      ELSE
         WRITE (IOUT,1000) IEPSIJ
         STOP
      ENDIF
C
C Finished with dissipation source terms
C
      RETURN
 1000 FORMAT (' ','Invalid epsilon_ij model ',I5,' in SRCEIJ.',
     1       /' ','RUN ABORTING')
      END






      SUBROUTINE SRCEIJW (IEPSIJ, NPTS, Q, PROPS, SN, EPSC, YNI, SRCW)
C
C Routine to calculate the wall normal dependent source terms for the 
C dissipation and epsilon equation in the Reynolds stress turbulence model.
C
C This routine is organized to calculate wall normal dependent source terms
C for epsilon_ij in the following sections:
C         Model                                        Wall Terms (YES/NO)
C      --------------------------------------------------------------------
C      0. Isotropic (High Reynolds Number Model)                    No
C      1. Lai-So model                                              YES
C      2. Fu-Launder-Tselepidakis Model                             No
C      3. Gatski-Speziale Algebraic Anisotropic Dissipation Rate    No
C
C IEPSIJ         : Model for epsilon_ij
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C PROPS          : Properties stored at cell centers
C                  PROPS(1) = RMU   molecular viscosity
C                  PROPS(2) = RMUT  turbulent eddy viscosity
C                  PROPS(3) = YPLUS Y+
C                  PROPS(4) = YNORML
C SN             : Wall normals
C EPSC           : Compressible dissipation
C YNI            : Inverse distance to wall
C SRCW           : Wall normal dependent terms for epsilon_ij tensor
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          PROPS (NPTS,NP),
     2          SN    (NPTS,3)
C
      DIMENSION EPSC  (NPTS),
     1          YNI   (NPTS),
     2          SRCW  (NPTS,6)
C
C One-Third
C Two-Thirds
C
      T13 = 1.E0 / 3.E0
      T23 = 2.E0 / 3.E0
C
C Add the source terms for the dissipation and the epsilon equation to SRC
C
C 1. Lai-So model
C    NASA Contractor Report 4369 modifications from Thomas Sommer
C
      IF (IEPSIJ .EQ. IELAIS) THEN
         DO 100 I = 1, NPTS
C
            RHO    = Q(I,1)
            P      = Q(I,5)
            TAUXX  = Q(I,6)
            TAUYY  = Q(I,7)
            TAUZZ  = Q(I,8)
            TAUXY  = Q(I,9)
            TAUXZ  = Q(I,10)
            TAUYZ  = Q(I,11)
            EPSLN  = Q(I,12)
C
            TKE    = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
            RMU    = PROPS(I,1)
            YPLUS  = PROPS(I,3)
C
C Damping
C
            RET    = RHO * TKE * TKE / (RMU * EPSLN) * RE / FSMACH
            FW1    = EXP( -(RET / 150.E0)**2 )
            IF (YPLUS .GT. 250.E0) FW1 = 0.E0
C
C Dissipation
C NOTE: PROBABLY ALSO NEED THE COMPRESSIBILITY CORRECTION MULTIPLYING EPSLN
C       IN THE FW1 TERM
C
            DISS   = T23 * RHO * (EPSLN + EPSC(I))
            TKLNKL =   TAUXX * SN(I,1) * SN(I,1) 
     1               + TAUYY * SN(I,2) * SN(I,2)
     2               + TAUZZ * SN(I,3) * SN(I,3)
     3               + 2.E0 * TAUXY * SN(I,1) * SN(I,2)
     4               + 2.E0 * TAUXZ * SN(I,1) * SN(I,3)
     5               + 2.E0 * TAUYZ * SN(I,2) * SN(I,3)
            DENOMI = 1.E0 / (1.E0 + 1.5E0 * TKLNKL / TKE)
            DISSXX = - FW1 * DISS +
     1               FW1*RHO*EPSLN/TKE*DENOMI * ( TAUXX  
     2+TAUXX*SN(I,1)*SN(I,1)+TAUXY*SN(I,2)*SN(I,1)+TAUXZ*SN(I,3)*SN(I,1)
     3+TAUXX*SN(I,1)*SN(I,1)+TAUXY*SN(I,2)*SN(I,1)+TAUXZ*SN(I,3)*SN(I,1)
     4+TKLNKL*SN(I,1)*SN(I,1) ) 
            DISSYY = - FW1 * DISS +
     1               FW1*RHO*EPSLN/TKE*DENOMI * ( TAUYY 
     2+TAUXY*SN(I,1)*SN(I,2)+TAUYY*SN(I,2)*SN(I,2)+TAUYZ*SN(I,3)*SN(I,2)
     3+TAUXY*SN(I,1)*SN(I,2)+TAUYY*SN(I,2)*SN(I,2)+TAUYZ*SN(I,3)*SN(I,2)
     4+TKLNKL*SN(I,2)*SN(I,2) ) 
            DISSZZ = - FW1 * DISS +
     1               FW1*RHO*EPSLN/TKE*DENOMI * ( TAUZZ
     2+TAUXZ*SN(I,1)*SN(I,3)+TAUYZ*SN(I,2)*SN(I,3)+TAUZZ*SN(I,3)*SN(I,3)
     3+TAUXZ*SN(I,1)*SN(I,3)+TAUYZ*SN(I,2)*SN(I,3)+TAUZZ*SN(I,3)*SN(I,3)
     4+TKLNKL*SN(I,3)*SN(I,3) ) 
            DISSXY = 
     1               FW1*RHO*EPSLN/TKE*DENOMI * ( TAUXY
     2+TAUXX*SN(I,1)*SN(I,2)+TAUXY*SN(I,2)*SN(I,2)+TAUXZ*SN(I,3)*SN(I,2)
     3+TAUXY*SN(I,1)*SN(I,1)+TAUYY*SN(I,2)*SN(I,1)+TAUYZ*SN(I,3)*SN(I,1)
     4+TKLNKL*SN(I,1)*SN(I,2) ) 
            DISSXZ = 
     1               FW1*RHO*EPSLN/TKE*DENOMI * ( TAUXZ
     2+TAUXX*SN(I,1)*SN(I,3)+TAUXY*SN(I,2)*SN(I,3)+TAUXZ*SN(I,3)*SN(I,3)
     3+TAUXZ*SN(I,1)*SN(I,1)+TAUYZ*SN(I,2)*SN(I,1)+TAUZZ*SN(I,3)*SN(I,1)
     4+TKLNKL*SN(I,1)*SN(I,3) ) 
            DISSYZ = 
     1               FW1*RHO*EPSLN/TKE*DENOMI * ( TAUYZ
     2+TAUXY*SN(I,1)*SN(I,3)+TAUYY*SN(I,2)*SN(I,3)+TAUYZ*SN(I,3)*SN(I,3)
     3+TAUXZ*SN(I,1)*SN(I,2)+TAUYZ*SN(I,2)*SN(I,2)+TAUZZ*SN(I,3)*SN(I,2)
     4+TKLNKL*SN(I,2)*SN(I,3) ) 
C
C Calculate Source Term
C
            SRCW(I,1) = SRCW(I,1) - YNI(I) * DISSXX
            SRCW(I,2) = SRCW(I,2) - YNI(I) * DISSYY
            SRCW(I,3) = SRCW(I,3) - YNI(I) * DISSZZ
            SRCW(I,4) = SRCW(I,4) - YNI(I) * DISSXY
            SRCW(I,5) = SRCW(I,5) - YNI(I) * DISSXZ
            SRCW(I,6) = SRCW(I,6) - YNI(I) * DISSYZ
  100    CONTINUE
C
      ENDIF
C
C Finished with wall surface normal dependent dissipation source terms
C
      RETURN
      END






      SUBROUTINE PSTRAN (IPSMDL, NPTS, Q, PROPS, SIJ, WIJ, PIIJ)
C
C Routine to calculate the Pressure-Strain Correlation for a Reynolds stress 
C turbulence model.
C
C This routine is organized into the following sections:
C      1. Gison-Launder Pressure-Strain Correlation
C         Launder-Reece-Rodi "Basic Model" Pressure-Strain Correlation.
C      2. Launder-Reece-Rodi Pressure Strain Correlation
C      3. Lai-So Pressure Strain Correlation
C      4. Speziale-Sarkar-Gatski Pressure Strain Correlation
C      5. Fu-Launder-Tselepidakis (FLT) Pressure Strain Correlation
C      6. So et al. 1995 SSG-based near-wall Pressure Strain Correlation
C
C IPSMDL         : Pressure-strain model
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C PROPS          : Properties stored at cell centers
C                  PROPS(1) = RMU   molecular viscosity
C                  PROPS(2) = RMUT  turbulent eddy viscosity
C                  PROPS(3) = YPLUS Y+
C                  PROPS(4) = YNORML
C SIJ            : Symmetric tensor du_i/dx_j + du_j/dx_i
C WIJ            : Anti-Symmetric tensor du_i/dx_j - du_j/dx_i
C PIIJ           : Pressure Strain Correlation
C                  On entry to this routine, PIIJ contains b_ij
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          PROPS (NPTS,NP)
C
C     +-------------+---------------+---------------+
C     +     SIJ     +       WIJ     +      PIIJ     +
C     + (i,j)  sub  +  (i,j)   sub  +  (i,j)   sub  +
C     +-------------+---------------+---------------+
C     +  x,x    1   +   y,x     1   +   x,x     1   +
C     +  y,y    2   +   z,x     2   +   y,y     2   +
C     +  z,z    3   +   z,y     3   +   z,z     3   +
C     +  x,y    4   +   x,y     4   +   x,y     4   +
C     +  x,z    5   +   x,z     5   +   x,z     5   +
C     +  y,z    6   +   y,z     6   +   y,z     6   +
C     +-------------+---------------+---------------+
C
      DIMENSION SIJ   (NPTS,6),
     1          WIJ   (NPTS,6),
     2          PIIJ  (NPTS,6)
C
C     One-Thirds
C     Two-Thirds
C
      T13 = 1.E0 / 3.E0
      T23 = 2.E0 / 3.E0
C
C     Calculate Pressure-Strain Correlation
C
      IF (IPSMDL .EQ. IPGL) THEN
C
C 1. Gison-Launder Pressure-Strain Correlation
C    Launder-Reece-Rodi "Basic Model" Pressure-Strain Correlation.
C    The form for this model is from Sarkar&Balakrishnan, ICASE 90-18
C
C     Coefficients
C
         CA = PSC2
         CB = PSC2
         CC = T23 * PSC2
         DO 100 I = 1, NPTS
C
            RHO   = Q(I,1)
            P     = Q(I,5)
            TAUXX = Q(I,6)
            TAUYY = Q(I,7)
            TAUZZ = Q(I,8)
            TAUXY = Q(I,9)
            TAUXZ = Q(I,10)
            TAUYZ = Q(I,11)
            EPSLN = Q(I,12)
C
            TKE   = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
C     Anisotropy tensor
C
            BXX    = PIIJ(I,1)
            BYY    = PIIJ(I,2)
            BZZ    = PIIJ(I,3)
            BXY    = PIIJ(I,4)
            BXZ    = PIIJ(I,5)
            BYZ    = PIIJ(I,6)
C
            SXX    = SIJ(I,1)
            SYY    = SIJ(I,2)
            SZZ    = SIJ(I,3)
            SXY    = SIJ(I,4)
            SXZ    = SIJ(I,5)
            SYZ    = SIJ(I,6)
C
            WYX    = WIJ(I,1)
            WZX    = WIJ(I,2)
            WZY    = WIJ(I,3)
            WXY    = WIJ(I,4)
            WXZ    = WIJ(I,5)
            WYZ    = WIJ(I,6)
C
C     S_kk = (S_xx + S_yy + S_zz) / 3
C
            SKK    = (SXX + SYY + SZZ) / 3.E0
C
C     b_lk * S_lk
C
            BLKSLK =       BXX*SXX + BYY*SYY + BZZ*SZZ
     1             + 2.E0*(BXY*SXY + BXZ*SXZ + BYZ*SYZ)
C
C     PI_ij = Phi*_ij,1 + Phi_ij,2 
C
            PIIJ(I,1) = - PSC1 * RHO * EPSLN * BXX
     1 + 2.E0 * RHO * TKE * 
     2   ( CA * (2.E0*(BXX*SXX + BXY*SXY + BXZ*SXZ) - T23*BLKSLK)
     3   + CB * (2.E0*(          BXY*WXY + BXZ*WXZ))
     4   + CC * (SXX - SKK) )
C
            PIIJ(I,2) = - PSC1 * RHO * EPSLN * BYY
     1 + 2.E0 * RHO * TKE * 
     2   ( CA * (2.E0*(BXY*SXY + BYY*SYY + BYZ*SYZ) - T23*BLKSLK)
     3   + CB * (2.E0*(BXY*WYX           + BYZ*WYZ))
     4   + CC * (SYY - SKK) )
C
            PIIJ(I,3) = - PSC1 * RHO * EPSLN * BZZ
     1 + 2.E0 * RHO * TKE * 
     2   ( CA * (2.E0*(BXZ*SXZ + BYZ*SYZ + BZZ*SZZ) - T23*BLKSLK)
     3   + CB * (2.E0*(BXZ*WZX + BYZ*WZY          ))
     4   + CC * (SZZ - SKK) )
C
            PIIJ(I,4) = - PSC1 * RHO * EPSLN * BXY 
     1 + 2.E0 * RHO * TKE * 
     2   (CA*(BXX*SXY + BXY*SYY + BXZ*SYZ + BXY*SXX + BYY*SXY + BYZ*SXZ)
     3   +CB*(BXX*WYX           + BXZ*WYZ           + BYY*WXY + BYZ*WXZ)
     4   +CC*(SXY) )
C
            PIIJ(I,5) = - PSC1 * RHO * EPSLN * BXZ 
     1 + 2.E0 * RHO * TKE * 
     2   (CA*(BXX*SXZ + BXY*SYZ + BXZ*SZZ + BXZ*SXX + BYZ*SXY + BZZ*SXZ)
     3   +CB*(BXX*WZX + BXY*WZY                     + BYZ*WXY + BZZ*WXZ)
     4   +CC*(SXZ) )
C
            PIIJ(I,6) = - PSC1 * RHO * EPSLN * BYZ 
     1 + 2.E0 * RHO * TKE * 
     2   (CA*(BXY*SXZ + BYY*SYZ + BYZ*SZZ + BXZ*SXY + BYZ*SYY + BZZ*SYZ)
     3   +CB*(BXY*WZX + BYY*WZY           + BXZ*WYX           + BZZ*WYZ)
     4   +CC*(SYZ) )
C
  100    CONTINUE
C
C 2. Launder-Reece-Rodi Pressure Strain Correlation
C    Cphi are coefficients casting the LRR form of the model into b,S,W.
C    Use PSCW = 1.e0 to turn wall reflection term ON,
C        PSCW = 0.e0 to turn wall reflection term OFF
C
      ELSE IF (IPSMDL .EQ. IPLRR1) THEN
         CPHI2A = ( 9.E0*PSC2 + 6.E0) / 11.E0
         CPHI2B = (10.E0 - 7.E0*PSC2) / 11.E0
         CPHI2C = 0.4E0
         DO 200 I = 1, NPTS
C
            RHO    = Q(I,1)
            P      = Q(I,5)
            TAUXX  = Q(I,6)
            TAUYY  = Q(I,7)
            TAUZZ  = Q(I,8)
            TAUXY  = Q(I,9)
            TAUXZ  = Q(I,10)
            TAUYZ  = Q(I,11)
            EPSLN  = Q(I,12)
C
            TKE    = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
            YNORML = ABS (PROPS(I,4))
C
C     Pressure Reflection term
C
            RLEN   = PSCW * TKE * SQRT(TKE) / (EPSLN * YNORML)
C
C     Anisotropy tensor
C
            BXX    = PIIJ(I,1)
            BYY    = PIIJ(I,2)
            BZZ    = PIIJ(I,3)
            BXY    = PIIJ(I,4)
            BXZ    = PIIJ(I,5)
            BYZ    = PIIJ(I,6)
C
            SXX    = SIJ(I,1)
            SYY    = SIJ(I,2)
            SZZ    = SIJ(I,3)
            SXY    = SIJ(I,4)
            SXZ    = SIJ(I,5)
            SYZ    = SIJ(I,6)
C
            WYX    = WIJ(I,1)
            WZX    = WIJ(I,2)
            WZY    = WIJ(I,3)
            WXY    = WIJ(I,4)
            WXZ    = WIJ(I,5)
            WYZ    = WIJ(I,6)
C
C     S_kk = (S_xx + S_yy + S_zz) / 3
C
            SKK    = (SXX + SYY + SZZ) / 3.E0
C
C     b_lk * S_lk
C
            BLKSLK =       BXX*SXX + BYY*SYY + BZZ*SZZ
     1             + 2.E0*(BXY*SXY + BXZ*SXZ + BYZ*SYZ)
C
C     Coefficients
C
            C1STAR = PSC1 - 2.E0 * 0.125E0 * RLEN
C
            CA     = CPHI2A 
            CB     = CPHI2B - 2.E0 * 0.015E0 * RLEN
            CC     = CPHI2C 
C
C     PI_ij = Phi_ij,1 + Phi_ij,2 
C
            PIIJ(I,1) = - C1STAR * RHO * EPSLN * BXX
     1 + 2.E0 * RHO * TKE * 
     2   ( CA * (2.E0*(BXX*SXX + BXY*SXY + BXZ*SXZ) - T23*BLKSLK)
     3   + CB * (2.E0*(          BXY*WXY + BXZ*WXZ))
     4   + CC * (SXX - SKK) )
C
            PIIJ(I,2) = - C1STAR * RHO * EPSLN * BYY
     1 + 2.E0 * RHO * TKE * 
     2   ( CA * (2.E0*(BXY*SXY + BYY*SYY + BYZ*SYZ) - T23*BLKSLK)
     3   + CB * (2.E0*(BXY*WYX           + BYZ*WYZ))
     4   + CC * (SYY - SKK) )
C
            PIIJ(I,3) = - C1STAR * RHO * EPSLN * BZZ
     1 + 2.E0 * RHO * TKE * 
     2   ( CA * (2.E0*(BXZ*SXZ + BYZ*SYZ + BZZ*SZZ) - T23*BLKSLK)
     3   + CB * (2.E0*(BXZ*WZX + BYZ*WZY          ))
     4   + CC * (SZZ - SKK) )
C
            PIIJ(I,4) = - C1STAR * RHO * EPSLN * BXY 
     1 + 2.E0 * RHO * TKE * 
     2   (CA*(BXX*SXY + BXY*SYY + BXZ*SYZ + BXY*SXX + BYY*SXY + BYZ*SXZ)
     3   +CB*(BXX*WYX           + BXZ*WYZ           + BYY*WXY + BYZ*WXZ)
     4   +CC*(SXY) )
C
            PIIJ(I,5) = - C1STAR * RHO * EPSLN * BXZ 
     1 + 2.E0 * RHO * TKE * 
     2   (CA*(BXX*SXZ + BXY*SYZ + BXZ*SZZ + BXZ*SXX + BYZ*SXY + BZZ*SXZ)
     3   +CB*(BXX*WZX + BXY*WZY                     + BYZ*WXY + BZZ*WXZ)
     4   +CC*(SXZ) )
C
            PIIJ(I,6) = - C1STAR * RHO * EPSLN * BYZ 
     1 + 2.E0 * RHO * TKE * 
     2   (CA*(BXY*SXZ + BYY*SYZ + BYZ*SZZ + BXZ*SXY + BYZ*SYY + BZZ*SYZ)
     3   +CB*(BXY*WZX + BYY*WZY           + BXZ*WYX           + BZZ*WYZ)
     4   +CC*(SYZ) )
C
  200    CONTINUE
C
C 3. Lai-So Pressure Strain Correlation
C    NASA Contractor Report 4369 (Based on LRR1)
C    Corrected form based on Zhang, So, Gatski, Speziale in the 
C    Near-Wall Turbulent Flows book later corrected from private communications
C    with Ron So.  Note that this is not strictly a pressure-strain model
C    but rather a model of pressure-strain and pressure-diffusion so that the
C    trace of the model is Non-Zero.
C    Modified based on private communication with Hou-Shun Zhang to include
C    the subtraction of the trace term from the S_ij term in the near-wall.
C
      ELSE IF (IPSMDL .EQ. IPLS) THEN
         CPHI2A = ( 9.E0*PSC2 + 6.E0) / 11.E0
         CPHI2B = (10.E0 - 7.E0*PSC2) / 11.E0
         CPHI2C = 0.4E0
         DO 300 I = 1, NPTS
C
            RHO    = Q(I,1)
            P      = Q(I,5)
            TAUXX  = Q(I,6)
            TAUYY  = Q(I,7)
            TAUZZ  = Q(I,8)
            TAUXY  = Q(I,9)
            TAUXZ  = Q(I,10)
            TAUYZ  = Q(I,11)
            EPSLN  = Q(I,12)
C
            TKE    = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
            RMU    = PROPS(I,1)
            YPLUS  = PROPS(I,3)
            YNORML = ABS (PROPS(I,4))
C
C     Anisotropy tensor
C
            BXX    = PIIJ(I,1)
            BYY    = PIIJ(I,2)
            BZZ    = PIIJ(I,3)
            BXY    = PIIJ(I,4)
            BXZ    = PIIJ(I,5)
            BYZ    = PIIJ(I,6)
C
            SXX    = SIJ(I,1)
            SYY    = SIJ(I,2)
            SZZ    = SIJ(I,3)
            SXY    = SIJ(I,4)
            SXZ    = SIJ(I,5)
            SYZ    = SIJ(I,6)
C
            WYX    = WIJ(I,1)
            WZX    = WIJ(I,2)
            WZY    = WIJ(I,3)
            WXY    = WIJ(I,4)
            WXZ    = WIJ(I,5)
            WYZ    = WIJ(I,6)
C
C     S_kk = (S_xx + S_yy + S_zz) / 3
C
            SKK    = (SXX + SYY + SZZ) / 3.E0
C
C     b_lk * S_lk
C
            BLKSLK =       BXX*SXX + BYY*SYY + BZZ*SZZ
     1             + 2.E0*(BXY*SXY + BXZ*SXZ + BYZ*SYZ)
C
C     Wall damping
C
            RET    = RHO * TKE * TKE / (RMU * EPSLN) * RE / FSMACH
            FW1    = EXP( -(RET / 150.E0)**2 )
            IF (YPLUS .GT. 250.E0) FW1 = 0.E0
            C1STAR = PSC1 * (1.E0 - FW1)
C
C     Coefficients
C
            CA     = CPHI2A - (PSALFA + PSBETA)*FW1
            CB     = CPHI2B - (PSALFA - PSBETA)*FW1
            CC     = CPHI2C - (T23*(PSALFA + PSBETA) - PSGAMA)*FW1
C
            CWALL  = 2.E0*PSCW*RHO*TKE*TKE*SQRT(TKE) / (EPSLN*YNORML)
C
C     PI_ij = Phi_ij,1 + Phi_ij,2 + Phi_ij,w + Phi_ij,R
C     Note: SN (wall normal) dependent portion of Phi_ij,w is now in PSTRANW
C
            PIIJ(I,1) = - C1STAR * RHO * EPSLN * BXX
     1 + 2.E0 * RHO * TKE * 
     2   ( CA * (2.E0*(BXX*SXX + BXY*SXY + BXZ*SXZ) - T23*BLKSLK)
     3   + CB * (2.E0*(          BXY*WXY + BXZ*WXZ))
     4   + CC * (SXX - SKK) )
     5 + CWALL * (SXX - SKK)
C
            PIIJ(I,2) = - C1STAR * RHO * EPSLN * BYY
     1 + 2.E0 * RHO * TKE * 
     2   ( CA * (2.E0*(BXY*SXY + BYY*SYY + BYZ*SYZ) - T23*BLKSLK)
     3   + CB * (2.E0*(BXY*WYX           + BYZ*WYZ))
     4   + CC * (SYY - SKK) )
     5 + CWALL * (SYY - SKK)
C
            PIIJ(I,3) = - C1STAR * RHO * EPSLN * BZZ
     1 + 2.E0 * RHO * TKE * 
     2   ( CA * (2.E0*(BXZ*SXZ + BYZ*SYZ + BZZ*SZZ) - T23*BLKSLK)
     3   + CB * (2.E0*(BXZ*WZX + BYZ*WZY          ))
     4   + CC * (SZZ - SKK) )
     5 + CWALL * (SZZ - SKK)
C
            PIIJ(I,4) = - C1STAR * RHO * EPSLN * BXY 
     1 + 2.E0 * RHO * TKE * 
     2   (CA*(BXX*SXY + BXY*SYY + BXZ*SYZ + BXY*SXX + BYY*SXY + BYZ*SXZ)
     3   +CB*(BXX*WYX           + BXZ*WYZ           + BYY*WXY + BYZ*WXZ)
     4   +CC*(SXY) )
     5 + CWALL * SXY
C
            PIIJ(I,5) = - C1STAR * RHO * EPSLN * BXZ 
     1 + 2.E0 * RHO * TKE * 
     2   (CA*(BXX*SXZ + BXY*SYZ + BXZ*SZZ + BXZ*SXX + BYZ*SXY + BZZ*SXZ)
     3   +CB*(BXX*WZX + BXY*WZY                     + BYZ*WXY + BZZ*WXZ)
     4   +CC*(SXZ) )
     5 + CWALL * SXZ
C
            PIIJ(I,6) = - C1STAR * RHO * EPSLN * BYZ 
     1 + 2.E0 * RHO * TKE * 
     2   (CA*(BXY*SXZ + BYY*SYZ + BYZ*SZZ + BXZ*SXY + BYZ*SYY + BZZ*SYZ)
     3   +CB*(BXY*WZX + BYY*WZY           + BXZ*WYX           + BZZ*WYZ)
     4   +CC*(SYZ) )
     5 + CWALL * SYZ
C
  300    CONTINUE
C
C 4. Speziale-Sarkar-Gatski Pressure Strain Correlation
C
      ELSE IF (IPSMDL .EQ. IPSSG) THEN
         CA0 =   3.40E0
         CA1 =   4.20E0
         CA3 = - 1.80E0
         CA4 =   1.25E0
         CA5 =   0.40E0
         DO 400 I = 1, NPTS
C
            RHO    = Q(I,1)
            P      = Q(I,5)
            TAUXX  = Q(I,6)
            TAUYY  = Q(I,7)
            TAUZZ  = Q(I,8)
            TAUXY  = Q(I,9)
            TAUXZ  = Q(I,10)
            TAUYZ  = Q(I,11)
            EPSLN  = Q(I,12)
C
            TKE    = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
            YNORML = ABS (PROPS(I,4))
C
C     Anisotropy tensor
C
            BXX    = PIIJ(I,1)
            BYY    = PIIJ(I,2)
            BZZ    = PIIJ(I,3)
            BXY    = PIIJ(I,4)
            BXZ    = PIIJ(I,5)
            BYZ    = PIIJ(I,6)
C
            SXX    = SIJ(I,1)
            SYY    = SIJ(I,2)
            SZZ    = SIJ(I,3)
            SXY    = SIJ(I,4)
            SXZ    = SIJ(I,5)
            SYZ    = SIJ(I,6)
C
            WYX    = WIJ(I,1)
            WZX    = WIJ(I,2)
            WZY    = WIJ(I,3)
            WXY    = WIJ(I,4)
            WXZ    = WIJ(I,5)
            WYZ    = WIJ(I,6)
C
C     S_kk = (S_xx + S_yy + S_zz) / 3
C
            SKK    = (SXX + SYY + SZZ) / 3.E0
C
C     b_lk * S_lk
C
            BLKSLK =       BXX*SXX + BYY*SYY + BZZ*SZZ
     1             + 2.E0*(BXY*SXY + BXZ*SXZ + BYZ*SYZ)
C
C     II = b_mn b_nm
C
            BMNBMN =       BXX*BXX + BYY*BYY + BZZ*BZZ
     1             + 2.E0*(BXY*BXY + BXZ*BXZ + BYZ*BYZ)
C
C     Coefficient for SSG
C
            CA2    =   0.80E0 - 1.3E0 * SQRT(BMNBMN)
C
C     P_k => Production of turbulent kinetic energy
C
            PK     =   - ( TAUXX*SXX + TAUYY*SYY + TAUZZ*SZZ
     1             + 2.E0*(TAUXY*SXY + TAUXZ*SXZ + TAUYZ*SYZ) )
C
C     PI_ij = Phi_ij,1 + Phi_ij,2 
C
            PIIJ(I,1) = - CA0 * RHO * EPSLN * BXX
     1 + CA1 * RHO * EPSLN * (BXX*BXX + BXY*BXY + BXZ*BXZ - T13*BMNBMN)
     2 + RHO * TKE * 
     3  ( CA4 * (2.E0*(BXX*SXX + BXY*SXY + BXZ*SXZ) - T23*BLKSLK)
     4  + CA5 * (2.E0*(          BXY*WXY + BXZ*WXZ))
     5  + CA2 * (SXX - SKK) )
     6 + CA3 * RHO * PK * BXX
C
            PIIJ(I,2) = - CA0 * RHO * EPSLN * BYY
     1 + CA1 * RHO * EPSLN * (BXY*BXY + BYY*BYY + BYZ*BYZ - T13*BMNBMN)
     2 + RHO * TKE * 
     3  ( CA4 * (2.E0*(BXY*SXY + BYY*SYY + BYZ*SYZ) - T23*BLKSLK)
     4  + CA5 * (2.E0*(BXY*WYX           + BYZ*WYZ))
     5  + CA2 * (SYY - SKK) )
     6 + CA3 * RHO * PK * BYY
C
            PIIJ(I,3) = - CA0 * RHO * EPSLN * BZZ
     1 + CA1 * RHO * EPSLN * (BXZ*BXZ + BYZ*BYZ + BZZ*BZZ - T13*BMNBMN)
     2 + RHO * TKE * 
     3  ( CA4 * (2.E0*(BXZ*SXZ + BYZ*SYZ + BZZ*SZZ) - T23*BLKSLK)
     4  + CA5 * (2.E0*(BXZ*WZX + BYZ*WZY          ))
     5  + CA2 * (SZZ - SKK) )
     6 + CA3 * RHO * PK * BZZ
C
            PIIJ(I,4) = - CA0 * RHO * EPSLN * BXY 
     1 + CA1 * RHO * EPSLN * (BXX*BXY + BXY*BYY + BXZ*BYZ)
     2 + RHO * TKE * 
     3  (CA4*(BXX*SXY + BXY*SYY + BXZ*SYZ + BXY*SXX + BYY*SXY + BYZ*SXZ)
     4  +CA5*(BXX*WYX           + BXZ*WYZ           + BYY*WXY + BYZ*WXZ)
     5  +CA2*(SXY) )
     6 + CA3 * RHO * PK * BXY
C
            PIIJ(I,5) = - CA0 * RHO * EPSLN * BXZ 
     1 + CA1 * RHO * EPSLN * (BXX*BXZ + BXY*BYZ + BXZ*BZZ)
     2 + RHO * TKE * 
     3  (CA4*(BXX*SXZ + BXY*SYZ + BXZ*SZZ + BXZ*SXX + BYZ*SXY + BZZ*SXZ)
     4  +CA5*(BXX*WZX + BXY*WZY                     + BYZ*WXY + BZZ*WXZ)
     5  +CA2*(SXZ) )
     6 + CA3 * RHO * PK * BXZ
C
            PIIJ(I,6) = - CA0 * RHO * EPSLN * BYZ 
     1 + CA1 * RHO * EPSLN * (BXY*BXZ + BYY*BYZ + BYZ*BZZ)
     2 + RHO * TKE * 
     3  (CA4*(BXY*SXZ + BYY*SYZ + BYZ*SZZ + BXZ*SXY + BYZ*SYY + BZZ*SYZ)
     4  +CA5*(BXY*WZX + BYY*WZY           + BXZ*WYX           + BZZ*WYZ)
     5  +CA2*(SYZ) )
     6 + CA3 * RHO * PK * BYZ
C
  400    CONTINUE
C
C 5. Fu-Launder-Tselepidakis (FLT) Pressure Strain Correlation
C
      ELSE IF (IPSMDL .EQ. IPFLT) THEN
         CA2 =   0.80E0
         CA3 =   1.20E0
         CA4 =   1.20E0
         CA5 =   1.73E0
         CA6 =   0.80E0
         CA7 =   0.80E0
         CA8 =  11.20E0
         DO 500 I = 1, NPTS
C
            RHO    = Q(I,1)
            P      = Q(I,5)
            TAUXX  = Q(I,6)
            TAUYY  = Q(I,7)
            TAUZZ  = Q(I,8)
            TAUXY  = Q(I,9)
            TAUXZ  = Q(I,10)
            TAUYZ  = Q(I,11)
            EPSLN  = Q(I,12)
C
            TKE    = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
C     Anisotropy tensor
C
            BXX    = PIIJ(I,1)
            BYY    = PIIJ(I,2)
            BZZ    = PIIJ(I,3)
            BXY    = PIIJ(I,4)
            BXZ    = PIIJ(I,5)
            BYZ    = PIIJ(I,6)
C
            SXX    = SIJ(I,1)
            SYY    = SIJ(I,2)
            SZZ    = SIJ(I,3)
            SXY    = SIJ(I,4)
            SXZ    = SIJ(I,5)
            SYZ    = SIJ(I,6)
C
            WYX    = WIJ(I,1)
            WZX    = WIJ(I,2)
            WZY    = WIJ(I,3)
            WXY    = WIJ(I,4)
            WXZ    = WIJ(I,5)
            WYZ    = WIJ(I,6)
C
C     S_kk = (S_xx + S_yy + S_zz) / 3
C
            SKK    = (SXX + SYY + SZZ) / 3.E0
C
C     b_lk * S_lk
C
            BLKSLK =       BXX*SXX + BYY*SYY + BZZ*SZZ
     1             + 2.E0*(BXY*SXY + BXZ*SXZ + BYZ*SYZ)
C
C     II = b_mn b_nm
C
            BMNBMN =       BXX*BXX + BYY*BYY + BZZ*BZZ
     1             + 2.E0*(BXY*BXY + BXZ*BXZ + BYZ*BYZ)
C
C     III = b_lm b_mn b_nl
C
            BIII   =      BXX*BXX*BXX + BYY*BYY*BYY + BZZ*BZZ*BZZ
     1             + 3.E0*BXY*BXY*(BXX + BYY)
     2             + 3.E0*BXZ*BXZ*(BXX + BZZ)
     3             + 3.E0*BYZ*BYZ*(BYY + BZZ)
     4             + 6.E0*BXY*BXZ*BYZ
C
C     Coefficients for FLT
C
            FSQRT  =   SQRT (1.E0 - 4.5E0*BMNBMN + 9.E0*BIII)
            CA0    =   60.E0 * BMNBMN * FSQRT
            CA1    = - 72.E0 * BMNBMN * FSQRT
C
C     P_k => Production of turbulent kinetic energy
C
            PK     =   - ( TAUXX*SXX + TAUYY*SYY + TAUZZ*SZZ
     1             + 2.E0*(TAUXY*SXY + TAUXZ*SXZ + TAUYZ*SYZ) )
C
            BW11   =   2.E0 * (BXY*WXY + BXZ*WXZ)
            BW22   =   2.E0 * (BYZ*WYZ - BXY*WXY)
            BW33   = - 2.E0 * (BXZ*WXZ + BYZ*WYZ)
            BW12   = (BYY-BXX)*WXY + BXZ*WYZ + BYZ*WXZ
            BW13   = (BZZ-BXX)*WXZ - BXY*WYZ + BYZ*WXY
            BW23   = (BZZ-BYY)*WYZ - BXY*WXZ - BXZ*WXY
C
C     PI_ij = Phi_ij,1 + Phi_ij,2 
C
            PIIJ(I,1) = - CA0 * RHO * EPSLN * BXX
     1 + CA1 * RHO * EPSLN * (BXX*BXX + BXY*BXY + BXZ*BXZ - T13*BMNBMN)
     2 + CA2 * RHO * TKE * (SXX - SKK) 
     3 + CA3 * RHO * PK  * BXX
     4 + RHO * TKE * 
     5  ( CA4 * (2.E0*(BXX*SXX + BXY*SXY + BXZ*SXZ) - T23*BLKSLK)
     6  + CA5 * (2.E0*(          BXY*WXY + BXZ*WXZ))
     7  + CA6 * ( BXX*2.E0*(BXX*SXX)
     8          + BYY*2.E0*(BXY*SXY)
     9          + BZZ*2.E0*(BXZ*SXZ)
     A          + BXY*2.E0*(BXX*SXY + BXY*SXX)
     B          + BXZ*2.E0*(BXX*SXZ + BXZ*SXX)
     C          + BYZ*2.E0*(BXY*SXZ + BXZ*SXY)
     D          - 2.E0 * (SXX*BXX*BXX + SYY*BXY*BXY + SZZ*BXZ*BXZ
     E                    + SXY*2.E0*(BXY*BXX)
     F                    + SXZ*2.E0*(BXZ*BXX)
     G                    + SYZ*2.E0*(BXZ*BXY)) )
     H  + CA7 * ( 
     I          + BYY*2.E0*(BXY*WXY)
     J          + BZZ*2.E0*(BXZ*WXZ)
     K          + BXY*2.E0*(BXX*WXY          )
     L          + BXZ*2.E0*(BXX*WXZ          )
     M          + BYZ*2.E0*(BXY*WXZ + BXZ*WXY) )
     N  + CA8 * ( BMNBMN
     O          *2.E0*(          BXY*WXY + BXZ*WXZ)
     P          +3.E0*( BW11*BXX*BXX + BW22*BXY*BXY + BW33*BXZ*BXZ
     Q                + BW12*(2.E0*BXX*BXY) + BW13*(2.E0*BXX*BXZ)
     R                + BW23*(2.E0*BXY*BXZ) ) )
     S  )
C
            PIIJ(I,2) = - CA0 * RHO * EPSLN * BYY
     1 + CA1 * RHO * EPSLN * (BXY*BXY + BYY*BYY + BYZ*BYZ - T13*BMNBMN)
     2 + CA2 * RHO * TKE * (SYY - SKK) 
     3 + CA3 * RHO * PK  * BYY
     4 + RHO * TKE * 
     5  ( CA4 * (2.E0*(BXY*SXY + BYY*SYY + BYZ*SYZ) - T23*BLKSLK)
     6  + CA5 * (2.E0*(BXY*WYX           + BYZ*WYZ))
     7  + CA6 * ( BXX*2.E0*(BXY*SXY)
     8          + BYY*2.E0*(BYY*SYY)
     9          + BZZ*2.E0*(BYZ*SYZ)
     A          + BXY*2.E0*(BXY*SYY + BYY*SXY)
     B          + BXZ*2.E0*(BXY*SYZ + BYZ*SXY)
     C          + BYZ*2.E0*(BYY*SYZ + BYZ*SYY)
     D          - 2.E0 * (SXX*BXY*BXY + SYY*BYY*BYY + SZZ*BYZ*BYZ
     E                    + SXY*2.E0*(BXY*BYY)
     F                    + SXZ*2.E0*(BXY*BYZ)
     G                    + SYZ*2.E0*(BYY*BYZ)) )
     H  + CA7 * ( BXX*2.E0*(BXY*WYX)
     I          
     J          + BZZ*2.E0*(BYZ*WYZ)
     K          + BXY*2.E0*(          BYY*WYX)
     L          + BXZ*2.E0*(BXY*WYZ + BYZ*WYX)
     M          + BYZ*2.E0*(BYY*WYZ          ) )
     N  + CA8 * ( BMNBMN
     O          *2.E0*(BXY*WYX           + BYZ*WYZ)
     P          +3.E0*( BW11*BXY*BXY + BW22*BYY*BYY + BW33*BYZ*BYZ
     Q                + BW12*(2.E0*BXY*BYY) + BW13*(2.E0*BXY*BYZ)
     R                + BW23*(2.E0*BYY*BYZ) ) )
     S  )
C
            PIIJ(I,3) = - CA0 * RHO * EPSLN * BZZ
     1 + CA1 * RHO * EPSLN * (BXZ*BXZ + BYZ*BYZ + BZZ*BZZ - T13*BMNBMN)
     2 + CA2 * RHO * TKE * (SZZ - SKK) 
     3 + CA3 * RHO * PK  * BZZ
     4 + RHO * TKE * 
     5  ( CA4 * (2.E0*(BXZ*SXZ + BYZ*SYZ + BZZ*SZZ) - T23*BLKSLK)
     6  + CA5 * (2.E0*(BXZ*WXZ + BYZ*WZY          ))
     7  + CA6 * ( BXX*2.E0*(BXZ*SXZ)
     8          + BYY*2.E0*(BYZ*SYZ)
     9          + BZZ*2.E0*(BZZ*SZZ)
     A          + BXY*2.E0*(BXZ*SYZ + BYZ*SXZ)
     B          + BXZ*2.E0*(BXZ*SZZ + BZZ*SXZ)
     C          + BYZ*2.E0*(BYZ*SZZ + BZZ*SYZ)
     D          - 2.E0 * (SXX*BXZ*BXZ + SYY*BYZ*BYZ + SZZ*BZZ*BZZ
     E                    + SXY*2.E0*(BXZ*BYZ)
     F                    + SXZ*2.E0*(BXZ*BZZ)
     G                    + SYZ*2.E0*(BYZ*BZZ)) )
     H  + CA7 * ( BXX*2.E0*(BXZ*WZX)
     I          + BYY*2.E0*(BYZ*WZY)
     J          
     K          + BXY*2.E0*(BXZ*WZY + BYZ*WZX)
     L          + BXZ*2.E0*(          BZZ*WZX)
     M          + BYZ*2.E0*(          BZZ*WZY) )
     N  + CA8 * ( BMNBMN
     O          *2.E0*(BXZ*WZX + BYZ*WZY          )
     P          +3.E0*( BW11*BXZ*BXZ + BW22*BYZ*BYZ + BW33*BZZ*BZZ
     Q                + BW12*(2.E0*BXZ*BYZ) + BW13*(2.E0*BXZ*BZZ)
     R                + BW23*(2.E0*BYZ*BZZ) ) )
     S  )
C
            PIIJ(I,4) = - CA0 * RHO * EPSLN * BXY 
     1 + CA1 * RHO * EPSLN * (BXX*BXY + BXY*BYY + BXZ*BYZ)
     2 + CA2 * RHO * TKE * SXY
     3 + CA3 * RHO * PK  * BXY
     4 + RHO * TKE * 
     5  (CA4*(BXX*SXY + BXY*SYY + BXZ*SYZ + BXY*SXX + BYY*SXY + BYZ*SXZ)
     6  +CA5*(BXX*WYX           + BXZ*WYZ           + BYY*WXY + BYZ*WXZ)
     7  + CA6 * ( BXX*(BXX*SXY + BXY*SXX)
     8          + BYY*(BXY*SYY + BYY*SXY)
     9          + BZZ*(BXZ*SYZ + BYZ*SXZ)
     A          + BXY*(BXX*SYY + BXY*SXY + BXY*SXY + BYY*SXX)
     B          + BXZ*(BXX*SYZ + BXZ*SXY + BXY*SXZ + BYZ*SXX)
     C          + BYZ*(BXY*SYZ + BXZ*SYY + BYY*SXZ + BYZ*SXY)
     D          - 2.E0 * (SXX*BXY*BXX + SYY*BYY*BXY + SZZ*BYZ*BXZ
     E                    + SXY*(BYY*BXX+BXY*BXY)
     F                    + SXZ*(BYZ*BXX+BXY*BXZ)
     G                    + SYZ*(BYZ*BXY+BYY*BXZ)) )
     H  + CA7 * ( BXX*(BXX*WYX          )
     I          + BYY*(          BYY*WXY)
     J          + BZZ*(BXZ*WYZ + BYZ*WXZ)
     K          + BXY*(          BXY*WYX + BXY*WXY          )
     L          + BXZ*(BXX*WYZ + BXZ*WYX + BXY*WXZ          )
     M          + BYZ*(BXY*WYZ +           BYY*WXZ + BYZ*WXY) )
     N  + CA8 * ( BMNBMN
     O      *(BXX*WYX           + BXZ*WYZ           + BYY*WXY + BYZ*WXZ)
     P          +3.E0*( BW11*BXX*BXY + BW22*BXY*BYY + BW33*BXZ*BYZ
     Q                + BW12*(BXX*BYY+BXY*BXY) + BW13*(BXX*BYZ+BXZ*BXY)
     R                + BW23*(BXY*BYZ+BXZ*BYY) ) )
     S  )
C
            PIIJ(I,5) = - CA0 * RHO * EPSLN * BXZ 
     1 + CA1 * RHO * EPSLN * (BXX*BXZ + BXY*BYZ + BXZ*BZZ)
     2 + CA2 * RHO * TKE * SXZ
     3 + CA3 * RHO * PK  * BXZ
     4 + RHO * TKE * 
     5  (CA4*(BXX*SXZ + BXY*SYZ + BXZ*SZZ + BXZ*SXX + BYZ*SXY + BZZ*SXZ)
     6  +CA5*(BXX*WZX + BXY*WZY                     + BYZ*WXY + BZZ*WXZ)
     7  + CA6 * ( BXX*(BXX*SXZ + BXZ*SXX)
     8          + BYY*(BXY*SYZ + BYZ*SXY)
     9          + BZZ*(BXZ*SZZ + BZZ*SXZ)
     A          + BXY*(BXX*SYZ + BXY*SXZ + BXZ*SXY + BYZ*SXX)
     B          + BXZ*(BXX*SZZ + BXZ*SXZ + BXZ*SXZ + BZZ*SXX)
     C          + BYZ*(BXY*SZZ + BXZ*SYZ + BYZ*SXZ + BZZ*SXY)
     D          - 2.E0 * (SXX*BXZ*BXX + SYY*BYZ*BXY + SZZ*BZZ*BXZ
     E                    + SXY*(BYZ*BXX+BXZ*BXY)
     F                    + SXZ*(BZZ*BXX+BXZ*BXZ)
     G                    + SYZ*(BZZ*BXY+BYZ*BXZ)) )
     H  + CA7 * ( BXX*(BXX*WZX          )
     I          + BYY*(BXY*WZY + BYZ*WXY)
     J          + BZZ*(          BZZ*WXZ)
     K          + BXY*(BXX*WZY + BXY*WZX + BXZ*WXY          )
     L          + BXZ*(          BXZ*WZX + BXZ*WXZ          )
     M          + BYZ*(          BXZ*WZY + BYZ*WXZ + BZZ*WXY) )
     N  + CA8 * ( BMNBMN
     O      *(BXX*WZX + BXY*WZY                     + BYZ*WXY + BZZ*WXZ)
     P          +3.E0*( BW11*BXX*BXZ + BW22*BXY*BYZ + BW33*BXZ*BZZ
     Q                + BW12*(BXX*BYZ+BXY*BXZ) + BW13*(BXX*BZZ+BXZ*BXZ)
     R                + BW23*(BXY*BZZ+BXZ*BYZ) ) )
     S  )
C
            PIIJ(I,6) = - CA0 * RHO * EPSLN * BYZ 
     1 + CA1 * RHO * EPSLN * (BXY*BXZ + BYY*BYZ + BYZ*BZZ)
     2 + CA2 * RHO * TKE * SYZ 
     3 + CA3 * RHO * PK  * BYZ
     4 + RHO * TKE * 
     5  (CA4*(BXY*SXZ + BYY*SYZ + BYZ*SZZ + BXZ*SXY + BYZ*SYY + BZZ*SYZ)
     6  +CA5*(BXY*WZX + BYY*WZY           + BXZ*WYX           + BZZ*WYZ)
     7  + CA6 * ( BXX*(BXY*SXZ + BXZ*SXY)
     8          + BYY*(BYY*SYZ + BYZ*SYY)
     9          + BZZ*(BYZ*SZZ + BZZ*SYZ)
     A          + BXY*(BXY*SYZ + BYY*SXZ + BXZ*SYY + BYZ*SXY)
     B          + BXZ*(BXY*SZZ + BYZ*SXZ + BXZ*SYZ + BZZ*SXY)
     C          + BYZ*(BYY*SZZ + BYZ*SYZ + BYZ*SYZ + BZZ*SYY)
     D          - 2.E0 * (SXX*BXZ*BXY + SYY*BYZ*BYY + SZZ*BZZ*BYZ
     E                    + SXY*(BYZ*BXY+BXZ*BYY)
     F                    + SXZ*(BZZ*BXY+BXZ*BYZ)
     G                    + SYZ*(BZZ*BYY+BYZ*BYZ)) )
     H  + CA7 * ( BXX*(BXY*WZX + BXZ*WYX)
     I          + BYY*(BYY*WZY          )
     J          + BZZ*(          BZZ*WYZ)
     K          + BXY*(BXY*WZY + BYY*WZX +           BYZ*WYX)
     L          + BXZ*(          BYZ*WZX + BXZ*WYZ + BZZ*WYX)
     M          + BYZ*(          BYZ*WZY + BYZ*WYZ          ) )
     N  + CA8 * ( BMNBMN
     O      *(BXY*WZX + BYY*WZY           + BXZ*WYX           + BZZ*WYZ)
     P          +3.E0*( BW11*BXY*BXZ + BW22*BYY*BYZ + BW33*BYZ*BZZ
     Q                + BW12*(BXY*BYZ+BYY*BXZ) + BW13*(BXY*BZZ+BYZ*BXZ)
     R                + BW23*(BYY*BZZ+BYZ*BYZ) ) )
     S  )
C
  500    CONTINUE
C
C 6. So et al. 1995 SSG-based near-wall Pressure Strain Correlation
C
      ELSE IF (IPSMDL .EQ. IPSO95) THEN
         CA0HR =   3.40E0
         CA1HR =   4.20E0
         CA3HR = - 1.80E0
         CA4HR =   1.25E0
         CA5HR =   0.40E0
         DO 600 I = 1, NPTS
C
            RHO    = Q(I,1)
            P      = Q(I,5)
            TAUXX  = Q(I,6)
            TAUYY  = Q(I,7)
            TAUZZ  = Q(I,8)
            TAUXY  = Q(I,9)
            TAUXZ  = Q(I,10)
            TAUYZ  = Q(I,11)
            EPSLN  = Q(I,12)
C
            TKE    = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
            RMU    = PROPS(I,1)
            YPLUS  = PROPS(I,3)
C
C     Anisotropy tensor
C
            BXX    = PIIJ(I,1)
            BYY    = PIIJ(I,2)
            BZZ    = PIIJ(I,3)
            BXY    = PIIJ(I,4)
            BXZ    = PIIJ(I,5)
            BYZ    = PIIJ(I,6)
C
            SXX    = SIJ(I,1)
            SYY    = SIJ(I,2)
            SZZ    = SIJ(I,3)
            SXY    = SIJ(I,4)
            SXZ    = SIJ(I,5)
            SYZ    = SIJ(I,6)
C
            WYX    = WIJ(I,1)
            WZX    = WIJ(I,2)
            WZY    = WIJ(I,3)
            WXY    = WIJ(I,4)
            WXZ    = WIJ(I,5)
            WYZ    = WIJ(I,6)
C
C     S_kk = (S_xx + S_yy + S_zz) / 3
C
            SKK    = (SXX + SYY + SZZ) / 3.E0
C
C     b_lk * S_lk
C
            BLKSLK =       BXX*SXX + BYY*SYY + BZZ*SZZ
     1             + 2.E0*(BXY*SXY + BXZ*SXZ + BYZ*SYZ)
C
C     II  = b_mn b_nm
C
            BMNBMN =       BXX*BXX + BYY*BYY + BZZ*BZZ
     1             + 2.E0*(BXY*BXY + BXZ*BXZ + BYZ*BYZ)
C
C     III = b_ij b_jk b_ki
C
            BIII   = BXX*BXX*BXX + BXX*BXY*BXY + BXX*BXZ*BXZ
     1             + BXY*BXY*BXX + BXY*BYY*BXY + BXY*BYZ*BXZ
     2             + BXZ*BXZ*BXX + BXZ*BYZ*BXY + BXZ*BZZ*BXZ
     3             + BXY*BXX*BXY + BXY*BXY*BYY + BXY*BXZ*BYZ
     4             + BYY*BXY*BXY + BYY*BYY*BYY + BYY*BYZ*BYZ
     5             + BYZ*BXZ*BXY + BYZ*BYZ*BYY + BYZ*BZZ*BYZ
     6             + BXZ*BXX*BXZ + BXZ*BXY*BYZ + BXZ*BXZ*BZZ
     7             + BYZ*BXY*BXZ + BYZ*BYY*BYZ + BYZ*BYZ*BZZ
     8             + BZZ*BXZ*BXZ + BZZ*BYZ*BYZ + BZZ*BZZ*BZZ
C
C     Coefficient for SSG
C
            CA2HR  =   0.80E0 - 1.3E0 * SQRT(BMNBMN)
C
C     Near-Wall damping
C
c-          AFLAT  = 1.E0 - 4.5E0*(BMNBMN - 2.E0*BIII)
            RET    = RHO * TKE * TKE / (RMU * EPSLN) * RE / FSMACH
c-          FW1    = EXP (-(AFLAT * RET / 57)**3)
            FW1    = EXP (-(RET / 175)**2)
C
C     Near-Wall coefficients from So et al. 1995
C
            CA0    = CA0HR - FW1 * CA0HR
            CA1    = CA1HR - FW1 * CA1HR
            CA2    = CA2HR - FW1 * 2.E0*(T23*(PSALFA+PSBETA)-PSGAMA)
            CA3    = CA3HR - FW1 * CA3HR
            CA4    = CA4HR - FW1 * 2.E0*(PSALFA+PSBETA)
            CA5    = CA5HR - FW1 * 2.E0*(PSALFA-PSBETA)
C
C     P_k => Production of turbulent kinetic energy
C
            PK     =   - ( TAUXX*SXX + TAUYY*SYY + TAUZZ*SZZ
     1             + 2.E0*(TAUXY*SXY + TAUXZ*SXZ + TAUYZ*SYZ) )
C
C     PI_ij = Phi_ij,1 + Phi_ij,2 
C
            PIIJ(I,1) = - CA0 * RHO * EPSLN * BXX
     1 + CA1 * RHO * EPSLN * (BXX*BXX + BXY*BXY + BXZ*BXZ - T13*BMNBMN)
     2 + RHO * TKE * 
     3  ( CA4 * (2.E0*(BXX*SXX + BXY*SXY + BXZ*SXZ) - T23*BLKSLK)
     4  + CA5 * (2.E0*(          BXY*WXY + BXZ*WXZ))
     5  + CA2 * (SXX - SKK) )
     6 + CA3 * RHO * PK * BXX
C
            PIIJ(I,2) = - CA0 * RHO * EPSLN * BYY
     1 + CA1 * RHO * EPSLN * (BXY*BXY + BYY*BYY + BYZ*BYZ - T13*BMNBMN)
     2 + RHO * TKE * 
     3  ( CA4 * (2.E0*(BXY*SXY + BYY*SYY + BYZ*SYZ) - T23*BLKSLK)
     4  + CA5 * (2.E0*(BXY*WYX           + BYZ*WYZ))
     5  + CA2 * (SYY - SKK) )
     6 + CA3 * RHO * PK * BYY
C
            PIIJ(I,3) = - CA0 * RHO * EPSLN * BZZ
     1 + CA1 * RHO * EPSLN * (BXZ*BXZ + BYZ*BYZ + BZZ*BZZ - T13*BMNBMN)
     2 + RHO * TKE * 
     3  ( CA4 * (2.E0*(BXZ*SXZ + BYZ*SYZ + BZZ*SZZ) - T23*BLKSLK)
     4  + CA5 * (2.E0*(BXZ*WZX + BYZ*WZY          ))
     5  + CA2 * (SZZ - SKK) )
     6 + CA3 * RHO * PK * BZZ
C
            PIIJ(I,4) = - CA0 * RHO * EPSLN * BXY 
     1 + CA1 * RHO * EPSLN * (BXX*BXY + BXY*BYY + BXZ*BYZ)
     2 + RHO * TKE * 
     3  (CA4*(BXX*SXY + BXY*SYY + BXZ*SYZ + BXY*SXX + BYY*SXY + BYZ*SXZ)
     4  +CA5*(BXX*WYX           + BXZ*WYZ           + BYY*WXY + BYZ*WXZ)
     5  +CA2*(SXY) )
     6 + CA3 * RHO * PK * BXY
C
            PIIJ(I,5) = - CA0 * RHO * EPSLN * BXZ 
     1 + CA1 * RHO * EPSLN * (BXX*BXZ + BXY*BYZ + BXZ*BZZ)
     2 + RHO * TKE * 
     3  (CA4*(BXX*SXZ + BXY*SYZ + BXZ*SZZ + BXZ*SXX + BYZ*SXY + BZZ*SXZ)
     4  +CA5*(BXX*WZX + BXY*WZY                     + BYZ*WXY + BZZ*WXZ)
     5  +CA2*(SXZ) )
     6 + CA3 * RHO * PK * BXZ
C
            PIIJ(I,6) = - CA0 * RHO * EPSLN * BYZ 
     1 + CA1 * RHO * EPSLN * (BXY*BXZ + BYY*BYZ + BYZ*BZZ)
     2 + RHO * TKE * 
     3  (CA4*(BXY*SXZ + BYY*SYZ + BYZ*SZZ + BXZ*SXY + BYZ*SYY + BZZ*SYZ)
     4  +CA5*(BXY*WZX + BYY*WZY           + BXZ*WYX           + BZZ*WYZ)
     5  +CA2*(SYZ) )
     6 + CA3 * RHO * PK * BYZ
C
  600    CONTINUE
C
C     Error
C
      ELSE
         WRITE (IOUT,1000) IPSMDL
         STOP
      ENDIF
C
C     Finished calculating pressure-strain correlation
C
      RETURN
 1000 FORMAT (' ','ERROR - RUN ABORTING',
     1       /' ','Incorrect Pressure-Strain Correlation ',
     2            'Specified in routine PSTRAN = ',I10)
      END






      SUBROUTINE PSTRANW (IPSMDL, NPTS, Q, PROPS, SN, YNI, SRCW)
C
C Routine to calculate the wall surface normal dependent Pressure-Strain 
C Correlation terms for a Reynolds stress turbulence model.
C
C This routine is organized to calculate wall normal dependent source terms
C for the Pressure-Strain in the following sections:
C         Model                                        Wall Terms (YES/NO)
C      --------------------------------------------------------------------
C      0. Isotropic (High Reynolds Number Model)                    No
C      1. Gison-Launder Pressure-Strain                             No
C         (Launder-Reece-Rodi "Basic Model" Pressure-Strain)
C      2. Launder-Reece-Rodi Pressure Strain                        No
C      3. Lai-So Pressure Strain                                    YES
C      4. Speziale-Sarkar-Gatski Pressure Strain                    No
C      5. Fu-Launder-Tselepidakis (FLT) Pressure Strain             No
C      6. So et al. 1995 Pressure Strain                            No
C
C IPSMDL         : Pressure-strain model
C NPTS           : Number of points to calculate source terms at
C Q              : Primitive variables at cell centers
C PROPS          : Properties stored at cell centers
C                  PROPS(1) = RMU   molecular viscosity
C                  PROPS(2) = RMUT  turbulent eddy viscosity
C                  PROPS(3) = YPLUS Y+
C                  PROPS(4) = YNORML
C SN             : Wall normals
C YNI            : Inverse distance to wall
C SRCW           : Wall normal dependent terms for Pressure-Strain
C
      include '../header/common.h'
C
      DIMENSION Q     (NPTS,NQ),
     1          PROPS (NPTS,NP),
     2          SN    (NPTS,3),
     3          YNI   (NPTS)
C
      DIMENSION SRCW  (NPTS,6)
C
C     One-Thirds
C     Two-Thirds
C
      T13 = 1.E0 / 3.E0
      T23 = 2.E0 / 3.E0
C
C Calculate wall normal dependent Pressure-Strain Correlation
C
C 3. Lai-So Pressure Strain Correlation
C    NASA Contractor Report 4369 (Based on LRR1)
C    Corrected form based on Zhang, So, Gatski, Speziale in the 
C    Near-Wall Turbulent Flows book later corrected from private communications
C    with Ron So.  Note that this is not strictly a pressure-strain model
C    but rather a model of pressure-strain and pressure-diffusion so that the
C    trace of the model is Non-Zero.
C    Modified based on private communication with Hou-Shun Zhang to include
C    the subtraction of the trace term from the S_ij term in the near-wall.
C
      IF (IPSMDL .EQ. IPLS) THEN
         DO 300 I = 1, NPTS
C
            RHO    = Q(I, 1)
            P      = Q(I, 5)
            TAUXX  = Q(I, 6)
            TAUYY  = Q(I, 7)
            TAUZZ  = Q(I, 8)
            TAUXY  = Q(I, 9)
            TAUXZ  = Q(I,10)
            TAUYZ  = Q(I,11)
            EPSLN  = Q(I,12)
C
            TKE    = 0.5E0 * (TAUXX + TAUYY + TAUZZ)
C
            RMU    = PROPS(I,1)
            YPLUS  = PROPS(I,3)
C
C     Wall damping
C
            RET    = RHO * TKE * TKE / (RMU * EPSLN) * RE / FSMACH
            FW1    = EXP( -(RET / 150.E0)**2 )
            IF (YPLUS .GT. 250.E0) FW1 = 0.E0
C
C     SRCW_ij = SRCW_ij + Phi_ij,w
C
            SRCW(I,1) = SRCW(I,1) - YNI(I) * FW1 * RHO * EPSLN / TKE * (
     1TAUXX*SN(I,1)*SN(I,1)+TAUXY*SN(I,2)*SN(I,1)+TAUXZ*SN(I,3)*SN(I,1)+
     2TAUXX*SN(I,1)*SN(I,1)+TAUXY*SN(I,2)*SN(I,1)+TAUXZ*SN(I,3)*SN(I,1))
C
            SRCW(I,2) = SRCW(I,2) - YNI(I) * FW1 * RHO * EPSLN / TKE * (
     1TAUXY*SN(I,1)*SN(I,2)+TAUYY*SN(I,2)*SN(I,2)+TAUYZ*SN(I,3)*SN(I,2)+
     2TAUXY*SN(I,1)*SN(I,2)+TAUYY*SN(I,2)*SN(I,2)+TAUYZ*SN(I,3)*SN(I,2))
C
            SRCW(I,3) = SRCW(I,3) - YNI(I) * FW1 * RHO * EPSLN / TKE * (
     1TAUXZ*SN(I,1)*SN(I,3)+TAUYZ*SN(I,2)*SN(I,3)+TAUZZ*SN(I,3)*SN(I,3)+
     2TAUXZ*SN(I,1)*SN(I,3)+TAUYZ*SN(I,2)*SN(I,3)+TAUZZ*SN(I,3)*SN(I,3))
C
            SRCW(I,4) = SRCW(I,4) - YNI(I) * FW1 * RHO * EPSLN / TKE * (
     1TAUXX*SN(I,1)*SN(I,2)+TAUXY*SN(I,2)*SN(I,2)+TAUXZ*SN(I,3)*SN(I,2)+
     2TAUXY*SN(I,1)*SN(I,1)+TAUYY*SN(I,2)*SN(I,1)+TAUYZ*SN(I,3)*SN(I,1))
C
            SRCW(I,5) = SRCW(I,5) - YNI(I) * FW1 * RHO * EPSLN / TKE * (
     1TAUXX*SN(I,1)*SN(I,3)+TAUXY*SN(I,2)*SN(I,3)+TAUXZ*SN(I,3)*SN(I,3)+
     2TAUXZ*SN(I,1)*SN(I,1)+TAUYZ*SN(I,2)*SN(I,1)+TAUZZ*SN(I,3)*SN(I,1))
C
            SRCW(I,6) = SRCW(I,6) - YNI(I) * FW1 * RHO * EPSLN / TKE * (
     1TAUXY*SN(I,1)*SN(I,3)+TAUYY*SN(I,2)*SN(I,3)+TAUYZ*SN(I,3)*SN(I,3)+
     2TAUXZ*SN(I,1)*SN(I,2)+TAUYZ*SN(I,2)*SN(I,2)+TAUZZ*SN(I,3)*SN(I,2))
C
  300    CONTINUE
C
      ENDIF
C
C Finished calculating wall surface normal dependent pressure-strain terms
C
      RETURN
      END





      SUBROUTINE RSDFEW (IEPSIJ, NPTS, QI, PROPSI, DQDX, DQDY, DQDZ, 
     1                   S, TAU, F)
C
C Routine to calculate the wall correction to the epsilon_ij tensor that
C depends on the diffusion term.  Modeled forms are:
C      1. So et al. (1995) model
C
C The flow variables, properties and derivatives of the flow variables are 
C provided to this routine from which it calculates the diffusive contribution
C to epsilon_wall_ij.
C
C This routine is organized as follows:
C      1. Calculate epsilon_wall_ij
C         a. So et al. (1995)
C
C IEPSIJ         : Epsilon_wall_ij model
C NPTS           : Number of interfaces to calculate viscous fluxes
C QI             : Primitive variables at interfaces
C PROPSI         : Properties at the interfaces
C                  PROPSI(1) = RMU
C                  PROPSI(2) = RMUT
C                  PROPSI(3) = Y+
C                  PROPSI(4) = YNORML
C DQDX           : X-Derivative of primitive variables at the interface
C DQDY           : Y-Derivative of primitive variables at the interface
C DQDZ           : Z-Derivative of primitive variables at the interface
C S              : Metrics
C TAU            : Reynold's stresses
C F              : Viscous fluxes at cell interfaces
C
      include '../header/common.h'
C
      DIMENSION QI    (NPTS,NQ),
     1          PROPSI(NPTS,NP),
     2          DQDX  (NPTS,NQ),
     3          DQDY  (NPTS,NQ),
     4          DQDZ  (NPTS,NQ),
     5          S     (NPTS,3),
     6          TAU   (NPTS,6),
     7          F     (NPTS,NF)
C
C     Two-Thirds
C
      T23 = 2.E0 / 3.E0
C
C
C 1.   Calculate epsilon_wall_ij from So et al. 1995
C
      IF (IEPSIJ .EQ. IESO95) THEN
         DO 100 I = 1, NPTS
C
            RMU    = PROPSI(I,1)
C
C     VRS : Coefficient for diffusion term for epsilon_wall_ij
C           Multiply by 1/2 as the diffusion term for the wall correction
C           to epsilon_ij has a 1/2 in front.
C
            VRS    = 0.5E0 * RMU * FSMACH / RE
C
C     d(K)/d(x_k) : 1/2 here since k = 1/2 Tau_ii
C
            DKDXK  = 0.5E0 * T23 * 
     1               ( S(I,1) * (DQDX(I,6) + DQDX(I,7) + DQDX(I,8))
     2               + S(I,2) * (DQDY(I,6) + DQDY(I,7) + DQDY(I,8))
     3               + S(I,3) * (DQDZ(I,6) + DQDZ(I,7) + DQDZ(I,8)) )
C
C      Calculate F(molecular diffusion)
C
            F(I, 6) = F(I, 6)
     1              + VRS * ( S(I,1)*DQDX(I, 6) + S(I,2)*DQDY(I, 6)
     2                      + S(I,3)*DQDZ(I, 6) - DKDXK )
            F(I, 7) = F(I, 7)
     1              + VRS * ( S(I,1)*DQDX(I, 7) + S(I,2)*DQDY(I, 7)
     2                      + S(I,3)*DQDZ(I, 7) - DKDXK )
            F(I, 8) = F(I, 8)
     1              + VRS * ( S(I,1)*DQDX(I, 8) + S(I,2)*DQDY(I, 8)
     2                      + S(I,3)*DQDZ(I, 8) - DKDXK )
            F(I, 9) = F(I, 9)
     1              + VRS * ( S(I,1)*DQDX(I, 9) + S(I,2)*DQDY(I, 9)
     2                      + S(I,3)*DQDZ(I, 9) )
            F(I,10) = F(I,10)
     1              + VRS * ( S(I,1)*DQDX(I,10) + S(I,2)*DQDY(I,10)
     2                      + S(I,3)*DQDZ(I,10) )
            F(I,11) = F(I,11)
     1              + VRS * ( S(I,1)*DQDX(I,11) + S(I,2)*DQDY(I,11)
     2                      + S(I,3)*DQDZ(I,11) )
  100    CONTINUE
C
C Error - Output error message and stop
C
      ELSE
         WRITE (IOUT,1000) IEPSIJ
         STOP
      ENDIF
C
C Finished with epsilon_wall_ij
C
      RETURN
 1000 FORMAT (' ','ERROR - RUN ABORTING',
     1       /' ','Invalid Epsilon_wall_ij in routine ',
     2            'RSDFEW  = ',I10)
      END






      SUBROUTINE RSMDIF (IDMMDL, NPTS, QI, PROPSI, DQDX, DQDY, DQDZ, 
     1                   S, TAU, F)
C
C Routine to calculate the Reynolds Stress equations' Molecular Diffusion.
C
C The flow variables, properties and derivatives of the flow variables are 
C provided to this routine from which it calculates the viscous fluxes.
C
C This routine is organized as follows:
C      1. Calculate Molecular Diffusion for Reynolds Stresses
C         a. Simple form              => Ron So
C         b. Complex form             => Charles Speziale
C
C IDMMDL         : Reynolds stress Molecular diffusion model
C NPTS           : Number of interfaces to calculate viscous fluxes
C QI             : Primitive variables at interfaces
C PROPSI         : Properties at the interfaces
C                  PROPSI(1) = RMU
C                  PROPSI(2) = RMUT
C                  PROPSI(3) = Y+
C                  PROPSI(4) = YNORML
C DQDX           : X-Derivative of primitive variables at the interface
C DQDY           : Y-Derivative of primitive variables at the interface
C DQDZ           : Z-Derivative of primitive variables at the interface
C S              : Metrics
C TAU            : Reynold's stresses
C F              : Viscous fluxes at cell interfaces
C
      include '../header/common.h'
C
      DIMENSION QI    (NPTS,NQ),
     1          PROPSI(NPTS,NP),
     2          DQDX  (NPTS,NQ),
     3          DQDY  (NPTS,NQ),
     4          DQDZ  (NPTS,NQ),
     5          S     (NPTS,3),
     6          TAU   (NPTS,6),
     7          F     (NPTS,NF)
C
C     Two-Thirds
C
      T23 = 2.E0 / 3.E0
C
C
C 1.   Calculate Molecular Diffusion
C
C 1-a. Simple form (So et al.)
C
      IF (IDMMDL .EQ. IDMSO) THEN
         DO 100 I = 1, NPTS
C
            RMU    = PROPSI(I,1)
C
C     VRS : Coefficient for molecular diffusion in Tau_ij
C
            VRS    = RMU * FSMACH / RE
C
C      Calculate F(molecular diffusion)
C
           F(I,6) = F(I, 6) 
     1            + VRS * ( S(I,1)*DQDX(I, 6) + S(I,2)*DQDY(I, 6)
     2                    + S(I,3)*DQDZ(I, 6) )
           F(I,7) = F(I, 7)
     1            + VRS * ( S(I,1)*DQDX(I, 7) + S(I,2)*DQDY(I, 7)
     2                    + S(I,3)*DQDZ(I, 7) )
           F(I,8) = F(I, 8)
     1            + VRS * ( S(I,1)*DQDX(I, 8) + S(I,2)*DQDY(I, 8)
     2                    + S(I,3)*DQDZ(I, 8) )
           F(I,9) = F(I, 9)
     1            + VRS * ( S(I,1)*DQDX(I, 9) + S(I,2)*DQDY(I, 9)
     2                    + S(I,3)*DQDZ(I, 9) )
           F(I,10)= F(I,10)
     1            + VRS * ( S(I,1)*DQDX(I,10) + S(I,2)*DQDY(I,10)
     2                    + S(I,3)*DQDZ(I,10) )
           F(I,11)= F(I,11)
     1            + VRS * ( S(I,1)*DQDX(I,11) + S(I,2)*DQDY(I,11)
     2                    + S(I,3)*DQDZ(I,11) )
  100    CONTINUE
C
C 1-b. Complex form (Speziale)
C
      ELSE IF (IDMMDL .EQ. IDMSPZ) THEN
         DO 200 I = 1, NPTS
C
            RMU    = PROPSI(I,1)
C
C      VRS : Coefficient for molecular diffusion in Tau_ij
C
            VRS    = RMU * FSMACH / RE
C
C      Calculate F(molecular diffusion)
C
           F(I,6) = F(I, 6) 
     1            + VRS * ( S(I,1)*(DQDX(I, 6)+DQDX(I, 6)+DQDX(I, 6))
     2                    + S(I,2)*(DQDY(I, 6)+DQDX(I, 9)+DQDX(I, 9))
     3                    + S(I,3)*(DQDZ(I, 6)+DQDX(I,10)+DQDX(I,10)) )
           F(I,7) = F(I, 7)
     1            + VRS * ( S(I,1)*(DQDX(I, 7)+DQDY(I, 9)+DQDY(I, 9))
     2                    + S(I,2)*(DQDY(I, 7)+DQDY(I, 7)+DQDY(I, 7))
     3                    + S(I,3)*(DQDZ(I, 7)+DQDY(I,11)+DQDY(I,11)) )
           F(I,8) = F(I, 8)
     1            + VRS * ( S(I,1)*(DQDX(I, 8)+DQDZ(I,10)+DQDZ(I,10))
     2                    + S(I,2)*(DQDY(I, 8)+DQDZ(I,11)+DQDZ(I,11))
     3                    + S(I,3)*(DQDZ(I, 8)+DQDZ(I, 8)+DQDZ(I, 8)) )
           F(I,9) = F(I, 9)
     1            + VRS * ( S(I,1)*(DQDX(I, 9)+DQDY(I, 6)+DQDX(I, 9))
     2                    + S(I,2)*(DQDY(I, 9)+DQDY(I, 9)+DQDX(I, 7))
     3                    + S(I,3)*(DQDZ(I, 9)+DQDY(I,10)+DQDX(I,11)) )
           F(I,10)= F(I,10)
     1            + VRS * ( S(I,1)*(DQDX(I,10)+DQDZ(I, 6)+DQDX(I,10))
     2                    + S(I,2)*(DQDY(I,10)+DQDZ(I, 9)+DQDX(I,11))
     3                    + S(I,3)*(DQDZ(I,10)+DQDZ(I,10)+DQDX(I, 8)) )
           F(I,11)= F(I,11)
     1            + VRS * ( S(I,1)*(DQDX(I,11)+DQDZ(I, 9)+DQDY(I,10))
     2                    + S(I,2)*(DQDY(I,11)+DQDZ(I, 7)+DQDY(I,11))
     3                    + S(I,3)*(DQDZ(I,11)+DQDZ(I,11)+DQDY(I, 8)) )
  200    CONTINUE
C
C Error - Output error message and stop
C
      ELSE
         WRITE (IOUT,1000) IDMMDL
         STOP
      ENDIF
C
C Finished with Molecular Diffusion Terms
C
      RETURN
 1000 FORMAT (' ','ERROR - RUN ABORTING',
     1       /' ','Invalid Molecular diffusion model in routine ',
     2            'RSMDIF  = ',I10)
      END






      SUBROUTINE RSTDIF (IDTMDL, NPTS, QI, PROPSI, DQDX, DQDY, DQDZ, 
     1                   S, TAU, F)
C
C Routine to calculate the Reynolds Stress equations' TURBULENT Diffusion.
C
C The flow variables, properties and derivatives of the flow variables are 
C provided to this routine from which it calculates the viscous fluxes.
C
C This routine is organized as follows:
C      1. Calculate TURBULENT Diffusion for Reynolds Stresses
C         a. DH         => Daly and Harlow 
C         b. HL72       => Hanjalic and Launder (1972)
C         c. MH73       => Mellor and Herring (1973)
C
C IDTMDL         : Reynolds stress Turbulent diffusion model
C NPTS           : Number of interfaces to calculate viscous fluxes
C QI             : Primitive variables at interfaces
C PROPSI         : Properties at the interfaces
C                  PROPSI(1) = RMU
C                  PROPSI(2) = RMUT
C                  PROPSI(3) = Y+
C                  PROPSI(4) = YNORML
C DQDX           : X-Derivative of primitive variables at the interface
C DQDY           : Y-Derivative of primitive variables at the interface
C DQDZ           : Z-Derivative of primitive variables at the interface
C S              : Metrics
C TAU            : Reynold's stresses
C F              : Viscous fluxes at cell interfaces
C
      include '../header/common.h'
C
      DIMENSION QI    (NPTS,NQ),
     1          PROPSI(NPTS,NP),
     2          DQDX  (NPTS,NQ),
     3          DQDY  (NPTS,NQ),
     4          DQDZ  (NPTS,NQ),
     5          S     (NPTS,3),
     6          TAU   (NPTS,6),
     7          F     (NPTS,NF)
C
C     Two-Thirds
C
      T23 = 2.E0 / 3.E0
C
C 1.   Calculate Turbulent Diffusion
C
C 1-a. Daly-Harlow
C
      IF (IDTMDL .EQ. IDTDH) THEN
         DO 100 I = 1, NPTS
C
            RHO    = QI(I,1)
            TKE    = 0.5E0 * (QI(I,6) + QI(I,7) + QI(I,8))
            EPSLN  = QI(I,12)
C
C      VDT : Coefficient for turbulent diffusion (DT_ij) in Tau_ij
C
            VDT    = CSUBS * RHO * TKE / EPSLN 
C
C      Calculate F(Turbulent diffusion)
C
           F(I, 6) = F(I, 6) + VDT * 
     1 ( S(I,1)*
     2   (QI(I, 6)*DQDX(I, 6)+QI(I, 9)*DQDY(I, 6)+QI(I,10)*DQDZ(I, 6))
     3 + S(I,2)*
     4   (QI(I, 9)*DQDX(I, 6)+QI(I, 7)*DQDY(I, 6)+QI(I,11)*DQDZ(I, 6))
     5 + S(I,3)*
     6   (QI(I,10)*DQDX(I, 6)+QI(I,11)*DQDY(I, 6)+QI(I, 8)*DQDZ(I, 6)) )
C
           F(I, 7) = F(I, 7) + VDT * 
     1 ( S(I,1)*
     2   (QI(I, 6)*DQDX(I, 7)+QI(I, 9)*DQDY(I, 7)+QI(I,10)*DQDZ(I, 7))
     3 + S(I,2)*
     4   (QI(I, 9)*DQDX(I, 7)+QI(I, 7)*DQDY(I, 7)+QI(I,11)*DQDZ(I, 7))
     5 + S(I,3)*
     6   (QI(I,10)*DQDX(I, 7)+QI(I,11)*DQDY(I, 7)+QI(I, 8)*DQDZ(I, 7)) )
C
           F(I, 8) = F(I, 8) + VDT * 
     1 ( S(I,1)*
     2   (QI(I, 6)*DQDX(I, 8)+QI(I, 9)*DQDY(I, 8)+QI(I,10)*DQDZ(I, 8))
     3 + S(I,2)*
     4   (QI(I, 9)*DQDX(I, 8)+QI(I, 7)*DQDY(I, 8)+QI(I,11)*DQDZ(I, 8))
     5 + S(I,3)*
     6   (QI(I,10)*DQDX(I, 8)+QI(I,11)*DQDY(I, 8)+QI(I, 8)*DQDZ(I, 8)) )
C
           F(I, 9) = F(I, 9) + VDT * 
     1 ( S(I,1)*
     2   (QI(I, 6)*DQDX(I, 9)+QI(I, 9)*DQDY(I, 9)+QI(I,10)*DQDZ(I, 9))
     3 + S(I,2)*
     4   (QI(I, 9)*DQDX(I, 9)+QI(I, 7)*DQDY(I, 9)+QI(I,11)*DQDZ(I, 9))
     5 + S(I,3)*
     6   (QI(I,10)*DQDX(I, 9)+QI(I,11)*DQDY(I, 9)+QI(I, 8)*DQDZ(I, 9)) )
C
           F(I,10) = F(I,10) + VDT * 
     1 ( S(I,1)*
     2   (QI(I, 6)*DQDX(I,10)+QI(I, 9)*DQDY(I,10)+QI(I,10)*DQDZ(I,10))
     3 + S(I,2)*
     4   (QI(I, 9)*DQDX(I,10)+QI(I, 7)*DQDY(I,10)+QI(I,11)*DQDZ(I,10))
     5 + S(I,3)*
     6   (QI(I,10)*DQDX(I,10)+QI(I,11)*DQDY(I,10)+QI(I, 8)*DQDZ(I,10)) )
C
           F(I,11) = F(I,11) + VDT * 
     1 ( S(I,1)*
     2   (QI(I, 6)*DQDX(I,11)+QI(I, 9)*DQDY(I,11)+QI(I,10)*DQDZ(I,11))
     3 + S(I,2)*
     4   (QI(I, 9)*DQDX(I,11)+QI(I, 7)*DQDY(I,11)+QI(I,11)*DQDZ(I,11))
     5 + S(I,3)*
     6   (QI(I,10)*DQDX(I,11)+QI(I,11)*DQDY(I,11)+QI(I, 8)*DQDZ(I,11)) )
  100    CONTINUE
C
C 1-b. HL72 - Hanjalic and Launder (1972)
C
      ELSE IF (IDTMDL .EQ. IDTHL) THEN
         DO 200 I = 1, NPTS
C
            RHO    = QI(I,1)
            TKE    = 0.5E0 * (QI(I,6) + QI(I,7) + QI(I,8))
            EPSLN  = QI(I,12)
C
C      VDT : Coefficient for turbulent diffusion (DT_ij) in Tau_ij
C
            VDT    = CSUBS * RHO * TKE / EPSLN 
C
C      Calculate F(Turbulent diffusion)
C
           F(I,6) = F(I, 6) + VDT * 
     1 ( S(I,1)*
     2   (QI(I, 6)*DQDX(I, 6)+QI(I, 9)*DQDY(I, 6)+QI(I,10)*DQDZ(I, 6)+
     3    QI(I, 6)*DQDX(I, 6)+QI(I, 9)*DQDY(I, 6)+QI(I,10)*DQDZ(I, 6)+
     4    QI(I, 6)*DQDX(I, 6)+QI(I, 9)*DQDY(I, 6)+QI(I,10)*DQDZ(I, 6))
     5 + S(I,2)*
     6   (QI(I, 9)*DQDX(I, 6)+QI(I, 7)*DQDY(I, 6)+QI(I,11)*DQDZ(I, 6)+
     7    QI(I, 6)*DQDX(I, 9)+QI(I, 9)*DQDY(I, 9)+QI(I,10)*DQDZ(I, 9)+
     8    QI(I, 6)*DQDX(I, 9)+QI(I, 9)*DQDY(I, 9)+QI(I,10)*DQDZ(I, 9))
     9 + S(I,3)*
     A   (QI(I,10)*DQDX(I, 6)+QI(I,11)*DQDY(I, 6)+QI(I, 8)*DQDZ(I, 6)+
     B    QI(I, 6)*DQDX(I,10)+QI(I, 9)*DQDY(I,10)+QI(I,10)*DQDZ(I,10)+
     C    QI(I, 6)*DQDX(I,10)+QI(I, 9)*DQDY(I,10)+QI(I,10)*DQDZ(I,10)) )
C
           F(I,7) = F(I, 7) + VDT *
     1 ( S(I,1)*
     2   (QI(I, 6)*DQDX(I, 7)+QI(I, 9)*DQDY(I, 7)+QI(I,10)*DQDZ(I, 7)+
     3    QI(I, 9)*DQDX(I, 9)+QI(I, 7)*DQDY(I, 9)+QI(I,11)*DQDZ(I, 9)+
     4    QI(I, 9)*DQDX(I, 9)+QI(I, 7)*DQDY(I, 9)+QI(I,11)*DQDZ(I, 9))
     5 + S(I,2)*
     6   (QI(I, 9)*DQDX(I, 7)+QI(I, 7)*DQDY(I, 7)+QI(I,11)*DQDZ(I, 7)+
     7    QI(I, 9)*DQDX(I, 7)+QI(I, 7)*DQDY(I, 7)+QI(I,11)*DQDZ(I, 7)+
     8    QI(I, 9)*DQDX(I, 7)+QI(I, 7)*DQDY(I, 7)+QI(I,11)*DQDZ(I, 7))
     9 + S(I,3)*
     A   (QI(I,10)*DQDX(I, 7)+QI(I,11)*DQDY(I, 7)+QI(I, 8)*DQDZ(I, 7)+
     B    QI(I, 9)*DQDX(I,11)+QI(I, 7)*DQDY(I,11)+QI(I,11)*DQDZ(I,11)+
     C    QI(I, 9)*DQDX(I,11)+QI(I, 7)*DQDY(I,11)+QI(I,11)*DQDZ(I,11)) )
C
           F(I,8) = F(I, 8) + VDT *
     1 ( S(I,1)*
     2   (QI(I, 6)*DQDX(I, 8)+QI(I, 9)*DQDY(I, 8)+QI(I,10)*DQDZ(I, 8)+
     3    QI(I,10)*DQDX(I,10)+QI(I,11)*DQDY(I,10)+QI(I, 8)*DQDZ(I,10)+
     4    QI(I,10)*DQDX(I,10)+QI(I,11)*DQDY(I,10)+QI(I, 8)*DQDZ(I,10))
     5 + S(I,2)*
     6   (QI(I, 9)*DQDX(I, 8)+QI(I, 7)*DQDY(I, 8)+QI(I,11)*DQDZ(I, 8)+
     7    QI(I,10)*DQDX(I,11)+QI(I,11)*DQDY(I,11)+QI(I, 8)*DQDZ(I,11)+
     8    QI(I,10)*DQDX(I,11)+QI(I,11)*DQDY(I,11)+QI(I, 8)*DQDZ(I,11))
     9 + S(I,3)*
     A   (QI(I,10)*DQDX(I, 8)+QI(I,11)*DQDY(I, 8)+QI(I, 8)*DQDZ(I, 8)+
     B    QI(I,10)*DQDX(I, 8)+QI(I,11)*DQDY(I, 8)+QI(I, 8)*DQDZ(I, 8)+
     C    QI(I,10)*DQDX(I, 8)+QI(I,11)*DQDY(I, 8)+QI(I, 8)*DQDZ(I, 8)) )
C
           F(I,9) = F(I, 9) + VDT *
     1 ( S(I,1)*
     2   (QI(I, 6)*DQDX(I, 9)+QI(I, 9)*DQDY(I, 9)+QI(I,10)*DQDZ(I, 9)+
     3    QI(I, 9)*DQDX(I, 6)+QI(I, 7)*DQDY(I, 6)+QI(I,11)*DQDZ(I, 6)+
     4    QI(I, 6)*DQDX(I, 9)+QI(I, 9)*DQDY(I, 9)+QI(I,10)*DQDZ(I, 9))
     5 + S(I,2)*
     6   (QI(I, 9)*DQDX(I, 9)+QI(I, 7)*DQDY(I, 9)+QI(I,11)*DQDZ(I, 9)+
     7    QI(I, 9)*DQDX(I, 9)+QI(I, 7)*DQDY(I, 9)+QI(I,11)*DQDZ(I, 9)+
     8    QI(I, 6)*DQDX(I, 7)+QI(I, 9)*DQDY(I, 7)+QI(I,10)*DQDZ(I, 7))
     9 + S(I,3)*
     A   (QI(I,10)*DQDX(I, 9)+QI(I,11)*DQDY(I, 9)+QI(I, 8)*DQDZ(I, 9)+
     B    QI(I, 9)*DQDX(I,10)+QI(I, 7)*DQDY(I,10)+QI(I,11)*DQDZ(I,10)+
     C    QI(I, 6)*DQDX(I,11)+QI(I, 9)*DQDY(I,11)+QI(I,10)*DQDZ(I,11)) )
C
           F(I,10)= F(I,10) + VDT *
     1 ( S(I,1)*
     2   (QI(I, 6)*DQDX(I,10)+QI(I, 9)*DQDY(I,10)+QI(I,10)*DQDZ(I,10)+
     3    QI(I,10)*DQDX(I, 6)+QI(I,11)*DQDY(I, 6)+QI(I, 8)*DQDZ(I, 6)+
     4    QI(I, 6)*DQDX(I,10)+QI(I, 9)*DQDY(I,10)+QI(I,10)*DQDZ(I,10))
     5 + S(I,2)*
     6   (QI(I, 9)*DQDX(I,10)+QI(I, 7)*DQDY(I,10)+QI(I,11)*DQDZ(I,10)+
     7    QI(I,10)*DQDX(I, 9)+QI(I,11)*DQDY(I, 9)+QI(I, 8)*DQDZ(I, 9)+
     8    QI(I, 6)*DQDX(I,11)+QI(I, 9)*DQDY(I,11)+QI(I,10)*DQDZ(I,11))
     9 + S(I,3)*
     A   (QI(I,10)*DQDX(I,10)+QI(I,11)*DQDY(I,10)+QI(I, 8)*DQDZ(I,10)+
     B    QI(I,10)*DQDX(I,10)+QI(I,11)*DQDY(I,10)+QI(I, 8)*DQDZ(I,10)+
     C    QI(I, 6)*DQDX(I, 8)+QI(I, 9)*DQDY(I, 8)+QI(I,10)*DQDZ(I, 8)) )
C
           F(I,11)= F(I,11) + VDT *
     1 ( S(I,1)*
     2   (QI(I, 6)*DQDX(I,11)+QI(I, 9)*DQDY(I,11)+QI(I,10)*DQDZ(I,11)+
     3    QI(I,10)*DQDX(I, 9)+QI(I,11)*DQDY(I, 9)+QI(I, 8)*DQDZ(I, 9)+
     4    QI(I, 9)*DQDX(I,10)+QI(I, 7)*DQDY(I,10)+QI(I,11)*DQDZ(I,10))
     5 + S(I,2)*
     6   (QI(I, 9)*DQDX(I,11)+QI(I, 7)*DQDY(I,11)+QI(I,11)*DQDZ(I,11)+
     7    QI(I,10)*DQDX(I, 7)+QI(I,11)*DQDY(I, 7)+QI(I, 8)*DQDZ(I, 7)+
     8    QI(I, 9)*DQDX(I,11)+QI(I, 7)*DQDY(I,11)+QI(I,11)*DQDZ(I,11))
     9 + S(I,3)*
     A   (QI(I,10)*DQDX(I,11)+QI(I,11)*DQDY(I,11)+QI(I, 8)*DQDZ(I,11)+
     B    QI(I,10)*DQDX(I,11)+QI(I,11)*DQDY(I,11)+QI(I, 8)*DQDZ(I,11)+
     C    QI(I, 9)*DQDX(I, 8)+QI(I, 7)*DQDY(I, 8)+QI(I,11)*DQDZ(I, 8)) )
  200    CONTINUE
C
C 1-c. Mellor and Herring (AIAA Journal, 1973, Vol. 11, pp. 590-599)
C      also recommended by Speziale as 'Isotropized' form 
C
      ELSE IF (IDTMDL .EQ. IDTMH) THEN
         DO 300 I = 1, NPTS
C
            RHO    = QI(I,1)
            TKE    = 0.5E0 * (QI(I,6) + QI(I,7) + QI(I,8))
            EPSLN  = QI(I,12)
C
C      VRS : Coefficient for turbulent diffusion in Tau_ij
C
            VRS    = CSUBS * RHO * TKE * TKE / EPSLN
C
C      Calculate F(Turbulent diffusion)
C
           F(I,6) = F(I, 6) 
     1            + VRS * ( S(I,1)*(DQDX(I, 6)+DQDX(I, 6)+DQDX(I, 6))
     2                    + S(I,2)*(DQDY(I, 6)+DQDX(I, 9)+DQDX(I, 9))
     3                    + S(I,3)*(DQDZ(I, 6)+DQDX(I,10)+DQDX(I,10)) )
           F(I,7) = F(I, 7)
     1            + VRS * ( S(I,1)*(DQDX(I, 7)+DQDY(I, 9)+DQDY(I, 9))
     2                    + S(I,2)*(DQDY(I, 7)+DQDY(I, 7)+DQDY(I, 7))
     3                    + S(I,3)*(DQDZ(I, 7)+DQDY(I,11)+DQDY(I,11)) )
           F(I,8) = F(I, 8)
     1            + VRS * ( S(I,1)*(DQDX(I, 8)+DQDZ(I,10)+DQDZ(I,10))
     2                    + S(I,2)*(DQDY(I, 8)+DQDZ(I,11)+DQDZ(I,11))
     3                    + S(I,3)*(DQDZ(I, 8)+DQDZ(I, 8)+DQDZ(I, 8)) )
           F(I,9) = F(I, 9)
     1            + VRS * ( S(I,1)*(DQDX(I, 9)+DQDY(I, 6)+DQDX(I, 9))
     2                    + S(I,2)*(DQDY(I, 9)+DQDY(I, 9)+DQDX(I, 7))
     3                    + S(I,3)*(DQDZ(I, 9)+DQDY(I,10)+DQDX(I,11)) )
           F(I,10)= F(I,10)
     1            + VRS * ( S(I,1)*(DQDX(I,10)+DQDZ(I, 6)+DQDX(I,10))
     2                    + S(I,2)*(DQDY(I,10)+DQDZ(I, 9)+DQDX(I,11))
     3                    + S(I,3)*(DQDZ(I,10)+DQDZ(I,10)+DQDX(I, 8)) )
           F(I,11)= F(I,11)
     1            + VRS * ( S(I,1)*(DQDX(I,11)+DQDZ(I, 9)+DQDY(I,10))
     2                    + S(I,2)*(DQDY(I,11)+DQDZ(I, 7)+DQDY(I,11))
     3                    + S(I,3)*(DQDZ(I,11)+DQDZ(I,11)+DQDY(I, 8)) )
  300    CONTINUE
C
C Error - Output error message and stop
C
      ELSE
         WRITE (IOUT,1000) IDTMDL
         STOP
      ENDIF
C
C Finished with Turbulent Diffusion Terms
C
      RETURN
 1000 FORMAT (' ','ERROR - RUN ABORTING',
     1       /' ','Invalid Turbulent diffusion model in routine ',
     2            'RSTDIF  = ',I10)
      END






      SUBROUTINE REALIZ (ITURB, IDIM, JDIM, KDIM, IPSTRT, IPEND, Q, 
     1                   NRELIZ)
C
C Routine to apply realizability to the Reynolds stresses based on 
C Schwartz's inequality (see Lumley and Panofsky, The Structure of 
C Atmospheric Turbulence).
C
C This routine is organized as follows
C     1. Schwartz's inequality for Tau_xy
C        a. Calculate magnitude of shear stress from Schwartz's inequality
C        b. Replace magnitude of shear stress if it exceeds above value
C        c. Output statistics
C     2. If 3D calculation then, Schwartz's inequality for Tau_xz and Tau_yz
C        a. Calculate magnitude of shear stress from Schwartz's inequality
C        b. Replace magnitude of shear stress if it exceeds above value
C        c. Output statistics
C
C ITURB          : Turbulence model
C IDIM,JDIM,KDIM : Dimensions of current block
C IPSTRT         : I Plane (cell center) to start calculations at 
C                  (Minimum Value: 2)
C IPEND          : I Plane (cell center) to end caclulations at
C                  (Maximum Value: IDIM)
C Q              : Primitive variables at cell centers
C NRELIZ         : Number of times violate realizability
C
      include '../header/common.h'
C
      DIMENSION Q     (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C
      DIMENSION NRELIZ(NQ)
C
C 1.   Schwartz's inequality for Tau_xy
C 1-a. Calculate magnitude of shear stress from Schwartz's inequality
C 1-b. Replace magnitude of shear stress if it exceeds above value
C
      DO 130 I = IPSTRT, IPEND
         DO 120 K = 2, KDIM
            DO 110 J = 2, JDIM
               TESTXY = SQRT (Q(J,K,I,6) * Q(J,K,I,7))
               IF (ABS(Q(J,K,I,9)) .GT. TESTXY) THEN
                  NRELIZ( 9) = NRELIZ( 9) + 1
                  Q(J,K,I,9) = SIGN (TESTXY, Q(J,K,I,9))
               ENDIF
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
C 2.   If 3D calculation then, Schwartz's inequality for Tau_xz and Tau_yz
C 2-a. Calculate magnitude of shear stress from Schwartz's inequality
C 2-b. Replace magnitude of shear stress if it exceeds above value
C
      IF (THREED) THEN
         DO 230 I = IPSTRT, IPEND
            DO 220 K = 2, KDIM
               DO 210 J = 2, JDIM
                  TESTXZ = SQRT (Q(J,K,I,6) * Q(J,K,I,8))
                  IF (ABS(Q(J,K,I,10)) .GT. TESTXZ) THEN
                     NRELIZ(10)  = NRELIZ(10) + 1
                     Q(J,K,I,10) = SIGN (TESTXZ, Q(J,K,I,10))
                  ENDIF
C
                  TESTYZ = SQRT (Q(J,K,I,7) * Q(J,K,I,8))
                  IF (ABS(Q(J,K,I,11)) .GT. TESTYZ) THEN
                     NRELIZ(11)  = NRELIZ(11) + 1
                     Q(J,K,I,11) = SIGN (TESTYZ, Q(J,K,I,11))
                  ENDIF
  210          CONTINUE
  220       CONTINUE
  230    CONTINUE
      ENDIF
C
C Finished with Realizability constraints
C
      RETURN
      END
