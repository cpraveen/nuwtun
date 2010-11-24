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



      COMMON /NOYES/  INO, IYES
C Global variables for NO or YES
C     INO    : No  (0)
C     IYES   : Yes (1)
C
      COMMON /CONST/  RMAX, RSMALL, RSMASQ
C Constants
C     RMAX   : Maximum real number - infinity for the current computer
C     RSMALL : Small real number to prevent division by zero, etc.
C     RSMASQ : Square of small number used as RSMALL but when need a square
C
      COMMON /UNITS/  INPUT, IOUT, IRDRST, IWRRST, IGRID, 
     1                IPLT3G, IPLT3Q, IPLT3F, IPLT3FN, IRDBC,
     2                ISTOPFL, IRESID, ICLCD, IRESULT
C Units to read/write from:
C     INPUT  : Unit to read input
C     IOUT   : Unit for standard output
C     IRDRST : Unit to read restart file from
C     IWRRST : Unit to write restart file to
C     IGRID  : Unit to read grid from
C     IPLT3G : Unit to write PLOT3D grid data to
C     IPLT3Q : Unit to write PLOT3D conserved variable data to
C     IPLT3F : Unit to write PLOT3D function file data to
C     IPLT3FN: Unit to write PLOT3D function names to
C     IRDBC  : Unit to read input profile for PROFILE bc
C     ISTOPFL: Unit to test for STOP file
C     IRESID : Unit to write flow residue for convergence monitoring
C     ICLCD  : Unit to write CL and CD convergence history
C
      COMMON /FILDAT/ IGTSEP, IGTP3D
C File type data
C     IGTSEP : Grid file type has separate files for each block
C     IGTP3D : Grid file type is single/multiblock Plot3D format
C
      COMMON /INPUT/  NQ, NF, NP, NRANK, NW, NMDL, 
     1                NTMLVL, NTIME, NTMTAU
C Storage variables
C     NQ     : Number of variables in the Q variable
C     NF     : Number of fluxes
C     NP     : Number of properties stored in PROPS
C     NRANK  : Rank of matrices to invert
C     NW     : Number of storage locations for work array
C     NMDL   : Number of model controls in ITURB
C     NTMLVL : Total number of additional time levels for unsteady calculation
C              NTMLVL = NTIME + NTMTAU
C     NTIME  : Number of additional time levels for unsteady calculation
C              = 1 for steady
C              = 2 for iterative, implicit unsteady
C     NTMTAU : Number of additional levels of storage for Q^(m) for pseudo
C              time stepping
C              = 1 for 2nd order time accurate with 2nd order pseudo time subit
C              = 0 for all other cases
C
C                                                         NQ    NF  NRANK
C     ---------------------------------------------------------------------
C     Perfect Gas (inviscid, laminar, algebraic turb.)     5     5      5
C     Perfect Gas - Two Equation Turbulence Model          7     7      7
C     Perfect Gas - Reynolds stress Turbulence Model      12    12     12
C
      COMMON /FLUID/  GAMMA, GAMM1, GSQM1
C Ratio of specific heats for perfect gas
C     GAMMA : Gamma
C     GAMM1 : Gamma - 1
C     GSQM1 : (Gamma**2) - 1
C
      COMMON /FLOW/   FSMACH, ALPHA, BETA, RE, PR, PRT, TINF, TWALL,
     1                AREARF, ALENRF
C Flow definition quantities
C     FSMACH : Free Stream Mach Number
C     ALPHA  : Angle of attack relative to the x-axis
C     BETA   : Side-Slip Angle
C     RE     : Reynolds Number per unit length = (rhoinf*unif/muinf)
C     PR     : Prandtl Number
C     PRT    : Turbulent Prandtl number
C     TINF   : Free Stream Temperature
C     TWALL  : Wall Temperature
C     AREARF : Wetted reference area used in calculation of coefficients
C     ALENRF : Reference length used in calculation of moment coefficients
C
      COMMON /QREF/   RHOINF, UINF, VINF, WINF, PREF, EINF, AINF, UREF
C Non-Dimensional Free Stream Flow Quantities
C     RHOINF : Density                    
C     UINF   : Velocity in the x-direction = u(inf) / ainf
C     VINF   : Velocity in the y-direction = v(inf) / ainf
C     WINF   : Velocity in the z-direction = w(inf) / ainf
C     PREF   : Pressure                    = p(inf) / (rhoinf * ainf**2)
C     EINF   : Total Energy                
C     AINF   : Speed of Sound
C     UREF
C
      COMMON /SHAT  / SHATX,  SHATY,  SHATZ
C Unit normal to use for collapsed cell faces
C     SHATX, SHATY, SHATZ: x,y,z components of unit normal
C     SHAT =  (1,1,1)/sqrt(3) for 3D
C             (1,1,0)/sqrt(2) for 2D
C
      COMMON /TRBREF/ TKEINF, EPSINF, OMEINF, TAUINF(6)
C Non-Dimensional Free Stream Turbulent Quantities
C     TKEINF : Turbulent kinetic energy    = 
C     EPSINF : Dissipation rate
C     OMEINF : Wilcox's omega variable (~epsilon/k)
C     TAUINF : Reynolds stresses (currently assumed isotropic at infinity)
C
      COMMON /QJET/   RHOJET, UJET, VJET, WJET, PJET, TJET,
     1                TKEJET, EPSJET, OMEJET, TAUJET(6)
C Jet Flow Quantities for JET boundary condition (Non-Dimensionalized
C by freestream quantities)
C     RHOJET : JET Density                     = rho(jet) / rhoinf
C                                              = GAMMA * PJET / TJET
C     UJET   : JET Velocity in the x-direction = u(jet) / ainf
C     VJET   : JET Velocity in the y-direction = v(jet) / ainf
C     WJET   : JET Velocity in the z-direction = w(jet) / ainf
C     PJET   : Jet Pressure                    = p(jet) / (rhoinf * ainf**2)
C     TJET   : Jet Temperature                 = T(jet) / Tinf
C
C     TKEJET : Jet turbulent kinetic energy
C     EPSJET : Jet dissipation rate
C     OMEJET :
C     TAUJET : Jet Reynolds stresses
C
      COMMON /QSUBBC/ PBAKBC, PTOTBC, TTOTBC
C Flow Quantities for SUBSONIC INFLOW and SUBSONIC OUTFLOW boundary conditions
C (non-dimensional variables)
C     PBAKBC : Back  pressure    for subsonic outflow
C     PTOTBC : Total pressure    for subsonic inflow
C     TTOTBC : Total temperature for subsonic inflow
C
      COMMON /BOUND/  ITAN,   IWALL,  IWALFN, IFAR,   IFAR2D, IEXT, 
     1                IFIX,   IJET,   ISYMXY, ISYMXZ, ISYMYZ, ISING, 
     2                IPERD,  IPROFL, IPROSB, IINFLO, IOUTFL, IAXISM, 
     3                IAXICL, IHSHR
C Boundary Condition Types
C     ITAN   : Tangency (Inviscid Wall)
C     IWALL  : Viscous Wall
C     IWALFN : Viscous Wall with Wall Functions
C     IFAR   : Farfield boundary
C     IFAR2D : 2D Farfield boundary condition with point vortex correction
C     IEXT   : Extrapolation boundary condition (first order)
C     IFIX   : Boundary condition fixed at freestream conditions
C     IJET   : Jet boundary condition (for jets or shear layers)
C     ISYMXY : Symmetry to the xy plane
C     ISYMXZ : Symmetry to the xz plane
C     ISYMYZ : Symmetry to the yz plane
C     ISING  : Singularity continuation
C     IPERD  : Periodic boundary condition in a coordinate direction
C     IPROFL : Profile boundary condition - input data from a file
C     IPROSB : Profile boundary condition for subsonic flows. Input data
C              from file and then use subsonic condition for pressure.
C     IINFLO : Subsonic inflow  boundary condition
C     IOUTFL : Subsonic outflow boundary condition
C     IAXISM : Axisymmetric symmetry   boundary condition
C     IAXICL : Axisymmetric centerline boundary condition
C     IHSHR  : Homogeneous shear periodic boundary condition
C
      COMMON /TIMEST/ CFL, DT, CFLFNL, ITUPDT, ITDFNL, TIME, TIMEDT
C Time Step Control
C     CFL    : CFL number for pseudo-time iteration
C     DT     : Time step  for pseudo-time iteration or non-iterative Runge-Kutta
C     CFLFNL : CFL number at iteration = ITDFNL
C     ITUPDT : Number of iterations at which to update the time step
C     ITDFNL : Final iteration at which CFL becomes CFLFNL
C     TIME   : Non-dimensional time (to be updated with solution when running
C              in time accurate mode)
C     TIMEDT : Time step for iterative, implicit time accurate
C
      COMMON /SCHEME/ IFROE, IFVL, IFHS,
     1                IFDS,  IFVS, IFCD,
     2                IFCHAR(3,1)
C Scheme Types
C     IFROE  : Roe's approximate Riemann solver flux evaluation 
C     IFVL   : van Leer's flux vector splitting flux evaluation
C     IFHS   : van Leer's flux vector splitting as modified by
C              Hanel and Schwane (AIAA 87-1105) flux evaluation
C     IFDS   : identifies a scheme as a Flux Difference Splitting scheme
C     IFVS   : identifies a scheme as a Flux Vector Splitting scheme
C     IFCD   : identifies a scheme as a Central Difference scheme
C     IFCHAR : identifies the characteristics of a scheme
C              IFCHAR(scheme,1) = IFDS if Flux Difference Splitting
C                               = IFVS if Flux Vector Splitting
C                               = IFCD if Central Difference
C
      COMMON /LMTNG/  ILNONE, ILSMTH, ILMNMD, ILVNKT
C Limiters:
C     ILNONE : Unlmiited scheme
C     ILSMTH : Smooth limiter
C     ILMNMD : MinMod limiter
C     ILVNKT : Venkatakrishnan's smooth limiter
C
      COMMON /NUMERX/ CNTRPY
C Numerical Constants
C     CNTRPY : Coefficient for entropy fix (used in froe)
C
      COMMON /DISSIP/ ADKP46
C Constants for Artificial Dissipation in Central Difference Scheme
C     ADKP46 : Coefficient for high order dissipation
C              2nd order scheme -> 4th difference dissipation
C              4th order scheme -> 6th difference dissipation (hyperviscosity)
C
      COMMON /MGCOMM/ CSMTH, CMXCHG
C Multigrid Control
C     CSMTH  : Coefficient for implicit smoothing in prolongation   (~ .25)
C     CMXCHG : Maximum relative change allowed for multigrid update (~ .25)
C
      COMMON /TRBMDL/ IEPSLN, ICMUST, ITQTAU, ITQDOT,
     1                ITNVSD, ITLMNR, ITBLMX, ITKE, ITKW, ITRS, ITLES,
     2                IFMHR,  IFMSAA, IFMZSG, IFMVD,
     3                IEEHR,  IEESAA, IEEZSG, IEEZSL, IEERNG, IEES95,
     4                IEEABD,
     5                IECNO,  IECSAR, IECRIS,
     6                IPDNO,  IPDSAR, IPDRIS,
     7                IADRNO, IADRGS,
     8                IASMBU, IASMGS,
     9                IPLRR1, IPLRR2, IPSSG, IPGL, IPLS, IPFLT, IPSO95,
     A                IEISO,  IELAIS, IEFLT, IEADGS, IESO95,
     B                IDMSO,  IDMSPZ,
     C                IDTDH,  IDTHL,  IDTMH,
     D                ISGSMG, ISGSFM,
     E                MUTALG, IFWALF, FOTURB, POSPRD
      LOGICAL         MUTALG, IFWALF, FOTURB, POSPRD
C Turbulence model data:
C     IEPSLN : Location in Q     array for epsilon (turb length scale, 7 or 12)
C     ICMUST : Location in PROPS array for CMUSTR (variable CMUSTR from ASM)
C     ITQTAU : Location in PROPS array for tau_wall  in wall functions 
C     ITQDOT : Location in PROPS array for qdot_wall in wall functions 
C              ITQTAU and ITQDOT are stored in unused locations of PROPS 
C              corresponding to y+ and y_normal as they are needed only at
C              block surfaces they use the 0 location which is not needed.
C Turbulence model types
C     ITNVSD : NO Turbulence model - Inviscid flow
C     ITLMNR : NO Turbulence model - Laminar flow
C     ITBLMX : Baldwin-Lomax algebraic eddy viscosity method
C     ITKE   : k-epsilon turbulence model
C     ITKW   : k-omega turbulence model
C     ITRS   : Reynolds stress turbulence model
C     ITLES  : Large-Eddy Simulation (LES)
C Fmu damping models: used in K-epsilon, Reynolds stress, and LES models
C                                                      K-E   Tau_ij    LES
C     IFMHR  : High Reynolds number form (Fmu = 1)     Yes    Yes      Yes
C     IFMSAA : Speziale-Abid-Anderson form             Yes    No       No
C     IFMZSG : Zhang-So-Gatski-Speziale form           Yes    Yes      No
C     IFMVD  : van Driest wall damping                 No     No       Yes
C Epsilon equation models:
C     IEEHR  : High Reynolds number form (no damping or additional terms)
C     IEESAA : Speziale-Abid-Anderson form
C     IEEZSG : Zhang-So-Gatski-Speziale form
C     IEEZSL : Zhang-So-Speziale-Lai form (AIAA J. Vol. 31, No. 1, Jan. 1993)
C     IEERNG : RNG form (Yakhot,Orszag,Thangam,Gatski,Speziale Phys.Fl.A 1992)
C     IEES95 : So et al. 1995 form
C     IEEABD : Abid's ASM (Abid,Morrison,Gatski,Speziale AIAA 95-0565)
C Compressible Dissipation models: (IEC)
C     IECNO  : No compressible dissipation corrections
C     IECSAR : Sarkar's compressible dissipation
C     IECRIS : Ristorcelli's compressible dissipation
C Pressure-Dilatation models: (IPD)
C     IPDNO  : No Pressure-Dilatation model
C     IPDSAR : Sarkar's Pressure-Dilatation model
C     IPDRIS : Ristorcelli's Pressure-Dilatation model
C Anisotropic Dissipation Rate Model
C     IADRNO : NO Anisotropic dissipation rate model (standard model)
C     IADRGS : Gatski-Speziale anisotropic dissipation rate model
C Algebraic Stress Models:
C     IASMBU : Boussinesq (standard eddy-viscosity)
C     IASMGS : Gatski-Speziale algebraic stress model (JFM Sept. 1993)
C Pressure-Strain Correlation Models:
C     IPLRR1 : Launder-Reece-Rodi model
C     IPLRR2 : Launder-Reece-Rodi model
C     IPSSG  : Speziale-Sarkar-Gatski model
C     IPGL   : Gibson-Launder model
C     IPLS   : Lai-So model
C     IPFLT  : Fu-Launder-Tselepidakis model
C     IPSO95 : So et al. 1995 model
C Epsilon_ij models:
C     IEISO  : Isotropic (high Reynolds number)
C     IELAIS : Lai-So model
C     IEFLT  : Fu-Launder-Tselepidakis model
C     IEADGS : Gatski-Speziale Anisotropic Dissipation Rate model for eps_ij
C     IESO95 : So et al. 1995 model
C Molecular Diffusion Models:
C     IDMSO  : So's model
C     IDMSPZ : Speziale's model
C Turbulent Diffusion Models:
C     IDTDH  : Daly-Harlow model
C     IDTHL  : Hanjalic-Launder 1972 model
C     IDTMH  : Mellor-Herring 1973 model 
C Sub-Grid Scale Stress Models:
C     ISGSMG : Smagorinsky model
C     ISGSFM : Structure Function model
C
C     MUTALG : Logical - update the eddy viscosity using algebraic model
C              this is used in the early iteration history to calculate
C              the eddy viscosity with the algebraic model while the
C              higher order methods adjust to the flow
C     IFWALF : Logical - use wall functions
C     FOTURB : Logical - use first order advection on turbulence equations
C              if this is true, else use higher order if false
C     POSPRD : Logical - TRUE  => Enforce positivity of turb. production 
C                        FALSE => Do not enforce positivity of production
C
      COMMON /TRBCON/ SIGK, SIGT2, SIGRHO, CPDFRO, CPDFT, CPVELK, 
     1                CMU, CEPSLN, CEPS1, CEPS2, CEPS3, CEPS4, CEPS5,
     2                A2KE, ALF1, ALF2, ALF3,
     3                GKW, BSTRKW, BKW, ROUGHK,
     4                CSUBS, PSC1, PSC2, PSCW, 
     5                PSALFA, PSBETA, PSGAMA,
     6                PRDLIM, PRDE,   PRDEM1
C Turbulence model constants 
C     SIGK   : (mu + sigk   * mut)   dk/dx(j)
C     SIGT2  : (mu + sigt2  * mut)  dt2/dx(j)
C     SIGRHO : (mu + sigrho * mut) drho/dx(j)
C     CPDFRO : Coefficient for density gradient terms in the
C              Pressure-Diffusion Correlation diffusion terms in K equation:
C              = 0   do not include terms
C              = 1/sigma_rho for the correct coefficient
C     CPDFT  : Coefficient for temperature gradient terms in the
C              Pressure-Diffusion Correlation diffusion terms in K equation:
C              = 0   do not include terms
C              = 1/sigma_T[PRT] for the correct coefficient
C     CPVELK : Coefficient for additional Pressure-Velocity Correlation
C              source terms in K equation:
C              = 0   do not include terms
C              = 1   include terms
C
C K-EPSILON MODEL->
C     CMU    : Cmu    - coefficient in definition of turbulent eddy viscosity
C     CEPS1  : Ceps_1 - coefficient on production term in epsilon equation
C              pe  = Ceps_1 * Pk * epsilon / k
C       where  Pk  = mut {1/2[dU(i)/dx(j)+dU(j)/dx(i))-2/3dU(k)/dx(k)]**2}
C                        - 2/3*rho*k*dU(k)/dx(k)
C     CEPS2  : Ceps_2 - coefficient on source term in epsilon equation
C              se  = Ceps_2 * rho * epsilon**2 / k
C     CEPS3-5: Additional constants for epsilon transport equation
C     A2KE   : Constant in wall damping of Speziale et al. AIAA 90-1481
C     ALF1KE : Constant for Sarkar's compressibility correction (JFM)
C              = 0   no compressibility correction
C              = 1.0 correction constant from Speziale&Sarkar ICASE 91-9
C     ALF2KE : Constant in Sarkar's pressure-dilatation model ICASE 91-42
C              = 0   no production term in pressure-dilatation model
C              = 0.4 constant from Sarkar ICASE 91-42
C     ALF3KE : Constant in Sarkar's pressure-dilatation model ICASE 91-42
C              = 0   no dissipation term in pressure-dilatation model
C              = 0.2 constant from Sarkar ICASE 91-42
C     SIGK   : 1 / Sigma_k
C     SIGT2  : 1 / Sigma_epsilon
C     SIGRHO : 1 / Sigma_rho
C     CPDIFP : 0 or 1/sigma_rho (see above)
C     CPDIFT : 0 or (1/PRT - 1/Sigma_rho) (see above)
C     CPVELK : 0 or 1 (see above)
C
C K-OMEGA MODEL-> Model from Wilcox, AIAA Journal, Vol. 26, No. 11
C     CMU    : Gamma* - coefficient in definition of turbulent eddy viscosity
C              mut = gamma* * rho * k / omega
C     GKW    : Gamma  - coefficient on omega equation
C              tw  = gamma * omega / k * tau(ij) * du(i)/dx(j)
C     BSTRKW : Beta*  - coefficient on source term in k equation 
C              sk  = beta* * rho * omega * k
C     BKW    : Beta   - coefficient on source term in omega equation
C              sw  = beta * rho * omega**2
C     ROUGHK : Kr - Grain height for roughness elements in boundary conditions
C     SIGK   : sigma*
C     SIGT2  : sigma
C     SIGRHO : = 0 (not included in this model)
C     CPDIFP : = 0 (not included in this model)
C     CPDIFT : = 0 (not included in this model)
C     CPVELK : = 0 (not included in this model)
C
C REYNOLD'S STRESS MODEL->
C     CSUBS  : Coefficient for triple velocity correlation
C   Pressure-Strain Correlation Coefficients:
C     PSC1   : Coefficient on rho epsilon b_ij term
C     PSC2   : Coefficient on (P_ij - P_kk delta_ij / 3) term
C     PSCW   : Coefficient for pressure-echo term in Lai-So model
C     PSALFA : Coefficient for Shima's Pressure-Strain Correlation
C     PSBETA : Coefficient for Shima's Pressure-Strain Correlation
C     PSGAMA : Coefficient for Shima's Pressure-Strain Correlation
C
C     PRDLIM : Limit on turbulent production
C              k-epsilon model : Prod = min (Prod, prdlim*rho*epsilon)
C     PRDE   : Production limiter for the epsilon/omega equation = (0, 1 ONLY)
C              = 1 then apply production limiter to the epsilon/omega equation
C              = 0 then do not apply production limiter to the epsilon/omega eq
C     PRDEM1 : (1 - PRDE) to avoid additional floating point operations
C
      COMMON /TRBFIX/ IFIXQ , 
     1                IFIXMN, IFIXAV
C Turbulence quantity positivity control
C                                                          -- Valid Values ---
C     IFIXQ  : Control of fix to Q variables (call FIXQ ): INO, IFIXMN, IFIXAV
C     IFIXMN : Minimum fix (set Q = Q_min)
C     IFIXAV : Average fix (set Q = average of neighbors)
C
      COMMON /SOLVER/ ISOLVR,
     1                IAF3F , ILGS, IRKN, IMARCH,
     2                IMPSRC,
     3                ISDIAG, ISBLOC,
     4                ITIMED,
     5                ISTDY , IUNSTD,
     6                ITMTYP,
     7                ITTS  , ITAUTS
C Solution procedure control
C     ISOLVR : Current solution scheme (= IAF3F, ILGS, IRKN, or IMARCH)
C     IAF3F  : Three factor spatially split approximate factorization
C     ILGS   : Line Gauss Seidel
C     IRKN   : N stage Runge-Kutta
C     IMARCH : Marching procedure
C     IMPSRC : Current implicit source term scheme (= ISDIAG or ISBLOC)
C     ISDIAG : Diagonalized implicit source terms
C     ISBLOC : Block        implicit source terms
C     ITIMED : Time dependent scheme               (= ISTDY  or IUNSTD)
C     ISTDY  : Steady   calculation - don't include dQ/dt time accurate
C     IUNSTD : Unsteady calculation - do    include dQ/dt time accurate
C     ITMTYP : Time dependendt scheme type         (= ITTS   or ITAUTS)
C     ITTS   : t-ts physical time subiteration
C     ITAUTS : tau-ts psuedo time subiteration
C
      COMMON /UNSTDY/ TDTHET, TDPHI, TDPHIP,
     1                TDPDTJ, TDPDEL
C Constants for unsteady calculation for dQ/dt term
C     t-ts: Physical Time Subiteration (Rumsey et al. Computers&Fluids, 1996)
C           dQ/dt = [(1 + phi) Q^(n+1) - (1 + 2 phi) Q^(n) + phi Q^(n-1)] / dt
C
C     tau-ts: Pseudo Time Subiteration (Rumsey et al. Computers&Fluids, 1996)
C           dQ/dt = [(1 + phi) Q^(n+1) - (1 + 2 phi) Q^(n) + phi Q^(n-1)] / dt
C                 + [phi' (Q^(m) - Q^(m-1))] / dtau
C
C     TDTHET : Theta for Time Dependent calculation (physical time term)
C     TDPHI  : Phi   for Time Dependent calculation (physical time term)
C     TDPHIP : Phi'  for Time Dependent calculation (pseudo   time term)
C     TDPDTJ : Coefficient of I/(DTJ)       term in implicit scheme
C     TDPDEL : Coefficient of I/(delta t J) term in implicit scheme
C +-------+-----------------+-------------+-------------------+---------------+
C |Scheme | Order           | Theta  Phi  |  Phi' store dQ^(m)| TDPDTJ TDPDEL |
C +-------+-----------------+-------------+-------------------+---------------+
C |t-ts   | 1st             |   1     0   |   -         No    | 1+Phi     0   |
C |       | 2nd             |   1    1/2  |   -         No    | 1+Phi     0   |
C |tau-ts | 2nd, 1st pseudo |   1    1/2  |   0         No    | 1+Phi'  1+Phi |
C |       | 2nd, 2nd pseudo |   1    1/2  |  1/2       Yes    | 1+Phi'  1+Phi |
C +-------+-----------------+-------------+-------------------+---------------+
C
      COMMON /CONTRL/ GLOBAL, UPDATJ, THNLYR, THREED, AXISYM,
     1                TWSPE,  SIUNIT, ISRCE,  FOURTH, INITTD
      LOGICAL         GLOBAL, UPDATJ, THNLYR, THREED, AXISYM,
     1                TWSPE,  SIUNIT, ISRCE,  FOURTH, INITTD
      COMMON /CONVAL/ AXIDTH
C
C Logical Control variables
C     GLOBAL :
C     UPDATJ : Control of updating of Jacobians in Marching scheme
C     THNLYR : Control of thin layer (.TRUE.) or full Navier-Stokes (.FALSE.)
C     THREED : 3-Dimensional run (.TRUE.)
C     AXISYM : Axi-symmetric run (.TRUE.)
C     TWSPE  : Wall temperature specified (.TRUE.) or adiabatic wall (.FALSE.)
C     SIUNIT : SI Unit input (.TRUE.) or English unit input (.FALSE.)
C     ISRCE  : Control of source term calculation
C              ISRCE = .TRUE.   calculate source terms
C                    = .FALSE.  don't calculate source terms
C     FOURTH :
C     INITTD : Initialize time dependent
C              INITTD = .TRUE.  read steady restart and initialize time depend.
C              INITTD = .FALSE. do not initialize time dependent. Read unsteady
C                               restart file if running unsteady calculation.
C Real    Control variables
C     AXIDTH : Angle of rotation to setup axisymmetric grid.
C
      COMMON /CALIND/ IKD (3,3), ICY (3,3)
C
C Data used to calculate indices
C     IKD    : Kronecker Delta array (IKD_ij = 1 if i=j, = 0 if i not = j)
C     ICY    : Cyclical array
C
C Data used to calculate moment
C Nodes AMOMAXIS(I,1,1:3) and AMOMAXIS(I,2,1:3) define I'th axis about which 
C moment is calculated
      CHARACTER*80 AXISNAME(100)
      COMMON /MOMDATA/ AMOMAXIS(100,3,3), NMOMENT, AXISNAME
