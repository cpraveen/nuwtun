C ======================================================================
C
C                            Main driver program
C
C ======================================================================
C
      program example
C
      implicit none
C
C     include the Ipopt return codes
C
      include 'IpReturnCodes.inc'
      include '../src-adj/header/refval.h'
C
C     Size of the problem (number of variables and equality constraints)
C
      integer     N,     M,     NELE_JAC,     NELE_HESS,      IDX_STY
      parameter  (N = 20, M = 1, NELE_JAC = 20, NELE_HESS = 20*20)
      parameter  (IDX_STY = 1 )
C
C     Space for multipliers and constraints
C
      double precision LAM(M)
      double precision G(M)
C
C     Vector of variables
C
      double precision X(N)
C
C     Vector of lower and upper bounds
C
      double precision X_L(N), X_U(N), Z_L(N), Z_U(N)
      double precision G_L(M), G_U(M)
C
C     Private data for evaluation routines
C     This could be used to pass double precision and integer arrays untouched
C     to the evaluation subroutines EVAL_*
C
      double precision DAT(2)
      integer IDAT(1)
C
C     Place for storing the Ipopt Problem Handle
C
c @BIT32FCOMMENT@C     for 32 bit platforms
c @BIT32FCOMMENT@      integer IPROBLEM
c @BIT32FCOMMENT@      integer IPCREATE
      integer IPROBLEM
      integer IPCREATE
c @BIT64FCOMMENT@C     for 64 bit platforms:
c @BIT64FCOMMENT@      integer*8 IPROBLEM
c @BIT64FCOMMENT@      integer*8 IPCREATE
C
      integer IERR
      integer IPSOLVE, IPADDSTROPTION
      integer IPADDNUMOPTION, IPADDINTOPTION
      integer IPOPENOUTPUTFILE
C
      double precision F
      integer i
      integer icftyp
C
C     The following are the Fortran routines for computing the model
C     functions and their derivatives - their code can be found furhter
C     down in this file.
C
      external EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS
C
C     Count number of function/grad evaluations
C
      integer nfeval, ngeval
      common/eval/nfeval, ngeval

      if(iargc().ne.1)then
         print*,'ipopt: Specify cost function type'
         stop
      endif

C     Set cost function type
      write(*,*)
      call sicftyp(icftyp)
      write(*,*)

      nfeval = 0
      ngeval = 0
C
C     Set initial point and bounds:
C
      do i=1,N
         X  (i) =  0.0d0
         X_L(i) = -0.01d0
         X_U(i) = +0.01d0
      enddo
C
C     Set bounds for the constraints
C
      G_L(1) = 0.0d0
      G_U(1) = 0.0d0
C
C     First create a handle for the Ipopt problem (and read the options
C     file)
C
      IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS,
     1     IDX_STY, EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS)
      if (IPROBLEM.eq.0) then
         write(*,*) 'Error creating an Ipopt Problem handle.'
         stop
      endif
C
C     Open an output file
C
      IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 5)
      if (IERR.ne.0 ) then
         write(*,*) 'Error opening the Ipopt output file.'
         goto 9000
      endif
C
C     Note: The following options are only examples, they might not be
C           suitable for your optimization problem.
C
C     Set a string option
C
      IERR = IPADDSTROPTION(IPROBLEM, 'mu_strategy', 'adaptive')
      if (IERR.ne.0 ) goto 9990
C
      IERR = IPADDSTROPTION(IPROBLEM, 'hessian_approximation',
     1                      'limited-memory')
      if (IERR.ne.0 ) goto 9990
C
C
C     Set an integer option
C
      IERR = IPADDINTOPTION(IPROBLEM, 'max_iter', 30)
      if (IERR.ne.0 ) goto 9990
C
C     Set a double precision option
C
      IERR = IPADDNUMOPTION(IPROBLEM, 'tol', 1.d-7)
      if (IERR.ne.0 ) goto 9990
C
C     As a simple example, we pass the constants in the constraints to
C     the EVAL_C routine via the "private" DAT array.
C
      DAT(1) = 0.d0
      DAT(2) = 0.d0
      IDAT(1)= icftyp
C
C     Call optimization routine
C
      IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IDAT, DAT)
C
C     Output:
C
      if( IERR.eq.IP_SOLVE_SUCCEEDED ) then
         write(*,*)
         write(*,*) 'The solution was found.'
         write(*,*)
         write(*,*) 'The final value of the objective function is ',F
         write(*,*)
         write(*,*) 'The optimal values of X are:'
         write(*,*)
         do i = 1, N
            write(*,*) 'X  (',i,') = ',X(i)
         enddo
         write(*,*)
         write(*,*) 'The multipliers for the lower bounds are:'
         write(*,*)
         do i = 1, N
            write(*,*) 'Z_L(',i,') = ',Z_L(i)
         enddo
         write(*,*)
         write(*,*) 'The multipliers for the upper bounds are:'
         write(*,*)
         do i = 1, N
            write(*,*) 'Z_U(',i,') = ',Z_U(i)
         enddo
         write(*,*)
         write(*,*) 'The multipliers for the equality constraints are:'
         write(*,*)
         do i = 1, M
            write(*,*) 'LAM(',i,') = ',LAM(i)
         enddo
         write(*,*)
      else
         write(*,*)
         write(*,*) 'An error occoured.'
         write(*,*) 'The error code is ',IERR
         write(*,*)
      endif
C
 9000 continue
C
C     Clean up
C
      call IPFREE(IPROBLEM)
      stop
C
 9990 continue
      write(*,*) 'Error setting an option'
      goto 9000
      end
C
C =============================================================================
C
C                    Computation of objective function
C
C =============================================================================
C
      subroutine EV_F(N, X, NEW_X, F, IDAT, DAT, IERR)
      implicit none
      include '../src-adj/header/refval.h'
      integer nfeval, ngeval
      common/eval/nfeval, ngeval
      integer N, NEW_X
      double precision F, X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR

      integer i, icftyp
      real    cl, cd, cost
      character*32 str

      print*,'in EV_F'
      print*,'Design variables:'
      write(*,*)(x(i),i=1,n)

c     write parameters to file
      open(10, file='hicks.in')
      write(10,'(e20.10)')(x(i),i=1,n)
      close(10)
      print*,'Deforming grid ...'
      call system("deform > def.log")

c     run flow solver
      call system("rm -f fort.19")
      print*,'Running flow solver ...'
      call system("$NUWTUN_HOME/src-flo/nuwtun_flo < flo.in > flo.log")

c     Read current values
      open(10, file='fort.19', status='old')
      read(10,*)
      read(10,*)
      read(10,*) str, cl
      read(10,*) str, cd
      close(10)

      print*,'cl =',cl
      print*,'cd =',cd

c     If first computation, set reference values
c     Create refval.dat file, needed by adjoint solver
      if(nfeval.eq.0)then
        clref = cl
        cdref = cd
        call system("cp fort.19 refval.dat")
      endif

      icftyp = IDAT(1)
      call costfun(icftyp, cl, cd, cost)
      f = cost

      nfeval = nfeval + 1

      IERR = 0
      return
      end
C
C =============================================================================
C
C                Computation of gradient of objective function
C
C =============================================================================
C
      subroutine EV_GRAD_F(N, X, NEW_X, GRAD, IDAT, DAT, IERR)
      implicit none
      integer nfeval, ngeval
      common/eval/nfeval, ngeval
      integer N, NEW_X
      double precision GRAD(N), X(N)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR

      integer i, j
      character*80 ct
      double precision F

      if(new_x)then
         print*,'EV_GRAD_F called for a new point'
         call EV_F(N, X, NEW_X, F, IDAT, DAT, IERR)
      endif

c     run adjoint solver
      call getarg(1, ct)
      print*,'Running adjoint flow solver ...'
      call system("nuwtun_adj "//ct//" < flo.in > adj.log")
      print*,'Running adjoint grid solver ...'
      call system("deform_adj > def_adj.log")

      open(10, file='gradient.dat', status='old')
      do i=1,n
         read(10,*) j, grad(i)
      enddo
      close(10)

      print*,'Gradient:'
      write(*,*)(grad(i),i=1,n)

      ngeval = ngeval + 1

      IERR = 0
      return
      end
C
C =============================================================================
C
C                     Computation of equality constraints
C
C =============================================================================
C
      subroutine EV_G(N, X, NEW_X, M, G, IDAT, DAT, IERR)
      implicit none
      integer N, NEW_X, M
      double precision G(*), X(*)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
      G(1) = 0.0d0
      IERR = 0
      return
      end
C
C =============================================================================
C
C                Computation of Jacobian of equality constraints
C
C =============================================================================
C
      subroutine EV_JAC_G(TASK, N, X, NEW_X, M, NZ, ACON, AVAR, A,
     1     IDAT, DAT, IERR)
      integer TASK, N, NEW_X, M, NZ
      double precision X(*), A(*)
      integer ACON(*), AVAR(*), I
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
C
      print*,'nz=',nz
      if(TASK.eq.0)then
         do i=1,n
            AVAR(I) = i
            ACON(I) = 1
         enddo
      else
         do i=1,n
            A(i) = 0.0d0
         enddo
      endif

      IERR = 0
      return
      end
C
C =============================================================================
C
C                Computation of Hessian of Lagrangian
C
C =============================================================================
C
      subroutine EV_HESS(TASK, N, X, NEW_X, OBJFACT, M, LAM, NEW_LAM,
     1     NNZH, IRNH, ICNH, HESS, IDAT, DAT, IERR)
      implicit none
      integer TASK, N, NEW_X, M, NEW_LAM, NNZH, i
      double precision X(*), OBJFACT, LAM(*), HESS(*)
      integer IRNH(*), ICNH(*)
      double precision DAT(*)
      integer IDAT(*)
      integer IERR
C
C     structure of Hessian:
C
      integer IRNH1(10), ICNH1(10)
      data  IRNH1 /1, 2, 2, 3, 3, 3, 4, 4, 4, 4/
      data  ICNH1 /1, 1, 2, 1, 2, 3, 1, 2, 3, 4/
      save  IRNH1, ICNH1

      print*,'in hessian'

      IERR = 1
      return
      end
