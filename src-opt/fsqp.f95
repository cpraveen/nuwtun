
MODULE database
   implicit none

   type db
      integer :: np
      integer :: nd
      real, dimension(100,50) :: x
      character(len=48) :: dir(100)
   end type db

   type(db) :: optdb

END MODULE database

program fsqp
   use database
   implicit none
   include '../src-adj/header/refval.h'

   integer :: nd
   real,allocatable :: x(:)

   integer :: i
   integer :: iwsize, nwsize, ndmax
   parameter(iwsize=100000)
   parameter(nwsize=100000)
   parameter(ndmax=20)
   integer :: nf, nineqn, nineq, neqn, neq, mode, iprint, &
              miter, inform, iw(iwsize)
   real    :: bigbnd, eps, epseqn, udelta, &
              xopt(ndmax), f(2), g(2), w(nwsize)
   external   objfun, confun, gobjfun, gconfun
   real :: bl(ndmax), bu(ndmax)
   character(len=56) :: dir

   call readinp(miter, neqn, nineqn, nd, bl, bu)

   call system('rm -rf workdir.*')

   ! Initialize database: empty to begin with
   optdb%np = 0
   optdb%nd = nd

   allocate( x(nd) )

   do i=1,nd
      x(i) = 0.0
   enddo

   nf     = 1
   !nineqn = 0     ! This is set in readinp
   nineq  = nineqn ! We have only non-linear inequality
   !neqn   = 1     ! This is set in readinp
   neq    = neqn   ! We have only non-linear equality
   mode   = 100
   iprint = 2
   !miter  = 15    ! This is set in readinp
   bigbnd = 1.0e20
   eps    = 1.0e-3
   epseqn = 1.0e-3
   udelta = 1.0e-8

   call FFSQP(nd,nf,nineqn,nineq,neqn,neq,mode,iprint, &
              miter,inform,bigbnd,eps,epseqn,udelta,bl,bu,x, &
              f,g,iw,iwsize,w,nwsize,objfun,confun,gobjfun, &
              gconfun)

   call finddir(x, dir)
   print*,'Best solution is in directory ',TRIM(dir)

end program fsqp

! Read input parameters from file
subroutine readinp(miter, neqn, nineqn, nd, bl, bu)
   implicit none
   integer :: miter, neqn, nineqn, nd
   real    :: bl(*), bu(*)

   integer :: fid, i, clcon

   fid = 10
   open(fid, file='fsqp.in', status='old')
   read(fid,*) miter
   read(fid,*) clcon
   read(fid,*) nd
   do i=1,nd
      read(fid,*) bl(i), bu(i)
   enddo
   close(fid)

   ! Print parameters to screen
   print*,'Max number of iterations =', miter
   if(clcon.eq.1)then
      print*,'Cl constraint is  Cl  = Cl0'
      nineqn = 0
      neqn   = 1
   else if(clcon.eq.2)then
      print*,'Cl constraint is  Cl >= Cl0'
      nineqn = 1
      neqn   = 0
   else
      print*,'Unknown value for Cl constraint. Valid options:'
      print*,' 1 : Cl  = Cl0'
      print*,' 2 : Cl >= Cl0'
      stop
   endif
   print*,'Lower/upper bounds:'
   do i=1,nd
      write(*,'(i5,2e20.10)') i, bl(i), bu(i)
   enddo

end subroutine readinp

! Objective function
subroutine objfun(nd, nf, x, f)
   use database
   implicit none
   include '../src-adj/header/refval.h'
   integer :: nd, nf
   real    :: x(*), f(*)

   character(len=56) :: dir
   real              :: cl, cd
   integer           :: i, icftyp

   icftyp = 2

   call finddir(x, dir)

   if(dir == '')then
      call add2db(x, dir)
      call system('cp -r templatedir '//TRIM(dir))
      open(10, file=TRIM(dir)//'/hicks.in')
      write(10,'(e20.10)')(x(i), i=1,nd)
      close(10)
      call system('cd ' // TRIM(dir) // ' && deform > def.log')
      call system('cd ' // TRIM(dir) // &
                  ' && $NUWTUN_HOME/src-flo/nuwtun_flo < flo.in > flo.log')
      call rdresu(dir, cl, cd)
      if(TRIM(dir) == 'workdir.1')then
         clref = cl
         cdref = cd
      endif
   else
      call rdresu(dir, cl, cd)
   endif

   f(1) = cd/cdref

end subroutine objfun

! Constraint function
subroutine confun(nd, nf, x, f)
   use database
   implicit none
   include '../src-adj/header/refval.h'
   integer :: nd, nf
   real    :: x(*), f(*)

   character(len=56) :: dir
   real              :: cl, cd
   integer           :: i, icftyp

   icftyp = 1

   call finddir(x, dir)

   if(dir == '')then
      call add2db(x, dir)
      call system('cp -r templatedir '//TRIM(dir))
      open(10, file=TRIM(dir)//'/hicks.in')
      write(10,'(e20.10)')(x(i), i=1,nd)
      close(10)
      call system('cd ' // TRIM(dir) // ' && deform > def.log')
      call system('cd ' // TRIM(dir) // &
                  ' && $NUWTUN_HOME/src-flo/nuwtun_flo < flo.in > flo.log')
      call rdresu(dir, cl, cd)
      if(TRIM(dir) == 'workdir.1')then
         clref = cl
         cdref = cd
      endif
   else
      call rdresu(dir, cl, cd)
   endif

   f(1) = cl/clref - 1.0

end subroutine confun

! Gradient of Objective function
subroutine gobjfun(nd, nf, x, gf, objfun)
   use database
   implicit none
   include '../src-adj/header/refval.h'
   integer :: nd, nf
   real    :: x(*), gf(nd,*)
   external   objfun

   character(len=56) :: dir
   integer           :: i
   real              :: f(2)

   call finddir(x, dir)

   if(dir == '')then
      call objfun(nd, nf, x, f(1))
      call finddir(x, dir)
   endif
   call system('cp workdir.1/fort.19 '//TRIM(dir)//'/refval.dat')
   call system('cd '//TRIM(DIR)//' && nuwtun_adj CD < flo.in > adj_cd.log')
   call system('cd '//TRIM(DIR)//' && deform_adj > def_adj.log')
   call rdgrad(dir, nd, gf(1,1))

end subroutine gobjfun

! Gradient of constraint function
subroutine gconfun(nd, nf, x, gf, confun)
   use database
   implicit none
   include '../src-adj/header/refval.h'
   integer :: nd, nf
   real    :: x(*), gf(nd,*)
   external   confun

   character(len=56) :: dir
   integer           :: i
   real              :: f(2)

   call finddir(x, dir)

   if(dir == '')then
      call confun(nd, nf, x, f(1))
      call finddir(x, dir)
   endif
   call system('cp workdir.1/fort.19 '//TRIM(dir)//'/refval.dat')
   call system('cd '//TRIM(DIR)//' && nuwtun_adj CL < flo.in > adj_cl.log')
   call system('cd '//TRIM(DIR)//' && deform_adj > def_adj.log')
   call rdgrad(dir, nd, gf(1,1))

end subroutine gconfun

! Read gradient from file
subroutine rdgrad(dir, nd, grad)
   implicit none
   integer :: nd
   character(len=56) :: dir
   real    :: grad(*)

   integer :: i, j

   open(30, file=TRIM(dir)//'/gradient.dat', status='old')
   do i=1,nd
      read(30,*) j, grad(i)
   enddo
   close(30)

end subroutine rdgrad

! Find if point x is already present in database
subroutine finddir(x, dir)
   use database
   implicit none
   real              :: x(*)
   character(len=56) :: dir

   integer           :: i, j, stat

   dir = ''
   do i=1,optdb%np
      stat = 1
      do j=1,optdb%nd
         if(optdb%x(i,j) /= x(j))then
            stat = 0
         endif
      enddo
      if(stat == 1)then
         dir = optdb%dir(i)
         return
      endif
   enddo

end subroutine finddir

! Find if point x is already present in database
subroutine add2db(x, dir)
   use database
   implicit none
   real              :: x(*)
   character(len=56) :: dir

   integer           :: i, j, dirnum

   i = optdb%np + 1
   do j=1,optdb%nd
      optdb%x(i,j) = x(j)
   enddo

   ! set dir name
   dirnum = i
   if(dirnum <= 9)then
      write(unit=dir, fmt='(A8,I1)') 'workdir.', dirnum
   elseif(dirnum <= 99)then
      write(unit=dir, fmt='(A8,I2)') 'workdir.', dirnum
   elseif(dirnum <= 999)then
      write(unit=dir, fmt='(A8,I3)') 'workdir.', dirnum
   else
      write(*,*)'Unable to set dir name'
      stop
   endif

   optdb%dir(i) = dir

   optdb%np = optdb%np + 1

end subroutine add2db

! Read cl and cd
subroutine rdresu(dir, cl, cd)
   implicit none
   character(len=56) :: dir
   real              :: cl, cd

   character(len=32) :: str

   open(20, file=TRIM(dir)//'/fort.19', status='old')
   read(20,*)
   read(20,*)
   read(20,*) str, cl
   read(20,*) str, cd
   close(20)

end subroutine rdresu
