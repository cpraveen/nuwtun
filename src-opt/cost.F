      program main
      include '../src-adj/header/refval.h'
      character*32 str

      if(iargc().ne.1)then
         print*,'cost: Specify cost function type'
         stop
      endif

c     Read reference values
      open(10, file='refval.dat', status='old')
      read(10,*)
      read(10,*)
      read(10,*) str, clref
      read(10,*) str, cdref
      close(10)

c     Read current values
      open(10, file='fort.19', status='old')
      read(10,*)
      read(10,*)
      read(10,*) str, cl
      read(10,*) str, cd
      close(10)

c     Set cost function type
      call sicftyp(icftyp)

      cost = 0.0
      call costfun(icftyp, cl, cd, apgint, cost)
      print*,'cd/cd0 =', cd/cdref
      print*,'cl/cl0 =', cl/clref
      print*,'cost = ', cost

      stop
      end
