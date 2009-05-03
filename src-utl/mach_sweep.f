C
C     Driver routine to run flow solver for different mach numbers
C
C     Specify min, max and increment for the mach number.
C     Cl and Cd will be written into mach.dat
C
      program mach_sweep
      implicit none
      character*64 flosolve, sdummy
      real      Cd, Cl
      integer   ifid, ifm, iter, i, imax
      real      rdummy
      real      mach_min, mach_max, dmach, mach

      print*,'Mach number: min, max, increment'
      read*, mach_min, mach_max, dmach

      if(mach_max.le.mach_min .or. dmach.le.0.0)then
         print*,'Error: give mach_max > mach_min'
         print*,'            dmach > 0.0'
         stop
      endif

      flosolve='$NUWTUN_HOME/src-flo/nuwtun_flo < flo.in > flo.log'
      ifid    = 10
      ifm     = 11

c     Check file exists
      open(ifid, file='flo.in.nomach', status='old')
      close(ifid)

      imax    = (mach_max - mach_min)/dmach
      open(ifm, file="mach.dat")

      do i = 0, imax

         mach = mach_min + dmach * i

         open(ifid, file='flo.in')
         write(ifid,*)'MACH    ', mach
         close(ifid)

         call system("cat flo.in.nomach >> flo.in")
         call system(flosolve)

         open(ifid,file='fort.19', status='old')
         read(ifid,*) sdummy, iter
         read(ifid,*) sdummy, rdummy, rdummy, rdummy, rdummy, rdummy
         read(ifid,*) sdummy, Cl
         read(ifid,*) sdummy, Cd
         close(ifid)
         print*,'mach, Cl, Cd =', mach, Cl, Cd
         write(ifm,*) mach, Cl, Cd

      enddo
      close(ifm)

      stop
      end program

