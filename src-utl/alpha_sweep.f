C
C     Driver routine to run flow solver for different ALPHA (AOA)
C
C     Specify min, max and increment for the AOA
C     Cl and Cd will be written into alpha.dat
C
      program alpha_sweep
      implicit none
      character*64 flosolve, sdummy
      real      Cd, Cl
      integer   ifid, ifm, iter, i, imax
      real      rdummy
      real      alpha_min, alpha_max, dalpha, alpha

      print*,'AOA: min, max, increment'
      read*, alpha_min, alpha_max, dalpha

      if(alpha_max.le.alpha_min .or. dalpha.le.0.0)then
         print*,'Error: give alpha_max > alpha_min'
         print*,'            dalpha > 0.0'
         stop
      endif

      flosolve='$NUWTUN_HOME/src-flo/nuwtun_flo < flo.in > flo.log'
      ifid    = 10
      ifm     = 11

c     Check file exists
      open(ifid, file='flo.in.noalpha', status='old')
      close(ifid)

      imax    = (alpha_max - alpha_min)/dalpha
      open(ifm, file="alpha.dat")

c     Loop over different angle of attack
      do i = 0, imax

         alpha = alpha_min + dalpha * i

c        Write angle of attack into file
         open(ifid, file='flo.in')
         write(ifid,*)'ALPHA    ', alpha
         close(ifid)
         call system("cat flo.in.noalpha >> flo.in")

c        Run CFD solver
         call system(flosolve)

c        Read solution
         open(ifid,file='fort.19', status='old')
         read(ifid,*) sdummy, iter
         read(ifid,*) sdummy, rdummy, rdummy, rdummy, rdummy, rdummy
         read(ifid,*) sdummy, Cl
         read(ifid,*) sdummy, Cd
         close(ifid)
         print*,'alpha, Cl, Cd =', alpha, Cl, Cd
         write(ifm,*) alpha, Cl, Cd

      enddo
      close(ifm)

      stop
      end program

