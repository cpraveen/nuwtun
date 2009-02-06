C
C     Program to compare AD and FD gradients
C
C     Number of design parameters and FD step size are hard-coded
C     nparam = Number of design parameters
C     dvar   = FD step size are hard-coded
C
      program validate
      implicit none
      include '../src-adj/header/refval.h'
      character*64 flosolve, adjsolve, deform, defadj, sdummy
      real cost(100), x(100), fd(100), dvar, cost0, Cd, Cl
      integer   i, j, icftyp, nparam
      integer   ifid, iter
      real      rdummy
      character*80 string

      if(iargc().ne.1)then
         print*,'steep: Specify cost function type'
         stop
      endif

c     Read and set cost function type
      call getarg(1, string)

c     Set cost function type
      write(*,*)
      call sicftyp(icftyp)
      write(*,*)

      flosolve='$NUWTUN_HOME/src-flo/nuwtun_flo < flo.in > flo.log'
      adjsolve=
     1 '$NUWTUN_HOME/src-adj/nuwtun_adj ' // TRIM(string) //
     2 ' < flo.in > adj.log'
      deform  ='$NUWTUN_HOME/src-grd/deform > def.log'
      defadj  ='$NUWTUN_HOME/src-grd/deform_adj > def_adj.log'
      
      nparam = 20
      dvar   = 1.0d-8

      print*,'Number of parameters =', nparam
      print*,'Perturbation         =', dvar

      ifid = 20
      open(ifid, file='hicks.in')
      write(ifid,'(e20.12)')(0.0, i=1,nparam)
      close(ifid)

      print*,'Running flow solver'
      call system(deform)
      call system(flosolve)

      open(ifid,file='fort.19', status='old')
      read(ifid,*) sdummy, iter
      read(ifid,*) sdummy, rdummy, rdummy, rdummy, rdummy, rdummy
      read(ifid,*) sdummy, Cl
      read(ifid,*) sdummy, Cd
      close(ifid)
      print*,'Cl, Cd=', Cl, Cd

c     Create file of reference values
      call system("cp fort.19 refval.dat")

      Clref = Cl
      Cdref = Cd

      call costfun(icftyp, cl, cd, cost0)

      print*,'Running adjoint flow solver'
      call system(adjsolve)
      print*,'Running adjoint grid'
      call system(defadj)

c     Perturb each variable in turn
      do i=1,nparam
         do j=1,nparam
            x(j) = 0.0d0
         enddo
         x(i) = dvar

         open(ifid, file='hicks.in')
         write(ifid,'(e20.12)')(x(j), j=1,nparam)
         close(ifid)

         print*
         print*,'Parameter =',i
         print*,'Running grid deformation ...'
         call system(deform)
         print*,'Running flow solver ...'
         call system(flosolve)

         open(ifid,file='fort.19', status='old')
         read(ifid,*) sdummy, iter
         read(ifid,*) sdummy, rdummy, rdummy, rdummy, rdummy, rdummy
         read(ifid,*) sdummy, Cl
         read(ifid,*) sdummy, Cd
         close(ifid)
         print*,'Cl, Cd=', Cl, Cd
         print*

         call costfun(icftyp, cl, cd, cost(i))
      enddo

c     One-sided FD
      do i=1,nparam
         fd(i) = (cost(i) - cost0)/dvar
      enddo

c     Save FD gradient into file
      print*,'FD gradients saved into fdgrad.dat'
      open(ifid,file='fdgrad.dat')
      do i=1,nparam
         write(ifid,*) i, fd(i)
      enddo
      close(ifid)

      print*,'AD gradients saved into gradient.dat'

      stop
      end program

