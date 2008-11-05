c This subroutine resets some parameters read from input file to modify
c it for the adjoint solver.
      subroutine rdadj(RESTRT, NCFS, NITS, FRDRST, FWRRST, ICFTYP)
      include '../header/refval.h'
      logical      RESTRT
      dimension    NITS(*)
      character*80 FRDRST, FWRRST

      if(iargc().ne.1)then
         print*,'RDADJ: Specify cost function type'
         stop
      endif

c     Set cost function type
      call sicftyp(icftyp)

c     Some cost functions require reference values
c     we read them from file refval.dat
      if(icftyp.eq.3)then
         call readrefval
      endif

      RESTRT = .TRUE.    ! So that we read the flow solution
      FRDRST = FWRRST    ! File to read flow solution
      NCFS   = 1         ! Only one grid level at present
      NITS(1)= 2*NITS(1) ! Double no. of iterations

      end

c     Read reference value for cl and cd
c     refval.dat has same format as fort.19 produced by nuwtun_flo
      subroutine readrefval
      include '../header/refval.h'
      character*32 str

      print*,'Reading reference values ...'
      open(30, file='refval.dat', status='old')
      read(30,*)
      read(30,*)
      read(30,*) str, clref
      read(30,*) str, cdref
      close(30)

      print*,'clref =', clref
      print*,'cdref =', cdref

      end
