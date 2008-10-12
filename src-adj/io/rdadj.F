c This subroutine resets some parameters read from input file to modify
c it for the adjoint solver.
      subroutine rdadj(RESTRT, NCFS, FRDRST, FWRRST, ICFTYP)
      logical      RESTRT
      character*80 FRDRST, FWRRST
      character*80 string

      if(iargc().ne.1)then
         print*,'RDADJ: Specify cost function type'
         stop
      endif
c     Read and set cost function type
      call getarg(1, string)
      if(string .eq. 'CL')then
         print*,'Cost function is lift coefficient'
         icftyp = 1
      else if(string .eq. 'CD')then
         print*,'Cost function is drag coefficient'
         icftyp = 2
      else
         print*,'RDADJ: cost function type not defined'
         stop
      endif

      RESTRT = .TRUE. ! So that we read the flow solution
      FRDRST = FWRRST ! File to read flow solution
      NCFS   = 1      ! Only one grid level at present

      end
