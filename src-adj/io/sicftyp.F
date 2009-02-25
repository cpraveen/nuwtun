      subroutine sicftyp(icftyp)
      character*80 string

c     Read and set cost function type
      call getarg(1, string)
      if(string .eq. 'CL')then
         print*,'Cost function is lift coefficient'
         icftyp = 1
      else if(string .eq. 'CD')then
         print*,'Cost function is drag coefficient'
         icftyp = 2
      else if(string .eq. 'CD/CL')then
         print*,'Cost function is drag/lift'
         icftyp = 3
      else if(string .eq. 'CDCLPENM')then
         print*,'Cost function is cd/cd0 + C * max(0,1-cl/cl0)'
         icftyp = 4
      else if(string .eq. 'CDCLPEN1')then
         print*,'Cost function is (cd/cd0) + C * abs(1-cl/cl0)'
         icftyp = 5
      else if(string .eq. 'CDCLPEN2')then
         print*,'Cost function is (cd/cd0)**2 + C * (1-cl/cl0)**2'
         icftyp = 6
      else
         print*,'SICFTYP: cost function type not defined'
         stop
      endif

      end
