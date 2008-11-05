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
      else if(string .eq. 'CDCLPEN1')then
         print*,'Cost function is cd with cl penalty'
         icftyp = 3
      else
         print*,'SICFTYP: cost function type not defined'
         stop
      endif

      end
