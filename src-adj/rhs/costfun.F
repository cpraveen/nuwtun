      subroutine costfun(icftyp, cl, cd, cost)
      include '../header/refval.h'
      integer icftyp
      real    cl, cd, cost

      if(icftyp.eq.3)then
         cost = cd/cdref + 1.0e4 * max(1.0 - cl/clref, 0.0)
      else
         print*,'costfun: unknown cost function type'
         stop
      endif

      end
