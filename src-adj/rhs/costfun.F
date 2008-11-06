      subroutine costfun(icftyp, cl, cd, cost)
      include '../header/refval.h'
      integer icftyp
      real    cl, cd, cost

      if(icftyp.eq.1)then
         cost = cl
      else if(icftyp.eq.2)then
         cost = cd
      else if(icftyp.eq.3)then
         cost = cd/cdref + 1.0e4 * max(1.0 - cl/clref, 0.0)
      else
         print*,'costfun: unknown cost function type'
         print*,'         icftyp =', icftyp
         stop
      endif

      end
