C
C     This file is part of NUWTUN, see <http://nuwtun.berlios.de>
C
      subroutine costfun(icftyp, cl, cd, apgint, cost)
C
C     Routine to compute cost function
C     icftyp : Cost function type, set in sicftyp.F
C     cl     : Lift coefficient
C     cd     : Drag coefficient
C     cost   : Cost function value
C
      include '../header/refval.h'
      integer icftyp
      real    cl, cd, apgint, cost

      if(icftyp.eq.1)then
         cost = cl/clref
      else if(icftyp.eq.2)then
         cost = cd/cdref
      else if(icftyp.eq.3)then
         cost = (cd/cdref)/(cl/clref)
      else if(icftyp.eq.4)then
         cost = cd/cdref + 1.0 * max(1.0 - cl/clref, 0.0)
      else if(icftyp.eq.5)then
         cost = cd/cdref + 1.0 * abs(1.0 - cl/clref)
      else if(icftyp.eq.6)then
         cost = (cd/cdref)**2 + 20.0 * (1.0 - cl/clref)**2
      else if(icftyp.eq.101)then
         cost = apgint/apref
      else if(icftyp.eq.102)then
         cost = apgint/apref + 1.0 * abs(1.0 - cl/clref)
      else
         print*,'costfun: unknown cost function type'
         print*,'         icftyp =', icftyp
         stop
      endif

      end
