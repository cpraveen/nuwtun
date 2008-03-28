c-----------------------------------------------------------------------------
c     Computes Sum(i=1,n) x(i)*H_i(t) where the H_i are Hicks-Henne
c     functions
c-----------------------------------------------------------------------------
      subroutine hickshenne_b(n, x, xb, hb, t)
      implicit none
      integer n
      real    x(*), xb(*), hb, t

      integer i
      real    h, HicksHenneFunc

      do i=1,n
         xb(i) = xb(i) + hb*HicksHenneFunc(n, i, t)
      enddo

      return
      end

