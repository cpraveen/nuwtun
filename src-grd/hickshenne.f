c-----------------------------------------------------------------------------
c     Computes Sum(i=1,n) x(i)*H_i(t) where the H_i are Hicks-Henne
c     functions
c-----------------------------------------------------------------------------
      real function hickshenne(n, x, t)
      implicit none
      integer n
      real    x(*), t

      integer i
      real    h, HicksHenneFunc

      h = 0.0
      do i=1,n
         h = h + x(i)*HicksHenneFunc(n, i, t)
      enddo

      HicksHenne = h

      return
      end

c-----------------------------------------------------------------------------
c     Compute i'th HicksHenne function
c-----------------------------------------------------------------------------
      real function HicksHenneFunc(n, i, t)
      implicit none
      integer n, i
      real    t

      real    xh, tmp1, tmp2, tmp3, tmp4, PI

      PI   = 4.0*atan(1.0)
      xh   = real(i)/(n+4.0)
      tmp1 = log(0.5)/log(xh)
      tmp2 = PI*t**tmp1
      tmp3 = sin(tmp2)
      tmp4 = tmp3**3

      HicksHenneFunc = tmp4

      return
      end
