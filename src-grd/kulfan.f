c-----------------------------------------------------------------------------
c     Computes Sum(i=0,degree) x(i)*B_i(t) where the B_i are Bernstein
c     polynomials of degree = n - 1
c-----------------------------------------------------------------------------
      real function kulfan(n, x, t)
      implicit none
      include 'kulfan.h'
      integer n
      real    x(*), t

      integer degree
      real    shapefun
      real    decas

      degree   = n - 1
      shapefun = decas(degree, x, t)

      kulfan   = shapefun * (t)**N1 * (1.0 - t)**N2 + t * yte

      end
c-----------------------------------------------------------------------------
