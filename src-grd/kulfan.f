c-----------------------------------------------------------------------------
c     Computes Sum(i=0,degree) x(i)*B_i(t) where the B_i are Bernstein
c     polynomials of degree = n - 1
c-----------------------------------------------------------------------------
      real function kulfan(n, x, t)
      implicit none
      integer n
      real    x(*), t

      integer degree
      real    shapefun
      real    decas
      real    yte

      yte = 0.0

      degree   = n - 1
      shapefun = decas(degree, x, t)

      kulfan   = shapefun * sqrt(t) * (1.0 - t) + t * yte

      end
c-----------------------------------------------------------------------------
