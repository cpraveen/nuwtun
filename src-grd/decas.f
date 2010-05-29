c 
c    uses  de Casteljau to compute one coordinate
c    value of a  Bezier curve. Has to be called 
c    for each coordinate  (x,y, and/or z) of a control polygon.
c    Input:   degree: degree of curve.
c             coeff:  array with coefficients of curve.
c             t:      parameter value.
c     Output: coordinate value.
      real function decas(degree, coeff, t)
      implicit none
      integer degree
      real    coeff(*), t

      integer r, i
      real    t1
c     an auxiliary array. Change dim. if too small
      real    coeffa(0:100)

      t1 = 1.0 - t;
      coeffa(0:degree) = coeff(1:degree+1)

      do r=1,degree
         do i=0,degree-r
            coeffa(i)= t1* coeffa(i)  +   t * coeffa(i+1)
         enddo
      enddo 

      decas = coeffa(0)

      end
