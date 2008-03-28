      subroutine testfun(idim, jdim, r, fun)
      implicit none
      integer idim, jdim
      real    r(idim,jdim,2), fun

      integer i, j

      do i=1,idim
         do j=1,jdim
c           fun = fun + r(i,j,1)**2 + r(i,j,2)**2
            fun = fun + r(i,j,1) * r(i,j,2)
         enddo
      enddo

      return
      end

      subroutine testgrad(idim, jdim, r, rb, fun)
      implicit none
      integer idim, jdim
      real    r(idim,jdim,2), fun
      real    rb(idim,jdim,2)

      integer i, j

      do i=1,idim
         do j=1,jdim
c           rb(i,j,1) = 2.0*r(i,j,1)
c           rb(i,j,2) = 2.0*r(i,j,2)
            rb(i,j,1) = r(i,j,2)
            rb(i,j,2) = r(i,j,1)
         enddo
      enddo

      return
      end
