      subroutine shepard(idim, jdim, r, nwp, rw, drw)
      implicit none
      include 'dim.h'
      integer idim, jdim, nwp
      real    r  (idim,jdim,NDIM),
     1        rw (NDIM,*),
     2        drw(NDIM,*)

      integer i, j, k
      real    dx, dy, dx1, dy1, dr1, den, sf

      do i=1,idim
         do j=1,jdim
            dx = 0.0
            dy = 0.0

            do k=1,nwp
               dx1 = r(i,j,1) - rw(1,k)
               dy1 = r(i,j,2) - rw(2,k)
               dr1 = sqrt(dx1*dx1 + dy1*dy1)
               if(dr1.lt.1.0e-15) then
                  dx = drw(1,k)
                  dy = drw(2,k)
                  goto 100
               endif
            enddo

            den = 0.0

            do k=1,nwp
               dx1 = r(i,j,1) - rw(1,k)
               dy1 = r(i,j,2) - rw(2,k)
               dr1 = sqrt(dx1*dx1 + dy1*dy1)
               sf  = 1.0/dr1/dr1/dr1
               dx  = dx + sf*drw(1,k)
               dy  = dy + sf*drw(2,k)
               den = den + sf
            enddo
            dx = dx/den
            dy = dy/den

100         continue
            r(i,j,1) = r(i,j,1) + dx
            r(i,j,2) = r(i,j,2) + dy
         enddo
      enddo

      return
      end
