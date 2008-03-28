      subroutine rbf_eval(idim, jdim, r, nwp, rw, wt)
      implicit none
      include 'dim.h'
      integer idim, jdim, nwp
      real    r  (idim,jdim,NDIM),
     1        rw (NDIM,*),
     2        wt (NDIM,*)

      integer i, j, k
      real    dx, dy, dx1, dy1, dr1, rbf

      do i=1,idim
         do j=1,jdim
            dx = 0.0
            dy = 0.0
            do k=1,nwp
               dx1 = r(i,j,1) - rw(1,k)
               dy1 = r(i,j,2) - rw(2,k)
               dr1 = sqrt(dx1*dx1 + dy1*dy1)
               if(dr1.eq.0.0)then
                  rbf = 0.0
               else
                  rbf = dr1*dr1*log(dr1)
               endif
               dx  = dx + wt(1,k)*rbf
               dy  = dy + wt(2,k)*rbf
            enddo
            dx       = dx + wt(1,nwp+1) + wt(1,nwp+2)*r(i,j,1) +
     1                 wt(1,nwp+3)*r(i,j,2)
            dy       = dy + wt(2,nwp+1) + wt(2,nwp+2)*r(i,j,1) +
     1                 wt(2,nwp+3)*r(i,j,2)
            r(i,j,1) = r(i,j,1) + dx
            r(i,j,2) = r(i,j,2) + dy
         enddo
      enddo

      return
      end
