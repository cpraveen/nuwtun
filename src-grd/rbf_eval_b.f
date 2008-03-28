      subroutine rbf_eval_b(idim, jdim, r, rb, nwp, rw, wt, wtb)
      implicit none
      include 'dim.h'
      integer idim, jdim, nwp
      real    r  (idim,jdim,NDIM),
     1        rw (NDIM,*),
     2        wt (NDIM,*)
      real    rb (idim,jdim,NDIM),
     1        wtb(NDIM,*)

      integer i, j, k
      real    dx, dy, dx1, dy1, dr1, rbf
      real    dxb, dyb

      do i=1,idim
         do j=1,jdim
            dxb = rb(i,j,1)
            dyb = rb(i,j,2)
            do k=1,nwp
               dx1 = r(i,j,1) - rw(1,k)
               dy1 = r(i,j,2) - rw(2,k)
               dr1 = sqrt(dx1*dx1 + dy1*dy1)
               if(dr1.eq.0.0)then
                  rbf = 0.0
               else
                  rbf = dr1*dr1*log(dr1)
               endif
               wtb(1,k) = wtb(1,k) + dxb * rbf
               wtb(2,k) = wtb(2,k) + dyb * rbf
            enddo
            wtb(1,nwp+1) = wtb(1,nwp+1) + dxb
            wtb(1,nwp+2) = wtb(1,nwp+2) + dxb*r(i,j,1)
            wtb(1,nwp+3) = wtb(1,nwp+3) + dxb*r(i,j,2)
            wtb(2,nwp+1) = wtb(2,nwp+1) + dyb
            wtb(2,nwp+2) = wtb(2,nwp+2) + dyb*r(i,j,1)
            wtb(2,nwp+3) = wtb(2,nwp+3) + dyb*r(i,j,2)
         enddo
      enddo

      return
      end
