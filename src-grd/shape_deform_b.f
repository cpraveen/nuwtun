      subroutine shape_deform_b(ibeg, jbeg, iend, jend, nhhp, 
     1                          idim, jdim, r, rw, drw, drwb,
     2                          nwp, xw, xwb, idx)
      implicit none
      include 'dim.h'
      integer ibeg, jbeg, iend, jend, nhhp, idim, jdim, nwp
      real    r   (idim,jdim,NDIM), 
     2        rw  (NDIM,*),
     3        drw (NDIM,*),
     4        xw  (*)
      real    drwb(NDIM,*),
     4        xwb (*)
      integer idx(*)

      integer i, j, nc, iw
      real    x1, y1, x2, y2, ln, xp, yp, xf, yf, l, t, dx, dy, ds, ex,
     1        ey, nx, ny, l1p, h, dh, dhb, A, B, C, HicksHenne

c     (x1,y1) = first point, (x2,y2) = second point on curve
      x1 = r(ibeg,jbeg,1)
      y1 = r(ibeg,jbeg,2)
      x2 = r(iend,jend,1)
      y2 = r(iend,jend,2)
      ln = sqrt( (x2-x1)**2 + (y2-y1)**2 )
      if(ln.eq.0.0 .and. nhhp.ne.0)then
         print*,'shape_deform: ln is <= 0'
         stop
      else if(ln.eq.0.0 .and. nhhp.eq.0)then
         ln = 1.0
      else if(ln.eq.0.0)then
         print*,'shape_deform: OOPS'
         print*,'ln =',ln
         print*,'nhhp =',nhhp
         stop
      endif
      nx = (x2-x1)/ln
      ny = (y2-y1)/ln

c     Equation of line joining (x1,y1) to (x2,y2) -> Ax + By + C = 0
      A  = (y2 - y1)
      B  =-(x2 - x1)
      C  = -(y2 - y1)*x1 + (x2 - x1)*y1

      nc          = 0
      xwb(1:nhhp) = 0.0

      do i=ibeg,iend
         do j=jbeg,jend
            xp        = r(i,j,1)
            yp        = r(i,j,2)

            if(nhhp.eq.0) goto 100

c           Find foot of perpendicular (xf,yf)
            h         = abs( A*xp + B*yp + C)/sqrt(A**2 + B**2)
            l1p       = sqrt( (x1-xp)**2 + (y1-yp)**2 )
            l         = l1p**2 - h**2
c           If l is close to zero it may be negative: clip it to zero
            if(abs(l).lt.1.0e-15) l = 0.0
            if(l.lt.0.0)then
               print*,'shape_deform: FATAL !!!'
               print*,'              l < 0'
               print*,'              l =',l
               stop
            endif
            l         = sqrt(l)
            t         = l/ln
            dh        = HicksHenne(nhhp, xw, t)

            xf        = x1 + nx*l
            yf        = y1 + ny*l
c           Unit vector from (xf,yf) to (xp,yp)
            dx        = xp - xf
            dy        = yp - yf
            ds        = sqrt(dx**2 + dy**2)
            if(ds.eq.0.0)then
               ex     = 0.0
               ey     = 0.0
            else
               ex     = dx/ds
               ey     = dy/ds
            endif

c           Reverse sweep

            nc        = nc + 1
            iw        = idx(nc)
            print*,nc,iw
            dhb       = drwb(1,iw)*ex + drwb(2,iw)*ey
            call HicksHenne_b(nhhp, xw, xwb, dhb, t)

100         continue
         enddo
      enddo

      return
      end
