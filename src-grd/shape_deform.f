c param_type = 1, Hicks-Henne
c              2, Kulfan
c              3, Clamped cubic spline
      subroutine shape_deform(ibeg, jbeg, iend, jend, nhhp, idim, jdim,
     1                        r, rw, drw, nwp, xw, idx, param_type)
      implicit none
      include 'dim.h'
      integer ibeg, jbeg, iend, jend, nhhp, idim, jdim, nwp
      real    r   (idim,jdim,NDIM), 
     2        rw  (NDIM,*),
     3        drw (NDIM,*),
     4        xw  (*)
      integer idx(*)
      integer param_type

      integer i, j, nc, iinc, jinc
      real    x1, y1, x2, y2, ln, xp, yp, xf, yf, l, t,
     1        nx, ny, h, dh, A, B, C, D
c     For cubic spline
      real    tau(nhhp+2), cc(4,nhhp+2)
      real    ppvalu
      real    HicksHenne
      real    kulfan

c     Set increment to +1 or -1
      iinc = 1
      jinc = 1
      if(iend.lt.ibeg) iinc = -1
      if(jend.lt.jbeg) jinc = -1

c     (x1,y1) = first point, (x2,y2) = second point on curve
      if(param_type .eq. 1)then
         x1 = r(ibeg,jbeg,1)
         y1 = r(ibeg,jbeg,2)
         x2 = r(iend,jend,1)
         y2 = r(iend,jend,2)
      elseif(param_type .eq. 2)then
         x1 = 1.0
         y1 = 0.0
         x2 = 0.0
         y2 = 0.0
      elseif(param_type .eq. 3)then
         x1 = r(ibeg,jbeg,1)
         y1 = r(ibeg,jbeg,2)
         x2 = r(iend,jend,1)
         y2 = r(iend,jend,2)
         cc(1,1) = 0.0
         tau(1) = 0.0
         do i=1,nhhp
            cc(1,i+1) = xw(i)
            tau(i+1)  = (i*1.0)/(nhhp+1)
         enddo
         cc(1,nhhp+2) = 0.0
         tau(nhhp+2)      = 1.0
c        zero slope at end points: clamped spline
         cc(2,1)      = 0.0
         cc(2,nhhp+2) = 0.0
         call cubspl(tau,cc,nhhp+2,1,1)
      endif
      print*,'   First  point =',x1,y1
      print*,'   Second point =',x2,y2

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
      nx = (y2-y1)/ln
      ny =-(x2-x1)/ln

      print*,'   Unit vector  =',nx,ny

c     Equation of line joining (x1,y1) to (x2,y2) -> Ax + By + C = 0
      A  = (y2 - y1)
      B  =-(x2 - x1)
      C  = -(y2 - y1)*x1 + (x2 - x1)*y1

c     Counter for idx array
      nc = 0

      do i=ibeg,iend,iinc
         do j=jbeg,jend,jinc
            nc        = nc + 1
            nwp       = nwp + 1
            xp        = r(i,j,1)
            yp        = r(i,j,2)
            rw(1,nwp) = xp 
            rw(2,nwp) = yp

            drw(1,nwp)= 0.0
            drw(2,nwp)= 0.0
            idx(nc)   = nwp
            if(nhhp.eq.0) goto 100

c           Equation of perpendicular line through (xp,yp)
c           -B x + A y + D = 0
            D = B * xp - A * yp
c           Solve for intersection point (xf,yf) = foot of perpendicular
            xf = -(A*C - B*D)/(A**2 + B**2)
            yf = -(B*C + A*D)/(A**2 + B**2)
            l  = sqrt( (xf-x1)**2 + (yf-y1)**2 )
            t         = l/ln
            if(param_type .eq. 1)then
               dh = HicksHenne(nhhp, xw, t)
               drw(1,nwp)= dh*nx
               drw(2,nwp)= dh*ny
            elseif(param_type .eq. 2)then
               t = xp
               h = kulfan(nhhp, xw, t)
               drw(1,nwp)= 0.0
               drw(2,nwp)= h - yp
            elseif(param_type .eq. 3)then
               dh = ppvalu(tau,cc,nhhp+1,4,t,0)
               drw(1,nwp)= dh*nx
               drw(2,nwp)= dh*ny
            endif

100         continue
         enddo
      enddo

      return
      end
