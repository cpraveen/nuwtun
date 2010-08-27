      subroutine shape_deform_b(ibeg, jbeg, iend, jend, nhhp, 
     1                          idim, jdim, r, rw, drw, drwb,
     2                          nwp, xw, xwb, idx, param_type)
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
      integer param_type

      integer i, j, nc, iw, iinc, jinc
      real    x1, y1, x2, y2, ln, xp, yp, xf, yf, l, t,
     1        nx, ny, dh, dhb, A, B, C, D, h
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
         print*,'shape_deform_b.f:  param_type=3 is not finished'
         stop
      endif
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

c     Equation of line joining (x1,y1) to (x2,y2) -> Ax + By + C = 0
      A  = (y2 - y1)
      B  =-(x2 - x1)
      C  = -(y2 - y1)*x1 + (x2 - x1)*y1

      nc          = 0
      xwb(1:nhhp) = 0.0

      do i=ibeg,iend,iinc
         do j=jbeg,jend,jinc
            xp        = r(i,j,1)
            yp        = r(i,j,2)

            if(nhhp.eq.0) goto 100

c           Equation of perpendicular line through (xp,yp)
c           -B x + A y + D = 0
            D = B * xp - A * yp
c           Solve for intersection point (xf,yf) = foot of perpendicular
            xf = -(A*C - B*D)/(A**2 + B**2)
            yf = -(B*C + A*D)/(A**2 + B**2)
            l  = sqrt( (xf-x1)**2 + (yf-y1)**2 )
            t         = l/ln
            nc        = nc + 1
            iw        = idx(nc)
            print*,nc,iw
            if(param_type .eq. 1)then
               dh        = HicksHenne(nhhp, xw, t)
c              Reverse sweep
               dhb       = drwb(1,iw)*nx + drwb(2,iw)*ny
               call HicksHenne_b(nhhp, xw, xwb, dhb, t)
            elseif(param_type .eq. 2)then
               t = xp
               h = kulfan(nhhp, xw, t)
c              Reverse sweep
               dhb = drwb(2,iw)
               call kulfan_b(nhhp, xw, xwb, t, dhb)
            endif

100         continue
         enddo
      enddo

      return
      end
