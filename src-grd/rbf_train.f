      subroutine rbf_train(nwp, rw, drw, wt)
      implicit none
      include 'dim.h'
      integer nwp
      real    rw  (NDIM,*),
     1        drw (NDIM,*),
     2        wt  (NDIM,*)

      integer i, j, n, ipvt(nwp+NDIM+1)
      real    dx, dy, ds, z(nwp+NDIM+1), rcond
      real    a(nwp+NDIM+1, nwp+NDIM+1), b(nwp+NDIM+1)

      print*,'Constructing RBF approximation'

      n      = nwp+NDIM+1
      a(:,:) = 0.0

c     construct rbf coefficient matrix
      do i=1,nwp
         a(i,i) = 0.0
         do j=i+1,nwp
            dx = rw(1,i) - rw(1,j)
            dy = rw(2,i) - rw(2,j)
            ds = sqrt(dx*dx + dy*dy)
            a(i,j) = ds*ds*log(ds)
         enddo
         a(i,nwp+1) = 1.0
         a(i,nwp+2) = rw(1,i)
         a(i,nwp+3) = rw(2,i)
      enddo

c     Copy lower diagonal part due to symmetry
      do i=1,nwp
         do j=1,i-1
            a(i,j) = a(j,i)
         enddo
      enddo

c     Extra equations due to linear polynomial
      do i=1,nwp
         a(nwp+1,i) = 1.0
         a(nwp+2,i) = rw(1,i)
         a(nwp+3,i) = rw(2,i)
      enddo

c     Perform LU Decomposition
      call dgeco(a, n, n, ipvt, rcond, z)

      write(*,'("    Condition number = ",e10.4)') 1.0/rcond

c     Solve for x-coordinate
      do i=1,nwp
         b(i) = drw(1,i)
      enddo
      b(nwp+1) = 0.0
      b(nwp+2) = 0.0
      b(nwp+3) = 0.0
      call dgesl(a, n, n, ipvt, b, 0)
      do i=1,nwp+3
         wt(1,i) = b(i)
      enddo

c     Solve for y-coordinate
      do i=1,nwp
         b(i) = drw(2,i)
      enddo
      b(nwp+1) = 0.0
      b(nwp+2) = 0.0
      b(nwp+3) = 0.0
      call dgesl(a, n, n, ipvt, b, 0)
      do i=1,nwp+3
         wt(2,i) = b(i)
      enddo

      return
      end
