      subroutine rbf_train_b(nwp, rw, drw, drwb, wt, wtb)
      implicit none
      include 'dim.h'
      integer nwp
      real    rw  (NDIM,*),
     1        drw (NDIM,*),
     2        wt  (NDIM,*)
      real    drwb(NDIM,*),
     2        wtb (NDIM,*)

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
      do i=1,nwp+3
         b(i) = wtb(1,i)
      enddo
      call dgesl(a, n, n, ipvt, b, 0)
      do i=1,nwp+3
         drwb(1,i) = b(i)
      enddo

c     Solve for y-coordinate
      do i=1,nwp+3
         b(i) = wtb(2,i)
      enddo
      call dgesl(a, n, n, ipvt, b, 0)
      do i=1,nwp+3
         drwb(2,i) = b(i)
      enddo

      return
      end
