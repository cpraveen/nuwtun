c-----------------------------------------------------------------------------
c     ADJOINT of Hicks-Henne based shape deformation parameterization
c     Deformation is interpolated to interior nodes using Thin-Plate Spline RBF 
c     interpolation
c
c     Input files:
c       shape.in  : Defines fixed and moving surfaces
c       hicks.in  : Hicks-Henne parameters
c       grid.0    : Base grid
c       rb.dat    : Gradient wrt grid coordinates, given by flow adjoint solver
c
c     Output files:
c-----------------------------------------------------------------------------
      program deform_adj
      implicit none
      include 'dim.h'
      include 'param.h'
      include 'kulfan.h'
      integer total_pts
      real,allocatable:: r(:), dr(:), rw(:,:), drw(:,:), wt(:,:)
      real,allocatable:: rb(:), drwb(:,:), wtb(:,:)
      integer,allocatable:: idx(:,:)

      integer i, j, ioffrt, iblk, npts, nwp, ir, di, dj, mn
      integer fid

c     Read input file defining moving and fixed surfaces
      call readin()

c     Total number of grid points to allocate memory
      npts = total_pts(nblks, idim, jdim)
      if(nblks.gt.nblkmax)then
         print*,'main: increase maxblks'
         print*,'      nblks   =',nblks
         print*,'      maxblks =',nblkmax
         stop
      endif
      allocate ( r (npts*NDIM) )
      allocate (rb (npts*NDIM) )
      allocate (dr (npts*NDIM) )

c     Create offset array for coordinates storage
      ioffrt = 0
      do i=1,nblks
         ioffr(i) = ioffrt
         ioffrt   = ioffrt + idim(i) * jdim(i)
      enddo
      call readgrd(nblks, idim, jdim, ioffr, r)

c     Count number of boundary points to allocate memory
      nwp = 0
      mn  = 0
      do i=1,nsurf
         di  = abs(iend(i) - ibeg(i)) + 1
         dj  = abs(jend(i) - jbeg(i)) + 1
         nwp = nwp + di*dj
         mn  = max( mn, di*dj )
         nspt(i) = di*dj
      enddo
      allocate( rw  (NDIM, nwp) )
      allocate( drw (NDIM, nwp) )
      allocate( wt  (NDIM, nwp+NDIM+1) )
      allocate( wtb (NDIM, nwp+NDIM+1) )
      allocate( drwb(NDIM, nwp) )
      allocate( idx (mn, nsurf) )

c     For mirror surfaces like in symmetric airfoil
      do i=1,nsurf
         nhhpo(i) = nhhp(i)
         if(nhhp(i).lt.0)then
            ir      =-nhhp(i)
            nhhp(i) = nhhp(ir)
            do j=1,nhhp(i)
               xw(j,i) = -xw(j,ir)
            enddo
         endif
      enddo

c     Loop over surfaces
c     Put boundary points in rw and deformation into drw
c     Count number nwp of boundary points: there may be duplicate points
      nwp = 0
      do i=1,nsurf
         print*,'Applying deformation to surface =',i
         iblk = blk(i)
         ir   = ioffr(iblk)*NDIM + 1
         call shape_deform(ibeg(i), jbeg(i), iend(i), jend(i), nhhp(i),
     1                     idim(iblk), jdim(iblk), r(ir),
     2                     rw, drw, nwp, xw(1,i), idx(1,i), param_type)
      enddo
      print*,'Number of wall points (inc. duplicates) =',nwp

c     Remove duplicate boundary points
      call rem_dup(nwp, rw, drw, nsurf, nspt, mn, idx)

c-------------------
c     Reverse sweep
c-------------------

      wtb(:,:)  = 0.0
      drwb(:,:) = 0.0

c     Read gradient file, output by adjoint solver
      call readgrad(nblks, idim, jdim, ioffr, rb)

c     Test gradient
c     do iblk=1,nblks
c        ir   = ioffr(iblk)*NDIM + 1
c        call testgrad(idim(iblk), jdim(iblk), r(ir), rb(ir))
c     enddo

c     Compute deformation for interior points
      do iblk=1,nblks
         ir   = ioffr(iblk)*NDIM + 1
         call rbf_eval_b(idim(iblk), jdim(iblk), r(ir), rb(ir),
     1                   nwp, rw, wt, wtb)
      enddo

c     Construct RBF interpolant
      call rbf_train_b(nwp, rw, drw, drwb, wt, wtb)

      do i=1,nwp
      print*,i,drwb(1,i),drwb(2,i)
      enddo

      do i=1,nsurf
         print*,'Applying deformation to surface =',i
         iblk = blk(i)
         ir   = ioffr(iblk)*NDIM + 1
         call shape_deform_b(ibeg(i), jbeg(i), iend(i), jend(i),
     1                     nhhp(i), idim(iblk), jdim(iblk), r(ir),
     2                     rw, drw, drwb, nwp, xw(1,i), xwb(1,i),
     3                     idx(1,i), param_type)
         do j=1,nhhp(i)
            print*,j,xwb(j,i)
         enddo
      enddo

c     Differentiation for symmetric case
      do i=1,nsurf
         if(nhhpo(i).lt.0)then
            ir      =-nhhpo(i)
            do j=1,nhhp(i)
               xwb(j,ir) = xwb(j,ir) - xwb(j,i)
               xwb(j,i)  = 0.0
            enddo
         endif
      enddo

c     Save gradients
      fid = 15
      mn  = 0 ! Counter for parameter number
      open(fid, file='gradient.dat')
      do i=1,nsurf
         do j=1,nhhpo(i)
            mn = mn + 1
            write(fid,*) mn, xwb(j,i)
         enddo
      enddo
      close(fid)

      stop
      end
