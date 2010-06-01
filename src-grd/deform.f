c-----------------------------------------------------------------------------
c     Hicks-Henne based shape deformation parameterization
c     Deformation is interpolated to interior nodes using Thin-Plate Spline RBF 
c     interpolation
c
c     Input files:
c       shape.in  : Defines fixed and moving surfaces
c       hicks.in  : Hicks-Henne parameters
c       grid.0    : Base grid
c
c     Output files:
c-----------------------------------------------------------------------------
      program deform
      implicit none
      include 'dim.h'
      include 'param.h'
      include 'kulfan.h'
      integer total_pts
      real,allocatable:: r(:), dr(:), rw(:,:), drw(:,:), wt(:,:)
      integer,allocatable:: idx(:,:)

      integer i, j, ioffrt, iblk, npts, nwp, ir, di, dj, mn, ifid
      real    fun

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

c     Construct RBF interpolant
      call rbf_train(nwp, rw, drw, wt)
c     do i=1,nwp
c        print*,i,drw(1,i),drw(2,i)
c     enddo
c     stop
c     do i=1,nwp+3
c        print*,i,wt(1,i),wt(2,i)
c     enddo

c     Copy initial grid into dr, only for visualization
      dr(:) = r(:)

c     Compute deformation for interior points
      print*,'Applying deformation to grid'
      do iblk=1,nblks
         ir   = ioffr(iblk)*NDIM + 1
         print*,'   Block =',iblk
         call rbf_eval(idim(iblk), jdim(iblk), r(ir), nwp, rw, wt)
c        call shepard(idim(iblk), jdim(iblk), r(ir), nwp, rw, drw)
      enddo

c     write out new grid
      ifid = 20
      open(ifid,file='grid.unf',form='unformatted')
      if(nblks.gt.1) write(ifid) nblks
      write(ifid)(idim(i),jdim(i),i=1,nblks)
      do i=1,nblks
         ir   = ioffr(i)*NDIM + 1
         call wrtgrd(ifid, idim(i), jdim(i), r(ir), dr(ir))
      enddo
      close(ifid)
      print*,'Deformed grid written into grid.unf'

c     Test function
c     fun = 0
c     do i=1,nblks
c        ir   = ioffr(iblk)*NDIM + 1
c        call testfun(idim(i), jdim(i), r(ir), fun)
c     enddo
c     print*,'Function =',fun

      stop
      end
