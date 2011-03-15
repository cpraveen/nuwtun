      program volume
      implicit none
      include 'volume.h'
      real, allocatable :: r(:)
      integer, allocatable :: idim(:), jdim(:), kdim(:), ioffr(:)
      integer nblks

      integer i, j, k, ioffrt, imem, ir, blk
      integer nvol, nsurf(nvolmax), iblk(nvolmax,nsurfmax), 
     1        imin(nvolmax,nsurfmax), jmin(nvolmax,nsurfmax),
     2        kmin(nvolmax,nsurfmax), imax(nvolmax,nsurfmax),
     3        jmax(nvolmax,nsurfmax), kmax(nvolmax,nsurfmax)
      real    totvol
      character grid_format*32, grid_file*32

      if(iargc() .ne. 2)then
         print*,'Usage: volume grid_format grid_file'
         print*,'       grid_format = plot3d or hyena'
         stop
      endif

      call getarg(1, grid_format)
      call getarg(2, grid_file)

c     read volume definition
      write(*,*)'Reading volume.in'
      open(10, file='volume.in', status='old')
      read(10,*) nvol
      if(nvol.gt.nvolmax)then
         print*,'volume: insufficient memory, increase nvolmax'
         stop
      endif
      do i=1,nvol
         read(10,*) nsurf(i)
         if(nsurf(i).gt.nsurfmax)then
            print*,'volume: insufficient memory, increase nsurfmax'
            stop
         endif
         do j=1,nsurf(i)
            read(10,*) iblk(i,j), imin(i,j), jmin(i,j), kmin(i,j),
     1                            imax(i,j), jmax(i,j), kmax(i,j)
         enddo
      enddo
      close(10)

c     read grid
      imem = 0
      nblks = 1 ! limited to single block at present
      write(*,*)'NOTE: WORKS ONLY FOR SINGLE BLOCK GRID'
      allocate( idim(nblks) )
      allocate( jdim(nblks) )
      allocate( kdim(nblks) )
      allocate( ioffr(nblks) )
      write(*,*)'Reading grid file ', grid_file
      open(10, file=grid_file, status='old')
      if(grid_format.eq.'hyena')then
         read(10,*) nblks
         if(nblks.ne.1) stop "Error: cannot read multi-block grid"
      endif
      do i=1,nblks
         read(10,*) idim(i), jdim(i), kdim(i)
         imem = imem + idim(i)*jdim(i)*kdim(i)
      enddo
      imem = 3*imem
      allocate( r(imem) )
c     create offset array
      ioffrt = 0
      do i=1,nblks
         ioffr(i) = ioffrt
         ioffrt   = ioffrt + idim(i)*jdim(i)*kdim(i)
      enddo
c     read each block of grid points
      do i=1,nblks
         ir = ioffr(i)*3 + 1
         if(grid_format.eq.'plot3d')then
            call rdp3d(10, idim(i), jdim(i), kdim(i), r(ir))
         else if(grid_format.eq.'hyena')then
            call rdhyena(10, idim(i), jdim(i), kdim(i), r(ir))
         else
            print*,'Error: Unknown grid format ', grid_format
            stop
         endif
      enddo
      close(10)


c     compute volumes
      totvol = 0.0
      call totalvolume(nvol,nsurf,iblk,ioffr,idim,jdim,kdim,r,
     1                 imin,jmin,kmin,imax,jmax,kmax,totvol)
      print*, 'VOL      ', totvol

      deallocate( idim )
      deallocate( jdim )
      deallocate( kdim )
      deallocate( ioffr)
      deallocate( r    )

      stop
      end

c-----------------------------------------------------------------------------
c     Read one block of 3-d plot3d grid
c-----------------------------------------------------------------------------
      subroutine rdp3d(fid, idim, jdim, kdim, r)
      implicit none
      include 'dim.h'
      integer fid, idim, jdim, kdim
      real    r(idim,jdim,kdim,3)

      integer i, j, k, l

      read(fid,*) ((((r(i,j,k,l),i=1,idim),j=1,jdim),k=1,kdim),l=1,3)

      return
      end
c-----------------------------------------------------------------------------
c     Read one block of 3-d hyena grid
c-----------------------------------------------------------------------------
      subroutine rdhyena(fid, idim, jdim, kdim, r)
      implicit none
      include 'dim.h'
      integer fid, idim, jdim, kdim
      real    r(idim,jdim,kdim,3)

      integer i, j, k, l

      read(fid,*) ((((r(i,j,k,l),l=1,3),i=1,idim),j=1,jdim),k=1,kdim)

      return
      end
c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------
      subroutine totalvolume(nvol,nsurf,iblk,ioffr,idim,jdim,kdim,r,
     1                       imin,jmin,kmin,imax,jmax,kmax,totvol)
      implicit none
      include 'volume.h'
      integer ioffr(*), idim(*), jdim(*), kdim(*)
      real    r(*)
      integer nvol, nsurf(nvolmax), iblk(nvolmax,nsurfmax), 
     1        imin(nvolmax,nsurfmax), jmin(nvolmax,nsurfmax),
     2        kmin(nvolmax,nsurfmax), imax(nvolmax,nsurfmax),
     3        jmax(nvolmax,nsurfmax), kmax(nvolmax,nsurfmax)
      real    totvol

      integer i, j, ir, blk
      real    vol(nvol)

      do i=1,nvol
         vol(i) = 0.0
         do j=1,nsurf(i)
            blk= iblk(i,j)
            ir = ioffr(blk)*3 + 1
            call surfvol(idim(blk),jdim(blk),kdim(blk),r(ir),
     1                   imin(i,j),jmin(i,j),kmin(i,j),
     2                   imax(i,j),jmax(i,j),kmax(i,j),vol(i))
         enddo
         totvol = totvol + vol(i)
      enddo

      return
      end

c-----------------------------------------------------------------------------
c     Computes surface integral
c-----------------------------------------------------------------------------
      subroutine surfvol(idim,jdim,kdim,r,imin,jmin,kmin,
     1                   imax,jmax,kmax,vol)
      implicit none
      integer idim, jdim, kdim, imin, jmin, kmin, imax, jmax, kmax
      real    r(idim,jdim,kdim,3), vol

      integer i, j, k
      real    x(3), y(3), z(3)


      if(imin.eq.imax)then
         i = imin
         do j=jmin,jmax-1
            do k=kmin,kmax-1
               x(1) = r(i,j,  k,  1)
               x(2) = r(i,j+1,k,  1)
               x(3) = r(i,j+1,k+1,1)
               y(1) = r(i,j,  k,  2)
               y(2) = r(i,j+1,k,  2)
               y(3) = r(i,j+1,k+1,2)
               z(1) = r(i,j,  k,  3)
               z(2) = r(i,j+1,k,  3)
               z(3) = r(i,j+1,k+1,3)
               call triang(x,y,z,vol)

               x(1) = r(i,j,  k,  1)
               x(2) = r(i,j+1,k+1,1)
               x(3) = r(i,j,  k+1,1)
               y(1) = r(i,j,  k,  2)
               y(2) = r(i,j+1,k+1,2)
               y(3) = r(i,j,  k+1,2)
               z(1) = r(i,j,  k,  3)
               z(2) = r(i,j+1,k+1,3)
               z(3) = r(i,j,  k+1,3)
               call triang(x,y,z,vol)
            enddo
         enddo
      endif

      if(jmin.eq.jmax)then
         j = jmin
         do i=imin,imax-1
            do k=kmin,kmax-1
               x(1) = r(i,  j,k,  1)
               x(2) = r(i+1,j,k+1,1)
               x(3) = r(i+1,j,k,  1)
               y(1) = r(i,  j,k,  2)
               y(2) = r(i+1,j,k+1,2)
               y(3) = r(i+1,j,k,  2)
               z(1) = r(i,  j,k,  3)
               z(2) = r(i+1,j,k+1,3)
               z(3) = r(i+1,j,k,  3)
               call triang(x,y,z,vol)

               x(1) = r(i,  j,k,  1)
               x(2) = r(i  ,j,k+1,1)
               x(3) = r(i+1,j,k+1,1)
               y(1) = r(i,  j,k,  2)
               y(2) = r(i  ,j,k+1,2)
               y(3) = r(i+1,j,k+1,2)
               z(1) = r(i,  j,k,  3)
               z(2) = r(i  ,j,k+1,3)
               z(3) = r(i+1,j,k+1,3)
               call triang(x,y,z,vol)
            enddo
         enddo
      endif

      if(kmin.eq.kmax)then
         stop "volume: kmin==kmax not finished"
         k = kmin
      endif

      end
c-----------------------------------------------------------------------------
c     Integrate over a 3-d triangle
c-----------------------------------------------------------------------------
      subroutine triang(x, y, z, vol)
      implicit none
      real x(*), y(*), z(*), vol

      real xa, ya, za
      real dx1, dy1, dz1
      real dx2, dy2, dz2
      real nx, ny, nz

      xa = (x(1) + x(2) + x(3))/3.0
      ya = (y(1) + y(2) + y(3))/3.0
      za = (z(1) + z(2) + z(3))/3.0

      dx1= x(2) - x(1)
      dy1= y(2) - y(1)
      dz1= z(2) - z(1)

      dx2= x(3) - x(1)
      dy2= y(3) - y(1)
      dz2= z(3) - z(1)

c     area normal vector
      nx = (dy1*dz2 - dy2*dz1)/2.0
      ny =-(dx1*dz2 - dx2*dz1)/2.0
      nz = (dx1*dy2 - dx2*dy1)/2.0

      vol = vol + (xa*nx + ya*ny + za*nz)/3.0

      end
