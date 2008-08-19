c-----------------------------------------------------------------------------
c     Read basic grid : 2-d only at present
c-----------------------------------------------------------------------------
      subroutine readgrd(nblks, idim, jdim, ioffr, r)
      implicit none
      include 'dim.h'
      integer nblks, idim(*), jdim(*), ioffr(*)
      real    r(*)

      integer fid, i, ir

      print*,'Reading base grid from file grid.0'

      fid = 10
      open(fid, file='grid.0')
      read(fid,*) nblks
      print*,'Number of blocks =',nblks
      do i=1,nblks
         read(fid,*) idim(i), jdim(i)
         print*,i,idim(i),jdim(i)
      enddo
      do i=1,nblks
         ir = ioffr(i)*NDIM + 1
         call rdp2d( fid, idim(i), jdim(i), r(ir) )
      enddo
      close(fid)

      return
      end
