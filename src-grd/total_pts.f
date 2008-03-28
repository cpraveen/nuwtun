c-----------------------------------------------------------------------------
c     Find total grid points in all blocks for allocating memory
c-----------------------------------------------------------------------------
      integer function total_pts(nblks, idim, jdim)
      implicit none
      integer nblks, idim(*), jdim(*)

      integer id, jd, fid, i

      total_pts = 0

      fid = 10
      open(fid, file='grid.0')
      read(fid,*) nblks
      do i=1,nblks
         read(fid,*) id, jd
         total_pts = total_pts + id*jd
         idim(i)   = id
         jdim(i)   = jd
      enddo
      close(fid)

      return
      end
