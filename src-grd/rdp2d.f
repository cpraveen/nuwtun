c-----------------------------------------------------------------------------
c     Read one block of 2-d grid
c-----------------------------------------------------------------------------
      subroutine rdp2d(fid, idim, jdim, r)
      implicit none
      include 'dim.h'
      integer fid, idim, jdim
      real    r(idim,jdim,NDIM)

      integer i, j, l

      read(fid,*) (((r(i,j,l),i=1,idim),j=1,jdim),l=1,NDIM)

      return
      end
c-----------------------------------------------------------------------------
