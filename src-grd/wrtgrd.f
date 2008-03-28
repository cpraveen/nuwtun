      subroutine wrtgrd(ifid, idim, jdim, r, dr)
      implicit none
      include 'dim.h'
      integer ifid, idim, jdim
      real    r (idim,jdim,NDIM)
      real   dr (idim,jdim,NDIM)

      integer i, j, l

      do i=1,idim
         do j=1,jdim
            write(50,*) r(i,j,1), dr(i,j,2)
            write(51,*) r(i,j,1), r(i,j,2)
         enddo
         write(50,*)
         write(51,*)
      enddo

         do j=1,jdim
      do i=1,idim
            write(50,*) r(i,j,1), dr(i,j,2)
            write(51,*) r(i,j,1), r(i,j,2)
         enddo
         write(50,*)
         write(51,*)
      enddo

c     Write deformed grid in unformatted style
      write(ifid)(((r(i,j,l),i=1,idim),j=1,jdim),l=1,NDIM)

      return
      end


