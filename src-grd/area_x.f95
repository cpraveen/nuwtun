!     Computes area enclosed by curves. These curves are defined from
!     the input file area.in by specifying the grid indices
      program area_x
      USE subroutines
      implicit none
      include 'area.h'
      real, pointer :: r(:) => NULL()
      real, pointer :: rb(:) => NULL()
      integer idim(nblkmax)
      integer jdim(nblkmax)
      integer ioffr(nblkmax)
      integer nblks

      integer i, j, k, ioffrt, imem, ir, blk, ierr
      integer narea, nsurf(narmax), iblk(narmax,nsurfmax),  &
              imin(narmax,nsurfmax), jmin(narmax,nsurfmax), &
              imax(narmax,nsurfmax), jmax(narmax,nsurfmax)
      real    totarea, totareab

      call common_area(nblks,narea,nsurf,iblk,imin,jmin,imax,jmax, &
                       idim,jdim,ioffr,r)

      allocate( rb(size(r)) )
      rb = 0.0

!     compute volumes
      totarea = 0.0
      totareab= 1.0
      call totalarea_b(narea,nsurf,iblk,ioffr,idim,jdim,r,rb, &
                       imin,jmin,imax,jmax,totarea,totareab)

      print*,'Writing rb.dat ...'
      open(20, file='rb.dat') 
      write(20,*) nblks
      do i=1,nblks
         write(20,*) idim(i), jdim(i)
      enddo
      do i=1,nblks
         ir = ioffr(i)*2 + 1
         call wtp2d(20,idim(i),jdim(i),rb(ir))
      enddo
      close(20)

      deallocate( r )
      deallocate( rb)

      stop
      end
