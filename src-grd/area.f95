!     Computes area enclosed by curves. These curves are defined from
!     the input file area.in by specifying the grid indices
      program area
      USE subroutines
      implicit none
      include 'area.h'
      real, pointer :: r(:) => NULL()
      integer idim(nblkmax)
      integer jdim(nblkmax)
      integer ioffr(nblkmax)
      integer nblks

      integer i, j, k, ioffrt, imem, ir, blk, ierr
      integer narea, nsurf(narmax), iblk(narmax,nsurfmax),  &
              imin(narmax,nsurfmax), jmin(narmax,nsurfmax), &
              imax(narmax,nsurfmax), jmax(narmax,nsurfmax)
      real    totarea

      call common_area(nblks,narea,nsurf,iblk,imin,jmin,imax,jmax, &
                       idim,jdim,ioffr,r)

!     compute volumes
      totarea = 0.0
      call totalarea(narea,nsurf,iblk,ioffr,idim,jdim,r, &
                       imin,jmin,imax,jmax,totarea)
      print*, 'AREA     ', totarea

      open(20, file='area.out')
      write(20,'(e24.15)') totarea
      close(20)

      deallocate( r    )

      stop
      end
