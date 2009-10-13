!     Computes area enclosed by curves. These curves are defined from
!     the input file area.in by specifying the grid indices
MODULE subroutines
      CONTAINS

      subroutine common_area(nblks,narea,nsurf,iblk,imin,jmin,imax,jmax, &
                             idim,jdim,ioffr,r)
      implicit none
      include 'area.h'
      integer nblks
      real,pointer :: r(:)
      integer idim(*), jdim(*), ioffr(*)
      integer narea, nsurf(narmax), iblk(narmax,nsurfmax),  &
              imin(narmax,nsurfmax), jmin(narmax,nsurfmax), &
              imax(narmax,nsurfmax), jmax(narmax,nsurfmax)
     
      integer i, j, k, ioffrt, imem, ir, blk, ierr
      
!     read area definition
      ierr = 0
      write(*,*)'Reading area.in'
      open(10, file='area.in', status='old')
      read(10,*) narea
      if(narea.gt.narmax)then
         print*,'area: insufficient memory, increase narmax'
         stop
      endif
      do i=1,narea
         read(10,*) nsurf(i)
         if(nsurf(i).gt.nsurfmax)then
            print*,'area: insufficient memory, increase nsurfmax'
            stop
         endif
         do j=1,nsurf(i)
            read(10,*) iblk(i,j), imin(i,j), jmin(i,j), &
                                  imax(i,j), jmax(i,j)
            write(*,*) iblk(i,j), imin(i,j), jmin(i,j), &
                                  imax(i,j), jmax(i,j)
            if(imin(i,j).eq.imax(i,j).and.jmin(i,j).eq.jmax(i,j))then
               print*,'area.in: imin=imax and jmin=jmax for surf=',j
               ierr = 1
            endif
         enddo
      enddo
      close(10)
      if(ierr.ne.0)then
         print*,'Errors detected in area.in'
         stop
      endif

!     read grid
      imem = 0
      nblks = 1 ! limited to single block at present
      write(*,*)'Reading grid.unf'
      open(10, file='grid.unf', status='old', form='unformatted')
      do i=1,nblks
         read(10) idim(i), jdim(i)
         imem = imem + idim(i)*jdim(i)
      enddo
      imem = 2*imem
      allocate( r(imem) )
!     create offset array
      ioffrt = 0
      do i=1,nblks
         ioffr(i) = ioffrt
         ioffrt   = ioffrt + idim(i)*jdim(i)
      enddo
!     read each block of grid points
      do i=1,nblks
         ir = ioffr(i)*2 + 1
         call rdp2d(10, idim(i), jdim(i), r(ir))
      enddo
      close(10)

      return
      end subroutine common_area

END MODULE subroutines
