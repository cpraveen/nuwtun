!-----------------------------------------------------------------------------
!     Read one block of 3-d grid
!-----------------------------------------------------------------------------
      subroutine rdp2d(fid, idim, jdim, r)
      implicit none
      integer fid, idim, jdim
      real    r(idim,jdim,2)

      integer i, j, k, l

      read(fid) (((r(i,j,l),i=1,idim),j=1,jdim),l=1,2)

      return
      end
!-----------------------------------------------------------------------------
!     Write one block of 2-d grid in ascii format
!-----------------------------------------------------------------------------
      subroutine wtp2d(fid, idim, jdim, r)
      implicit none
      integer fid, idim, jdim
      real    r(idim,jdim,2)

      integer i, j, k, l

      write(fid,*) (((r(i,j,l),i=1,idim),j=1,jdim),l=1,2)

      return
      end
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
      subroutine totalarea(narea,nsurf,iblk,ioffr,idim,jdim,r, &
                             imin,jmin,imax,jmax,totarea)
      implicit none
      include 'area.h'
      integer ioffr(*), idim(*), jdim(*)
      real    r(*)
      integer narea, nsurf(narmax), iblk(narmax,nsurfmax), &
              imin(narmax,nsurfmax), jmin(narmax,nsurfmax),&
              imax(narmax,nsurfmax), jmax(narmax,nsurfmax)
      real    totarea

      integer i, j, ir, blk
      real    area

!$AD II-LOOP
      do i=1,narea
         area = 0.0
!$AD II-LOOP
         do j=1,nsurf(i)
            blk= iblk(i,j)
            ir = ioffr(blk)*2 + 1
            call surfarea(idim(blk),jdim(blk),r(ir), &
                         imin(i,j),jmin(i,j), &
                         imax(i,j),jmax(i,j),area)
         enddo
         totarea = totarea + area
      enddo

      return
      end
!-----------------------------------------------------------------------------
!     Computes surface integral
!-----------------------------------------------------------------------------
      subroutine surfarea(idim,jdim,r,imin,jmin,imax,jmax,area)
      implicit none
      integer idim, jdim, imin, jmin, imax, jmax
      real    r(idim,jdim,2), area

      integer i, j, k
      real    dx, dy, xc, yc

      if(imin.eq.imax)then

         i = imin
!$AD II-LOOP
         do j=jmin,jmax-1
            dx = r(i,j+1,1) - r(i,j,1)
            dy = r(i,j+1,2) - r(i,j,2)
            xc = 0.5*(r(i,j+1,1) + r(i,j,1))
            yc = 0.5*(r(i,j+1,2) + r(i,j,2))
            area = area + 0.5*(-xc*dy + yc*dx)
         enddo

      elseif(jmin.eq.jmax)then

         j = jmin
!$AD II-LOOP
         do i=imin,imax-1
            dx = r(i+1,j,1) - r(i,j,1)
            dy = r(i+1,j,2) - r(i,j,2)
            xc = 0.5*(r(i+1,j,1) + r(i,j,1))
            yc = 0.5*(r(i+1,j,2) + r(i,j,2))
            area = area + 0.5*(-xc*dy + yc*dx)
         enddo

      endif

      end
