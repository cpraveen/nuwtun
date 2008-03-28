      subroutine rem_dup(nwp, rw, drw, nsurf, nspt, lidx, idx)
      implicit none
      include 'dim.h'
      integer nwp, nsurf, nspt(*), lidx, idx(lidx,*)
      real    rw(NDIM,*), drw(NDIM,*)

c     Temporary storage
      real    rwt(NDIM,nwp), drwt(NDIM,nwp)

      integer i, j, nwpt, m, n, inew
      integer newid(nwp)
      logical isdup

      rwt (1:NDIM, 1:nwp) = rw (1:NDIM, 1:nwp)
      drwt(1:NDIM, 1:nwp) = drw(1:NDIM, 1:nwp)
      nwpt                = nwp

      nwp = 0
      do i=1,nwpt
         isdup = .false.
         do j=1,i-1
            if( rwt(1,i).eq.rwt(1,j) .and. 
     1          rwt(2,i).eq.rwt(2,j) .and. (.not.isdup)) then
               isdup = .true.
               inew  = j
            endif
         enddo
         if(.not.isdup)then
            nwp        = nwp + 1
            rw (1,nwp) = rwt (1,i)
            rw (2,nwp) = rwt (2,i)
            drw(1,nwp) = drwt(1,i)
            drw(2,nwp) = drwt(2,i)
            newid(i)   = nwp
         else
            newid(i)   = newid(inew)
         endif
      enddo

c     Now put new index location
      do m=1,nsurf
         do n=1,nspt(m)
            i = idx(n,m)
            idx(n,m) = newid(i)
         enddo
      enddo

      if(nwp.ne.nwpt) then
         print*,'Removed duplicate points'
         print*,'   Initial points =',nwpt
         print*,'   Final   points =',nwp
         print*,'   Removed points =',nwpt-nwp
      endif

c     do i=1,nsurf
c        print*,'------------------'
c        print*,'Surface =',i
c        do j=1,nspt(i)
c           print*,j,idx(j,i)
c        enddo
c     enddo
c     stop

      return
      end
