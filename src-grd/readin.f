c-----------------------------------------------------------------------------
c     Read input file describing surfaces
c-----------------------------------------------------------------------------
      subroutine readin()
      implicit none
      include 'param.h'

      integer fid, i, j

      fid = 10

      print*,'Reading surfaces from shape.in'
      open(fid, file='shape.in', status='old')
      read(fid,*) nsurf
      print*,'Number of surfaces =',nsurf
      do i=1,nsurf
         read(fid,*) blk(i), ibeg(i), jbeg(i), iend(i), jend(i), nhhp(i)
         write(*,10) blk(i), ibeg(i), jbeg(i), iend(i), jend(i), nhhp(i)
      enddo
      close(fid)

      print*,'Reading Hicks-Henne parameters from hicks.in'
      open(fid, file='hicks.in', status='old')
      do i=1,nsurf
         do j=1,nhhp(i)
            read(fid,*) xw(j,i)
         enddo
      enddo
      close(fid)

10    format(6I6)

      return
      end
