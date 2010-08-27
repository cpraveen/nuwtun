c-----------------------------------------------------------------------------
c     Read input file describing surfaces
c-----------------------------------------------------------------------------
      subroutine readin()
      implicit none
      include 'param.h'
      include 'kulfan.h'

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
      read(fid,*) param_type
      if(param_type.eq.2)then ! Some kulfan parameters
         read(fid,*) N1
         read(fid,*) N2
         read(fid,*) yte
         write(*,'(" N1, N2, yte =",2e12.4,e24.14)') N1, N2, yte
      endif
      close(fid)

      if(param_type.eq.1)then
         print*,'Reading HICKS-HENNE parameters from shape.dat'
      elseif(param_type.eq.2)then
         print*,'Reading KULFAN parameters from shape.dat'
      elseif(param_type.eq.3)then
         print*,'Reading CLAMPED CUBIC SPLINE parameters from shape.dat'
      else
         print*,'Unknown parameterization type'
         stop
      endif
      open(fid, file='shape.dat', status='old')

      do i=1,nsurf
         do j=1,nhhp(i)
            read(fid,*) xw(j,i)
         enddo
      enddo
      close(fid)

10    format(6I6)

      return
      end
