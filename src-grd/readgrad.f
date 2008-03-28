c-----------------------------------------------------------------------------
c     Read gradient output by adjoint solver
c-----------------------------------------------------------------------------
      subroutine readgrad(nblks, idim, jdim, ioffr, rb)
      implicit none
      include 'dim.h'
      integer nblks, idim(*), jdim(*), ioffr(*)
      real    rb(*)

      integer fid, i, ir, nb, id, jd

      print*,'Reading gradient file'

      fid = 10
      open(fid, file='rb.dat')
      read(fid,*) nb
      if(nb.ne.nblks)then
         print*,'Number of blocks is incorrect'
         stop
      endif
      do i=1,nblks
         read(fid,*) id, jd
         if(id.ne.idim(i) .or. jd.ne.jdim(i))then
            print*,'Dimensions dont agree for block =',i
            stop
         endif
      enddo
      do i=1,nblks
         ir = ioffr(i)*NDIM + 1
         call rdp2d( fid, idim(i), jdim(i), rb(ir) )
      enddo
      close(fid)

      return
      end
