c     scale a 3-d plot3d file
c     divides all coordinates by specified scaling factor
      program plot3_scale

      integer i
      integer j
      integer k
      integer m
      integer n
      integer nblks
      real    lscale
 
      integer, allocatable :: ni(:), nj(:), nk(:)
      real, allocatable :: x(:,:,:), y(:,:,:), z(:,:,:), q(:,:,:,:)

      character*64 gfile, ofile

      if(iargc().ne.1)then
         print*,'Specify one argument'
         stop
      endif

      print*,'Scaling factor'
      read*,lscale

      call getarg(1, gfile)
      ofile = trim(gfile)//'.scaled'
      print*,'Writing scaled grid into ', ofile

      open ( unit=7, form='formatted', file=gfile, status='old')
      open ( unit=8, form='formatted', file=ofile)

      read(7,*) nblks
      write(8,*) nblks
      allocate(ni(nblks), nj(nblks), nk(nblks))

      do i=1,nblks
         read(7,*) ni(i), nj(i), nk(i)
         write(8,*) ni(i), nj(i), nk(i)
      enddo

      do m=1,nblks
         allocate( x(ni(m),nj(m),nk(m)) )
         allocate( y(ni(m),nj(m),nk(m)) )
         allocate( z(ni(m),nj(m),nk(m)) )
         read(7,*) 
     &    ((( x(i,j,k), i=1,ni(m)), j=1,nj(m)), k=1,nk(m)),
     &    ((( y(i,j,k), i=1,ni(m)), j=1,nj(m)), k=1,nk(m)),
     &    ((( z(i,j,k), i=1,ni(m)), j=1,nj(m)), k=1,nk(m))

         x = x/lscale
         y = y/lscale
         z = z/lscale

         write(8,'(e24.14)') 
     &    ((( x(i,j,k), i=1,ni(m)), j=1,nj(m)), k=1,nk(m)),
     &    ((( y(i,j,k), i=1,ni(m)), j=1,nj(m)), k=1,nk(m)),
     &    ((( z(i,j,k), i=1,ni(m)), j=1,nj(m)), k=1,nk(m))

         deallocate(x)
         deallocate(y)
         deallocate(z)
      enddo

      close(7)


      stop
      end
