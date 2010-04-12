!! Program to generate cp files
!! First use "extract" to separate cp at different sections
!! Then run ext2cp
!! Example:
!!   ext2cp x 1 2 3 4
!!   This assumes x-direction to be allong chord and processes 
!!   fort.1 fort.2 fort.3 fort.4 files to generate fort.1.dat etc.
program extract
   implicit none
   integer :: fid
   character(LEN=80) :: str
   integer :: IOstatus
   integer :: i, j, k, np, narg, nfile
   integer, parameter :: npmax=1000
   real(kind=kind(1.0d0)) :: x(npmax), y(npmax), z(npmax), cp(npmax)
   real(kind=kind(1.0d0)) :: x1, y1, z1, cp1
   real(kind=kind(1.0d0)) :: xmin, xmax, chord
   character(len=56) :: dir, f
   integer :: iargc

   narg = iargc()
   if(narg < 2)then
      print*,'Must specify atleast two arguments'
      stop
   endif

   call getarg(1, dir)
   print*,'Direction = ',TRIM(dir)

   nfile = narg - 1
   print*,'Processing file = ', nfile

   fid = 20

   do i=1,nfile
      print*,'----------------------------------------'
      call getarg(i+1, f)
      print*,'Processing file ', 'fort.'//TRIM(f)
      open(fid, file='fort.'//TRIM(f), status='old')

      np = 0

      do
         read(fid,*,IOSTAT=IOstatus) x1, y1, z1, cp1
         if(IOstatus > 0)then
            print*,'Error reading'
         else if(IOstatus < 0)then
            close(fid)
            goto 123
         else
            np = np + 1
            x(np) = x1
            y(np) = y1
            z(np) = z1
            cp(np) = cp1
         endif
      enddo
123   print*,'Number of points = ', np
      close(fid)
      xmin = +1.0e20
      xmax = -1.0e20
      do j=1,np
        if(dir == 'x')then
           xmin = min(xmin, x(j))
           xmax = max(xmin, x(j))
        elseif(dir == 'y')then
           xmin = min(xmin, y(j))
           xmax = max(xmin, y(j))
        elseif(dir == 'z')then
           xmin = min(xmin, z(j))
           xmax = max(xmin, z(j))
        else
           print*,'Unknown direction = ', dir
        endif
      enddo
      chord = xmax - xmin
      write(*,'(" xmin, xmax, chord =",3e15.5)') xmin, xmax, chord
      open(fid,file='fort.'//TRIM(f)//'.dat')
      do j=1,np
         if(dir == 'x')then
            write(fid,'(2e20.10)') (x(j)-xmin)/chord, cp(j)
         elseif(dir == 'y')then
            write(fid,'(2e20.10)') (y(j)-xmin)/chord, cp(j)
         elseif(dir == 'z')then
            write(fid,'(2e20.10)') (z(j)-xmin)/chord, cp(j)
         endif
      enddo
      close(fid)
   enddo

end
