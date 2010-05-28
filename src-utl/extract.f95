!! Program to extract PRINT data from nuwtun log
!! Currently extracts only one set of data
program extract
   implicit none
   integer :: fid
   character(LEN=80) :: str
   integer :: IOstatus
   integer :: i, j, k, oid
   real(kind=kind(1.0d0)) :: x, y, z, rho, u, v, w, p, T
   real(kind=kind(1.0d0)) :: mach, gamma, cp

   gamma = 1.4d0

   str = ''
   fid = 20
   open(fid, file='flo.log', status='old')

   do while(str.ne.'MACH')
      read(fid,*) str
   enddo
   backspace(fid)
   read(fid,*) str, mach
   print*,'#Mach number =', mach
   str = ''
   do while(str .ne. 'WRTP3D:')
      read(fid,*) str
   enddo
   if( str == '')then
      print*,'Error: did not find WRTP3D:'
      stop
   endif

   oid = 100

   do
      read(fid,*,IOSTAT=IOstatus) i, j, k, x, y, z, rho, u, v, w, p, T
      if(IOstatus > 0)then
         read(fid,*)
         read(fid,*)
         read(fid,*)
         read(fid,*)
         oid = oid + 1
         write(*,'("Writing file fort.",i3)') oid
      else if(IOstatus < 0)then
         close(fid)
         stop
      else
         cp = 2.0d0*(1.0d0 - p)/mach**2/gamma ! actually -Cp
         write(oid,'(7e20.12)')x,y,z,cp,u,v,w
      endif
   enddo

end
