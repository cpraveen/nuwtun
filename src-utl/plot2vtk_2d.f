C     Read a single block 2-d formatted plot3d file and write a vtk 
C     structured grid file. Reads g-file and q-file.
C     Check gfile, qfile for the input files
      program plot2vtk

      integer i
      integer j
      integer m
      integer n
      integer ni
      integer nj
      integer nblks
 
      real mach   ! freestream Mach number
      real alpha  ! freestream angle-of-attack
      real reyn   ! freestream Reynolds number
      real time   ! time

      real, allocatable :: x(:,:), y(:,:), q(:,:,:)

      character*64 basename
      character*64 gfile, qfile, ofile

      integer ivtk
      real    gamma, d, du, dv, e, p, ma
      parameter(gamma=1.4)

      if(iargc().ne.1)then
         print*,'Specify one argument'
         stop
      endif

      call getarg(1, basename)
      gfile = trim(basename)//'.g.fmt'
      qfile = trim(basename)//'.q.fmt'
      ofile = trim(basename)//'.vtk'

      open ( unit=7, form='formatted', file=gfile, status='old')
      open ( unit=8, form='formatted', file=qfile, status='old' )

      read(7,*) nblks
      if(nblks.ne.1)then
         print*,'Cannot read multi-block file'
         print*,'Number of blocks =', nblks
         stop
      endif
      read(7,*) ni, nj
      allocate( x(ni,nj) )
      allocate( y(ni,nj) )
      allocate( q(ni,nj,4) )
      read(7,*) 
     &    (( x(i,j), i=1,ni), j=1,nj),
     &    (( y(i,j), i=1,ni), j=1,nj)

      read(8,*) nblks
      if(nblks.ne.1)then
         print*,'Cannot read multi-block file'
         print*,'Number of blocks =', nblks
         stop
      endif
      read(8,*) ni, nj
      read(8,*) mach, alpha, reyn, time
      read(8,*) ((( q(i,j,n), i=1,ni), j=1,nj), n=1,4)

      ivtk=10
      write(*,*)'Writing vtk file ', ofile
      open(unit=ivtk, file=ofile)
      write(ivtk,'("# vtk DataFile Version 2.0")')
      write(ivtk,'("NACA0012")')
      write(ivtk,'("ASCII")')
      write(ivtk,'("DATASET STRUCTURED_GRID")')
      write(ivtk,'("DIMENSIONS ",i10,i10,i10)') ni, nj, 1
      write(ivtk,'("POINTS ",i10," float")') ni*nj
      do j=1,nj
         do i=1,ni
            write(ivtk,*) x(i,j), y(i,j), 0.0
         enddo
      enddo
      write(ivtk,'("POINT_DATA",i10)') ni*nj
c     Write density
      write(ivtk,'("SCALARS density float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do j=1,nj
         do i=1,ni
            write(ivtk,*) q(i,j,1)
         enddo
      enddo
c     Write pressure
      write(ivtk,'("SCALARS pressure float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do j=1,nj
         do i=1,ni
            d = q(i,j,1)
            du= q(i,j,2)
            dv= q(i,j,3)
            e = q(i,j,4)
            q2= (du**2 + dv**2)/d/d
            p = (gamma-1.0)*(e - 0.5*d*q2)
            write(ivtk,*) p
         enddo
      enddo
c     Write mach
      write(ivtk,'("SCALARS mach float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do j=1,nj
         do i=1,ni
            d = q(i,j,1)
            du= q(i,j,2)
            dv= q(i,j,3)
            e = q(i,j,4)
            q2= (du**2 + dv**2)/d/d
            p = (gamma-1.0)*(e - 0.5*d*q2)
            ma= sqrt(q2)/sqrt(gamma*p/d)
            write(ivtk,*) ma
         enddo
      enddo
c     Write velocity vector
      write(ivtk,'("VECTORS  Velocity  float")')
      do j=1,nj
         do i=1,ni
            write(ivtk,*) q(i,j,2)/q(i,j,1), q(i,j,3)/q(i,j,1), 0.0
         enddo
      enddo

      close(ivtk)

      deallocate(x)
      deallocate(y)
      deallocate(q)

      stop
      end
