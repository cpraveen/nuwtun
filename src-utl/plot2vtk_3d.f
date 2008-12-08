C     Read a single block 3-d formatted plot3d file and write a vtk 
C     structured grid file. Reads g-file and q-file.
C     Check gfile, qfile for the input files
      program plot3vtk

      integer i
      integer j
      integer k
      integer m
      integer n
      integer ni
      integer nj
      integer nk
      integer nblks
 
      real mach   ! freestream Mach number
      real alpha  ! freestream angle-of-attack
      real reyn   ! freestream Reynolds number
      real time   ! time

      real, allocatable :: x(:,:,:), y(:,:,:), z(:,:,:), q(:,:,:,:)

      character*64 basename
      character*64 gfile, qfile, ofile

      integer ivtk
      real    gamma, d, du, dv, dw, e, p, ma
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
      read(7,*) ni, nj, nk
      allocate( x(ni,nj,nk) )
      allocate( y(ni,nj,nk) )
      allocate( z(ni,nj,nk) )
      allocate( q(ni,nj,nk,5) )
      read(7,*) 
     &    ((( x(i,j,k), i=1,ni), j=1,nj), k=1,nk),
     &    ((( y(i,j,k), i=1,ni), j=1,nj), k=1,nk),
     &    ((( z(i,j,k), i=1,ni), j=1,nj), k=1,nk)
      close(7)

      read(8,*) nblks
      if(nblks.ne.1)then
         print*,'Cannot read multi-block file'
         print*,'Number of blocks =', nblks
         stop
      endif
      read(8,*) ni, nj, nk
      read(8,*) mach, alpha, reyn, time
      read(8,*) (((( q(i,j,k,n), i=1,ni), j=1,nj), k=1,nk), n=1,5)
      close(8)

      ivtk=10
      write(*,*)'Writing vtk file ', ofile
      open(unit=ivtk, file=ofile)
      write(ivtk,'("# vtk DataFile Version 2.0")')
      write(ivtk,'("NACA0012")')
      write(ivtk,'("ASCII")')
      write(ivtk,'("DATASET STRUCTURED_GRID")')
      write(ivtk,'("DIMENSIONS ",i10,i10,i10)') ni, nj, nk
      write(ivtk,'("POINTS ",i10," float")') ni*nj*nk
      do k=1,nk
         do j=1,nj
            do i=1,ni
               write(ivtk,*) x(i,j,k), y(i,j,k), z(i,j,k)
            enddo
         enddo
      enddo
      write(ivtk,'("POINT_DATA",i10)') ni*nj*nk
c     Write density
      write(ivtk,'("SCALARS density float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do k=1,nk
         do j=1,nj
            do i=1,ni
               write(ivtk,*) q(i,j,k,1)
            enddo
         enddo
      enddo
c     Write pressure
      write(ivtk,'("SCALARS pressure float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do k=1,nk
         do j=1,nj
            do i=1,ni
               d = q(i,j,k,1)
               du= q(i,j,k,2)
               dv= q(i,j,k,3)
               dw= q(i,j,k,4)
               e = q(i,j,k,5)
               q2= (du**2 + dv**2 + dw**2)/d/d
               p = (gamma-1.0)*(e - 0.5*d*q2)
               write(ivtk,*) p
            enddo
         enddo
      enddo
c     Write mach
      write(ivtk,'("SCALARS mach float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do k=1,nk
         do j=1,nj
            do i=1,ni
               d = q(i,j,k,1)
               du= q(i,j,k,2)
               dv= q(i,j,k,3)
               dw= q(i,j,k,4)
               e = q(i,j,k,5)
               q2= (du**2 + dv**2 + dw**2)/d/d
               p = (gamma-1.0)*(e - 0.5*d*q2)
               ma= sqrt(q2)/sqrt(gamma*p/d)
               write(ivtk,*) ma
            enddo
         enddo
      enddo
      close(ivtk)

      deallocate(x)
      deallocate(y)
      deallocate(z)
      deallocate(q)

      stop
      end
