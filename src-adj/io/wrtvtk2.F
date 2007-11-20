C     Read a single block 2-d formatted plot3d file and write a vtk 
C     structured grid file. Reads g-file and q-file.
C     Check imax, jmax dimensions
C     Check gfile, qfile for the input files
      subroutine wrtvtk2(idim, jdim, kdim, rc, q, nq)
      dimension q(0:jdim+2, 0:kdim+2, 0:idim+2, nq),
     1          rc(jdim+1, kdim+1, idim+1, 3)
      character*32 ofile

      ofile = 'out.vtk'

      ibeg = 1
      iend = idim+1
      jbeg = 1
      jend = jdim+1
      ni   = iend - ibeg + 1
      nj   = jend - jbeg + 1

      ivtk = 49
      write(*,*)'Writing vtk file ', ofile
      open(unit=ivtk, file=ofile)
      write(ivtk,'("# vtk DataFile Version 2.0")')
      write(ivtk,'("NACA0012")')
      write(ivtk,'("ASCII")')
      write(ivtk,'("DATASET STRUCTURED_GRID")')
      write(ivtk,'("DIMENSIONS ",i10,i10,i10)') ni, nj, 1
      write(ivtk,'("POINTS ",i10," float")') ni*nj
      k = 2
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) rc(j,k,i,1), rc(j,k,i,2), 0.0
         enddo
      enddo
      write(ivtk,'("POINT_DATA",i10)') ni*nj
c     Write density
      write(ivtk,'("SCALARS density float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) q(j,k,i,1)
         enddo
      enddo
c     Write u velocity
      write(ivtk,'("SCALARS xmomentum float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) q(j,k,i,2)
         enddo
      enddo
c     Write v velocity
      write(ivtk,'("SCALARS ymomentum float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) q(j,k,i,3)
         enddo
      enddo
c     Write pressure
      write(ivtk,'("SCALARS energy float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) q(j,k,i,5)
         enddo
      enddo
      close(ivtk)

      stop
      end
