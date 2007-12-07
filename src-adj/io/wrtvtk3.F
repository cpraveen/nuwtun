C     Read a single block 2-d formatted plot3d file and write a vtk 
C     structured grid file. Reads g-file and q-file.
C     Check imax, jmax dimensions
C     Check gfile, qfile for the input files
      subroutine wrtvtk3(idim, jdim, kdim, r, rb, sb)
      dimension r (0:jdim+1, 0:kdim+1, 0:idim+1, 3),
     1          rb(0:jdim+1, 0:kdim+1, 0:idim+1, 3),
     1          sb(0:jdim+1, 0:kdim+1, 0:idim+1, 4, 3)
      character*32 ofile

      ofile = 'rb.vtk'

      ibeg = 1
      iend = idim
      jbeg = 1
      jend = jdim
      kbeg = 1
      kend = kdim
      ni   = iend - ibeg + 1
      nj   = jend - jbeg + 1
      nk   = kend - kbeg + 1

      ivtk = 49
      write(*,*)'Writing vtk file ', ofile

      open(unit=ivtk, file=ofile)
      write(ivtk,'("# vtk DataFile Version 2.0")')
      write(ivtk,'("NACA0012")')
      write(ivtk,'("ASCII")')
      write(ivtk,'("DATASET STRUCTURED_GRID")')
      write(ivtk,'("DIMENSIONS ",i10,i10,i10)') ni, nj, nk
      write(ivtk,'("POINTS ",i10," float")') ni*nj*nk

      do k=kbeg,kend
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) r(j,k,i,1), r(j,k,i,2), r(j,k,i,3)
         enddo
      enddo
      enddo
      write(ivtk,'("POINT_DATA",i10)') ni*nj*nk

c     Write density
      write(ivtk,'("SCALARS J_x float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do k=kbeg,kend
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) rb(j,k,i,1)
         enddo
      enddo
      enddo

c     Write u velocity
      write(ivtk,'("SCALARS J_y float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do k=kbeg,kend
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) rb(j,k,i,2)
         enddo
      enddo
      enddo

c     Write u velocity
      write(ivtk,'("SCALARS J_si1 float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do k=kbeg,kend
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) sb(j,k,i,1,2)
         enddo
      enddo
      enddo

c     Write u velocity
      write(ivtk,'("SCALARS J_si2 float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do k=kbeg,kend
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) sb(j,k,i,2,2)
         enddo
      enddo
      enddo

c     Write u velocity
      write(ivtk,'("SCALARS J_si3 float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do k=kbeg,kend
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) sb(j,k,i,3,2)
         enddo
      enddo
      enddo

c     Write u velocity
      write(ivtk,'("SCALARS J_si4 float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do k=kbeg,kend
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) sb(j,k,i,4,2)
         enddo
      enddo
      enddo

      close(ivtk)

      return
      end
