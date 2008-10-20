C     Read a single block 2-d formatted plot3d file and write a vtk 
C     structured grid file. Reads g-file and q-file.
C     Check imax, jmax dimensions
C     Check gfile, qfile for the input files
      subroutine wrtrb (ifid, idim, jdim, kdim, rb)
      dimension rb(0:jdim+1, 0:kdim+1, 0:idim+1, 3)

      print*,'Writing gradient into rb.dat'
      if(kdim.eq.2)then
         fact = 2.0
         write(ifid,*) (((fact*rb(j,1,i,l),i=1,idim),j=1,jdim),l=1,2)
      else
         fact = 1.0
         write(ifid,*) ((((fact*rb(j,k,i,l),i=1,idim),j=1,jdim),
     1                     k=1,kdim), l=1,3)
      endif

      return
      end
