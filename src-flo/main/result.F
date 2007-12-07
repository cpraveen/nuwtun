      SUBROUTINE RESULT(ITER, R2NORM, CL, CD)
      include '../header/common.h'
      DIMENSION R2NORM(*)

      write(iresult,10) ITER
      write(iresult,11) (R2NORM(I),I=1,NF)
      write(iresult,12) CL
      write(iresult,13) CD

   10 format('Iterations ', I8)
   11 format('Residue    ', 12E12.4)
   12 format('CL         ', E28.20)
   13 format('CD         ', E28.20)

      return
      end
