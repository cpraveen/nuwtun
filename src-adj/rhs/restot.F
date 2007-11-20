      SUBROUTINE RESTOT(IDIM, JDIM, KDIM, QB, QB0)
C 
C     Add QB0 to QB and return result in QB 
C 
      include '../../src-flo/header/common.h'
      DIMENSION QB (0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ),
     1          QB0(0:JDIM+2,0:KDIM+2,0:IDIM+2,NQ)
C 
      QB(:,:,:,:) = QB(:,:,:,:) + QB0(:,:,:,:)
      RETURN
      END
