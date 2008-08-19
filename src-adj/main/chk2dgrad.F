C This subroutine makes two checks for 2-d gradients. 
C 1. x and y gradients for k=1 and k=2 plane must be equal. 
C 2. Z gradient for k=1 and k=2 must sum to zero.
      subroutine chk2dgrad (idim, jdim, kdim, rb)
      dimension rb(0:jdim+1, 0:kdim+1, 0:idim+1, 3)

      do i=1,idim
         do j=1,jdim
            if(rb(j,1,i,1).ne.rb(j,2,i,1) .or.
     +         rb(j,1,i,2).ne.rb(j,2,i,2)) then
               print*,'k=1 and k=2 gradient differs:',i,j
            endif
            if(rb(j,1,i,3)+rb(j,2,i,3) .ne. 0.0) then
               print*,'z-gradient is not zero:',i,j
            endif
         enddo
      enddo

      return
      end
