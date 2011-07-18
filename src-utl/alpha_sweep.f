C
C     Driver routine to run flow solver for different ALPHA (AOA)
C
C     Specify min, max and increment for the AOA
C     Cl and Cd will be written into alpha.dat
C
      program alpha_sweep
      implicit none
      include 'mpif.h'
      character*64 flosolve, sdummy
      character*32, allocatable :: dir(:)
      integer, allocatable :: proc(:)
      integer   ifid, ifm, iter, i, imax, j, k
      integer   rank, ierr, nproc, n, remainder, n_elem
      integer   stat
      logical   has_C_m
      real, allocatable ::  alpha(:), C_l(:), C_d(:), ap_grad(:), C_m(:)
      real      Cd, Cl
      real      rdummy
      real      alpha_min, alpha_max, dalpha
      logical ex 
      call MPI_init(ierr)
      call MPI_comm_size(MPI_comm_world, nproc, ierr)
      call MPI_comm_rank(MPI_comm_world, rank, ierr)

c     Read alpha from input file
      open(unit =10,file = 'alpha.in', status='old')
      read (10,*) alpha_min, alpha_max, dalpha
      close(10)

c     check values of alpha read from input file
      if(rank==0 .and. (alpha_max<=alpha_min .or. dalpha<=0.0))then
         print*,'Error: give alpha_max > alpha_min'
         print*,'            dalpha > 0.0'
         call MPI_finalize(ierr)
         stop
      endif

      flosolve='$NUWTUN_HOME/src-flo/nuwtun_flo < flo.in > flo.log'
      ifid    = 10
      ifm     = 11

c     Check file exists
      if(rank==0)then
         open(ifid, file='template/flo.in.noalpha', status='old')
         close(ifid)
      endif

      imax    = (alpha_max - alpha_min)/dalpha
      n = (imax+1)/nproc
      remainder = MOD((imax+1),nproc)

      allocate (alpha(0:imax), proc(0:imax), C_l(0:imax))
      allocate (C_d(0:imax), dir(0:imax))
      allocate (ap_grad(0:imax))
      allocate (C_m(0:imax))

c     Loop over different angle of attack
      do i = 0, imax
         alpha(i) = alpha_min + dalpha * i
      end do

c     Create partition table for distributing alpha
      do i = 1, n
          do j = 0,nproc-1
              if (i.eq.1) then
                  k = j
              else    
                  k =(i-1)* nproc+j
              end if
              proc(k) = j
          end do
      end do
 
      if(remainder.ne.0) then
          i = 0
          do j = (n*nproc),(n*nproc)+remainder-1
              k = j
              proc(k) = i
              i=i+1
          end do 
      end if

c     Create directories for running cfd
      do i = 0,imax
          write(dir(i),'(a,i3.3,a)'),'RUN',i
      end do

c     Create individual directory for each processor
      do i = 0,imax
          if(rank == proc(i)) then
              inquire (file = trim(dir(i)), exist = ex)
              if (ex .eqv. .TRUE.) then
                  call system("rm -r "//trim(dir(i)))
              end if
              call system("mkdir "//trim(dir(i)))
              call system("cp template/* "//trim(dir(i)))
          end if
      end do

      do i = 0, imax
          if (proc(i).eq.rank) then
           
c        Write angle of attack into file
              open(ifid, file=trim(dir(i))//'/flo.in')
              write(ifid,*)'ALPHA    ', alpha(i)
              close(ifid)
              call system("cd "//trim(dir(i))//" && "//
     +                    "cat flo.in.noalpha >> flo.in")

c        Run CFD solver
              call system("cd "//trim(dir(i))//" && "//flosolve)

c        Read solution
              open(ifid,file=trim(dir(i))//'/fort.19', status='old')
              read(ifid,*) sdummy, iter
              read(ifid,*) sdummy, rdummy, rdummy, rdummy, rdummy
     +         ,rdummy
              read(ifid,*) sdummy, C_l(i)
              read(ifid,*) sdummy, C_d(i)
              read(ifid,*) sdummy, ap_grad(i)
              read(ifid,*,IOSTAT=stat) sdummy, C_m(i)
              if(stat.gt.0)then
                 print*,'Check file, something wrong'
                 stop
              else if(stat.lt.0)then
                 has_C_m = .false.
              else
                 has_C_m = .true.
              endif
              close(ifid)
          end if
      end do

      do i = 0, imax
          call MPI_Bcast(C_d(i), 1, MPI_real, proc(i), MPI_comm_world
     +      ,ierr)
          call MPI_Bcast(C_l(i), 1, MPI_real, proc(i), MPI_comm_world
     +      ,ierr)
          if(has_C_m) then
             call MPI_Bcast(C_m(i), 1, MPI_real, proc(i), MPI_comm_world
     +         ,ierr)
          endif
      end do
      
      
c     Save Cl/Cd into file
      if (rank==0) then
         open(ifm, file="alpha.dat")
         do i = 0,imax
            if(has_C_m)then
               write(ifm,*) alpha(i), C_l(i), C_d(i), C_m(i)
            else
               write(ifm,*) alpha(i), C_l(i), C_d(i)
            endif
         end do
         close(ifm)
      end if


      call MPI_finalize(ierr)

      stop
      end program
