# Include parameters from file
include ../makefile.in

# Flags for .f95 files
F95FLAGS=

# Extra compiler flag for F95 with ifort
ifeq ($(FC),ifort)
	F95FLAGS = -free -Tf
endif

# Executables to generate
TARGETS=extract ext2cp plot2vtk_2d plot2vtk_3d ad_vs_fd mach_sweep \
		  alpha_sweep mglvls plot3d_scale

all: $(TARGETS)

ad_vs_fd: ad_vs_fd.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o
	$(FC) $(LDFLAGS) -o ad_vs_fd ad_vs_fd.o \
		                 ../src-adj/rhs/costfun.o \
							  ../src-adj/io/sicftyp.o

extract: extract.f95 
	$(FC) $(LDFLAGS) $(F95FLAGS) extract.f95 -o extract

ext2cp: ext2cp.f95 
	$(FC) $(LDFLAGS) $(F95FLAGS) ext2cp.f95 -o ext2cp

mglvls: mglvls.c
	$(CC)  -o mglvls mglvls.c

plot3d_scale: plot3d_scale.f
	$(FC)  $(FFLAGS) -o plot3d_scale plot3d_scale.f

alpha_sweep: alpha_sweep.f
	$(MPIFC)  -o alpha_sweep alpha_sweep.f

# Remove object, library and executable files
clean:
	$(RM) -f $(TARGETS)
	$(RM) -f *.o

cleanall:
	make clean
