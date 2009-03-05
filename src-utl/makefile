# Include parameters from file
include ../makefile.in

# Flags for .f95 files
F95FLAGS=

# Extra compiler flag for F95 with ifort
ifeq ($(FC),ifort)
	F95FLAGS = -free -Tf
endif

# Executables to generate
TARGETS=extract ext2cp plot2vtk_2d plot2vtk_3d ad_vs_fd

all: $(TARGETS)

ad_vs_fd: ad_vs_fd.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o
	$(FC) -o ad_vs_fd ad_vs_fd.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o

extract: extract.f95 
	$(FC) $(F95FLAGS) extract.f95 -o extract

ext2cp: ext2cp.f95 
	$(FC) $(F95FLAGS) ext2cp.f95 -o ext2cp

# Remove object, library and executable files
clean:
cleanall:
	$(RM) -f $(TARGETS)
	$(RM) -f *.o
