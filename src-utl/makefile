# Include parameters from file
include ../makefile.in

F95FLAGS=
ifeq ($(FC),ifort)
	F95FLAGS = -free -Tf
endif

TARGETS=extract plot2vtk_2d plot2vtk_3d
all: $(TARGETS)

extract: extract.f95 
	$(FC) $(F95FLAGS) extract.f95 -o extract

# Remove object, library and executable files
clean:
	$(RM) -f $(TARGETS)
	$(RM) -f *.o
