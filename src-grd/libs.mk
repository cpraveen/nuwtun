# This is executed inside the sub-directories, hence the relative path ../
include ../makefile.in

SRCf = $(wildcard *.f)
SRCF = $(wildcard *.F)
OBJS = $(patsubst %.f,%.o,$(SRCf)) $(patsubst %.F,%.o,$(SRCF))

#
all: .lib.a

.lib.a: $(OBJS)
	@echo "Building  .lib.a"
	@$(AR) $(OPTAR) .lib.a $(OBJS)

%.o: %.F
	@echo "Compiling " $<
	@$(FC) $(FFLAGS) -c $<

%.o: %.f
	@echo "Compiling " $<
	@$(FC) $(FFLAGS) -c $<
