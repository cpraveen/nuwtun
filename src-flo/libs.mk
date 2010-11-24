#This will be executed inside the sub-directories
include ../makefile.in

SRCF = $(wildcard *.F)
HDRS = $(wildcard ../header/*.h)
OBJS = $(patsubst %.F,%.o,$(SRCF))

#
all: .lib.a

.lib.a: $(OBJS) $(HDRS)
	@echo "Building  .lib.a"
	@$(AR) $(OPTAR) .lib.a $(OBJS)

%.o: %.F $(HDRS)
	@echo "Compiling " $<
	@$(FC) $(FFLAGS) -c $<
