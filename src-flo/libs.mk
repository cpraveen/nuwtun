include ../makefile.in

SRCF = $(wildcard *.F)
OBJS = $(patsubst %.F,%.o,$(SRCF))

#
all: .lib.a

.lib.a: $(OBJS)
	@echo "Building  .lib.a"
	@$(AR) $(OPTAR) .lib.a $(OBJS)

%.o: %.F
	@echo "Compiling " $<
	@$(FC) $(FFLAGS) -c $<
