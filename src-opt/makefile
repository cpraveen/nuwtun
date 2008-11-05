# Include parameters from file
include ../makefile.in

all: costfun

costfun: cost.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o
	$(FC) -o cost cost.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o

# Remove object, library and executable files
clean:
	@for dir in $(DIRS); do cd $$dir; $(MAKE) clean; cd ..; done
	$(RM) -f $(OBJS)
	$(RM) -f $(PGM)
