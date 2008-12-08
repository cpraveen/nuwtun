# Include parameters from file
include ../makefile.in

all: cost steep

cost: cost.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o
	$(FC) -o cost cost.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o

steep: steep.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o
	$(FC) -o steep steep.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o

# Remove object, library and executable files
clean:
	$(RM) -f cost
	$(RM) -f *.o
