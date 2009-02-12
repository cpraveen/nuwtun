# Include parameters from file
include ../makefile.in

#FFSQPLIB=/opt/fsqp/lib/libffsqp.a
FFSQPFLAGS= $(FSQP_DIR)/lib/libffsqp.a
IPOPTFLAGS= -L$(IPOPT_DIR)/lib -lipopt -lblas -llapack

all: cost steep ipopt fsqp

cost: cost.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o
	$(FC) -o cost cost.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o

steep: steep.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o
	$(FC) -o steep steep.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o

ipopt: ipopt.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o
	$(FC) -o ipopt ipopt.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o \
		$(IPOPTFLAGS)

ipopt.o: ipopt.f
	$(FC) $(FFLAGS) -I$(IPOPT_DIR)/include/coin -c $<

fsqp: fsqp.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o
	$(FC) -o fsqp fsqp.o ../src-adj/rhs/costfun.o ../src-adj/io/sicftyp.o \
		$(FFSQPFLAGS)

fsqp.o: fsqp.f95
	$(FC) -c $(FFLAGS) -free -Tf $<

# Remove object, library and executable files
clean:
	$(RM) -f cost
	$(RM) -f steep
	$(RM) -f ipopt
	$(RM) -f fsqp
	$(RM) -f *.o
	$(RM) -f *.mod
