# Include parameters from file
include makefile.in

# Name of the executable
PGM  = nuwtun_adj

# List of directories containing sources, except main
DIRS = bc rhs solve io

# Libraries to create in each of the DIRS
LIBS = $(foreach dir,$(DIRS),$(addsuffix /.lib.a,$(dir)))

# Some object files
OBJS = main/main.o     \
       main/cpavar.o   \
       main/fvres.o    \
       main/intqa1.o   \
       main/adddpt.o   \
       main/adBuffer.o \
       main/adStack.o

# Library files from src-flo
FLOLIB   = $(SRCFLO)/rhs/.lib.a      \
           $(SRCFLO)/bc/.lib.a       \
           $(SRCFLO)/io/.lib.a       \
           $(SRCFLO)/turb/.lib.a     \
           $(SRCFLO)/implicit/.lib.a

# This is executed for every make
all:
	@cd $(SRCFLO); $(MAKE) libs
	@for dir in $(DIRS); do cd $$dir && $(MAKE); cd .. ; done
	@$(MAKE) nuwtun_adj

nuwtun_adj: $(OBJS) $(LIBS) $(FLOLIB)
	$(FC) -o $(PGM) $(OBJS) $(LIBS) $(FLOLIB)

main/%.o: main/%.F
	@echo "Compiling " $<; cd main; $(FC) $(FFLAGS) -c $*.F

main/adBuffer.o: $(TAPENADE_HOME)/lib/adBuffer.f
	@echo "Compiling " $<; cd main; $(FC) $(FFLAGS) -c $<

main/adStack.o: $(TAPENADE_HOME)/lib/adStack.c
	@echo "Compiling " $<; cd main; $(CC) -c $<

# Remove object, library and executable files
clean:
	@for dir in $(DIRS); do cd $$dir; $(MAKE) clean; cd ..; done
	$(RM) -f $(OBJS)
	$(RM) -f $(PGM)

# Remove everything, including AD generated files
cleanall:
	@for dir in $(DIRS); do cd $$dir; $(MAKE) cleanall; cd ..; done
	$(RM) -f $(OBJS)
	$(RM) -f $(PGM)

# Print some information
info:
	@echo "DIRS =" $(DIRS)
	@echo "LIBS =" $(LIBS)
	@echo "OBJS =" $(OBJS)
