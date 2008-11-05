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
       main/chk2dgrad.o\
       main/cpavar.o   \
       main/fvres.o    \
       main/intqa1.o   \
       main/adddpt.o   \
       main/sumgrad.o  \
       main/adBuffer.o \
       main/adStack.o

HDIF = rhs/fluxi_q.F \
       rhs/fluxi_x.F \
       rhs/fluxj_q.F \
       rhs/fluxj_x.F \
       rhs/fluxk_q.F \
       rhs/fluxk_x.F \
       rhs/muscli_q.F \
       rhs/musclj_q.F \
       rhs/musclk_q.F \
       rhs/resid_q.F \
       rhs/resid_x.F \
       rhs/fds_q.F \
       rhs/fds_x.F \
       bc/bc_q.F \
       bc/bcfar_q.F \
       bc/bcsym_q.F  \
       bc/bctan_q.F  \
       bc/bcext_q.F \
       bc/cut_q.F \

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

nuwtun_adj: $(OBJS) $(LIBS) $(FLOLIB) $(HDIF)
	$(FC) -o $(PGM) $(OBJS) $(LIBS) $(FLOLIB)

normal:
	cd main; $(FC) $(FFLAGS) -c main.F
	@make nuwtun_adj

dptest:
	cd main; $(FC) -D_DPTEST $(FFLAGS) -c main.F
	@make nuwtun_adj

clalpha:
	cd main; $(FC) -D_CL_ALPHA $(FFLAGS) -c main.F
	@make nuwtun_adj

main/%.o: main/%.F
	@echo "Compiling " $<; cd main; $(FC) $(FFLAGS) -c $*.F

main/adBuffer.o: $(TAPENADE_HOME)/ADFirstAidKit/adBuffer.f
	@echo "Compiling " $<; cd main; $(FC) -c $<

main/adStack.o: $(TAPENADE_HOME)/ADFirstAidKit/adStack.c
	@echo "Compiling " $<; cd main; $(CC) -c $<

# Hand differentiated files

rhs/cvtc2p_q.F : rhs/cvtc2p.F
rhs/fluxi_q.F  : $(SRCFLO)/rhs/fluxi.F
rhs/fluxi_x.F  : $(SRCFLO)/rhs/fluxi.F
rhs/fluxj_q.F  : $(SRCFLO)/rhs/fluxj.F
rhs/fluxj_x.F  : $(SRCFLO)/rhs/fluxj.F
rhs/fluxk_q.F  : $(SRCFLO)/rhs/fluxk.F
rhs/fluxk_x.F  : $(SRCFLO)/rhs/fluxk.F
rhs/muscli_q.F : $(SRCFLO)/rhs/muscli.F
rhs/musclj_q.F : $(SRCFLO)/rhs/musclj.F
rhs/musclk_q.F : $(SRCFLO)/rhs/musclk.F
rhs/fds_q.F    : $(SRCFLO)/rhs/fds.F
rhs/fds_x.F    : $(SRCFLO)/rhs/fds.F
rhs/limitr_q.F : $(SRCFLO)/rhs/limitr.F
bc/bc_q.F      : $(SRCFLO)/bc/bc.F
bc/bcfar_q.F   : $(SRCFLO)/bc/bcfar.F
bc/bcsym_q.F   : $(SRCFLO)/bc/bcsym.F
bc/bctan_q.F   : $(SRCFLO)/bc/bctan.F
bc/bcext_q.F   : $(SRCFLO)/bc/bcext.F
bc/cut_q.F     : $(SRCFLO)/bc/cut.F
	@echo "$< seems to have changed"
	@false

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
