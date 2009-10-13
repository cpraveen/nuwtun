include $(NUWTUN_HOME)/makefile.in

# Path to linpack and blas libraries, for example
# LIBPATH = -L/usr/local/lib
LIBPATH= -L$(HOME)/usr/local/lib
LIBPATH=

DEFORM = deform
DEFADJ = deform_adj

# Object files for deform
ODEF = deform.o       \
       hickshenne.o   \
       rbf_eval.o     \
       rbf_train.o    \
       rdp2d.o        \
       readgrd.o      \
       readin.o       \
       rem_dup.o      \
       shape_deform.o \
       shepard.o \
       test.o         \
       total_pts.o    \
       wrtgrd.o


# Object files for deform_adj
OADJ = deform_adj.o \
       rbf_eval_b.o \
       rbf_train_b.o \
       rdp2d.o        \
       readin.o \
       readgrd.o \
       readgrad.o \
       shape_deform.o \
       shape_deform_b.o \
       rem_dup.o \
       test.o \
       total_pts.o \
       hickshenne.o \
       hickshenne_b.o \

all: $(DEFORM) $(DEFADJ) volume area area_x

deform: $(ODEF)
	@echo "Building   deform"
	@$(FC) -o deform $(ODEF) $(LIBPATH) -llinpack -lblas

deform_adj: $(OADJ)
	@echo "Building   deform_adj"
	@$(FC) -o deform_adj $(OADJ) $(LIBPATH) -llinpack -lblas

volume: volume.o
	@echo "Building   volume"
	@$(FC) -o volume volume.o $(LIBPATH)

# Area computation routine
area: area.o common_area.o routines_area.o
	@echo "Building   area"
	@$(FC) -o area area.o common_area.o routines_area.o $(LIBPATH)

area.o: common_area.o

# Area derivative computation routine
OAREA = area_x.o common_area.o routines_area.o totalarea_b.o adBuffer.o \
		  adStack.o
area_x: $(OAREA)
	@echo "Building   area_x"
	@$(FC) -o area_x $(OAREA) $(LIBPATH)

totalarea_b.o: routines_area.f95
	$(TPN) -backward -vars    "r totarea" \
		              -outvars "r totarea" \
						  -head    totalarea \
						  -o       totalarea \
						           routines_area.f95
	@$(FC) -c $(FFLAGS) $(IFORTFLAGS) totalarea_b.f95

# Generate *.o from *.f
%.o: %.f
	@echo "Compiling " $<
	@$(FC) $(FFLAGS) -c $<

# Generate *.o from *.f95
%.o: %.f95
	@echo "Compiling " $<
	@$(FC) -c $(FFLAGS) $(IFORTFLAGS) $<

adBuffer.o: $(TAPENADE_HOME)/ADFirstAidKit/adBuffer.f
	@echo "Compiling " $<; $(FC) -c $<

adStack.o: $(TAPENADE_HOME)/ADFirstAidKit/adStack.c
	@echo "Compiling " $<; $(CC) -c $<

clean:
cleanall:
	$(RM) -f $(ODEF) $(OADJ) $(DEFORM) $(DEFADJ)
	$(RM) -f $(OAREA) *.mod *.msg *~
