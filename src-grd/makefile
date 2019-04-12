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
       shepard.o      \
       test.o         \
       total_pts.o    \
       wrtgrd.o       \
       kulfan.o       \
       decas.o        \
       cubspl.o       \
       interv.o       \
       ppvalu.o       \
       dgeco.o        \
       dgefa.o        \
       dgesl.o


# Object files for deform_adj
OADJ = deform_adj.o     \
       rbf_eval_b.o     \
       rbf_train_b.o    \
       rdp2d.o          \
       readin.o         \
       readgrd.o        \
       readgrad.o       \
       shape_deform.o   \
       shape_deform_b.o \
       rem_dup.o        \
       test.o           \
       total_pts.o      \
       hickshenne.o     \
       hickshenne_b.o   \
       kulfan.o         \
       kulfan_b.o       \
       decas.o          \
       cubspl.o         \
       interv.o         \
       ppvalu.o         \
       dgeco.o          \
       dgefa.o          \
       dgesl.o          \
       adStack.o        \
       adBuffer.o

TARGETS = $(DEFORM) $(DEFADJ) volume area area_x

all: $(TARGETS)

# Compile linpack files
dgeco.o: dgeco.f
	$(FC) -O3 -c dgeco.f

dgefa.o: dgefa.f
	$(FC) -O3 -c dgefa.f

dgesl.o: dgesl.f
	$(FC) -O3 -c dgesl.f

deform: $(ODEF)
	@echo "Building   deform"
	@$(FC) $(LDFLAGS) -o deform $(ODEF) $(LIBPATH) -lblas

deform_adj: $(OADJ)
	@echo "Building   deform_adj"
	@$(FC) $(LDFLAGS) -o deform_adj $(OADJ) $(LIBPATH) -lblas

volume: volume.o
	@echo "Building   volume"
	@$(FC) $(LDFLAGS) -o volume volume.o $(LIBPATH)

# Area computation routine
area: area.o common_area.o routines_area.o
	@echo "Building   area"
	@$(FC) $(LDFLAGS) -o area area.o common_area.o routines_area.o $(LIBPATH)

area.o: common_area.o

# Area derivative computation routine
OAREA = area_x.o common_area.o routines_area.o totalarea_b.o adBuffer.o \
		  adStack.o
area_x: $(OAREA)
	@echo "Building   area_x"
	@$(FC) $(LDFLAGS) -o area_x $(OAREA) $(LIBPATH)

totalarea_b.o: routines_area.f95
	$(TPN) -backward -vars    "r totarea" \
		              -outvars "r totarea" \
						  -head    totalarea \
						  -o       totalarea \
						           routines_area.f95
	@$(FC) -c $(FFLAGS) $(IFORTFLAGS) totalarea_b.f95

kulfan_b.o: kulfan.f decas.f
	$(TPN) -backward -vars    "x kulfan" \
                    -outvars "x kulfan" \
                    -head    kulfan    \
                    -o       kulfan    \
                             kulfan.f decas.f
	@$(FC) -c $(FFLAGS) kulfan_b.f

kulfan_b.o: kulfan.f decas.f
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
	$(RM) -f *.o
	$(RM) -f *.mod *.msg *~
	$(RM) -f $(TARGETS)

cleanall:
	make clean
