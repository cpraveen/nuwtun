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

all: $(DEFORM) $(DEFADJ)

deform: $(ODEF)
	@echo "Building   deform"
	@$(FC) -o deform $(ODEF) $(LIBPATH) -llinpack -lblas

deform_adj: $(OADJ)
	@echo "Building   deform_adj"
	@$(FC) -o deform_adj $(OADJ) $(LIBPATH) -llinpack -lblas

%.o: %.f
	@echo "Compiling " $<
	@$(FC) $(FFLAGS) -c $<

clean:
	$(RM) -f $(ODEF) $(OADJ) $(DEFORM) $(DEFADJ)
