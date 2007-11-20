include ../makefile.in

SRCP = $(SRCFLO)/rhs
SRCD = $(SRCADJ)/rhs

TPNBQ    = $(TPN) -b -difffuncname "_bq" -O $(SRCD)

#Backward diff routines wrt state Q
BQ = fluroe_bq.f \
	  rmmlim_bq.f \
	  rnolim_bq.f \
	  smilim_bq.f \
	  smtlim_bq.f \
	  frcint_bq.f \
	  frcel_bq.f

all: $(BQ)
	@$(MKLIB)

##############################################################################
# backward differentiation wrt state Q
##############################################################################
fluroe_bq.f: $(SRCP)/fluroe.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBQ) \
		-head    fluroe     \
		-vars    "F QM QP"  \
		-outvars "F QM QP"  \
		$<
#-----------------------------------------------------------------------------
rmmlim_bq.f: $(SRCP)/limfun.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBQ)       \
		-head    rmmlim           \
		-vars    "DQP1 DQP2"      \
		-outvars "DQP1 DQP2 DQM1" \
		$<
#-----------------------------------------------------------------------------
rnolim_bq.f: $(SRCP)/limfun.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBQ)       \
		-head    rnolim           \
		-vars    "DQP1 DQP2"      \
		-outvars "DQP1 DQP2 DQM1" \
		$<
#-----------------------------------------------------------------------------
smilim_bq.f: $(SRCP)/limfun.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBQ)       \
		-head    smilim           \
		-vars    "DQP1 DQP2"      \
		-outvars "DQP1 DQP2 DQM1" \
		$<
#-----------------------------------------------------------------------------
smtlim_bq.f: $(SRCP)/limfun.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBQ)       \
		-head    smtlim           \
		-vars    "DQP1 DQP2"      \
		-outvars "DQP1 DQP2 DQM1" \
		$<
#-----------------------------------------------------------------------------
frcint_bq.f: $(SRCP)/force.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBQ) \
		-head    frcint     \
		-vars    "Q CL CD"  \
		-outvars "Q CL CD"  \
		$<
#-----------------------------------------------------------------------------
frcel_bq.f: $(SRCP)/force.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBQ) \
		-head    frcel      \
		-vars    "Q FORCE"  \
		-outvars "Q FORCE"  \
		$<
##############################################################################
clean:
	$(RM) -f .lib.a *.o
#Remove everything
cleanall:
	$(RM) -f .lib.a *.o $(BQ) *.msg *~
#Remove everything and recompile
new:
	$(MAKE) cleanall
	$(MAKE) all

info:
	@echo "FC     =" $(FC)
	@echo "CFLAGS =" $(CFLAGS)
	@echo "FFLAGS =" $(FFLAGS)
	@echo "PGM    =" $(PGM)
