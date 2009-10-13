include ../makefile.in

SRCP = $(SRCFLO)/rhs
SRCD = $(SRCADJ)/rhs

TPNBQ    = $(TPN) -b -difffuncname "_bq" -O $(SRCD)
TPNBX    = $(TPN) -b -difffuncname "_bx" -O $(SRCD)

#Backward diff routines wrt state Q
BQ = fluroe_bq.f \
     fluroe_bx.f \
	  rmmlim_bq.f \
	  rnolim_bq.f \
	  smilim_bq.f \
	  smtlim_bq.f \
	  frcint_bq.f \
	  frcint_bx.f \
	  frcel_bq.f  \
	  frcel_bx.f  \
	  apgrad_bq.f \
	  apgrad_bx.f \
	  metbc_bx.f  \
	  metric_bx.f \
	  grdcut_bx.f \
	  gdxtnd_bx.f \
	  costfun_bq.f

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
fluroe_bx.f: $(SRCP)/fluroe.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBX) \
		-head    fluroe     \
		-vars    "F S"      \
		-outvars "F S"      \
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
frcint_bx.f: $(SRCP)/force.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBX)  \
		-head    frcint      \
		-vars    "Q S CL CD" \
		-outvars "Q S CL CD" \
		$<
#-----------------------------------------------------------------------------
frcel_bq.f: $(SRCP)/force.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBQ) \
		-head    frcel      \
		-vars    "Q FORCE"  \
		-outvars "Q FORCE"  \
		$<
#-----------------------------------------------------------------------------
frcel_bx.f: $(SRCP)/force.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBX)  \
		-head    frcel       \
		-vars    "Q S FORCE" \
		-outvars "Q S FORCE" \
		$<
#-----------------------------------------------------------------------------
apgrad_bq.f: $(SRCP)/apgrad.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBQ) \
		-head    apgrad     \
		-vars    "Q APGINT"  \
		-outvars "Q APGINT"  \
		$<
#-----------------------------------------------------------------------------
apgrad_bx.f: $(SRCP)/apgrad.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBX)  \
		-head    apgrad      \
		-vars    "Q S APGINT" \
		-outvars "Q S APGINT" \
		$<
#-----------------------------------------------------------------------------
gdxtnd_bx.f: $(SRCP)/gdxtnd.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBX)  \
		-head    gdxtnd      \
		-vars    "R"         \
		-outvars "R"         \
		$<
#-----------------------------------------------------------------------------
grdcut_bx.f: $(SRCP)/grdcut.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBX)  \
		-head    grdcut      \
		-vars    "R1 R2"     \
		-outvars "R1 R2"     \
		$<
#-----------------------------------------------------------------------------
metric_bx.f: $(SRCP)/metric.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBX)  \
		-head    metric      \
		-vars    "R S"       \
		-outvars "R S"       \
		$<
#-----------------------------------------------------------------------------
metbc_bx.f: $(SRCP)/metbc.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBX)  \
		-head    metbc       \
		-vars    "S"         \
		-outvars "S"         \
		$<
#-----------------------------------------------------------------------------
costfun_bq.f: $(SRCD)/costfun.F
	$(RM) -f $@
	cd $(SRCD) && $(TPNBQ)  \
		-head    costfun     \
		-vars    "cl cd apgint cost"\
		-outvars "cl cd apgint cost"\
		$<
##############################################################################
clean:
	$(RM) -f .lib.a *.o *.msg *~
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
