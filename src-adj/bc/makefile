include ../makefile.in

SRCP = $(SRCFLO)/bc
SRCD = $(SRCADJ)/bc

TPNBQ    = $(TPN) -b -difffuncname "_bq" -O $(SRCD)
TPNBX    = $(TPN) -b -difffuncname "_bx" -O $(SRCD)

#Backward diff routines wrt state Q
BQ = bctinv_bq.f \
	  bcfinv_bq.f

all: $(BQ)
	@$(MKLIB)

##############################################################################
# backward differentiation wrt state Q
##############################################################################
bctinv_bq.f: $(SRCP)/bctan.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBQ) \
		-head    bctinv     \
		-vars    "Q S"      \
		-outvars "Q S"      \
		$<
#-----------------------------------------------------------------------------
bcfinv_bq.f: $(SRCP)/bcfar.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBQ) \
		-head    bcfinv     \
		-vars    "Q S"      \
		-outvars "Q S"      \
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
