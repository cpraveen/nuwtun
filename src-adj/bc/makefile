include ../makefile.in

SRCP = $(SRCFLO)/bc
SRCD = $(SRCADJ)/bc

TPNBQ    = $(TPN) -b -difffuncname "_bq" -O $(SRCD)

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
		-vars    "Q"        \
		-outvars "Q"        \
		$<
#-----------------------------------------------------------------------------
bcfinv_bq.f: $(SRCP)/bcfar.F
	$(RM) -f $@
	cd $(SRCP) && $(TPNBQ) \
		-head    bcfinv     \
		-vars    "Q"        \
		-outvars "Q"        \
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
