include ../makefile.in

all:
	@$(MKLIB)

##############################################################################
clean:
	$(RM) -f .lib.a *.o
#Remove everything
cleanall:
	$(RM) -f .lib.a *.o
#Remove everything and recompile
new:
	$(MAKE) cleanall
	$(MAKE) all
