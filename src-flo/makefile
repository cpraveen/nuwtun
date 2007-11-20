include makefile.in

#Executable program
PGM  = nuwtun_flo

#Object files
OBJS = main/main.o

#Libraries
LIBS = rhs/.lib.a      \
       bc/.lib.a       \
       march/.lib.a    \
       explicit/.lib.a \
       implicit/.lib.a \
       io/.lib.a       \
       turb/.lib.a

all: nuwtun_flo

#Build the executable
nuwtun_flo: $(OBJS) $(LIBS)
	$(FC) -o $(PGM) $(OBJS) $(LIBS)

#Just build the libraries
libs: $(LIBS)
	@echo "All libraries are updated"

##############################################################################
#Rule to make individual .lib.a files
##############################################################################
rhs/.lib.a: $(wildcard rhs/*.F)
	@cd rhs; $(MKLIB)
#-----------------------------------------------------------------------------
bc/.lib.a: $(wildcard bc/*.F)
	@cd bc; $(MKLIB)
#-----------------------------------------------------------------------------
explicit/.lib.a: $(wildcard explicit/*.F)
	@cd explicit; $(MKLIB)
#-----------------------------------------------------------------------------
implicit/.lib.a: $(wildcard implicit/*.F)
	@cd implicit; $(MKLIB)
#-----------------------------------------------------------------------------
io/.lib.a: $(wildcard io/*.F)
	@cd io; $(MKLIB)
#-----------------------------------------------------------------------------
turb/.lib.a: $(wildcard turb/*.F)
	@cd turb; $(MKLIB)
#-----------------------------------------------------------------------------
march/.lib.a: $(wildcard march/*.F)
	@cd march; $(MKLIB)
#-----------------------------------------------------------------------------
main/%.o: main/%.F
	@echo "Compiling " main/$*.F; cd main; $(FC) $(FFLAGS) -c $*.F

##############################################################################
# Generate call tree using ftnchek
##############################################################################
tree:
	$(RM) -rf $(CTDIR)
	mkdir -p $(CTDIR)
	cp --parents `find . -name "*.F" -o -name "*.h"` $(CTDIR)
	cd $(CTDIR) && ftnchek -nochek -mkhtml -call=tree */*.F > ftnchek.log
	@echo "To see the call-tree, use a web-browser, for example"
	@echo "  firefox $(CTDIR)/CallTree.html"

##############################################################################
clean:
	$(RM) -f $(LIBS) $(PGM)
#Remove everything
cleanall:
	$(RM) -f $(LIBS) $(OBJS) $(PGM) */*.o
#Remove everything and recompile
new:
	$(MAKE) cleanall
	$(MAKE) all

info:
	@echo "FC     =" $(FC)
	@echo "CFLAGS =" $(CFLAGS)
	@echo "FFLAGS =" $(FFLAGS)
	@echo "PGM    =" $(PGM)
