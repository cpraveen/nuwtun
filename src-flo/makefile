include makefile.in

#Executable program
PGM  = nuwtun_flo

#Object files
MAINF = $(wildcard main/*.F)
OBJS  = $(patsubst %.F,%.o,$(MAINF))
HDRS  = $(wildcard header/*.h)

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
	$(FC) $(LDFLAGS) -o $(PGM) $(OBJS) $(LIBS)

#Just build the libraries
libs: $(LIBS)
	@echo "All libraries are updated"

##############################################################################
#Rule to make individual .lib.a files
##############################################################################
rhs/.lib.a: $(wildcard rhs/*.F) $(HDRS)
	@cd rhs; $(MKLIB)
#-----------------------------------------------------------------------------
bc/.lib.a: $(wildcard bc/*.F) $(HDRS)
	@cd bc; $(MKLIB)
#-----------------------------------------------------------------------------
explicit/.lib.a: $(wildcard explicit/*.F) $(HDRS)
	@cd explicit; $(MKLIB)
#-----------------------------------------------------------------------------
implicit/.lib.a: $(wildcard implicit/*.F) $(HDRS)
	@cd implicit; $(MKLIB)
#-----------------------------------------------------------------------------
io/.lib.a: $(wildcard io/*.F) $(HDRS)
	@cd io; $(MKLIB)
#-----------------------------------------------------------------------------
turb/.lib.a: $(wildcard turb/*.F) $(HDRS)
	@cd turb; $(MKLIB)
#-----------------------------------------------------------------------------
march/.lib.a: $(wildcard march/*.F) $(HDRS)
	@cd march; $(MKLIB)
#-----------------------------------------------------------------------------
main/%.o: main/%.F $(HDRS)
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
	$(RM) -f $(LIBS) $(OBJS) */*.o $(PGM)
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
	@echo "MAINO  =" $(MAINO)
