#List of files/directories to make tgz file
FILES=nuwtun/src-flo     \
		nuwtun/src-adj     \
		nuwtun/src-grd     \
		nuwtun/src-opt     \
		nuwtun/src-utl     \
		nuwtun/docs        \
		nuwtun/makefile    \
		nuwtun/makefile.in \
		nuwtun/README

#Things to exclude in the above list
EXCLU=nuwtun/docs/tree

#Name of tgz file
TGZFILE=nuwtun.tgz

all: flo adj grd opt utl

flo:
	@cd src-flo && make

adj:
	@cd src-adj && make

grd:
	@cd src-grd && make

opt:
	@cd src-opt && make

utl:
	@cd src-utl && make

tgz:
	cd src-flo && make cleanall
	cd src-adj && make clean
	cd .. && tar zcvf $(TGZFILE) $(FILES) --exclude "$(EXCLU)" --exclude ".svn"
	@echo "TGZ file is ../$(TGZFILE)"

clean:
	cd src-flo && make clean
	cd src-adj && make clean
	cd src-grd && make clean
	cd src-opt && make clean
	cd src-utl && make clean

cleanall:
	cd src-flo && make cleanall
	cd src-adj && make cleanall
	cd src-grd && make cleanall
	cd src-opt && make cleanall
	cd src-utl && make cleanall
