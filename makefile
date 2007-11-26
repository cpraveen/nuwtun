#List of files/directories to make tgz file
FILES=nuwtun/src-flo     \
		nuwtun/src-adj     \
		nuwtun/docs        \
		nuwtun/makefile    \
		nuwtun/makefile.in \
		nuwtun/README

#Things to exclude in the above list
EXCLU=nuwtun/docs/tree

#Name of tgz file
TGZFILE=nuwtun.tgz

all: flo

flo:
	@cd src-flo && make

tgz:
	cd src-flo && make cleanall
	cd src-adj && make clean
	cd .. && tar zcvf $(TGZFILE) $(FILES) --exclude "$(EXCLU)" --exclude ".svn"
	@echo "TGZ file is ../$(TGZFILE)"
