#List of files/directories to make tgz file
FILES=nuwtun/src-flo  \
		nuwtun/docs     \
		nuwtun/makefile \
		nuwtun/README

#Things to exclude in the above list
EXCLU=nuwtun/docs/tree

#Name of tgz file
TGZFILE=nuwtun.tgz

all: flo tgz

flo:
	cd src-flo && make

tgz:
	cd src-flo && make cleanall
	cd .. && tar zcvf $(TGZFILE) $(FILES) --exclude "$(EXCLU)"
	@echo "TGZ file is ../$(TGZFILE)"
