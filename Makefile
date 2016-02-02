all: makefiles
	cd src && $(MAKE)

clean: 
	cd src && $(MAKE) cleanall
	rm -f src/Makefile

makefiles:
	@if [ ! -f src/Makefile ]; then \
	cd src && opp_makemake -O out -o koisim;\
	fi
