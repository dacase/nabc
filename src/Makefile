#   top-level Makefile for nab

include ../config.h

install: $(INSTALLTYPE)

serial:
	@cd leap && $(MAKE) tleap
	@cd antechamber && $(MAKE) ac
	@cd sff && $(MAKE) libsff
	@cd nabc && $(MAKE) libnabc
	@cd nab && $(MAKE) install

parallel:

clean:
	-(cd leap && $(MAKE) clean )
	-(cd antechamber && $(MAKE) clean )
	-(cd sff && $(MAKE) clean )
	-(cd nabc && $(MAKE) clean )

uninstall: clean
	-(cd $(BINDIR) && /bin/rm -f teLeap tleap)
	/bin/rm -f $(INCDIR)/* ../config.h

distclean: clean
	-(cd .. && /bin/rm -rf $(BINDIR) $(LIBDIR) $(INCDIR) logs config.h)

