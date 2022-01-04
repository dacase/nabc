#   top-level Makefile for nab

include ../config.h

install: $(INSTALLTYPE)

serial:
	@cd ucpp-1.3 && $(MAKE) ucpp
	@cd sff && $(MAKE) libsff
	@cd nabc && $(MAKE) libnabc
	@cd nab && $(MAKE) install

parallel:

clean:
	-(cd ucpp-1.3 && $(MAKE) clean )
	-(cd sff && $(MAKE) clean )
	-(cd nabc && $(MAKE) clean )
	-(cd nab && $(MAKE) clean )

uninstall: clean
	/bin/rm -f $(INCDIR)/* ../config.h

distclean: clean
	-(cd .. && /bin/rm -rf $(BINDIR) $(LIBDIR) $(INCDIR) logs config.h)

