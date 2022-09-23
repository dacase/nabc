#   top-level Makefile for nabc

include ../config.h

install:
	@cd cifparse && $(MAKE) $(CIFTARGET)
	@cd sff && $(MAKE) libsff
	@cd nabc && $(MAKE) libnabc
	@cd arpack && $(MAKE) libarpack

clean:
	-(cd cifparse && $(MAKE) clean )
	-(cd sff && $(MAKE) clean )
	-(cd nabc && $(MAKE) clean )
	-(cd arpack && $(MAKE) clean )

uninstall: clean
	/bin/rm -f $(INCDIR)/* ../config.h

distclean: clean
	-(cd .. && /bin/rm -rf $(BINDIR) $(LIBDIR) $(INCDIR) logs config.h)

