#   top-level Makefile for nabc

include ../config.h

install:
	@cd cifparse && $(MAKE) $(CIFTARGET)
	@cd sff && $(MAKE) libsff
	@cd nabc && $(MAKE) libnabc
	@cd arpack && $(MAKE) libarpack

nab: install
	@cd ucpp-1.3 && make ucpp
	@cd nab && make
	@cd nss & make

clean:
	-(cd cifparse && $(MAKE) clean )
	-(cd sff && $(MAKE) clean )
	-(cd nabc && $(MAKE) clean )
	-(cd arpack && $(MAKE) clean )
	-(cd ucpp-1.3 && $(MAKE) clean )
	-(cd nab && $(MAKE) clean )
	-(cd nss && $(MAKE) clean )

uninstall: clean
	/bin/rm -f $(INCDIR)/* ../config.h

distclean: clean
	-(cd .. && /bin/rm -rf $(BINDIR) $(LIBDIR) $(INCDIR) logs config.h)

