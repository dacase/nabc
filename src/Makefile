#   top-level Makefile for nabc

include ../config.h

install:
	@cd cifparse && $(MAKE) libcifparse
	@cd sff && $(MAKE) libsff
	@cd nabc && $(MAKE) libnabc
	@cd arpack && $(MAKE) libarpack
	@cd lapack && $(MAKE) liblapack
	@cd blas && $(MAKE) libblas

nab: install
	@cd ucpp-1.3 && make ucpp
	@cd nab && make install
	@cd nss & make install

clean:
	-(cd cifparse && $(MAKE) clean )
	-(cd sff && $(MAKE) clean )
	-(cd nabc && $(MAKE) clean )
	-(cd arpack && $(MAKE) clean )
	-(cd lapack && $(MAKE) clean )
	-(cd blas && $(MAKE) clean )
	-(cd ucpp-1.3 && $(MAKE) clean )
	-(cd nab && $(MAKE) clean )
	-(cd nss && $(MAKE) clean )
	-(cd nss && $(MAKE) clean )
	-(cd netcdf-c-4.7.2 && $(MAKE) clean )
	-(cd netcdf-fortran-4.5.2 && $(MAKE) clean )

uninstall: clean
	/bin/rm -f $(INCDIR)/* ../config.h

distclean: clean
	-(cd netcdf-c-4.7.2 && $(MAKE) distclean )
	-(cd netcdf-fortran-4.5.2 && $(MAKE) distclean )
	-(cd .. && /bin/rm -rf $(BINDIR) $(LIBDIR) $(INCDIR) logs share config.h)

