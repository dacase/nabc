#  just redirect things to lower-level Makefiles

install::
	cd src && make install

nab::
	cd src && make nab

shifts::
	cd src/shifts/src && make install

test::
	cd test && make test

clean::
	cd src && make clean

uninstall::
	cd src && make uninstall

distclean::
	cd src && make distclean
