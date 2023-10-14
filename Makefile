#  just redirect things to lower-level Makefiles

install::
	cd src && make install

nab::
	cd src && make nab

test::
	cd test && make test

clean::
	cd src && make clean

uninstall::
	cd src && make uninstall

distclean::
	cd src && make distclean
