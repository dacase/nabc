test:
	@./run_tests.sh

test.all: test.tleap test.antechamber

test.tleap:
	-(cd tleap/OL3 && ./Run.OL3.sh)
	-(cd tleap/BondAugmentations && ./Run.1p7e)

test.antechamber:
	-(cd antechamber/sustiva && Run.sustiva)
	-(cd antechamber/sustiva && Run.sustiva_g2)

clean::
	/bin/rm -f TEST_FAILURES.diff
