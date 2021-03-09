# Makefile for antechamber

include ../../config.h

ACDOCTOROBJS = acdoctor.o eprintf.o
AM1BCCOBJS = am1bcc.o eprintf.o
ANTECHAMBEROBJS = antechamber.o eprintf.o
ATOMTYPEOBJS = atomtype.o eprintf.o
BONDTYPEOBJS = bondtype.o eprintf.o ncsu-penalties.o
CHARMMGENOBJS = charmmgen.o eprintf.o
CRDGROWOBJS = crdgrow.o eprintf.o
DATABASEOBJS = database.o eprintf.o
DELPHIGENOBJS = delphigen.o
ESPGENOBJS = espgen.o eprintf.o
MATCHOBJS = eprintf.o match.o
MATCH_ATOMNAMEOBJS = match_atomname.o eprintf.o
PARMCALOBJS = eprintf.o parmcal.o
PARMCHK2OBJS = parmchk2.o eprintf.o
PREPGENOBJS = eprintf.o prepgen.o
RESIDUEGENOBJS = residuegen.o eprintf.o
RESPGENOBJS = respgen.o eprintf.o
TRANSLATEOBJS = translate.o eprintf.o
TOP2FFOBJS = top2ff.o eprintf.o
TOP2MOL2OBJS = top2mol2.o eprintf.o

# Dispatched programs are wrapper scripts of programs; these scripts
# automatically set MSANDERHOME; the programs themselves live in 
# $(BINDIR)/x/
DISPATCHEDPROGS = \
	$(BINDIR)/x/am1bcc$(SFX) \
	$(BINDIR)/x/antechamber$(SFX) \
	$(BINDIR)/x/atomtype$(SFX) \
	$(BINDIR)/x/bondtype$(SFX) \
	$(BINDIR)/x/espgen$(SFX) \
	$(BINDIR)/x/parmchk2$(SFX) \
	$(BINDIR)/x/prepgen$(SFX) \
	$(BINDIR)/x/residuegen$(SFX) \
	$(BINDIR)/x/respgen$(SFX) \

LOCALPROGS = \
	am1bcc$(SFX) \
	antechamber$(SFX) \
	atomtype$(SFX) \
	bondtype$(SFX) \
	espgen$(SFX) \
	parmchk2$(SFX) \
	prepgen$(SFX) \
	residuegen$(SFX) \
	respgen$(SFX) \

PROGS = \
	$(BINDIR)/am1bcc$(SFX) \
	$(BINDIR)/antechamber$(SFX) \
	$(BINDIR)/atomtype$(SFX) \
	$(BINDIR)/bondtype$(SFX) \
	$(BINDIR)/espgen$(SFX) \
	$(BINDIR)/parmchk2$(SFX) \
	$(BINDIR)/prepgen$(SFX) \
	$(BINDIR)/residuegen$(SFX) \
	$(BINDIR)/respgen$(SFX) \

TESTS = test-penalties test-int-part test-combination

#############################################################################

all: $(PROGS)

dispatch:
	SFX=$(SFX) BINDIR=$(BINDIR) sh patch.sh

tests: $(TESTS)

.PHONY: ac
ac: $(DISPATCHEDPROGS) dispatch

clean:
	/bin/rm -f *.o $(LOCALPROGS) $(TESTS)

uninstall:
	/bin/rm -f $(PROGS) $(DISPATCHEDPROGS)

#############################################################################

$(LIBDIR)/libcifparse.a:
	cd ../cifparse && $(MAKE) install

.c.o:
	@echo "[ANTECHAMBER]  CC $<"
	$(VB)$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

#_____________________________________________________________________________
# compilation rules.

$(BINDIR)/x/acdoctor$(SFX): $(ACDOCTOROBJS) $(LIBDIR)/libcifparse.a
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) $(ACDOCTOROBJS) \
	    -L$(LIBDIR) -lcifparse $(LM) -o $@

$(BINDIR)/x/am1bcc$(SFX): $(AM1BCCOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) $(AM1BCCOBJS) \
	    $(LM) -o $@

$(BINDIR)/x/antechamber$(SFX): $(ANTECHAMBEROBJS) $(LIBDIR)/libcifparse.a
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) $(ANTECHAMBEROBJS) \
	    -L$(LIBDIR) -lcifparse $(LM) -o $@

$(BINDIR)/x/atomtype$(SFX): $(ATOMTYPEOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(ATOMTYPEOBJS) $(LM) -o $@

$(BINDIR)/x/bondtype$(SFX): $(BONDTYPEOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(BONDTYPEOBJS) $(LM) -o $@

$(BINDIR)/x/charmmgen$(SFX): $(CHARMMGENOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(CHARMMGENOBJS) $(LM) -o $@

$(BINDIR)/x/crdgrow$(SFX): $(CRDGROWOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(CRDGROWOBJS) $(LM) -o $@

$(BINDIR)/x/database$(SFX): $(DATABASEOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(DATABASEOBJS) $(LM) -o $@

$(BINDIR)/x/delphigen$(SFX): $(DELPHIGENOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(DELPHIGENOBJS) -o $@

$(BINDIR)/x/espgen$(SFX): $(ESPGENOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(ESPGENOBJS) $(LM) -o $@

$(BINDIR)/x/match$(SFX): $(MATCHOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(MATCHOBJS) $(LM) -o $@

$(BINDIR)/x/match_atomname$(SFX): $(MATCH_ATOMNAMEOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(MATCH_ATOMNAMEOBJS) $(LM) -o $@

$(BINDIR)/x/parmcal$(SFX): $(PARMCALOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(PARMCALOBJS) $(LM) -o $@

$(BINDIR)/x/parmchk2$(SFX): $(PARMCHK2OBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(PARMCHK2OBJS) $(LM) -o $@

$(BINDIR)/x/parmjoin$(SFX): parmjoin.o
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) parmjoin.o -o $@

$(BINDIR)/x/prepgen$(SFX): $(PREPGENOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(PREPGENOBJS) $(LM) -o $@

$(BINDIR)/x/residuegen$(SFX): $(RESIDUEGENOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(RESIDUEGENOBJS) $(LM) -o $@

$(BINDIR)/x/respgen$(SFX): $(RESPGENOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(RESPGENOBJS) $(LM) -o $@

$(BINDIR)/x/translate$(SFX): $(TRANSLATEOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(TRANSLATEOBJS) $(LM) -o $@

$(BINDIR)/x/top2ff$(SFX): $(TOP2FFOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(TOP2FFOBJS) $(LM) -o $@

$(BINDIR)/x/top2mol2$(SFX): $(TOP2MOL2OBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
	    $(TOP2MOL2OBJS) $(LM) -o $@

#_____________________________________________________________________________
# Programs

$(BINDIR)/acdoctor$(SFX): acdoctor$(SFX)
	$(VB)$(MV) acdoctor$(SFX) $@

$(BINDIR)/am1bcc$(SFX): am1bcc$(SFX)
	$(VB)$(MV) am1bcc$(SFX) $@

$(BINDIR)/antechamber$(SFX): antechamber$(SFX)
	$(VB)$(MV) antechamber$(SFX) $@

$(BINDIR)/atomtype$(SFX): atomtype$(SFX)
	$(VB)$(MV) atomtype$(SFX) $@

$(BINDIR)/bondtype$(SFX): bondtype$(SFX)
	$(VB)$(MV) bondtype$(SFX) $@

$(BINDIR)/charmmgen$(SFX): charmmgen$(SFX)
	$(VB)$(MV) charmmgen$(SFX) $@

$(BINDIR)/crdgrow$(SFX): crdgrow$(SFX)
	$(VB)$(MV) crdgrow$(SFX) $@

$(BINDIR)/database$(SFX): database$(SFX)
	$(VB)$(MV) database$(SFX) $@

$(BINDIR)/delphigen$(SFX): delphigen$(SFX)
	$(VB)$(MV) delphigen$(SFX) $@

$(BINDIR)/espgen$(SFX): espgen$(SFX)
	$(VB)$(MV) espgen$(SFX) $@

$(BINDIR)/match$(SFX): match$(SFX)
	$(VB)$(MV) match$(SFX) $@

$(BINDIR)/match_atomname$(SFX): match_atomname$(SFX)
	$(VB)$(MV) match_atomname$(SFX) $@

$(BINDIR)/parmcal$(SFX): parmcal$(SFX)
	$(VB)$(MV) parmcal$(SFX) $@

$(BINDIR)/parmchk2$(SFX): parmchk2$(SFX)
	$(VB)$(MV) parmchk2$(SFX) $@

$(BINDIR)/parmjoin$(SFX): parmjoin$(SFX)
	$(VB)$(MV) parmjoin$(SFX) $@

$(BINDIR)/prepgen$(SFX): prepgen$(SFX)
	$(VB)$(MV) prepgen$(SFX) $@

$(BINDIR)/residuegen$(SFX): residuegen$(SFX)
	$(VB)$(MV) residuegen$(SFX) $@

$(BINDIR)/respgen$(SFX): respgen$(SFX)
	$(VB)$(MV) respgen$(SFX) $@

$(BINDIR)/translate$(SFX): translate$(SFX)
	$(VB)$(MV) translate$(SFX) $@

$(BINDIR)/top2ff$(SFX): top2ff$(SFX)
	$(VB)$(MV) top2ff$(SFX) $@

$(BINDIR)/top2mol2$(SFX): top2mol2$(SFX)
	$(VB)$(MV) top2mol2$(SFX) $@

#_____________________________________________________________________________
# Tests

test:
	cd ../../test/antechamber; make clean; make; make check

test-penalties: ncsu-penalties.o test-penalties.o
	$(CC) $(CFLAGS) $(AMBERCFLAGS) -o test-penalties $(LDFLAGS) $(AMBERLDFLAGS) ncsu-penalties.o test-penalties.o

test-int-part: ncsu-penalties.h ncsu-penalties.c
	$(CC)  $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o test-int-part ncsu-penalties.c $(LDFLAGS) $(AMBERLDFLAGS) \
   -DENABLE_INT_PART_MAIN=1

test-combination: ncsu-penalties.h ncsu-penalties.c
	$(CC) $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o test-combination ncsu-penalties.c $(LDFLAGS) $(AMBERLDFLAGS) \
   -DENABLE_COMBINATION_MAIN=1

#############################################################################
# DO NOT DELETE

#_____________________________________________________________________________
# Dependencies, from ./makedepend:

acdoctor.o: common.h
acdoctor.o: define.h
acdoctor.o: atom.h
acdoctor.o: utility.c
acdoctor.o: common.c
acdoctor.o: rotate.c
acdoctor.o: ac.c
acdoctor.o: charmm.c
acdoctor.o: mol2.c
acdoctor.o: mmcif.c
acdoctor.o: mopcrt.c
acdoctor.o: divcrt.c
acdoctor.o: mopint.c
acdoctor.o: mopout.c
acdoctor.o: divout.c
acdoctor.o: sqmcrt.c
acdoctor.o: sqmout.c
acdoctor.o: gesp.c
acdoctor.o: gcrt.c
acdoctor.o: gzmat.c
acdoctor.o: gout.c
acdoctor.o: gamess.c
acdoctor.o: pdb.c
acdoctor.o: csd.c
acdoctor.o: mdl.c
acdoctor.o: alc.c
acdoctor.o: hin.c
acdoctor.o: prep.c
acdoctor.o: rst.c
acdoctor.o: jzmat.c
acdoctor.o: jcrt.c
acdoctor.o: jout.c
acdoctor.o: checkmolecule.c
acdoctor.o: fileformat.c
am1bcc.o: common.h
am1bcc.o: define.h
am1bcc.o: atom.h
am1bcc.o: utility.c
am1bcc.o: common.c
am1bcc.o: rotate.c
am1bcc.o: ac.c
am1bcc.o: pdb.c
antechamber.o: define.h
antechamber.o: atom.h
antechamber.o: eprintf.h
antechamber.o: utility.c
antechamber.o: common.c
antechamber.o: equatom.c
antechamber.o: ac.c
antechamber.o: charmm.c
antechamber.o: mol2.c
antechamber.o: mmcif.c
antechamber.o: mopcrt.c
antechamber.o: divcrt.c
antechamber.o: sqmcrt.c
antechamber.o: sqmout.c
antechamber.o: mopint.c
antechamber.o: mopout.c
antechamber.o: divout.c
antechamber.o: gcrt.c
antechamber.o: gzmat.c
antechamber.o: gout.c
antechamber.o: gamess.c
antechamber.o: rotate.c
antechamber.o: pdb.c
antechamber.o: csd.c
antechamber.o: mdl.c
antechamber.o: alc.c
antechamber.o: hin.c
antechamber.o: prep.c
antechamber.o: rst.c
antechamber.o: jzmat.c
antechamber.o: jcrt.c
antechamber.o: jout.c
antechamber.o: gesp.c
antechamber.o: charge.c
antechamber.o: fileformat.c
antechamber.o: checkmolecule.c
atomtype.o: common.h
atomtype.o: define.h
atomtype.o: atom.h
atomtype.o: utility.c
atomtype.o: common.c
atomtype.o: ring.c
atomtype.o: rotate.c
atomtype.o: ac.c
atomtype.o: mol2.c
bondtype.o: ncsu-penalties.h
bondtype.o: common.h
bondtype.o: define.h
bondtype.o: atom.h
bondtype.o: utility.c
bondtype.o: common.c
bondtype.o: ring.c
bondtype.o: ac.c
bondtype.o: mol2.c
charmmgen.o: common.h
charmmgen.o: define.h
charmmgen.o: atom.h
charmmgen.o: utility.c
charmmgen.o: common.c
charmmgen.o: rotate.c
charmmgen.o: ac.c
charmmgen.o: mol2.c
charmmgen.o: prep.c
checkmolecule.o: define.h
checkmolecule.o: atom.h
common.o: define.h
crdgrow.o: common.h
crdgrow.o: define.h
crdgrow.o: atom.h
crdgrow.o: utility.c
crdgrow.o: rotate.c
crdgrow.o: common.c
crdgrow.o: pdb.c
crdgrow.o: ac.c
crdgrow.o: prep.c
database.o: define.h
database.o: utility.c
delphigen.o: common.h
delphigen.o: define.h
delphigen.o: atom.h
delphigen.o: common.c
delphigen.o: ac.c
delphigen.o: pdb.c
eprintf.o: eprintf.h
equatom.o: rotate.c
equatom2.o: rotate.c
errormsg.o: errormsg.h
espgen.o: common.h
espgen.o: define.h
espgen.o: atom.h
fileformat.o: define.h
fileformat.o: atom.h
fileformat.o: checkmolecule.c
match.o: common.h
match.o: define.h
match.o: atom.h
match.o: utility.c
match.o: common.c
match.o: ring.c
match.o: ac.c
match.o: pdb.c
match.o: mol2.c
match.o: mdl.c
match.o: rst.c
match_atomname.o: common.h
match_atomname.o: define.h
match_atomname.o: atom.h
match_atomname.o: utility.c
match_atomname.o: common.c
match_atomname.o: equatom2.c
match_atomname.o: pdb.c
match_atomname.o: ac.c
match_atomname.o: prep.c
match_atomname.o: mol2.c
mmcif.o: ../cifparse/cifparse.h
ncsu-penalties.o: ncsu-penalties.h
parmcal.o: common.h
parmcal.o: define.h
parmcal.o: atom.h
parmcal.o: utility.c
parmcal.o: common.c
parmchk2.o: common.h
parmchk2.o: define.h
parmchk2.o: atom.h
parmchk2.o: utility.c
parmchk2.o: common.c
parmchk2.o: rotate.c
parmchk2.o: ac.c
parmchk2.o: mol2.c
parmchk2.o: prep.c
parmjoin.o: common.h
parmjoin.o: define.h
parmjoin.o: utility.c
prepgen.o: common.h
prepgen.o: define.h
prepgen.o: atom.h
prepgen.o: utility.c
prepgen.o: common.c
prepgen.o: rotate.c
prepgen.o: ring.c
prepgen.o: ac.c
residuegen.o: common.h
residuegen.o: define.h
residuegen.o: atom.h
residuegen.o: utility.c
residuegen.o: common.c
residuegen.o: ac.c
respgen.o: common.h
respgen.o: define.h
respgen.o: atom.h
respgen.o: utility.c
respgen.o: common.c
respgen.o: ac.c
respgen.o: equatom.c
rotate.o: common.h
smiles.o: define.h
smiles.o: eprintf.h
test-penalties.o: ncsu-penalties.h
translate.o: common.h
translate.o: define.h
translate.o: atom.h
translate.o: utility.c
translate.o: common.c
translate.o: rotate.c
translate.o: ac.c
translate.o: pdb.c
translate.o: mol2.c
translate.o: prep.c
translate.o: lsfit.c
utility.o: define.h
