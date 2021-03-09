include ../../config.h

SUFFIXES:
SUFFIXES: .c .o

test:  test_mask test_xmin test_md test_lmod test_nmode

test_xmin:  txmin
	time ./txmin xfin.parm7 xfin.rst7 foo.rst7
	@/bin/rm -f foo.rst7 txmin txmin.o

test_lmod:  tlmod
	./tlmod trpcage.parm7 trpcage.rst7
	@/bin/rm -f conflib.dat tlmod tlmod.o

test_md:  tmd
	time ./tmd xfin.parm7 xfin.md1.x
	@/bin/rm -f xfin.md2.x

test_md.MPI:  tmd.MPI
	time mpirun -np 4 ./tmd.MPI xfin.parm7 xfin.md1.x

test_nmode:  tnmode
	./tnmode crambin.parm7 crambin.rst7 > nmode.out
	@/bin/rm -f vecs

test_mask:  tmask
	./tmask xfin.parm7 xfin.md1.x ":1-3"

test_hcp1:  thcp1
	./thcp1

test_hcp2:  thcp2
	./thcp2

test_nabc:  tnabc
	./tnabc 2.pdb 2.rst7
	@/bin/rm -f tnabc tnabc.o

test_netcdf:  tnetcdf
	./tnetcdf

test_cif:  tcif
	./tcif
	/bin/rm tcif tcif.o

LIBS = -L$(LIBDIR) -lsff $(FLIBSF) $(NETCDFLIBF)  $(LM) -lgfortran -w

#---------------------------------------------------------------------------
txmin:  txmin.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o txmin txmin.o $(LIBS)

tlmod:  tlmod.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o tlmod tlmod.o $(LIBS) \
                $(BASEDIR)/src/msander/dsarpack.o -lgfortran -w

tmd:  tmd.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS)

tmd.MPI:  tmd.MPI.o
	@echo "[NABC] CC $@"
	$(VB)mpicc -o $@ $<  $(LIBS)

tnmode:  tnmode.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS)

tmask:  tmask.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS)

thcp1:  thcp1.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS)

thcp2:  thcp2.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS)

tnabc:  tnabc.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  -L$(LIBDIR)  -lnabc -lcifparse $(LIBS)

tcif:  tcif.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  -L$(LIBDIR)  -lnabc -lcifparse $(LIBS)

tnetcdf:  tnetcdf.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS)

#---------------------------------------------------------------------------
#  non-tests:

fd:  fd.c
	@echo "[NABC] CC $<"
	$(VB)$(CC) $(CFLAGS) -I$(INCDIR) -o $@ $< -L$(LIBDIR) \
        -lnabc -lcifparse  $(LIBS)

wch:  wch.c
	@echo "[NABC] CC $<"
	$(VB)$(CC) $(CFLAGS) -I$(INCDIR) -o $@ $< -L$(LIBDIR) \
        -lnabc -lcifparse $(LIBS)

circle:  circle.c
	@echo "[NABC] CC $<"
	$(VB)$(CC) $(CFLAGS) -I$(INCDIR) -o $@ $< -L$(LIBDIR) \
        -lnabc -lcifparse $(LIBS)

loop:  loop.c
	@echo "[NABC] CC $<"
	$(VB)$(CC) $(CFLAGS) -I$(INCDIR) -o $@ $< -L$(LIBDIR) \
        -lnabc -lcifparse $(LIBS)

#---------------------------------------------------------------------------
.c.o:
	@echo "[NABC] CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) -I$(INCDIR) -o $@ $<

tmd.MPI.o: tmd.c
	@echo "[NABC] CC $<"
	$(VB)mpicc -DMPI -c $(COPTFLAGS) $(CFLAGS) -I$(INCDIR) -o $@ $<

clean:
	/bin/rm -f *.o txmin tmd tmd.MPI tnabc tlmod tnmode tmask thcp1 thcp2 \
       trpcageTest.nc
