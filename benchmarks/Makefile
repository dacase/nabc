include ../config.h

SUFFIXES:
SUFFIXES: .c .o

LIBS = -L$(LIBDIR) $(FLIBS) $(NETCDFLIBF)  $(LM) 

#---------------------------------------------------------------------------
bench_md:  bench_md.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS)

bench_md.MPI:  bench_md.MPI.o
	@echo "[NABC] CC $@"
	$(VB)mpicc -o $@ $<  $(LIBS)

.c.o:
	@echo "[NABC] CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) -I$(INCDIR) -o $@ $<

bench_md.MPI.o: bench_md.c
	@echo "[NABC] CC $<"
	$(VB)mpicc -DMPI -c $(COPTFLAGS) $(CFLAGS) -I$(INCDIR) -o $@ $<

