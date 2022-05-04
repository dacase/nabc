include ../config.h

SUFFIXES:
SUFFIXES: .c .o

LIBS = -L$(LIBDIR) $(FLIBS) $(NETCDFLIBF)  $(LM) 

#---------------------------------------------------------------------------
bench_md:  bench_md.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS)

.c.o:
	@echo "[NABC] CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) -I$(INCDIR) -o $@ $<

