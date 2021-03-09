include ../../config.h

.SUFFIXES:
.SUFFIXES:  .nab .c .o .f

%.o: %.nab $(BINDIR)/$(MPI)nab$(SFX) $(BINDIR)/$(MPI)nab2c$(SFX)
	@echo "[NAB]  NAB $<"
	$(VB)$(WINE) $(BINDIR)/$(MPI)nab$(SFX) -c $(NABFLAGS) $< -o $@
	@rm -f $*.c

NAB2COBJS=	\
	nab2c.o	\
	cgen.o	\
	checkexpr.o	\
	dumpnode.o\
	errormsg.o	\
	fixexpr.o	\
	node.o	\
	symbol.o	\
	y.tab.o

LIBNABOBJS=	\
	gsub.o	\
	hashutil.o	\
	traceback.o \
    time.o \
	tm_malloc.o	


#===========================================================================

install: $(BINDIR)/$(MPI)nab$(SFX) $(BINDIR)/$(MPI)nab2c$(SFX)
	@rm -f $(INCDIR)/nabcode.h $(INCDIR)/defreal.h $(INCDIR)/nab.h
	@cp nabcode.h $(INCDIR)
	@cp nabtypes.h $(INCDIR)
	@cp defreal.h $(INCDIR)
	@cp nab.h $(INCDIR)
#	@cp nab_netcdf.h $(INCDIR)
#	@cp xmin_opt.h lmod_opt.h $(INCDIR)
	@mv $(BINDIR)/$(MPI)nab$(SFX) $(BINDIR)/x/
	@cat nab_template | sed "s/REPLACE_ME/$(MPI)nab$(SFX)/g" > $(BINDIR)/$(MPI)nab$(SFX)
	@chmod +x $(BINDIR)/$(MPI)nab$(SFX)

#  Force rebuilding of nab, since there is no easy to to see if FLIBS has
#  changed.
$(BINDIR)/$(MPI)nab$(SFX)::
	@echo "[NAB]  CC $@"
	$(VB)$(CC) -DCC='"$(CC)"' -DCPP='"$(CPP)"' -DFLIBS='"$(FLIBS)"' \
		$(CNOOPTFLAGS) $(CFLAGS) \
		$(LDFLAGS) -o $(BINDIR)/$(MPI)nab$(SFX) nab.c

$(BINDIR)/$(MPI)nab2c$(SFX):	$(NAB2COBJS)
	@echo "[NAB]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
		-o $(BINDIR)/$(MPI)nab2c$(SFX) $(NAB2COBJS)

#===========================================================================

y.tab.h: y.tab.c

y.tab.c:  nabgrm.y lex.yy.c
	@echo "[NAB]  BISON $<"
	$(VB)$(YACC) -d nabgrm.y

y.tab.o: y.tab.c
	@echo "[NAB]  CC $<"
	$(VB)$(CC) -c -D$(LEX) $(CFLAGS) $(AMBERCFLAGS) y.tab.c

lex.yy.c: nablex.l
	@echo "[NAB]  FLEX $<"
	$(VB)$(LEX) nablex.l

lex.dgo.c: dg_options.l
	@echo "[NAB]  FLEX $<"
	$(VB)$(LEX) dg_options.l

checkexpr.c: y.tab.o
	(cd semantics; $(MAKE))
	mv semantics/checkexpr.c .

.c.o:
	@echo "[NAB]  CC $<"
	$(VB)$(CC) -c $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

#===========================================================================

clean:
	-rm -f *.o nab$(SFX) nab2c$(SFX) checkexpr.c
	-rm -f y.tab.c y.tab.h lex.dgo.c lex.yy.c
	(cd semantics && $(MAKE) clean )

uninstall: clean
	-rm -f \
	$(BINDIR)/to_be_dispatched/nab$(SFX) \
	$(BINDIR)/to_be_dispatched/mpinab$(SFX) \
	$(BINDIR)/nab$(SFX) \
        $(BINDIR)/nab2c$(SFX) \
	$(BINDIR)/mpinab$(SFX) \
        $(BINDIR)/mpinab2c$(SFX) \
	$(INCDIR)/nabcode.h \
	$(INCDIR)/nabtypes.h \
	$(INCDIR)/defreal.h \
        $(INCDIR)/nab.h \
	$(INCDIR)/xmin_opt.h \
	$(INCDIR)/lmod_opt.h \

#===========================================================================

depend: $(NAB2COBJS:.o=.c)
	@echo "[NAB2C] make depend"
	$(VB)$(CC) -MM $^ > $@

include depend

