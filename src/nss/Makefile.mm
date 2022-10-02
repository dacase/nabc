###
#
# Makefile for mm
#
###

AVS_LIBS=	$(AVS_PATH)/lib
AVS_INC=	$(AVS_PATH)/include
NAB_INC=	$(NABHOME)/include
NAB_LIBS=	$(NABHOME)/lib/$(ARCH)
SYM_LIBS=	-L$(SYM_HOME) -lsym -L$(NAB_LIBS) -lnab
BASELIBS=	-lgeom -lutil -lm $(LASTLIBS)
CFLAGS=		-I$(AVS_INC) -I$(NAB_INC) $(AOPTCFLAGS) $(LOCAL_CFLAGS) $(G)
CFLOWLIBS=	-L$(AVS_LIBS) -lflow_c $(BASELIBS)

mm:	mm.o
	$(CC) $(CFLAGS) -o mm mm.o $(SYM_LIBS) $(CFLOWLIBS)
