NABCHOME=/home/case/nabc
NABFLAGS=""

SUBS= shifts.o shifts_emp.o ring_current.o  get_banis_info.o \
		plane.o read_obs_shifts.o ax_susc.o\
		sbcoil.o swap_shifts.o get_pep_info.o \
		get_sugar_info.o get_chi_info.o el.o \
		get_dist_ang.o write_headers.o get_vdwa.o \
		shifts_qdb.o

.SUFFIXES:
.SUFFIXES: .nab .o

.nab.o: $<
	$(NABCHOME)/bin/nab -c $(NABFLAGS) $<
	@/bin/rm -f $*.c

install:  shifts

shifts: $(SUBS)
	$(NABCHOME)/bin/nab $(NABFLAGS) -o ../bin/shifts $(SUBS) 

shifts_qdb.o:		shifts_qdb.h
shifts_emp.o:       shifts_emp.h constants.h map.h
ring_current.o:     ringinfo.h get_ring_info.h
read_obs_shifts.o:  omap.h
sbcoil.o:           omap.h

clean:
	/bin/rm -f *.o

