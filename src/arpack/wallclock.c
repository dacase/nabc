#include <sys/time.h>
#include <unistd.h>

void arsecond_( double *wallc ){

#  ifdef NO_DETAILED_TIMINGS
	*wallc = 0.0;
#  else
	struct timeval tv;

	long rsec;

	gettimeofday( &tv, NULL);

	rsec = tv.tv_sec - 959800000;
	*wallc = (double)rsec + (double)(tv.tv_usec)/1000000.;
#  endif

	return;

}

