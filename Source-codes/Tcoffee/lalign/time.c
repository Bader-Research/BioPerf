#include <stdio.h>
#ifndef MSDOS
#define TIMES
#ifndef VMS
#include <sys/types.h>
#include <sys/times.h>
#else
struct tms { long tms_utime, tms_stime, tms_cutime, tms_cstime; }
#endif /* VMS */
#else /* MSDOS */
#undef TIMES
#endif

#ifndef HZ
#define HZ 100
#endif

#ifndef linux 
long stime()		/* gets the time as an int */
{
#ifndef TIMES
	long time(), tt;
	return time(&tt)*HZ;
#else
	struct tms tt;
	times(&tt);
	return (long)tt.tms_utime;
#endif
	}
#endif
ptime(fd,time)
	FILE *fd; long time;
{
	long dtime, mtime, htime, stime;

	stime = HZ;
	mtime = stime*60;
	htime = mtime*60;
	dtime = htime*24;

	if (time < 0) time = time + dtime;
	fprintf(fd,"%2ld:%02ld:%02ld",
		time/htime,(time%htime)/mtime,(time%mtime)/stime);
	}

