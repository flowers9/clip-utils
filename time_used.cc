#include "time_used.h"
#include <limits.h>	/* ULONG_MAX */
#include <sys/times.h>	/* struct tms, times() */
#include <sys/types.h>	/* clock_t */
#include <unistd.h>	/* _SC_CLK_TCK, sysconf() */

static struct tms start;
static clock_t start_real;

/* start a timing loop */

void start_time() {
	start_real = times(&start);
}

/* return user time since start_time(), in seconds */

double used_time() {
	struct tms stop;
	times(&stop);
	return (double)(stop.tms_utime - start.tms_utime) / sysconf(_SC_CLK_TCK);
}

/* return real time since start_time(), in seconds */

double elapsed_time() {
	struct tms stop;
	return static_cast<double>(times(&stop) - start_real) / sysconf(_SC_CLK_TCK);
}

/* find out how fast machine currently is - return loops per 0.1 second */

double timing_loop() {
	unsigned long i;
	start_time();
	for (i = 0; i != ULONG_MAX && elapsed_time() < 0.1; ++i) { }
	return i / elapsed_time();
}
