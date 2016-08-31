#include "timer.h"
#include<iostream>
#include <unistd.h>

#ifndef WINDOWS
#include <sys/time.h>
#include <sys/resource.h>
#endif


// extern "C" int getrusage(__rusage_who, struct rusage *rusage);



/*
 *  The virtual time of day and the real time of day are calculated and
 *  stored for future use.  The future use consists of subtracting these
 *  values from similar values obtained at a later time to allow the user
 *  to get the amount of time used by the backtracking routine.
 */

#ifdef WINDOWS
void start_timers() {};
void stop_timers(Timer type) {};

#else
static struct rusage res;
static struct timeval tp;
static Time virtual_utime, virtual_stime;
Time virtual_ulapse, virtual_slapse;
static Time real_time;
Time real_lapse;


void start_timers()
{
	getrusage(RUSAGE_SELF, &res);
	virtual_utime = (Time) res.ru_utime.tv_sec +
	                (Time) res.ru_utime.tv_usec / 1000000.0;
	virtual_stime = (Time) res.ru_stime.tv_sec +
	                (Time) res.ru_stime.tv_usec / 1000000.0;

	gettimeofday(&tp, NULL);
	real_time = (Time) tp.tv_sec +
	            (Time) tp.tv_usec / 1000000.0;
}


/*
 *  Stop the stopwatch and return the time used in seconds (either
 *  REAL or VIRTUAL time, depending on ``type'').
 */
void stop_timers(Timer type)
{
	if (type == REAL) {
		gettimeofday(&tp, NULL);
		real_lapse = (Time) tp.tv_sec +
		             (Time) tp.tv_usec / 1000000.0
		             - real_time;
	} else {
		getrusage(RUSAGE_SELF, &res);
		virtual_ulapse = (Time) res.ru_utime.tv_sec +
		                 (Time) res.ru_utime.tv_usec / 1000000.0
		                 - virtual_utime;
		virtual_slapse = (Time) res.ru_stime.tv_sec +
		                 (Time) res.ru_stime.tv_usec / 1000000.0
		                 - virtual_stime;
	}
}
#endif
