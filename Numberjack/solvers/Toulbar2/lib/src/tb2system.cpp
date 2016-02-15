/*
 * ****** System dependent functions.
 */

#include "tb2types.hpp"
#include "tb2system.hpp"

/* --------------------------------------------------------------------
// Timer management functions
// -------------------------------------------------------------------- */
#ifdef LINUX 
#include <unistd.h> 
#include <sys/time.h>
#include <sys/times.h>

double cpuTime()
{
    static struct tms buf;

    times(&buf);
    double res = ((double) (buf.tms_utime+buf.tms_stime+buf.tms_cutime+buf.tms_cstime)) / ((double) sysconf(_SC_CLK_TCK));
    return (res>0)?res:0;
}

void timeOut(int sig)
{
    if(ToulBar2::verbose>=0) cout << endl << "Time limit expired... Aborting..." << endl;
    if (ToulBar2::timeOut) ToulBar2::timeOut();
    else exit(0);
}

static struct itimerval thetimer = { {0, 0}, {0, 0} };

/* set a timer (in seconds) */
void timer(int t)
{
    ToulBar2::interrupted = false;
    signal(SIGVTALRM, timeOut);
    thetimer.it_interval.tv_sec = 0;
    thetimer.it_interval.tv_usec = 0;
    thetimer.it_value.tv_sec = t;
    thetimer.it_value.tv_usec = 0;
    setitimer(ITIMER_VIRTUAL, &thetimer, NULL);
}

/* stop the current timer */
void timerStop()
{
    thetimer.it_value.tv_sec = 0;
    thetimer.it_value.tv_usec = 0;
    setitimer(ITIMER_VIRTUAL, &thetimer, NULL);
    ToulBar2::interrupted = false;
}

#else
double cpuTime()
{
    return (double) (clock() / CLOCKS_PER_SEC);
}
void timer(int t) {}
void timerStop() {}
#endif

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

