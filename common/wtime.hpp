#include "argo.hpp"

#include <omp.h>
#include <sys/time.h>

void wtime(double *t){
	static int sec = -1;
	struct timeval tv;
	gettimeofday(&tv, 0);
	if (sec < 0) sec = tv.tv_sec;
	*t = (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;
}

typedef struct {
	int locks{0};
	int barrs{0};
	double locktime{0.0};
	double barrtime{0.0};
} lock_barr_t;

lock_barr_t argo_stats;
double lock_t1, lock_t2;
double barr_t1, barr_t2;

/**
 * @note: this will only work if the application is written for
 * only one thread per node to capture the global lock.
 */
static inline __attribute__((always_inline))
void argo_lock(argo::globallock::global_tas_lock *lock) {
	wtime(&lock_t1);
	lock->lock();
}

static inline __attribute__((always_inline))
void argo_unlock(argo::globallock::global_tas_lock *lock) {
	lock->unlock();
	wtime(&lock_t2);
	argo_stats.locktime += lock_t2 - lock_t1;
	argo_stats.locks++;
}

/**
 * @note: we need to overload this function and not supply a
 * default argument to make sure that it is inlined.
 */ 
static inline __attribute__((always_inline))
void argo_barrier() {
	wtime(&barr_t1);
	argo::barrier();
	wtime(&barr_t2);
	argo_stats.barrtime += barr_t2 - barr_t1;
	argo_stats.barrs++;
}

static inline __attribute__((always_inline))
void argo_barrier(int nthreads) {
	#pragma omp master
		wtime(&barr_t1);
	argo::barrier(nthreads);
	#pragma omp master
	{
		wtime(&barr_t2);
		argo_stats.barrtime += barr_t2 - barr_t1;
		argo_stats.barrs++;
	}
}

static inline
void print_argo_stats() {
	printf("#####################STATISTICS#########################\n");
	printf("Argo locks : %d, barriers : %d\n",
		argo_stats.locks, argo_stats.barrs);
	printf("Argo locktime : %.3lf sec., barriertime : %.3lf sec.\n",
		argo_stats.locktime, argo_stats.barrtime);
	printf("########################################################\n\n");
}
