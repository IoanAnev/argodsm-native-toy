/* An ArgoDSM@OpenMP implementation of the matrix-vector multiplication.
 * kernel: y += A*x.
 *
 * It receives as input the dimensions ([M, N]) of the problem. We can
 * control the number of times that we will execute the matvec kernel
 * with the ITER argument.
 *
 * We initialize the vectors with prefixed values which we can later check to
 * ensure the correctness of the computation.
 *
 * ArgoDSM/OpenMP version written by Ioannis Anevlavis - Eta Scale AB
 */

#include "argo.hpp"

#include <omp.h>
#include <iostream>

int workrank, numtasks, nthreads;

void distribute(size_t& beg, size_t& end, const size_t& loop_size, const size_t& beg_offset, const size_t& less_equal){
	size_t chunk = loop_size / numtasks;
	beg = workrank * chunk + ((workrank == 0) ? beg_offset : less_equal);
	end = (workrank != numtasks - 1) ? workrank * chunk + chunk : loop_size;
}

void matvec(const size_t& beg, const size_t& end, double *A, size_t N, double *x, double *y)
{
	#pragma omp for nowait schedule(static)
	for (size_t i = beg; i < end; ++i) {
		double res = 0.0;
		for (size_t j = 0; j < N; ++j) {
			res += A[i * N + j] * x[j];
		}
		
		y[i] += res;
	}
}

void init(const size_t& beg, const size_t& end, double *vec, double value)
{
	#pragma omp for schedule(static)
	for (size_t i = beg; i < end; ++i) {
		vec[i] = value;
	}
}

void check_result(size_t M, double *A, size_t N, double *x, double *y,
		size_t ITER)
{
	double *y_serial = new double[M];

	#pragma omp parallel
	{
		init(0, M, y_serial, 0);
		
		for (size_t iter = 0; iter < ITER; ++iter) {
			matvec(0, M, A, N, x, y_serial);
		}
	}
	
	for (size_t i = 0; i < M; ++i) {
		if (y_serial[i] != y[i]) {
			printf("FAILED\n");
			delete[] y_serial;
			return;
		}
	}
	
	printf("SUCCESS\n");
	delete[] y_serial;
}

void usage()
{
	fprintf(stderr, "usage: matvec_strong M N ITER [CHECK]\n");
	return;
}

int main(int argc, char *argv[])
{
	argo::init(10*1024*1024*1024UL);

	workrank = argo::node_id();
	numtasks = argo::number_of_nodes();

	#pragma omp parallel
	{
		#if defined(_OPENMP)
			#pragma omp master
			nthreads = omp_get_num_threads();
		#endif /* _OPENMP */
	}

	size_t M, N, ITER;
	double *A, *x, *y;
	int check = false;
	struct timespec tp_start, tp_end;

	if (argc != 4 && argc != 5) {
		if (workrank == 0) usage();
		return -1;
	}
	
	M = atoi(argv[1]);
	N = atoi(argv[2]);
	ITER = atoi(argv[3]);
	
	if (argc == 5) {
		check = atoi(argv[4]);
	}
	
	A = argo::conew_array<double>(M * N);
	y = argo::conew_array<double>(M);
	x = argo::conew_array<double>(N);


	size_t begM, endM, begN, endN;
	distribute(begM, endM, M, 0, 0);
	distribute(begN, endN, N, 0, 0);

	#pragma omp parallel
	{
		init(begM, endM, y, 0);
		init(begN, endN, x, 1);
		
		#pragma omp for schedule(static)
		for (size_t i = begM; i < endM; ++i) {
			init(0, N, &A[i * N], 2);
		}
	}
	argo::barrier();
	
	clock_gettime(CLOCK_MONOTONIC, &tp_start);
	
	#pragma omp parallel
	{
		for (size_t iter = 0; iter < ITER; ++iter) {
			matvec(begM, endM, A, N, x, y);
		}
	}
	argo::barrier();

	clock_gettime(CLOCK_MONOTONIC, &tp_end);
	
	if (workrank == 0) {
		if (check) {
			check_result(M, A, N, x, y, ITER);
		}
			
		double time_msec = (tp_end.tv_sec - tp_start.tv_sec) * 1e3
			+ ((double)(tp_end.tv_nsec - tp_start.tv_nsec) * 1e-6);
		
		double mflops =
			ITER 			/* 'ITER' times of kernel FLOPS */
			* 3 * M * N 		/* 3 operations for every element of A */
			/ (time_msec / 1000.0) 	/* time in seconds */
			/ 1e6; 			/* convert to Mega */

		printf("M:%zu N:%zu ITER:%zu NR_PROCS:%d CPUS:%d TIME_MSEC:%.2lf MFLOPS:%.2lf\n",
			M, N, ITER, numtasks, nthreads,
			time_msec, mflops);
	}
	
	argo::codelete_array(A);
	argo::codelete_array(y);
	argo::codelete_array(x);

	argo::finalize();
	
	return 0;
}
