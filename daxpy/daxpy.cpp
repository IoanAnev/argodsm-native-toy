/* An ArgoDSM@OpenMP implementation of the daxpy BLAS operation.
 *
 * It receives as input the dimension N of the vectors. We can
 * control the number of times that the daxpy kernel will execute
 * with the ITER argument.
 *
 * We initialize the vectors with prefixed values which we can later
 * check to ensure correctness of the computations
 *
 * ArgoDSM/OpenMP version written by Ioannis Anevlavis - Eta Scale AB
 */

#include "argo.hpp"

#include <thread>
#include <vector>
#include <iostream>

int workrank, numtasks, nthreads;

void distribute(size_t& beg, size_t& end, const size_t& loop_size, const size_t& beg_offset, const size_t& less_equal){
	size_t chunk = loop_size / numtasks;
	beg = workrank * chunk + ((workrank == 0) ? beg_offset : less_equal);
	end = (workrank != numtasks - 1) ? workrank * chunk + chunk : loop_size;
}

void daxpy(const size_t& beg, const size_t& end, double *x, double alpha, double *y)
{
	for (size_t i = beg; i < end; ++i) {
		y[i] += alpha * x[i];
	}
}

void init(const size_t& beg, const size_t& end, double *vector, double value)
{	
	for (size_t i = beg; i < end; ++i) {
		vector[i] = value;
	}
}

void check_result(size_t N, double *x, double alpha, double *y, size_t ITER)
{
	double *y_serial = new double[N];

	init(0, N, y_serial, 0);

	for (size_t iter = 0; iter < ITER; ++iter) {
		daxpy(0, N, x, alpha, y_serial);
	}

	for (size_t i = 0; i < N; ++i) {
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
	fprintf(stderr, "usage: daxpy_strong N ITER [CHECK]\n");
	return;
}

int main(int argc, char *argv[])
{
	argo::init(256*1024*1024UL,
             256*1024*1024UL);

	workrank = argo::node_id();
	numtasks = argo::number_of_nodes();

	nthreads = std::thread::hardware_concurrency();
	std::vector<std::thread> threads(nthreads);

	size_t N, ITER;
	double alpha = 3.14, *x, *y;
	bool check = false;
	struct timespec tp_start, tp_end;
	
	if (argc != 3 && argc != 4) {
		if (workrank == 0) usage();
		return -1;
	}
	
	N = atol(argv[1]);
	ITER = atol(argv[2]);
	
	if (argc == 4) {
		check = atoi(argv[3]);
	}
	
	x = argo::conew_array<double>(N);
	y = argo::conew_array<double>(N);
	
	size_t beg, end;
	distribute(beg, end, N, 0, 0);
	
	for (int i = 0; i < threads.size(); ++i) {
		threads[i] = std::thread([&,i] {
			int chunk = (end - beg) / nthreads;
			int beg_i = beg+i*chunk;
			int end_i = beg_i+chunk;

			init(beg_i, end_i, y, 0);
			init(beg_i, end_i, x, 42);
		});
	}

	for (auto& t : threads)
		t.join();
	argo::barrier();
	
	clock_gettime(CLOCK_MONOTONIC, &tp_start);

	for (int i = 0; i < threads.size(); ++i) {
		threads[i] = std::thread([&,i] {
			int chunk = (end - beg) / nthreads;
			int beg_i = beg+i*chunk;
			int end_i = beg_i+chunk;

			for (size_t iter = 0; iter < ITER; ++iter) {
				daxpy(beg_i, end_i, x, alpha, y);
			}
		});
	}

	for (auto& t : threads)
		t.join();
	argo::barrier();
	
	clock_gettime(CLOCK_MONOTONIC, &tp_end);
	
	if (workrank == 0) {
		if (check) {
			check_result(N, x, alpha, y, ITER);
		}

		double time_msec = (tp_end.tv_sec - tp_start.tv_sec) * 1e3
			+ ((tp_end.tv_nsec - tp_start.tv_nsec) * 1e-6);
		
		double mflops =
			ITER 			/* 'ITER' times of kernel FLOPS */
			* 3 * N 		/* 3 operations for every vector element */
			/ (time_msec / 1000.0) 	/* time in seconds */
			/ 1e6;			/* convert to Mega */
		
		printf("N:%zu ITER:%zu NR_PROCS:%d CPUS:%d TIME_MSEC:%.2lf MFLOPS:%.2lf\n",
			N, ITER, numtasks, nthreads,
			time_msec, mflops);
	}

	argo::codelete_array(x);
	argo::codelete_array(y);

	argo::finalize();

	return 0;
}
