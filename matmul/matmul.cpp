/* An ArgoDSM@OpenMP implementation of matrix multiplication.
 *
 * It receives as input the dimension N and constructs three NxN matrices
 * (+1 for verification). We can enable verification with the -v argument.
 *
 * We initialize the matrices with prefixed values which we can later
 * check to ensure correctness of the computations.
 *
 * ArgoDSM/OpenMP version written by Ioannis Anevlavis - Eta Scale AB
 */

#include "argo.hpp"
#include "../common/wtime.hpp"

#include <omp.h>
#include <string>
#include <iostream>
#include <unistd.h>
#include <sys/time.h>

int size;
int impl;
int verify;
int nthreads;
int workrank;
int numtasks;

double *mat_a;
double *mat_b;
double *mat_c;
double *mat_r;

static void info();
static double get_time();
static void matmul_opt();
static void matmul_ref();
static int verify_result();
static void init_matrices();
static void run_multiply(const int &);
static void usage(std::ostream &, const char *);
static void distribute(int& beg, int& end, const int& loop_size, 
		const int& beg_offset, const int& less_equal);

#define at(x, y) ((x) * (size) + (y))

int main(int argc, char *argv[])
{
        argo::init(10*1024*1024*1024UL);

        int c;
        int errexit = 0;
        extern char *optarg;
        extern int optind, optopt, opterr;

        while ((c = getopt(argc, argv, "i:s:vh")) != -1) {
                switch (c) {
			case 'i':
				impl = atoi(optarg);
				break;
			case 's':
				size = atoi(optarg);
				break;
			case 'v':
				verify = 1;
				break;
			case 'h':
				usage(std::cout, argv[0]);
				exit(0);
				break;
			case ':':
				std::cerr << argv[0] << ": option -" << (char)optopt << " requires an operand"
					<< std::endl;
				errexit = 1;
				break;
			case '?':
				std::cerr << argv[0] << ": illegal option -- " << (char)optopt
					<< std::endl;
				errexit = 1;
				break;
			default:
				abort();
                }
        }

        if (errexit) {
                if (workrank == 0) usage(std::cerr, argv[0]);
                exit(2);
        }

        workrank = argo::node_id();
        numtasks = argo::number_of_nodes();

        #pragma omp parallel
        {
                #if defined(_OPENMP)
                	#pragma omp master
                        nthreads = omp_get_num_threads();
                #endif /* _OPENMP */
        }

        info();

        mat_a = argo::conew_array<double>(size * size);
        mat_b = argo::conew_array<double>(size * size);
        mat_c = argo::conew_array<double>(size * size);
        mat_r = argo::conew_array<double>(size * size);

        #pragma omp parallel
        {
                init_matrices();
                run_multiply(verify);
        }

        argo::codelete_array(mat_a);
        argo::codelete_array(mat_b);
        argo::codelete_array(mat_c);
        argo::codelete_array(mat_r);

        argo::finalize();

        return 0;
}

void init_matrices()
{
	int beg, end;
        int i, j;

	distribute(beg, end, size, 0, 0);

        #pragma omp for schedule(static)
        for (i = beg; i < end; i++) {
                for (j = 0; j < size; j++) {
                        mat_c[at(i,j)] = 0.0;
                        mat_r[at(i,j)] = 0.0;
                        mat_a[at(i,j)] = ((i + j) & 0x0F) * 0x1P-4;
                        mat_b[at(i,j)] = (((i << 1) + (j >> 1)) & 0x0F) * 0x1P-4;
                }
        }

        argo_barrier(nthreads);
}

void matmul_ref()
{
        int i, j, k;

        for (j = 0; j < size; j++) {
                for (i = 0; i < size; i++) {
                        for (k = 0; k < size; k++) {
                                mat_r[at(i,j)] += mat_a[at(i,k)] * mat_b[at(k,j)];
                        }
                }
        }
}

void matmul_opt()
{
	int beg, end;
        int i, j, k;
        double temp;

	distribute(beg, end, size, 0, 0);

        if      (impl == 0) {
                #pragma omp for schedule(static)
                for (j = beg; j < end; j++) {
                        for (i = 0; i < size; i++) {
                                for (k = 0; k < size; k++) {
                                        mat_c[at(i,j)] += mat_a[at(i,k)] * mat_b[at(k,j)];
                                }
                        }
                }
        }
        else if (impl == 1) {
                #pragma omp for schedule(static)
                for (i = beg; i < end; i++) {
                        for (j = 0; j < size; j++) {
                                for (k = 0; k < size; k++) {
                                        mat_c[at(i,j)] += mat_a[at(i,k)] * mat_b[at(k,j)];
                                }
                        }
                }
        }
        else if (impl == 2) {
                #pragma omp for schedule(static)
                for (i = beg; i < end; i++) {
                        for (k = 0; k < size; k++) {
                                temp = mat_a[at(i,k)];
                                for (j = 0; j < size; j++) {
                                        mat_c[at(i,j)] += temp * mat_b[at(k,j)];
                                }
                        }
                }
        }

        argo_barrier(nthreads);
}

int verify_result()
{
        int i, j;
        double e_sum = 0.0;

        for (i = 0; i < size; i++) {
                for (j = 0; j < size; j++) {
                        e_sum += mat_c[at(i,j)] < mat_r[at(i,j)] ?
                                 mat_r[at(i,j)] - mat_c[at(i,j)] :
                                 mat_c[at(i,j)] - mat_r[at(i,j)];
                }
        }

        return e_sum < 1E-6;
}

static void run_multiply(const int &verify)
{
        double time_start, time_stop;

        #pragma omp master
        time_start = get_time();
        matmul_opt();
        #pragma omp master
        time_stop = get_time();
        
        #pragma omp master
        if (workrank == 0) {
                std::cout.precision(4);
                std::cout << "Time: " << time_stop - time_start
                          << std::endl;
        }

        #pragma omp master
        if (workrank == 0) {
                if (verify) {
                        std::cout << "Verifying solution... ";
                        
                        time_start = get_time();
                        matmul_ref();
                        time_stop = get_time();

                        if (verify_result())
                                std::cout << "OK"
                                          << std::endl;
                        else
                                std::cout << "MISMATCH"
                                          << std::endl;

                        std::cout << "Reference runtime: " << time_stop - time_start
                                  << std::endl;
                }
		print_argo_stats();
        }
}

double get_time()
{
        struct timeval tv;

        if (gettimeofday(&tv, NULL)) {
                std::cerr << "gettimeofday failed. Aborting."
                          << std::endl;
                abort();
        }
        
        return tv.tv_sec + tv.tv_usec * 1E-6;
}

void usage(std::ostream &os, const char *argv0)
{
        if (workrank == 0) {
                os << "Usage: " << argv0 << " [OPTION]...\n"
                << "\n"
                << "Options:\n"
                << "\t-i\tSelect implementation <0:jik, 1:ijk, 2:ikj>\n"
                << "\t-s\tSize of matrices <N>\n"
                << "\t-v\tVerify solution\n"
                << "\t-h\tDisplay usage"
                << std::endl;
        }
}

void info()
{
        if (workrank == 0) {
                const std::string sverif = (verify == 0) ? "OFF" : "ON";
                const std::string sorder = (impl == 0) ? "jik" :
                                           (impl == 1) ? "ijk" :
                                           (impl == 2) ? "ikj" : "inv";

                std::cout << "MatMul: "           << size << "x" << size
                          << ", implementation: " << sorder
                          << ", verification: "   << sverif
                          << ", numtasks: "       << numtasks
                          << ", nthreads: "       << nthreads
                          << std::endl;
        }
}

void distribute(int& beg, int& end, const int& loop_size, 
		const int& beg_offset, const int& less_equal)
{
	int chunk = loop_size / numtasks;
	beg = workrank * chunk + ((workrank == 0) ? beg_offset : less_equal);
	end = (workrank != numtasks - 1) ? workrank * chunk + chunk : loop_size;
}
