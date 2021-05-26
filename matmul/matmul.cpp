/**
 * @file
 * @brief An ArgoDSM@OpenMP implementation of matrix multiplication.
 * @copyright ArgoDSM@OpenMP version written by Ioannis Anevlavis - Eta Scale AB
 */

#include "argo.hpp"

#include <omp.h>
#include <string>
#include <iostream>
#include <unistd.h>
#include <sys/time.h>

/* arrays' size */
int NSIZE;
/* size of tile */
int BSIZE;

/* global vars  */
int nthreads;
int workrank;
int numtasks;

/* working mtxs */
double* mat_a;
double* mat_b;
double* mat_c;
double* mat_r;

void info(const int&,
	  const int&);
void usage(std::ostream&,
	   const char*);
void init_matrices();
void run_multiply(const int&,
		  const int&);
void matmul_opt(const int&);
void matmul_ref(const int&);
void distribute(int&,
		int&,
		const int&, 
		const int&,
		const int&);
double get_time();
int verify_result();

#define at(x, y) ((x) * (NSIZE) + (y))

int
main(int argc, char *argv[])
{
        argo::init(500*1024*1024UL);

        int c;
	int verify{0};
        int errexit{0};
	int transpose{0};
        extern char *optarg;
        extern int optind, optopt, opterr;

	workrank = argo::node_id();
        numtasks = argo::number_of_nodes();

        #pragma omp parallel
        {
                #if defined(_OPENMP)
                	#pragma omp master
                        nthreads = omp_get_num_threads();
                #endif /* _OPENMP */
        }

        while ((c = getopt(argc, argv, "s:b:vth")) != -1) {
                switch (c) {
			case 's':
				NSIZE = atoi(optarg);
				break;
			case 'b':
				BSIZE = atoi(optarg);
				break;
			case 'v':
				verify = 1;
				break;
			case 't':
				transpose = 1;
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

	if (workrank == 0) {
		/* Check 0 */
		if (errexit) {
			usage(std::cerr, argv[0]);
			exit(errexit);
		} else {
			info(verify, transpose);
		}
		/* Check 1 */
		if (NSIZE % BSIZE) {
			std::cerr << "Error: The block size needs to "
				  << "divide the matrix dimensions."
				  << std::endl;
			exit(EXIT_FAILURE);
		}
	}
        
        mat_a = argo::conew_array<double>(NSIZE * NSIZE);
        mat_b = argo::conew_array<double>(NSIZE * NSIZE);
        mat_c = argo::conew_array<double>(NSIZE * NSIZE);
        mat_r = argo::conew_array<double>(NSIZE * NSIZE);

        #pragma omp parallel
        {
                init_matrices();
                run_multiply(verify, transpose);
        }

        argo::codelete_array(mat_a);
        argo::codelete_array(mat_b);
        argo::codelete_array(mat_c);
        argo::codelete_array(mat_r);

        argo::finalize();

        return 0;
}

void
init_block(const int& i,
	   const int& j)
{
	for (int ii = i; ii < i+BSIZE; ii++)
		for (int jj = j; jj < j+BSIZE; jj++) {
			mat_c[at(ii,jj)] = 0.0;
			mat_r[at(ii,jj)] = 0.0;
			mat_a[at(ii,jj)] = ((ii + jj) & 0x0F) * 0x1P-4;
			mat_b[at(ii,jj)] = (((ii << 1) + (jj >> 1)) & 0x0F) * 0x1P-4;
		}
}

void
init_matrices()
{
	int beg, end;
	distribute(beg, end, NSIZE, 0, 0);

        #pragma omp for schedule(static)
        for (int i = beg; i < end; i += BSIZE)
                for (int j = 0; j < NSIZE; j += BSIZE)
			init_block(i, j);

        argo::barrier(nthreads);
}

void
multiply_block(const int& i,
	       const int& j,
	       const int& k,
	       const int& transpose)
{
	for (int ii = i; ii < i+BSIZE; ii++)
		for (int jj = j; jj < j+BSIZE; jj++)
			for (int kk = k; kk < k+BSIZE; kk++)
if (!transpose)
				mat_c[at(ii,jj)] += mat_a[at(ii,kk)] * mat_b[at(kk,jj)];
else
				mat_c[at(ii,jj)] += mat_a[at(ii,kk)] * mat_b[at(jj,kk)];
}

void
matmul_opt(const int& transpose)
{
	int beg, end;
	distribute(beg, end, NSIZE, 0, 0);

	#pragma omp for schedule(static)
	for (int i = beg; i < end; i += BSIZE)
		for (int j = 0; j < NSIZE; j += BSIZE)
			for (int k = 0; k < NSIZE; k += BSIZE)
				multiply_block(i, j, k, transpose);

        argo::barrier(nthreads);
}

void
matmul_ref(const int& transpose)
{
        for (int i = 0; i < NSIZE; i++)
                for (int j = 0; j < NSIZE; j++)
                        for (int k = 0; k < NSIZE; k++)
if (!transpose)
                                mat_r[at(i,j)] += mat_a[at(i,k)] * mat_b[at(k,j)];
else
				mat_r[at(i,j)] += mat_a[at(i,k)] * mat_b[at(j,k)];
}

int
verify_result()
{
	double e_sum{0.0};

        for (int i = 0; i < NSIZE; i++) {
                for (int j = 0; j < NSIZE; j++) {
			e_sum += (mat_c[at(i,j)] < mat_r[at(i,j)])
				? mat_r[at(i,j)] - mat_c[at(i,j)]
				: mat_c[at(i,j)] - mat_r[at(i,j)];
                }
        }

        return e_sum < 1E-6;
}

void
run_multiply(const int& verify,
	     const int& transpose)
{
        double time_start, time_stop;

        #pragma omp master
        time_start = get_time();
        matmul_opt(transpose);
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
                        matmul_ref(transpose);
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
        }
}

double
get_time()
{
        struct timeval tv;

        if (gettimeofday(&tv, NULL)) {
                std::cerr << "gettimeofday failed. Aborting."
                          << std::endl;
                abort();
        }
        
        return tv.tv_sec + tv.tv_usec * 1E-6;
}

void
usage(std::ostream &os,
      const char *argv0)
{
	os << "Usage: " << argv0 << " [OPTION]...\n"
	   << "\n"
	   << "Options:\n"
	   << "\t-s\tSize of matrices <N>\n"
	   << "\t-v\tVerify solution\n"
	   << "\t-h\tDisplay usage"
	   << std::endl;
}

void
info(const int& verify,
     const int& transpose)
{
	const std::string sverif = (verify == 0)    ? "OFF" : "ON";
	const std::string strans = (transpose == 0) ? "OFF" : "ON";

	std::cout << "MatMul: "          << NSIZE << "x" << NSIZE
		  << ", block size: "    << BSIZE
		  << ", verification: "  << sverif
		  << ", transposition: " << strans
		  << ", numtasks: "      << numtasks
		  << ", nthreads: "      << nthreads
		  << std::endl;
}

void
distribute(int& beg,
	   int& end,
	   const int& loop_size, 
	   const int& beg_offset,
	   const int& less_equal)
{
	int chunk = loop_size / numtasks;
	beg = workrank * chunk + ((workrank == 0) ? beg_offset : less_equal);
	end = (workrank != numtasks - 1) ? workrank * chunk + chunk : loop_size;
}
