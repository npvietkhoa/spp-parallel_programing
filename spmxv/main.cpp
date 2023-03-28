// C header
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// utility header
#include "utils/ooo_cmdline.h"

void spmxv(ooo_options *tOptions, ooo_input *tInput)
{
	int iNumRepetitions = tOptions->iNumRepetitions; // set with -r <numrep>
	bool bCCNuma = tOptions->bOptimizeCcNuma; // set with -c

	// setup data structures
	// your task: decide on appropriate data types and, if necessary, implement conversion
	double *y = (double*) malloc(sizeof(double) * (tInput->stNumRows)); // result
	double *Aval = (double*) malloc(sizeof(double) * (tInput->stNumNonzeros)); // values
	int *Acol = (int*) malloc(sizeof(int) * (tInput->stNumNonzeros)); // column indices
	int *Arow = (int*) malloc(sizeof(int) * (tInput->stNumRows + 1)); // begin of each row
	double *x = (double*) malloc(sizeof(double) * (tInput->stNumRows)); // RHS

	// allocate helper data
	timespan *timings = (timespan*) malloc(sizeof(timespan) * iNumRepetitions);
	double t1, t2;

	int i, rep;

	// helper data to manually calculate the "schedule"
	int nthreads = tOptions->iNumThreads;
	int *rowStart = (int*) malloc(sizeof(int) * (nthreads));
	int *rowEnd   = (int*) malloc(sizeof(int) * (nthreads)); //<- consider this variable as the **start** of the next batch of rows, analog to how Arow works
	int nnzSum = 0, currThread = 0, nnzPerThread = ((tInput->stNumNonzeros) / nthreads);

	// manually calculate the "schedule"
	cout << endl;
	cout << "Calculating work distribution..." << endl;
	
	// actual allocation
	rowStart[currThread] = 0;
	for (i = 0; i < tInput->stNumRows; i++) {
		nnzSum += tInput->row[i + 1] - tInput->row[i];
		// cout << "nnzSum at i=" << i << ": " << nnzSum << endl;

		if (nnzSum >= nnzPerThread) {

			rowEnd[currThread] = i;
			currThread++;
			rowStart[currThread] = i;
			nnzSum = 0;
		}
	}
	rowEnd[currThread] = tInput->stNumRows;

	// print the work allocation for debugging purposes:
	cout << "Work allocation for " << nthreads << " Threads - ideal: " 
		<< nnzPerThread << " Elements per thread" << endl;

	for (i = 0; i < nthreads; i++) {
		cout << "Thread " << i << "\t " << rowStart[i] << " to \t " << rowEnd[i] 
			<< "\t (" << tInput->row[rowEnd[i]] - tInput->row[rowStart[i]] << " Elements)" << endl;
	}
	cout << endl;

	cout << "Initialization..." << endl;
	// initialize data structures
	#pragma omp parallel
	{
		int thread = omp_get_thread_num();

		for (i = rowStart[thread]; i < rowEnd[thread]; i++)
		{
			Arow[i] = tInput->row[i];
			y[i] = 0.0;
			x[i] = 1;
		}
	
		for (i = rowStart[thread]; i < rowEnd[thread]; i++)
		{
			unsigned int rowbeg = Arow[i];
			unsigned int rowend = Arow[i+1];
			unsigned int nz;
			for (nz = rowbeg; nz < rowend; nz++)
			{
				Aval[nz] = tInput->val[nz];
				Acol[nz] = tInput->col[nz];
			}
		}
	}
	Arow[tInput->stNumRows] = tInput->stNumNonzeros;

	cout << "Kernel..." << endl;

	// take the time: start
	t1 = omp_get_wtime();

	// run kernel iNumRepetitions times
	#pragma omp parallel
	for (rep = 0; rep < iNumRepetitions; rep++)
	{
		timings[rep].dBegin = omp_get_wtime();
		// cout << "Iteration " << rep << endl;
		// yAx-kernel

		// iterate over all rows
		// #pragma omp parallel 
		// {
			int thread = omp_get_thread_num();

			for (int i = rowStart[thread]; i < rowEnd[thread]; i++) {
				// reset result vector element
				y[i] = 0;

				// do maths, iterate over column elements of A and RHS vector
				// #pragma omp parallel for reduction(+:y[i]) num_threads(4)
				for (int j = Arow[i]; j < Arow[i + 1]; j++) {

					// y(i) += A(j) + x(j)
					y[i] += Aval[j];
				}
			}

			// #pragma omp barrier
		// }
	}

	t2 = omp_get_wtime();

	// error check
	print_error_check(y, tInput);

	// process_results
	print_performance_results(tOptions, t1, t2, timings, tInput);

	// cleanup
	free(y);
	free(Aval);
	free(Acol);
	free(Arow);
	free(x);
	free(timings);
	free(rowStart);
	free(rowEnd);
}


int main(int argc, char* argv[])
{
	// parse command line
	ooo_options tOptions;
	if (! parseCmdLine(&tOptions, argc, argv))
	{
		return EXIT_FAILURE;
	}

	// load filename
	ooo_input tInput;
	if (! loadInputFile_4SMXV(&tOptions, &tInput))
	{
		return EXIT_FAILURE;
	}

	// SpMXV-Kernel
	spmxv(&tOptions, &tInput);

	return EXIT_SUCCESS;
}
