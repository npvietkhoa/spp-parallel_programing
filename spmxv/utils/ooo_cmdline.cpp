#include "ooo_cmdline.h"
#include <math.h>

void print_performance_results(ooo_options *tOptions, double t1, double t2, timespan *timings, ooo_input *tInput)
{
	// compute metrics
	double dTotalExperimentTime = t2 - t1;
	double dMinTime = numeric_limits<double>::max();
	double dMaxTime = numeric_limits<double>::min();
	double dMeanTime = 0.0;
	for (int i = 0; i < tOptions->iNumRepetitions; i++)
	{
		double d = timings[i].dEnd - timings[i].dBegin;
		dMinTime = std::min<double>(dMinTime, d);
		dMaxTime = std::max<double>(dMaxTime, d);
		dMeanTime += d;
	}
	dMeanTime /= tOptions->iNumRepetitions;

	double temp = tInput->stNumNonzeros / 1000000. * 2.;
	double mTotalFlops = tOptions->iNumRepetitions / dTotalExperimentTime * temp;
	double mMinFlops = temp / dMaxTime;
	double mMaxFlops = temp / dMinTime;
	double mMeanFlops = temp / dMeanTime;

	// print results
	std::cout << std::endl;
	std::cout << "Configuration              " << std::endl;
	std::cout << "Number of Threads:         " << tOptions->iNumThreads << std::endl;
	std::cout << "Number of Repetitions:     " << tOptions->iNumRepetitions << std::endl;
	std::cout << "Support for cc-NUMA arch.: " << tOptions->bOptimizeCcNuma << std::endl;
	std::cout << "Input filename:            " << tOptions->strFilename << std::endl;
	std::cout << std::endl;
	std::cout << "Time measurements          " << std::endl;
	std::cout << "Total experiment time:     " << dTotalExperimentTime << std::endl;
	std::cout << "Minimum kernel time:       " << dMinTime << std::endl;
	std::cout << "Maximum kernel time:       " << dMaxTime << std::endl;
	std::cout << "Arithm. Mean kernel time:  " << dMeanTime << std::endl;
	std::cout << std::endl;
	std::cout << "Performance results        " << std::endl;
	std::cout << "Total MFlops/s:            " << mTotalFlops << std::endl;
	std::cout << "Minimum MFlops/s:          " << mMinFlops << std::endl;
	std::cout << "Maximum MFlops/s:          " << mMaxFlops << std::endl;
	std::cout << "Arithm. Mean MFlops/s:     " << mMeanFlops << std::endl;
	std::cout << std::endl;
}

bool parseCmdLine(ooo_options *options, int argc, char* argv[])
{
	AnyOption *opt = new AnyOption();

	// set usage help text
	opt->addUsage("");
	opt->addUsage("Usage:");
	opt->addUsage("");
	opt->addUsage(" -h  --help:                 Prints this help text.");
	opt->addUsage(" -t  --threads num:          Number of threads to be used.");
	opt->addUsage(" -c  --ccnuma:               Optimize for cc-NUMA architecture (default: off).");
	opt->addUsage(" -f  --filename name:        Use input file (default: none).");
	opt->addUsage(" -r  --repetitions num:      Number of repetitions (default: 10).");
	opt->addUsage("");
	opt->addUsage("");

	// set options and flags
	opt->setFlag("help", 'h');
	opt->setOption("threads", 't');
	opt->setFlag("ccnuma", 'c');
	opt->setOption("filename", 'f');
	opt->setOption("repetitions", 'r');

	// process commandline
	opt->processCommandArgs(argc, argv);
	if (! opt->hasOptions())
	{
		opt->printUsage();
		delete opt;
		return false;
	}

	// get the values
	if (opt->getFlag("help") || opt->getFlag('h'))
	{
		opt->printUsage();
		delete opt;
		return false;
	}
	if (opt->getValue("threads") != NULL || opt->getValue('t') != NULL)
	{
		char *strNumThreads = opt->getValue('t');
		options->iNumThreads = atoi(strNumThreads);
	}
	else
	{
		opt->printUsage();
		delete opt;
		return false;
	}
	if (opt->getFlag("ccnuma") || opt->getFlag('c'))
	{
		options->bOptimizeCcNuma = true;
	}
	else
	{
		options->bOptimizeCcNuma = false;
	}
	if (opt->getValue("filename") != NULL || opt->getValue('f') != NULL)
	{
		options->strFilename = opt->getValue('f');
	}
	else
	{
		std::cerr << "ERROR: input filename required for this experiment!" << std::endl;
		std::cout << std::endl;
		delete opt;
		return false;
	}
	if (opt->getValue("repetitions") != NULL || opt->getValue('r') != NULL)
	{
		char *strNumRepetitions = opt->getValue('r');
		int iNumRepetitions = atoi(strNumRepetitions);
		if (iNumRepetitions <= 0)
		{
			opt->printUsage();
			delete opt;
			return false;
		}
		options->iNumRepetitions = iNumRepetitions;
	}
	else
	{
		options->iNumRepetitions = 10;
	}

	// done
	delete opt;
	return true;
}

bool loadInputFile_4SMXV(ooo_options *tOptions, ooo_input *tInput)
{
	// open file
	std::ifstream fdat(tOptions->strFilename.c_str(), std::fstream::binary);
	if (! fdat)
	{
		std::cerr << "ERROR: could not open input filename: " << tOptions->strFilename << std::endl;
		std::cout << std::endl;
		return false;
	}

	int iDummy;
	load_drops_matlab_matrix<double, int>(tOptions->strFilename.c_str(),
		tInput->row, tInput->col, tInput->val, tInput->stNumRows, iDummy,
		tInput->stNumNonzeros);

	return true;
}

// checks if the MXV is done correctly by computing the sum of all elements in the result vector
void print_error_check(double *result, ooo_input *tInput)
{
	printf("\nCorrectness check\n");
	
	// determine correct result for input matrix
	double correctResult, errorMargin;
	if (tInput->stNumRows == 59319 && tInput->stNumNonzeros == 853747) // Small matrix
	{
		correctResult = 229.02006;
		errorMargin = 1e-5;
	}
	else if (tInput->stNumRows == 493039 && tInput->stNumNonzeros == 7246747) // Large matrix
	{
		correctResult = 469.02944;
		errorMargin = 1e-5;
	}
	else if (tInput->stNumRows == 1585478 && tInput->stNumNonzeros == 7660826) // G3_circuit matrix
	{
		correctResult = 691071695.44214;
		errorMargin = 1e-3;
	}
	else if (tInput->stNumRows == 1391349 && tInput->stNumNonzeros == 64531701) // Serena matrix
	{
		correctResult = 9735110469227511.20605;
		errorMargin = 1000;
	}
	else
	{
		printf("Unknown Matrix, result correctness can not be determined!\n");
		return;
	}
	
	// get the computed result
	double sumResult = 0.0;
	for (int i = 0; i < tInput->stNumRows; i++) {
		sumResult += result[i];
	}
	
	// check if the computed result is correct
	double diff = fabs(sumResult - correctResult);
	if (diff < errorMargin)
	{
		printf("Success, correct result: %.5f +- %f\n", sumResult, errorMargin);
	}
	else
	{
		printf("Error, difference to correct result: %.5f\n", diff);
	}

}

