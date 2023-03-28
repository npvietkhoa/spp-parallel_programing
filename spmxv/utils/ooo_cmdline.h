
#ifndef INC_OOOCMDLINE_H
#define INC_OOOCMDLINE_H

// C++ header
#include <string>
#include <limits>
#include <ios>

// utility header
#include "anyoption.h"

struct ooo_options
{
	int			iNumThreads;					// number of threads to be used
	bool		bOptimizeCcNuma;				// optimize for cc-NUMA architecture
	int			iNumRepetitions;				// number of repetitions
	string		strFilename;					// use input file
};

struct timespan
{
	double dBegin;
	double dEnd;
};

struct ooo_input
{
	int				stNumNonzeros;
	int				stNumRows;
	double			*val;
	int				*col;
	int				*row;
	double			*x;
};

void print_performance_results(ooo_options *tOptions, double t1, double t2, timespan *timings, ooo_input *tInput);
bool parseCmdLine(ooo_options *options, int argc, char* argv[]);
bool loadInputFile_4SMXV(ooo_options *tOptions, ooo_input *tInput);
void print_error_check(double *result, ooo_input *tInput);

/**
 * Loads sparse crs matrix in DROPS human readable format
 * @param strFilename Filename of matrix
 */
template<typename T, typename indexT>
void load_drops_matlab_matrix( std::string strFilename_LHS, indexT* & row, indexT* & col, T* & val, indexT & stNumRows, indexT & stNumCols, indexT & stNumNonzeros)
{
	indexT stDimension = 0;

    std::ifstream fdat( strFilename_LHS.c_str() );
    if (!fdat)
    {
    	std::cout << "\n";
        std::cout << "ERROR : could not open input file!\n";
        std::cout << "\n\n";
        throw new ios_base::failure("File Input Error!");
	}

    std::cout << "reading matrix in matlab format from " << strFilename_LHS << std::endl;

    // Read the leading comment: % rows 'x' columns nonzeros "nonzeros\n"
    while ( fdat && fdat.get() != '%');
    fdat >> stNumRows;
    while ( fdat && fdat.get() != 'x');
    fdat >> stNumCols >> stNumNonzeros;
    while ( fdat && fdat.get() != '\n');
    if (!fdat)
    {
        std::cout << "SparseMatBaseCL operator>>: Missing \"% rows cols nz\" comment.\n" << std::endl;
		throw new ios_base::failure("File Input Error!");
    }

    val = new T[stNumNonzeros];
    col = new indexT[stNumNonzeros];
    row = new indexT[stNumRows+1];
    row[stNumRows] = stNumNonzeros;

    row[0] = 0;

    size_t nz, cr;
    cr = 0;
    for ( nz= 0; nz < stNumNonzeros; ++nz )
    {
        indexT r, c;
        T v;
        fdat >> r >> c >> v;
        if ( !fdat )
		{
            std::cout << "SparseMatBaseCL operator>>: Stream corrupt.\n" << std::endl;
			throw new ios_base::failure("File Input Error!");
		}

        if ( r-1 != cr )
        	row[ ++cr ] = nz;
        col[ nz ] = c-1;
        val[ nz ] = v; // save nonzero
    }

    stDimension = stNumCols;
    fdat.close();
}

#endif
