#include "depends.h"
#ifdef HAVE_GSL_HEADER

#include "commonFunctions.h"
#include "computeLinearRidge.h"
#include "ReadInData.h"
#include "ridgeRegressionFunctions.h"

/* Header file for thin data functions in the regression package */

gsl_vector_int * readThinFile(char * thinfilename,
			      char ** SNPNAMES,
			      int thinning_distance,
			      int NINDIV,
			      int NSNPS,
			      int * nThinnedSnps,
			      int verbose);

int readSNPsThinAndComputePCs(char * genofilename,
			      gsl_vector_int * thin,
			      GSL_TYPE(matrix) * Z,
			      GSL_TYPE(matrix) * thinnedGenotypes,
			      GSL_TYPE(vector) * D2,
			      int * howManyK);

#endif

