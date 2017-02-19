#include "depends.h"
#ifdef HAVE_GSL_HEADER
#include "ridgeRegressionFunctions.h"

#if _CUDA_
#include "cudaOnlyFunctions.h"
#endif

/* Function prototype - read in phenotypes for linear regression */
GSL_TYPE(vector) * readLinearPhenotypes(char * phenotypefilename,
					int NINDIV);

		 /* Function prototype - compute Kr based on a, Z, y and r */
		 int computeLinearKr(GSL_TYPE(vector) * a,
				     GSL_TYPE(matrix) * Z,
				     GSL_TYPE(vector) * y,
				     GSL_TYPE(vector) * D2,
				     int r,
				     PREC * kr,
				     PREC * DofF);

#endif

