#include "depends.h"
#ifdef HAVE_GSL_HEADER

/* Function prototype - read in phenotypes for logistic regression */
gsl_vector_int * readLogisticPhenotypes(char * phenotypefilename,
	int NINDIV);

/* Function prototype - returning to original scale - component-wise computation of regression coefficients */
int returnToOriginalScaleLogistic(GSL_TYPE(vector) * betaOut,
				  GSL_TYPE(vector) * Bridge,
				  GSL_TYPE(vector) * means,
				  GSL_TYPE(vector) * scales,
				  int intercept_flag);
#endif

