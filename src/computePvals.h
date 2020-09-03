/* 
   header file for computePvals.h 
   Contains prototypes for all the computePvals functions we could need
   i.e. the approx ones and the permutation ones for 
   different sorts of models
*/

/* includes */
#include "depends.h"
#ifdef HAVE_GSL_HEADER

#include "ReadInData.h"
#include "ridgeRegressionFunctions.h"
#include "coordinateDescent.h"


/* float version of cumulative distribution function */
float my_ugaussian_function(float x);

/* Compute Approx Ps - Linear */
int computeApproxPsLinear(GSL_TYPE(vector) * B,
			  GSL_TYPE(vector) * y,
			  GSL_TYPE(matrix) * U,
			  GSL_TYPE(vector) * D,
			  GSL_TYPE(vector) * D2,
			  GSL_TYPE(matrix) * V,
			  PREC k,
			  GSL_TYPE(vector) * approxPs);

/* Compute Approx Ps - Generalized Linear*/
int computeApproxPsGeneralizedLinear(GSL_TYPE(vector) * beta,
				     GSL_TYPE(matrix) * predictors,
				     GSL_TYPE(vector) * y,
				     GSL_TYPE(vector) * shrinkage,
				     int intercept_flag,
				     GSL_TYPE(vector) * approxPs);

/* Compute ApproxPs - Logistic */
int computeApproxPsLogistic(GSL_TYPE(vector) * B,
			    GSL_TYPE(matrix) * X, 
			    GSL_TYPE(vector) * shrinkage, 
			    int intercept_flag,
			    GSL_TYPE(vector) * approxPs);


/* Compute PermPs - all models */
int computePermPs(GSL_TYPE(vector) * permPs,
		  GSL_TYPE(matrix) * pred,
		  GSL_TYPE(vector) * pheno_linear,
		  gsl_vector_int * pheno_logistic,
		  GSL_TYPE(vector) * Bridge,
		  PREC lambda,
		  GSL_TYPE(vector) * shrinkage,
		  int NPERM,
		  int SEED,
		  int intercept_flag,
		  char *model);

#endif
