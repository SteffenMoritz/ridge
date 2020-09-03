
/* Header file for ridgeRegressionFunctions.c */
#include "depends.h"
#ifdef HAVE_GSL_HEADER

/* includes */
#if _CUDA_
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <cula_blas.h>
#include <cula.h>
#include "cudaOnlyFunctions.h"
#endif

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <math.h>

#include "computeLinearRidge.h"

#include "ReadInData.h"

/* Compute linear ridge coefficeints */
int computeLinearRidge(GSL_TYPE(vector) * ahat, 
		       GSL_TYPE(vector) * B, 
		       GSL_TYPE(vector) * D2, 
		       GSL_TYPE(matrix) * V, 
		       PREC lambda);

/* linear generalized ridge regression */
GSL_TYPE(vector) * computeLinearGeneralizedRidge(GSL_TYPE(vector) * beta, 
					   GSL_TYPE(matrix) * pred, 
					   GSL_TYPE(vector) * pheno, 
					   GSL_TYPE(vector) * shrinkage, 
					   int intercept_flag);



int computeLogisticRidge(GSL_TYPE(vector) * beta,
				GSL_TYPE(matrix) * pred,
				GSL_TYPE(vector) * pheno,
				GSL_TYPE(vector) * shrinkage,
				int intercept_flat,
				int Doff_flat,
				PREC * DofF);

PREC objectiveFunction(GSL_TYPE(vector) * beta,
				GSL_TYPE(matrix) * X,
				GSL_TYPE(vector) * pheno,
				GSL_TYPE(vector) * shrinkage,
				int intercept_flat);

int updateBeta(GSL_TYPE(vector) * beta,
			GSL_TYPE(matrix) * X,
			GSL_TYPE(vector) * pheno,
			GSL_TYPE(matrix) * kI,
			int intercept_flag,
			int DofF_flag,
			GSL_TYPE(matrix) * invtXWX_return,
			GSL_TYPE(matrix) * W_return);

int getProb(GSL_TYPE(vector) * p, GSL_TYPE(vector) * XB);

/* NB my_gsl_solve only works for matrices of type double 
  due to the linalg functions only having been written for this type */
int my_gsl_solve(gsl_matrix * X,
	gsl_matrix * solvedX);

int compute_XB_and_p(GSL_TYPE(matrix) * X,
	GSL_TYPE(vector) * B,
	GSL_TYPE(vector) * XB,
	GSL_TYPE(vector) * p);

int chooseHowManyK(GSL_TYPE(vector) * D);

int returnToOriginalScaleLinear(GSL_TYPE(vector) * betaOut, 
				GSL_TYPE(vector) * Bridge, 
				GSL_TYPE(vector) * means,
				GSL_TYPE(vector) * scales,
				PREC y_mean,
				int intercept_flag);

/* return to original scale - generalized linear ridge regression */
int returnToOriginalScaleGenLinear(GSL_TYPE(vector) * Bridge, 
					    GSL_TYPE(vector) * betaOut, 
					    GSL_TYPE(matrix) * pred, 
					    GSL_TYPE(vector) * pheno, 
					    GSL_TYPE(vector) * scales, 
					    int intercept_flag);



/* Convert trans */
char setTrans(CBLAS_TRANSPOSE_t Trans);

/* prepare for linear ridge 
	used when we are going to call linear rr on a range of different shrinkage parameters */
int prepareForLinearRidge(GSL_TYPE(matrix) * X, 
			  GSL_TYPE(vector) * y, 
			  GSL_TYPE(matrix) * U, 
			  GSL_TYPE(matrix) * V, 
			  GSL_TYPE(vector) * D, 
			  GSL_TYPE(vector) * D2, 
			  GSL_TYPE(matrix) * Z,
			  GSL_TYPE(vector) * ahat);


/* compute DofF for ridge model */
int computeDofF(GSL_TYPE(vector) * D2,
		PREC Kr,
		PREC * DofF);

/* compute (A + BDC')^(-1) */ 
int invert_sum_of_matrices(const PREC Ainv,
			   const GSL_TYPE(matrix) * B,
			   const GSL_TYPE(vector) * Dinv,
			   const GSL_TYPE(matrix) * tC,
			   GSL_TYPE(matrix) * out);


/* ordinary regression */
gsl_vector * my_gsl_linear_fit(gsl_matrix * X, 
			       gsl_vector * y, 
			       int NROW, 
			       int NCOL);

#endif

