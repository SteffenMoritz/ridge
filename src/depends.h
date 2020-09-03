/* Pre-requisites for the package regression */
#include <stdlib.h>
#include "config.h"

#ifdef HAVE_GSL_HEADER

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_vector.h>

#include "R.h"

#if _CUDA_
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <cula_blas.h>
#include <cula.h>
#endif

#if _CUDA_
	#define PREC float
	#define PREC_EPS FLT_EPSILON
	#define GSL_TYPE(type) gsl_ ## type ## _float
	#define	GSL_FUNCTION(type,name) gsl_ ## type ## _float_ ## name
	#define GSL_STATS_FUNCTION(name) gsl_stats ## _float_ ## name
	#define MATHS_FUNCTION(name) name ## f
	#define SVD_FUNCTION svdAnyMatCuda
	#define BLAS_FUNCTION(name) my_cula_s ## name
	#define MY_FUNCTION(name) my_cula_ ## name
	#define UGAUSSIAN_FUNCTION my_ugaussian_function
	#define GSL_BLAS_FUNCTION(name) gsl_blas_s ## name
	#define PREPARE_FUNCTION(name) prepare ## name ## ForCoordinateDescentCuda
	#define PREC_DIFF -4
#else
	#define PREC double
	#define PREC_EPS DBL_EPSILON
	#define GSL_TYPE(type) gsl_ ## type
	#define	GSL_FUNCTION(type,name) gsl_ ## type ## _ ## name
	#define GSL_STATS_FUNCTION(name) gsl_stats ## _ ## name
	#define MATHS_FUNCTION(name) name
	#define SVD_FUNCTION svdAnyMat
	#define BLAS_FUNCTION(name) gsl_blas_d ## name
	#define MY_FUNCTION(name) my_gsl_ ## name
	#define UGAUSSIAN_FUNCTION gsl_cdf_ugaussian_P
	#define GSL_BLAS_FUNCTION(name) gsl_blas_d ## name
	#define PREPARE_FUNCTION(name) prepare ## name ## ForCoordinateDescent
	#define PREC_DIFF -6
#endif

#endif
