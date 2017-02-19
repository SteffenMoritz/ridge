/* computeLinearRidge.h is the header file for the source file computeLinearRidge.c */

/* This file defines the function computeLinearRidge, computeLinearGeneralizedRidge, computeLinGenExpRidge */

/* includes */
#include "depends.h"
#ifdef HAVE_GSL_HEADER

/* SVD of any matrix */
int svdAnyMat(gsl_matrix * X, 
		 gsl_matrix * U, 
		 gsl_matrix * V, 
		 gsl_vector * D);

/* Prepare the lambdas */
int prepareLambdas(gsl_vector * y, 
		   gsl_matrix * U, 
		   gsl_vector * D2, 
		   gsl_vector * lambdaVeckHKB, 
		   char * skhkbfilename, 
		   char * sklwfilename, 
		   gsl_vector * lambdaVeckLW, 
		   int randomized, 
		   int s);

/* compute pvals linear ridge regression */
void computeLinearPvalsApprox(gsl_matrix * X, 
			      gsl_vector * B, 
			      gsl_vector * y, 
			      double lambda, 
			      char * pvalsfilename);

/* compute pvals linear ridge regression - Malo's method */
void computeLinearPvaslMalo(gsl_matrix * X, 
			    gsl_vector * B, 
			    gsl_vector * y, 
			    double lambda, 
			    char * pvalsfilename);

/* compute pvals linear ridge regression - our method */
void computeLinearPvalsOurs(gsl_matrix * X, 
			    gsl_vector * B, 
			    gsl_vector * y, 
			    double lambda, 
			    char * pvalsfilename);

#endif

