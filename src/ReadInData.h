#include "depends.h"
#ifdef HAVE_GSL_HEADER
/* ReadInData.h is the header file for the source file ReadInData.c */

/*This file defines the functions writeGenotypes and writePhenotypes*/

/* scale the ys */
int scaley(GSL_TYPE(vector) * y, PREC * y_mean);

/*get number of individuals by reading phenotypes file*/
int getNROW(FILE *fp);

/*function prototypes - writing*/
void writePhenotypes(int NINDIV);

/* Append char to string */
int appendToString(char * ptr, int *currentpos, char charToAppend);

/* function prototype - prepareShrinkage */
/* The vector shrinkage is intercept_flag + NSNPS + NCOVAR long */
/* Not used in R package */
/* int prepareShrinkage(char * model, */
/* 		     PREC * lambda,  */
/* 		     PREC * lambda_c,  */
/* 		     char * lambdafilename,  */
/* 		     char * lambdacovarfilename,  */
/* 		     int NSNPS,  */
/* 		     int NCOVAR,  */
/* 		     int intercept_flag,  */
/* 		     GSL_TYPE(vector) * shrinkage,  */
/* 		     int * automaticK, */
/* 		     int * singleK); */

/* function prototype - concatenate two vectors */
gsl_vector * concatenateTwoVectors(gsl_vector * result,
	gsl_vector * vec1,
	gsl_vector * vec2);

/* function prototype - print matrix */
void printMatrix(GSL_TYPE(matrix) * mat);

/* function prototype - print matrix */
void printMatrixTen(GSL_TYPE(matrix) * mat);

/* function prototype - print matrix */
void printIntMatrix(gsl_matrix_int * mat);

/* function prototype - print matrix */
void printIntMatrixTen(gsl_matrix_int * mat);

/* function prototype - print vector */
void printVector(GSL_TYPE(vector) * Vec);

/* function prototype - print vector of ints */
void printIntVector(gsl_vector_int * Vec);

/* function prototype - print vector */
void printVectorTen(GSL_TYPE(vector) * Vec);

/* function prototype - print vector of ints - first ten elements */
void printIntVectorTen(gsl_vector_int * Vec);

/* function prototype - compute sum of vector */
int sumVector(GSL_TYPE(vector) * vector, PREC * sum);

/* function prototype - compute sum of vector - double precision - for coordinateDescent */
int sumVectorDouble(gsl_vector * vector, double * sum);

/* check Model */
int checkModel(char * model);

/* function prototype - prepare matrix of predictors */
GSL_TYPE(matrix) * preparePredictors(int NINDIV, 
				     int NSNPS, 
				     char ** SNPnames,
				     int NCOVAR, 
				     char ** COVARnames,
				     char * genofilename, 
				     char * covarfilename, 
				     int intercept_flag, 
				     int standardize_flag, 
				     int standardize_c_flag, 
				     GSL_TYPE(vector) * means, 
				     GSL_TYPE(vector) * scales,
				     int automaticK);

/* function prototype - prepare matrix of genotypes */
/* Not used in R package */
/* int prepareGenotypes(gsl_matrix * genotypes, int NINDIV, int NSNPS, char * genofilename); */

/* Get header row of file */
char ** getHeaderRow(char * filename, int *N);

/* Check file exists and can be opened for reading */
void checkFileForReading(char * filename);

/* Function prototype - get data from file into matrix 
   excluding header row */
GSL_TYPE(matrix) * getDataWithoutHeaderRow(char * filename, int NROW, int NCOL);

/* Check for invariant SNPs and covariates, standardize */
int checkForInvAndStandardize(GSL_TYPE(matrix) * mat, 
			      int START, 
			      int END, 
			      int standardize_flag,
			      int corr_form_flag,
			      GSL_TYPE(vector) * means, 
			      GSL_TYPE(vector) * sds,
			      char ** names);

int sumIntVec(gsl_vector_int * vec);

int writeOut(int intercept_flag,
	     int NSNPS,
	     int NCOVAR,
	     char ** SNPnames,
	     char ** COVARnames,
	     char * betafilename,
	     GSL_TYPE(vector) * betaOut);

int printBeta(char * name, PREC beta, FILE * file);

/* Functions to safely free matrices and vectors */

int safelyFreeVector(GSL_TYPE(vector) * vec);

int safelyFreeMatrix(GSL_TYPE(matrix) * mat);

/* Print opening blurb */
/* Function not used in R package */
// void printOpening(void,);

/* scaleX not used in R package */
/* int scaleX(gsl_matrix * X, int n2, int NSNPS); */
#endif
