
#include "depends.h"
#ifdef HAVE_GSL_HEADER
#include "ReadInData.h" // Required for getNROW in readCoefficients

gsl_matrix_int * readShortGenotypes(char * genofilename,
	int NINDIV,
	int NSNPS);

GSL_TYPE(matrix) * readGenotypes(char * genofilename,
	int NINDIV,
	int NSNPS);

GSL_TYPE(vector) * readCoefficients(char * betafilename,
	int * intercept_flag,
	PREC * intercept_coefficient);

int getGenotypeInfo(gsl_matrix_int * genotypes,
	int standardize_flag,
	int corr_form_flag,
	GSL_TYPE(vector) * means,
	GSL_TYPE(vector) * scales,
	char ** names);

int convert_int_vector(const gsl_vector_int * src, GSL_TYPE(vector)
 * dest);

int checkOperationType(PREC lambda,
		       PREC lambda_c,
		       char * lambdafilename,
		       char * lambdacovarfilename,
		       char * approxfilename,
		       int howManyK,
		       int individualK,
		       int * automaticK,
		       int * singleK,
		       int predict_flag);

int checkGenotypes(gsl_matrix_int * mat);

int checkForInvariantPredictors(char * genofilename, 
				int NINDIV);

#endif
