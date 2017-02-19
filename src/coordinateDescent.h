#include "depends.h"
#ifdef HAVE_GSL_HEADER

#include "commonFunctions.h"
#include "ReadInData.h"
#include "ridgeRegressionFunctions.h"

int coordinateDescentLogisticGenotypes(gsl_matrix_int * genotypes,
		gsl_vector_int * phenotypes,
		int intercept_flag,
		int standardize_flag,
		int unpenalized,
		GSL_TYPE(vector) * tau_vector,
		GSL_TYPE(vector) * means,
		GSL_TYPE(vector) * scales,
		GSL_TYPE(vector) * B,
		PREC epsilon);


int coordinateDescentLogistic(GSL_TYPE(vector) * B,
			      GSL_TYPE(matrix) * X,
			      gsl_vector_int * y,
			      GSL_TYPE(vector) * tau_vector,
			      int intercept_flag,
			      int unpenalized,
			      PREC epsilon);

PREC Fr(PREC r, PREC delta);

PREC computeUpdate(GSL_TYPE(vector) * X_column,
	GSL_TYPE(vector) * y,
	GSL_TYPE(vector) * rvector,
	PREC B_element,
	PREC delta,
	PREC tau,
	int unpen_flag);

int convergenceCheckLogistic(GSL_TYPE(vector) * deltar,
	GSL_TYPE(vector) * rvector,
	PREC epsilon);

int preparePhenotypesForCoordinateDescent(GSL_TYPE(vector) * y_cd, const GSL_TYPE(vector) * y);

int prepareShrinkageForCoordinateDescent(GSL_TYPE(vector) * shrinkage_cd, const GSL_TYPE(vector) * shrinkage);

PREC computeDofFLogistic(GSL_TYPE(matrix) * X,
			 GSL_TYPE(vector) * beta,
			 PREC k);

int coordinateDescentLinearFloat(GSL_TYPE(matrix) * Z,
				 GSL_TYPE(vector) * y,
				 GSL_TYPE(vector) * a,
				 PREC epsilon);

int updateYtilde(GSL_TYPE(vector) * ytilde,
		 GSL_TYPE(matrix) * Z,
		 GSL_TYPE(vector) * B,
		 int j);

int updateBetaLinear(GSL_TYPE(vector) * Bpen,
		     GSL_TYPE(matrix) * Z,
		     GSL_TYPE(vector) * y,
		     GSL_TYPE(vector) * ytilde,
		     int j,
		     PREC penalty);

int convergenceCheckLinear(GSL_TYPE(vector) * Bold,
			   GSL_TYPE(vector) * Bpen,
			   GSL_TYPE(vector) * B,
			   PREC epsilon);

int updateBetaLinearGenotypes(GSL_TYPE(vector) * Bpen,
			      gsl_matrix_int * X,
			      GSL_TYPE(vector) * means,
			      GSL_TYPE(vector) * scales,
			      GSL_TYPE(vector) * y,
			      GSL_TYPE(vector) * ytilde,
			      int j,
			      PREC penalty);

int updateYtildeGenotypes(GSL_TYPE(vector) * ytilde,
			  gsl_matrix_int * X,
			  GSL_TYPE(vector) * means,
			  GSL_TYPE(vector) * scales,
			  GSL_TYPE(vector) * B,
			  int j);

GSL_TYPE(vector) * getScaledColOfX(gsl_matrix_int * X,
				   GSL_TYPE(vector) * means,
				   GSL_TYPE(vector) * scales,
				   int i);


int coordinateDescentLinearGenotypes(gsl_matrix_int * X,
				     GSL_TYPE(vector) * y,
				     int intercept_flag,
				     int standardize_flag,
				     PREC lambda,
				     GSL_TYPE(vector) * means,
				     GSL_TYPE(vector) * scales,
				     GSL_TYPE(vector) * Bout,
				     PREC epsilon);
		 
int get_prev_variant_col(int current_pos, 
			 int number_of_columns);
#endif

