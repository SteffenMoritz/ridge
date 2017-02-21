#include "coordinateDescent.h"
#ifdef HAVE_GSL_HEADER

int coordinateDescentLogisticGenotypes(gsl_matrix_int * X,
		gsl_vector_int * y,
		int intercept_flag,
		int standardize_flag,
		int unpenalized,
		GSL_TYPE(vector) * tau_vector,
		GSL_TYPE(vector) * means,
		GSL_TYPE(vector) * scales,
		GSL_TYPE(vector) * B,
		PREC epsilon)
{
  int j = 0;
  int not_converged = 1;
  int unpen_flag = 0;
  PREC deltav = 0.0;
  int n = X->size1;
  int p = X->size2 + intercept_flag;
  GSL_TYPE(vector) * change_in_linear_scores = GSL_FUNCTION(vector,calloc)(n);
  GSL_TYPE(vector) * delta = GSL_FUNCTION(vector,alloc)(p);
  GSL_TYPE(vector) * deltar = GSL_FUNCTION(vector,calloc)(n);
  GSL_TYPE(vector) * rvector = GSL_FUNCTION(vector,calloc)(n);
  GSL_TYPE(vector) * X_column_prec = GSL_FUNCTION(vector,calloc)(n);
  GSL_TYPE(vector) * y_prec = GSL_FUNCTION(vector,alloc)(n);
  GSL_TYPE(vector) * Xdoty = GSL_FUNCTION(vector,alloc)(n);
  convert_int_vector(y, y_prec);
  // Initial values for delta
  GSL_FUNCTION(vector,set_all)(delta, 1.0);
  // Vector view to store X_column
  gsl_vector_int_view X_column;
  // PREC to store element of B
  PREC B_element = 0.0;
  // PREC to store current deltaB
  PREC deltaB_element = 0.0;
  // PREC to store element of tau
  PREC tau_element = 0.0;
  // PREC to store element of delta
  PREC deltaj = 0.0;
  // PREC to temporarily store mean and sd
  PREC mean = 0.0;
  PREC sd = 0.0;
  while(not_converged)
    {
      // Set deltar vector to all zeros
      GSL_FUNCTION(vector,set_all)(deltar, 0);
      for(j = 0; j < p; j++)
	{
	  if((j == 0 && intercept_flag) || unpenalized)
	    {
	      unpen_flag = 1;
	    } else
	    {
	      unpen_flag = 0;
	    }
	  // First part: set X_column_prec
	  // If we are not at the intercept 
	  if(!(j==0 && intercept_flag))
	    {
	      // Extract column of X
	      X_column = gsl_matrix_int_column(X, j - intercept_flag);
	      // Convert to PREC
	      convert_int_vector(&X_column.vector, X_column_prec);
	      // Scale it
	      if(intercept_flag)
		{
		  mean = GSL_FUNCTION(vector,get)(means, j - intercept_flag);
		  GSL_FUNCTION(vector, add_constant)(X_column_prec, -1*mean);
		}
	      sd = GSL_FUNCTION(vector,get)(scales, j - intercept_flag);
	      GSL_FUNCTION(vector, scale)(X_column_prec, 1 / sd);
	      // Else if we are at the intercept 
	    } else if (j == 0 && intercept_flag) {
	    GSL_FUNCTION(vector,set_all)(X_column_prec, 1.0);
	  }
	  // Next part
	  // 
	  // Compute Xdoty
	  GSL_FUNCTION(vector,memcpy)(Xdoty, X_column_prec);
	  GSL_FUNCTION(vector,mul)(Xdoty, y_prec);
	  // Extract element of B
	  B_element = GSL_FUNCTION(vector,get)(B,j);
	  // Extract element of tau
	  tau_element = GSL_FUNCTION(vector,get)(tau_vector, j);
	  // Extract element of delta
	  deltaj = GSL_FUNCTION(vector,get)(delta,j);
	  // Compute tentative step deltav
	  deltav = computeUpdate(X_column_prec,
				 y_prec,
				 rvector,
				 B_element,
				 deltaj,
				 tau_element,
				 unpen_flag);
	  deltaB_element = GSL_MIN( GSL_MAX(deltav, -1.0), deltaj);
	  GSL_FUNCTION(vector,memcpy)(change_in_linear_scores, Xdoty);
	  GSL_FUNCTION(vector,scale)(change_in_linear_scores, deltaB_element);
	  GSL_FUNCTION(vector,add)(deltar, change_in_linear_scores);
	  GSL_FUNCTION(vector,add)(rvector, change_in_linear_scores);
	  GSL_FUNCTION(vector,set)(B, j, B_element + deltaB_element);
	  GSL_FUNCTION(vector,set)(delta, j, GSL_MAX(2 * deltaB_element, deltaj / 2));
	}
      not_converged = convergenceCheckLogistic(deltar, rvector, epsilon);
    }
  GSL_FUNCTION(vector,free)(change_in_linear_scores);
  GSL_FUNCTION(vector,free)(X_column_prec);
  GSL_FUNCTION(vector,free)(y_prec);
  GSL_FUNCTION(vector,free)(deltar);
  GSL_FUNCTION(vector,free)(rvector);
  GSL_FUNCTION(vector,free)(delta);
  GSL_FUNCTION(vector,free)(Xdoty);
  return 0;
} // Ends function

int coordinateDescentLogistic(GSL_TYPE(vector) * B,
			      GSL_TYPE(matrix) * X,
			      gsl_vector_int * y,
			      GSL_TYPE(vector) * tau_vector,
			      int intercept_flag,
			      int unpenalized,
			      PREC epsilon)
{
	int not_converged = 1; // Slightly strange logical construction - so that the while loop continues until convergence
	int j = 0; // Counter for iterating over predictors
	int unpen_flag = 0; // Flag to indicate when not to penalize (used for the first predictor if intercept_flag = 1)
	PREC deltav = 0.0; // tentative step
	// y should already be in -1, 1 format before calling coordinateDescentLogistic
	// and tau_vector should already be in 1 / 2 * k format
	int n = X->size1;
	int p = X->size2;
	GSL_TYPE(vector) * change_in_linear_scores = GSL_FUNCTION(vector,calloc)(n);
	GSL_TYPE(vector) * delta = GSL_FUNCTION(vector,alloc)(p);
	GSL_TYPE(vector) * deltar = GSL_FUNCTION(vector,calloc)(n);
	GSL_TYPE(vector) * rvector = GSL_FUNCTION(vector,calloc)(n);
	GSL_TYPE(vector) * y_prec = GSL_FUNCTION(vector,alloc)(n);
	convert_int_vector(y, y_prec);
	// Initial values for delta
	GSL_FUNCTION(vector,set_all)(delta, 1.0);
	// Vector view to store X_column
	GSL_FUNCTION(vector,view) X_column;
	// PREC to store element of B
	PREC B_element = 0.0;
	// PREC to store current deltaB
	PREC deltaB_element = 0.0;
	// PREC to store element of tau
	PREC tau_element = 0.0;
	// PREC to store element of delta
	PREC deltaj = 0.0;
	// Make column-wise product of X[,j] and y for all columns of X
	GSL_TYPE(matrix) * Xdoty = GSL_FUNCTION(matrix,calloc)(n,p);
	GSL_FUNCTION(vector,view) Xdoty_column;
	for(j = 0; j < p; j++)
	{
		X_column = GSL_FUNCTION(matrix,column)(X, j);
		Xdoty_column = GSL_FUNCTION(matrix,column)(Xdoty, j);
		GSL_FUNCTION(vector,memcpy)(&Xdoty_column.vector, y_prec);
		GSL_FUNCTION(vector,mul)(&Xdoty_column.vector, &X_column.vector);
	}
	while(not_converged)
	{
	  // Set deltar vector to all zeros
	  GSL_FUNCTION(vector,set_all)(deltar, 0);
	  for(j = 0; j < p; j++)
	    {
	      if((j == 0 && intercept_flag) || unpenalized) {
		unpen_flag = 1;
	      } else {
		unpen_flag = 0;
	      }
	      // Extract column of X
	      X_column = GSL_FUNCTION(matrix,column)(X, j);
	      // Extract column of Xdoty
	      Xdoty_column = GSL_FUNCTION(matrix,column)(Xdoty, j);
	      // Extract element of B
	      B_element = GSL_FUNCTION(vector,get)(B,j);
	      // Extract element of tau
	      tau_element = GSL_FUNCTION(vector,get)(tau_vector, j);
	      // Extract element of delta
	      deltaj = GSL_FUNCTION(vector,get)(delta,j);
	      // Compute tentative step deltav
	      deltav = computeUpdate(&X_column.vector,
				     y_prec,
				     rvector,
				     B_element,
				     deltaj,
				     tau_element,
				     unpen_flag);
	      deltaB_element = GSL_MIN( GSL_MAX(deltav, -1.0), deltaj);
	      GSL_FUNCTION(vector,memcpy)(change_in_linear_scores, &Xdoty_column.vector);
	      GSL_FUNCTION(vector,scale)(change_in_linear_scores, deltaB_element);
	      GSL_FUNCTION(vector,add)(deltar, change_in_linear_scores);
	      GSL_FUNCTION(vector,add)(rvector, change_in_linear_scores);
	      GSL_FUNCTION(vector,set)(B, j, B_element + deltaB_element);
	      GSL_FUNCTION(vector,set)(delta, j, GSL_MAX(2 * deltaB_element, deltaj / 2));
	    }
	  not_converged = convergenceCheckLogistic(deltar, rvector, epsilon);
	}
	GSL_FUNCTION(vector,free)(change_in_linear_scores);
	GSL_FUNCTION(vector,free)(deltar);
	GSL_FUNCTION(vector,free)(rvector);
	GSL_FUNCTION(vector,free)(delta);
	GSL_FUNCTION(matrix,free)(Xdoty);
	GSL_FUNCTION(vector,free)(y_prec);
	return 0;
}

PREC Fr(PREC r, PREC delta)
{
  PREC res = 0.0;
  PREC criterion = 0.25;
  PREC absr = MATHS_FUNCTION(fabs)(r);
  if(absr < delta)
    {
      res = criterion;
    } else {
    res = 1 / (2 + MATHS_FUNCTION(exp)(absr - delta) + MATHS_FUNCTION(exp)(delta - absr));
  }
  return res;
}

PREC computeUpdate(GSL_TYPE(vector) * X_column,
	GSL_TYPE(vector) * y,
	GSL_TYPE(vector) * rvector,
	PREC B_element,
	PREC deltaj,
	PREC tau,
	int unpen_flag)
{
	int i = 0;
	int n = X_column->size;
	PREC tmp = 0.0;
	PREC X_element = 0.0;
	PREC deltav = 0.0; // Return value
	PREC one_plus_exp_rvector_element = 0.0;
	GSL_TYPE(vector) * numerator = GSL_FUNCTION(vector,calloc)(n);
	PREC numerator_sum  = 0.0;
	PREC denominator_sum = 0.0;
	GSL_TYPE(vector) * denominator = GSL_FUNCTION(vector,calloc)(n);
	for(i = 0; i < n; i++)
	{
		X_element = GSL_FUNCTION(vector,get)(X_column, i);
		one_plus_exp_rvector_element = 1 + MATHS_FUNCTION(exp)(GSL_FUNCTION(vector,get)(rvector,i));
		tmp = X_element * GSL_FUNCTION(vector,get)(y, i) / one_plus_exp_rvector_element;
		GSL_FUNCTION(vector,set)(numerator, i, tmp);
		tmp = MATHS_FUNCTION(pow)(X_element, 2) * Fr( GSL_FUNCTION(vector,get)(rvector, i) , deltaj * MATHS_FUNCTION(fabs)(X_element) ) ;
		GSL_FUNCTION(vector,set)(denominator, i, tmp);
	}
	sumVector(numerator,&numerator_sum);
	sumVector(denominator, &denominator_sum);
	if(!unpen_flag)
	{
		numerator_sum = numerator_sum - B_element / tau;
		denominator_sum = denominator_sum + 1 / tau;
	}
	deltav = numerator_sum / denominator_sum;
	GSL_FUNCTION(vector,free)(numerator);
	GSL_FUNCTION(vector,free)(denominator);
	return deltav;
}

int convergenceCheckLogistic(GSL_TYPE(vector) * deltar,
	GSL_TYPE(vector) * rvector,
	PREC epsilon)
{
	int not_converged = 1; // Not converged is true (bit of a strange logical setup here but follows with rest of algorithm)
	PREC test = 0.0;
	sumVector(deltar, &test);
	test = MATHS_FUNCTION(fabs)(test);
	PREC sumRvector = 0.0;
	sumVector(rvector, &sumRvector);
	sumRvector = MATHS_FUNCTION(fabs)(sumRvector);
	test = test / (1 + sumRvector);
	if(test <= epsilon)
	{
		not_converged = 0; // i.e. convergence has been reached
	}
	return not_converged;
}


int preparePhenotypesForCoordinateDescent(GSL_TYPE(vector) * y_cd, const GSL_TYPE(vector) * y)
{
  // y <- y - (1 - y) = y <- y - 1 + y = y <- 2y - 1
  GSL_FUNCTION(vector,memcpy)(y_cd, y);
  GSL_FUNCTION(vector,scale)(y_cd, 2);
  GSL_FUNCTION(vector,add_constant)(y_cd, -1);
  return 0;
}

int prepareShrinkageForCoordinateDescent(GSL_TYPE(vector) * shrinkage_cd, const GSL_TYPE(vector) * shrinkage)
{
  // shrinkage_cd = 1 / 2 * shrinkage
  GSL_FUNCTION(vector,set_all)(shrinkage_cd, 0.5);
  GSL_FUNCTION(vector,div)(shrinkage_cd, shrinkage);
  return 0;
}

PREC computeDofFLogistic(GSL_TYPE(matrix) * X,
			 GSL_TYPE(vector) * beta,
			 PREC k)
{
  const size_t NINDIV = X->size1;
  const size_t NPRED = X->size2;
  PREC DofF = 0.0;
  GSL_TYPE(vector) * XB = GSL_FUNCTION(vector,calloc)(NINDIV);
  GSL_TYPE(vector) * p = GSL_FUNCTION(vector,calloc)(NINDIV);
  compute_XB_and_p(X, beta, XB, p);
  GSL_FUNCTION(vector,free)(XB);
  GSL_TYPE(vector) * negP = GSL_FUNCTION(vector,alloc)(NINDIV);
  GSL_TYPE(vector) * PnegP = GSL_FUNCTION(vector,alloc)(NINDIV);
  GSL_FUNCTION(vector,memcpy)(negP, p);
  GSL_FUNCTION(vector,scale)(negP, -1.0);
  GSL_FUNCTION(vector,add_constant)(negP, 1.0);
  GSL_FUNCTION(vector,memcpy)(PnegP, p);
  GSL_FUNCTION(vector,mul)(PnegP, negP);
  GSL_FUNCTION(vector,free)(negP);
  GSL_FUNCTION(vector,free)(p);
  // make the matrix tX
  GSL_TYPE(matrix) * tX = GSL_FUNCTION(matrix,alloc)(NPRED,NINDIV);
  GSL_FUNCTION(matrix,transpose_memcpy)(tX,X);
  // make the constant Ainv
  PREC Ainv = 1 / (2 * k);
  // make the vector Dinv
  GSL_TYPE(vector) * Dinv = GSL_FUNCTION(vector,calloc)(NINDIV);
  GSL_FUNCTION(vector,set_all)(Dinv, 1.0);
  GSL_FUNCTION(vector,div)(Dinv,PnegP);
  // Matrix to store the inverted matrix
  GSL_TYPE(matrix) * inv = GSL_FUNCTION(matrix,calloc)(NPRED,NPRED);
  invert_sum_of_matrices(Ainv,
			 tX,
			 Dinv,
			 X,
			 inv);
  GSL_FUNCTION(matrix,free)(tX);
  GSL_TYPE(matrix) * W = GSL_FUNCTION(matrix,calloc)(NINDIV,NINDIV);
  GSL_FUNCTION(matrix,set_identity)(W);
  GSL_FUNCTION(vector,view) diagW = GSL_FUNCTION(matrix,diagonal)(W);
  GSL_FUNCTION(vector,mul)(&diagW.vector, PnegP);
  GSL_FUNCTION(vector,free)(PnegP);
  GSL_FUNCTION(vector,free)(Dinv);
  // Make tXW
  GSL_TYPE(matrix) * WX = GSL_FUNCTION(matrix,calloc)(W->size1,X->size2);
  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, 1.0, W, X, 0.0, WX);
  GSL_FUNCTION(matrix,free)(W);
  // Make XWinvtXWX
  GSL_TYPE(matrix) * WXinvtXWX = GSL_FUNCTION(matrix,calloc)(WX->size1,inv->size2);
  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, 1.0, WX, inv, 0.0, WXinvtXWX);
  GSL_FUNCTION(matrix,free)(WX);
  GSL_FUNCTION(matrix,free)(inv);
  // Make H
  GSL_TYPE(matrix) * H = GSL_FUNCTION(matrix,calloc)(NINDIV,NINDIV);
  BLAS_FUNCTION(gemm)(CblasNoTrans,CblasTrans,1.0,WXinvtXWX,X,0.0,H);
  GSL_FUNCTION(matrix,free)(WXinvtXWX);
  // Make HH'
  GSL_TYPE(matrix) * HH = GSL_FUNCTION(matrix,calloc)(NINDIV,NINDIV);
  BLAS_FUNCTION(gemm)(CblasNoTrans,CblasTrans,1.0,H,H,0.0,HH);
  GSL_FUNCTION(matrix,free)(H);
  // Get diagonal
  GSL_FUNCTION(vector,view) diagHH = GSL_FUNCTION(matrix,diagonal)(HH);
  sumVector(&diagHH.vector,&DofF);
  GSL_FUNCTION(matrix,free)(HH);
  return DofF;
}

int coordinateDescentLinearFloat(GSL_TYPE(matrix) * Z,
				 GSL_TYPE(vector) * y,
				 GSL_TYPE(vector) * a,
				 PREC epsilon)
{
  /* j is a counter */
  int j = 0;
  /* Z are already scaled */
  /* y are already scaled */
  /* get dimensions */
  int n = Z->size1;
  int p = Z->size2;
  /* Counter to iterate over predictors */
  /* Initial estimates of B */
  GSL_TYPE(vector) * B = GSL_FUNCTION(vector, calloc)(p);
  GSL_TYPE(vector) * Bpen = GSL_FUNCTION(vector, calloc)(p);
  GSL_TYPE(vector) * Bold = GSL_FUNCTION(vector, calloc)(p);
  GSL_TYPE(vector) * ytilde = GSL_FUNCTION(vector, calloc)(n);
  /* Flag to indicate we have not converged yet */
  int not_converged = 1;
  while(not_converged)
    {
      GSL_FUNCTION(vector, set_all)(Bpen, 0);
      GSL_FUNCTION(vector, memcpy)(Bold, B);
      for(j = 0; j < p; j++)
	{
	  /* Compute ytilde(j) using existing B: */
	  updateYtilde(ytilde,
		       Z,
		       B,
		       j);
	  /* Compute penalized Bpen[j] */
	  updateBetaLinear(Bpen,
			   Z,
			   y,
			   ytilde,
			   j,
			   0);
	  /* Update B[j] */
	  GSL_FUNCTION(vector, set)(B, j, GSL_FUNCTION(vector, get)(Bpen, j));
	}
      /* Check for convergence */
      /* Also puts Bpen in B */
      not_converged = convergenceCheckLinear(Bold, Bpen, B, epsilon);
    }
  GSL_FUNCTION(vector, memcpy)(a, B);
  GSL_FUNCTION(vector, free)(B);
  GSL_FUNCTION(vector, free)(Bpen);
  GSL_FUNCTION(vector, free)(Bold);
  GSL_FUNCTION(vector, free)(ytilde);
  return 0;
}

int updateYtilde(GSL_TYPE(vector) * ytilde,
		 GSL_TYPE(matrix) * Z,
		 GSL_TYPE(vector) * B,
		 int j)
{
  int p =  B->size;
  int n = ytilde->size;
  GSL_FUNCTION(vector, set_all)(ytilde, 0);
  int l = 0;
  PREC B_element;
  GSL_TYPE(vector) * X_column = GSL_FUNCTION(vector, calloc)(n);
  for(l = 0; l < p; l++)
    {
      if(l != j)
	{
	  /* Extract column of X */
	  GSL_FUNCTION(matrix, get_col)(X_column, Z, l);
	  /* Extract element of B */
	  B_element = GSL_FUNCTION(vector, get)(B, l);
	  /* Multiply it by B */
	  GSL_FUNCTION(vector, scale)(X_column, B_element);
	  /* Add it to ytilde */
	  GSL_FUNCTION(vector,add)(ytilde, X_column);
	}
    }
  GSL_FUNCTION(vector, free)(X_column);
  return 0;
}

int updateBetaLinear(GSL_TYPE(vector) * Bpen,
		     GSL_TYPE(matrix) * Z,
		     GSL_TYPE(vector) * y,
		     GSL_TYPE(vector) * ytilde,
		     int j,
		     PREC penalty)
{
  /* Number of individuals */
  int n = y->size;
  /* Beta element */
  PREC B_element = 0;
  /* Compute partial residual for fitting Bj */
  GSL_TYPE(vector) * partial_residual = GSL_FUNCTION(vector, calloc)(n);
  GSL_FUNCTION(vector, memcpy)(partial_residual, y);
  GSL_FUNCTION(vector, scale)(ytilde, -1);
  GSL_FUNCTION(vector, add)(partial_residual, ytilde);
  /* Extract column of X */
  GSL_TYPE(vector) * X_column = GSL_FUNCTION(vector, calloc)(n);
  GSL_FUNCTION(matrix, get_col)(X_column, Z, j);
  /* When working on short genotypes, the part for scaling the short genotypes goes here */
  /**/
  /**/
  /* Multiply X_column by partial_residual, store the result in B_element */
  GSL_BLAS_FUNCTION(dot)(X_column, partial_residual, &B_element);
  /* If shrinking, do it here */
  B_element = B_element / (1 + penalty);
  /* Put B_element into Bpen */
  GSL_FUNCTION(vector, set)(Bpen, j, B_element);
  GSL_FUNCTION(vector, free)(partial_residual);
  GSL_FUNCTION(vector, free)(X_column);
  return 0;
}

int convergenceCheckLinear(GSL_TYPE(vector) * Bold,
			   GSL_TYPE(vector) * Bpen,
			   GSL_TYPE(vector) * B,
			   PREC epsilon)
{
  /* Flag to return */
  int not_converged = 1;
  /* Number of predictors */
  int p = B->size;
  /* Counter for iterating */
  int j = 0;
  /* Temporary PREC for element of B_diff */
  PREC B_diff_element = 0;
  /* Allocate vector of differences between betas */
  GSL_TYPE(vector) * B_diff = GSL_FUNCTION(vector, calloc)(p);
  /* Copy Bold to it */
  GSL_FUNCTION(vector, memcpy)(B_diff, Bold);
  /* Subtract Bpen */
  GSL_FUNCTION(vector, sub)(B_diff, Bpen);
  /* Vector to keep track of whether B has converged or not */
  gsl_vector_int * track_convergence = gsl_vector_int_calloc(p);
  /* Check each element of B_diff */
  for(j = 0; j < p; j++)
    {
      B_diff_element = GSL_FUNCTION(vector, get)(B_diff, j);
      B_diff_element = MATHS_FUNCTION(fabs)(B_diff_element);
      if(B_diff_element > epsilon)
	{
	  gsl_vector_int_set(track_convergence, j, 1);
	}
    }
  /* sum track_convergence */
  int number_not_converged = 0;
  number_not_converged = sumIntVec(track_convergence);
  /* Free the track_convergence vector */
  gsl_vector_int_free(track_convergence);
  /* Free B_diff vector */
  GSL_FUNCTION(vector,free)(B_diff);
  if(number_not_converged == 0)
    {
      not_converged = 0;
    }
  /* Put Bpen in B */
  GSL_FUNCTION(vector, memcpy)(B, Bpen);
  return not_converged;
}

int coordinateDescentLinearGenotypes(gsl_matrix_int * X,
				     GSL_TYPE(vector) * y,
				     int intercept_flag,
				     int standardize_flag,
				     PREC lambda,
				     GSL_TYPE(vector) * means,
				     GSL_TYPE(vector) * scales,
				     GSL_TYPE(vector) * Bout,
				     PREC epsilon)
{
  /* j is a counter */
  int j = 0;
  /* y are already scaled */
  /* get dimensions */
  int n = X->size1;
  int p = X->size2;
  /* Counter to iterate over predictors */
  /* Initial estimates of B */
  GSL_TYPE(vector) * B = GSL_FUNCTION(vector, calloc)(p);
  GSL_TYPE(vector) * Bpen = GSL_FUNCTION(vector, calloc)(p);
  GSL_TYPE(vector) * Bold = GSL_FUNCTION(vector, calloc)(p);
  GSL_TYPE(vector) * ytilde = GSL_FUNCTION(vector, calloc)(n);
  /* Flag to indicate we have not converged yet */
  int not_converged = 1;
  while(not_converged)
    {
      GSL_FUNCTION(vector, set_all)(Bpen, 0);
      GSL_FUNCTION(vector, memcpy)(Bold, B);
      for(j = 0; j < p; j++)
	{
	  /* Compute ytilde(j) using existing B: */
	  updateYtildeGenotypes(ytilde,
				X,
				means,
				scales,
				B,
				j);
	  /* Compute penalized Bpen[j] */
	  updateBetaLinearGenotypes(Bpen,
				    X,
				    means,
				    scales,
				    y,
				    ytilde,
				    j,
				    lambda);
	  /* Update B[j] */
	  GSL_FUNCTION(vector, set)(B, j, GSL_FUNCTION(vector, get)(Bpen, j));
	}
      /* Check for convergence */
      /* Also puts Bpen in B */
      not_converged = convergenceCheckLinear(Bold, Bpen, B, epsilon);
    }
  GSL_FUNCTION(vector, memcpy)(Bout, B);
  GSL_FUNCTION(vector, free)(B);
  GSL_FUNCTION(vector, free)(Bpen);
  GSL_FUNCTION(vector, free)(Bold);
  GSL_FUNCTION(vector, free)(ytilde);
  return 0;
}

GSL_TYPE(vector) * getScaledColOfX(gsl_matrix_int * X,
				   GSL_TYPE(vector) * means,
				   GSL_TYPE(vector) * scales,
				   int i)
{
  /* Vector to return */
  GSL_TYPE(vector) * scaledX = GSL_FUNCTION(vector, calloc)(X->size1);
  /* Vector view of the column of X that we wish to scale and return */
  gsl_vector_int_view Xcolumn = gsl_matrix_int_column(X, i);
  /* Temporary variable to store the mean */
  PREC mean = 0;
  /* Temporary variable to store the scale */
  PREC scale = 0;
  /* Get the mean */
  mean = GSL_FUNCTION(vector, get)(means, i);
  /* Get the scale */
  scale = GSL_FUNCTION(vector, get)(scales, i);
  /* Copy the vector view into the vector to return */
  convert_int_vector(&Xcolumn.vector, scaledX);
  /* Subtract the mean */
  GSL_FUNCTION(vector, add_constant)(scaledX, -1 * mean);
  /* Divide by the scale */
  GSL_FUNCTION(vector, scale)(scaledX, 1 / scale);
  /* Return the scaled vector */
  return scaledX;
 }

int updateBetaLinearGenotypes(GSL_TYPE(vector) * Bpen,
			      gsl_matrix_int * X,
			      GSL_TYPE(vector) * means,
			      GSL_TYPE(vector) * scales,
			      GSL_TYPE(vector) * y,
			      GSL_TYPE(vector) * ytilde,
			      int j,
			      PREC penalty)
{
  /* Number of individuals */
  int n = y->size;
  /* Beta element */
  PREC B_element = 0;
  /* Compute partial residual for fitting Bj */
  GSL_TYPE(vector) * partial_residual = GSL_FUNCTION(vector, calloc)(n);
  GSL_FUNCTION(vector, memcpy)(partial_residual, y);
  GSL_FUNCTION(vector, sub)(partial_residual, ytilde);
  /* Extract column of X and scale it */
  GSL_TYPE(vector) * X_column = getScaledColOfX(X, 
						means, 
						scales,
						j);
  /* Multiply X_column by partial_residual, store the result in B_element */
  GSL_BLAS_FUNCTION(dot)(X_column, partial_residual, &B_element);
  /* If shrinking, do it here */
  B_element = B_element / (1 + penalty);
  /* Put B_element into Bpen */
  GSL_FUNCTION(vector, set)(Bpen, j, B_element);
  GSL_FUNCTION(vector, free)(partial_residual);
  GSL_FUNCTION(vector, free)(X_column);
  return 0;
}

/* This version updates only the parts of ytilde that we need to update */
int updateYtildeGenotypes(GSL_TYPE(vector) * ytilde,
			  gsl_matrix_int * X,
			  GSL_TYPE(vector) * means,
			  GSL_TYPE(vector) * scales,
			  GSL_TYPE(vector) * B,
			  int j)
{
  int colToAdd = 0;
  PREC B_element = 0.0;
  /* Extract column of X, scaled */
  GSL_TYPE(vector) * scaledX = getScaledColOfX(X, 
					       means,
					       scales,
					       j);
  /* Extract element of B */
  B_element = GSL_FUNCTION(vector, get)(B, j);
  /* Multiply scaledX by B_element */
  GSL_FUNCTION(vector, scale)(scaledX, B_element);
  /* Subtract it from ytilde */
  GSL_FUNCTION(vector, sub)(ytilde, scaledX);
  /* Free the vector scaledX */
  GSL_FUNCTION(vector, free)(scaledX);
  /* Get the column of which we should add the scaledX * B_element */
  colToAdd = get_prev_variant_col(j, means->size);
  /* Extract the relevant column of X */
  scaledX = getScaledColOfX(X,
			    means,
			    scales,
			    colToAdd);
  /* Extract the relevant element of B */
  B_element = GSL_FUNCTION(vector, get)(B, colToAdd);
  /* Multiply scaledX by B_element */
  GSL_FUNCTION(vector, scale)(scaledX, B_element);
  /* Add to ytilde */
  GSL_FUNCTION(vector, add)(ytilde, scaledX);
  GSL_FUNCTION(vector, free)(scaledX);
  return 0;
}

int get_prev_variant_col(int current_pos, int number_of_columns)
{
  /* Function to get the previous variant column */
  int new_pos = 0;
  if(current_pos == 0)
    {
      new_pos = number_of_columns - 1;
    } else {
    new_pos = current_pos - 1;
  }
  return new_pos;
}
#endif

typedef int make_iso_compilers_happy;

