/* ridgeRegressionFunctions.c */

/* includes */
#include "ridgeRegressionFunctions.h"
#ifdef HAVE_GSL_HEADER

/* Compute Linear Ridge */
int computeLinearRidge(GSL_TYPE(vector) * ahat, 
		       GSL_TYPE(vector) * B, 
		       GSL_TYPE(vector) * D2, 
		       GSL_TYPE(matrix) * V, 
		       PREC lambda)
{
  /*
     In this function we want to compute
     aridge = ahat * (D2 / (D2 + lambda)
     then return to the original orientation
  */
  GSL_TYPE(vector) * aridge = GSL_FUNCTION(vector,alloc)(ahat->size);
  GSL_TYPE(vector) * vec_to_multiply_ahat_by = GSL_FUNCTION(vector,calloc)(D2->size);
  GSL_TYPE(vector) * denominator = GSL_FUNCTION(vector,calloc)(D2->size);
  GSL_FUNCTION(vector,memcpy)(vec_to_multiply_ahat_by, D2);
  GSL_FUNCTION(vector,memcpy)(denominator, D2);
  GSL_FUNCTION(vector,add_constant)(denominator, lambda);
  GSL_FUNCTION(vector,div)(vec_to_multiply_ahat_by, denominator);
  GSL_FUNCTION(vector,free)(denominator);
  GSL_FUNCTION(vector,memcpy)(aridge, ahat);
  GSL_FUNCTION(vector,mul)(aridge, vec_to_multiply_ahat_by);
  BLAS_FUNCTION(gemv)(CblasNoTrans, 1.0, V, aridge, 0.0, B);
  GSL_FUNCTION(vector,free)(aridge);
  GSL_FUNCTION(vector,free)(vec_to_multiply_ahat_by);
  return(0);
}

/* computeLinearGeneralizedRidge  */
GSL_TYPE(vector) * computeLinearGeneralizedRidge(GSL_TYPE(vector) * beta, 
					   GSL_TYPE(matrix) * pred, 
					   GSL_TYPE(vector) * pheno, 
					   GSL_TYPE(vector) * shrinkage, 
					   int intercept_flag)
{
  int NINDIV = pred->size1;
  int NCOVAR = pred->size2;
  int i;
  GSL_TYPE(matrix) * X;
  GSL_TYPE(vector) * shrinkage_to_use;
  if(intercept_flag)
    {
      NCOVAR = NCOVAR - 1;
      X = GSL_FUNCTION(matrix,alloc)(NINDIV, NCOVAR);
      shrinkage_to_use = GSL_FUNCTION(vector,calloc)(NCOVAR);
      // make a submatrix as a matrix view
      GSL_FUNCTION(matrix,view) x = GSL_FUNCTION(matrix,submatrix)(pred, 0, 1, NINDIV, NCOVAR);
      GSL_FUNCTION(matrix,memcpy)(X, &x.matrix);
      // make a subvector as a vector view
      GSL_FUNCTION(vector,view) shrinkage_vector_view = GSL_FUNCTION(vector,subvector)(shrinkage, 1, NCOVAR);
      GSL_FUNCTION(vector,memcpy)(shrinkage_to_use, &shrinkage_vector_view.vector);
    }
  else
    {
      // X is an exact copy of pred
      X = GSL_FUNCTION(matrix,alloc)(NINDIV, NCOVAR);
      GSL_FUNCTION(matrix,memcpy)(X, pred);
      // shrinkage_to_use is an exact copy of shrinkage
      shrinkage_to_use = GSL_FUNCTION(vector,alloc)(NCOVAR);
      GSL_FUNCTION(vector,memcpy)(shrinkage_to_use, shrinkage);
    }
  // Compute X'X
  GSL_TYPE(matrix) * XX = GSL_FUNCTION(matrix,alloc)(NCOVAR, NCOVAR);
  BLAS_FUNCTION(gemm)(CblasTrans, CblasNoTrans, 1.0, X, X, 0.0, XX);
  // Compute X'X + shrinkage
  GSL_TYPE(matrix) * K = GSL_FUNCTION(matrix,calloc)(shrinkage_to_use->size, shrinkage_to_use->size);
  for(i = 0; i < shrinkage_to_use->size; ++i)
    {
      GSL_FUNCTION(matrix,set)(K, i, i, GSL_FUNCTION(vector,get)(shrinkage_to_use, i));
    }
  GSL_FUNCTION(matrix,add)(XX, K);
  // Compute solve(X'X + shrinkage)
  GSL_TYPE(matrix) * invXXK = GSL_FUNCTION(matrix,alloc)(XX->size1, XX->size2);
  MY_FUNCTION(solve)(XX, invXXK);
  // Compute (X'X + shrinkage)^(-1) X'
  GSL_TYPE(matrix) * invXXKXprime = GSL_FUNCTION(matrix,alloc)(X->size2, X->size1);
  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasTrans, 1.0, invXXK, X, 0.0, invXXKXprime);
  // Compute (X'X + shrinkage)^(-1) X'y
  BLAS_FUNCTION(gemv)(CblasNoTrans, 1.0, invXXKXprime, pheno, 0.0, beta);
  return beta;
}

int computeLogisticRidge(GSL_TYPE(vector) * beta,
				GSL_TYPE(matrix) * pred,
				GSL_TYPE(vector) * pheno,
				GSL_TYPE(vector) * shrinkage,
				int intercept_flag,
				int DofF_flag,
				PREC * DofF)
{
  // pred has dimensions (1) NINDIV (2) intercept_flag + NSNPS + NCOVAR
  int NINDIV = pred->size1;
  int NPRED = pred->size2;
  int i = 0;
  // Not allocated until we need the DofF part
  GSL_TYPE(matrix) * invtXWX = NULL;
  GSL_TYPE(matrix) * W_out = NULL;
  if(DofF_flag)
    {
      invtXWX = GSL_FUNCTION(matrix,alloc)(NPRED, NPRED);
      W_out = GSL_FUNCTION(matrix,alloc)(NINDIV, NINDIV);
    }
  // Start Loop here - we already have initial estimates of beta
  // PREC to store objective function
  PREC objOld = 0.0, objNew = 0.0;
  // First estimates of beta are in the vector in the input function 
  // We need an initial objective function
  // Depending how the inverting function goes we could rewrite this as a vector
  // Such that updateBeta takes a vector kI not a matrix
  	GSL_TYPE(matrix) * kI = GSL_FUNCTION(matrix,calloc)(NPRED, NPRED);
  	for(i = intercept_flag; i < shrinkage->size; ++i)
	  {
	    GSL_FUNCTION(matrix,set)(kI, i, i, 2*GSL_FUNCTION(vector,get)(shrinkage, i));
	  }
   	objOld = objectiveFunction(beta, pred, pheno, shrinkage, intercept_flag);
	updateBeta(beta, pred, pheno, kI, intercept_flag, DofF_flag, invtXWX, W_out);
  	objNew = objectiveFunction(beta, pred, pheno, shrinkage, intercept_flag);
  	while(MATHS_FUNCTION(fabs)(objOld - objNew) > pow(10,PREC_DIFF))
	  {
	    objOld = objNew;
	    updateBeta(beta, pred, pheno, kI, intercept_flag, DofF_flag, invtXWX, W_out);
	    objNew=0.0;
	    objNew = objectiveFunction(beta, pred, pheno, shrinkage, intercept_flag);
	  }
	// compute DofF
	if(DofF_flag)
	  {
      /* // Need to get W, X, invtXWX <- the latter already includes kI */
      /* // X = pred */
      /* // Compute WX */
      GSL_TYPE(matrix) * WX = GSL_FUNCTION(matrix,alloc)(W_out->size1, pred->size2);
      BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, 1.0, W_out, pred, 0.0, WX);
      GSL_FUNCTION(matrix,free)(W_out);
      GSL_TYPE(matrix) * WXinvtXWX = GSL_FUNCTION(matrix,alloc)(WX->size1, invtXWX->size2);
      BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, 1.0, WX, invtXWX, 0.0, WXinvtXWX);
      GSL_FUNCTION(matrix,free)(WX);
      GSL_FUNCTION(matrix,free)(invtXWX);
      GSL_TYPE(matrix) * hat = GSL_FUNCTION(matrix,alloc)(NINDIV, NINDIV);
      GSL_TYPE(matrix) * hat2 = GSL_FUNCTION(matrix,alloc)(NINDIV, NINDIV);
      BLAS_FUNCTION(gemm)(CblasNoTrans, CblasTrans, 1.0, WXinvtXWX, pred, 0.0, hat);
      GSL_FUNCTION(matrix,free)(WXinvtXWX);
      // Compute hat2
      BLAS_FUNCTION(gemm)(CblasNoTrans, CblasTrans, 1.0, hat, hat, 0.0, hat2);
      GSL_FUNCTION(matrix,free)(hat);
      GSL_FUNCTION(vector,view) diag = GSL_FUNCTION(matrix,diagonal)(hat2);
      sumVector(&diag.vector, DofF);
      GSL_FUNCTION(matrix,free)(hat2);
    }
  GSL_FUNCTION(matrix,free)(kI);
   return 0;
}

PREC objectiveFunction(GSL_TYPE(vector) * beta,
				GSL_TYPE(matrix) *X,
				GSL_TYPE(vector) * pheno,
				GSL_TYPE(vector) * shrinkage,
				int intercept_flag)
{
	int i;
  	PREC obj = 0.0;
  	const size_t NINDIV = X->size1;
  	GSL_TYPE(vector) * XB = GSL_FUNCTION(vector,calloc)(NINDIV);
  	GSL_TYPE(vector) * p = GSL_FUNCTION(vector,calloc)(NINDIV);
  	compute_XB_and_p(X, beta, XB, p);
  	GSL_FUNCTION(vector,free)(XB);
  	GSL_TYPE(vector) * top = GSL_FUNCTION(vector,calloc)(NINDIV);
  	for(i = 0; i < NINDIV; ++i)
    	{
      		if (GSL_FUNCTION(vector,get)(pheno, i) == 0.0)
  		{
			GSL_FUNCTION(vector,set)(top, i, MATHS_FUNCTION(log)(1 - GSL_FUNCTION(vector,get)(p, i)));
  		} else if (GSL_FUNCTION(vector,get)(pheno, i) == 1.0)
  		{
			GSL_FUNCTION(vector,set)(top, i, MATHS_FUNCTION(log)(GSL_FUNCTION(vector,get)(p,i)));
		}
	}
	GSL_FUNCTION(vector,free)(p);
	sumVector(top, &obj);
	GSL_FUNCTION(vector,free)(top);
	PREC res = 0.0;
	for(i = intercept_flag; i < shrinkage->size; ++i)
	{
		res = res + GSL_FUNCTION(vector,get)(shrinkage, i) * GSL_FUNCTION(vector,get)(beta,
		i) *
		GSL_FUNCTION(vector,get)(beta, i);
	}
	obj = obj - res;
	return obj;
}

int updateBeta(GSL_TYPE(vector) * beta,
			GSL_TYPE(matrix) * X,
			GSL_TYPE(vector) * pheno,
			GSL_TYPE(matrix) * kI,
			int intercept_flag,
			int DofF_flag,
			GSL_TYPE(matrix) * invtXWX_return,
			GSL_TYPE(matrix) * W_return)
{
	int i = 0;
	const size_t NINDIV = X->size1;
	const size_t NPRED = X->size2;
	GSL_TYPE(vector) * XB = GSL_FUNCTION(vector,calloc)(NINDIV);
	GSL_TYPE(vector) * p = GSL_FUNCTION(vector,calloc)(NINDIV);
	compute_XB_and_p(X, beta, XB, p);
  	GSL_TYPE(vector) * negP = GSL_FUNCTION(vector,alloc)(NINDIV);
  	GSL_TYPE(vector) * PnegP = GSL_FUNCTION(vector,alloc)(NINDIV);
  	GSL_FUNCTION(vector,memcpy)(negP, p);
  	GSL_FUNCTION(vector,scale)(negP, -1.0);
  	GSL_FUNCTION(vector,add_constant)(negP, 1.0);
  	GSL_FUNCTION(vector,memcpy)(PnegP, p);
  	GSL_FUNCTION(vector,mul)(PnegP, negP);
  	GSL_FUNCTION(vector,free)(negP);
  	//
  	GSL_TYPE(matrix) * W_in = GSL_FUNCTION(matrix,calloc)(NINDIV, NINDIV);
  	for(i = 0; i < NINDIV; ++i)
    	{
      		GSL_FUNCTION(matrix,set)(W_in, i, i, GSL_FUNCTION(vector,get)(PnegP, i));
    	}
	GSL_TYPE(vector) * WZ = GSL_FUNCTION(vector,alloc)(NINDIV);
	for(i = 0; i < NINDIV; ++i)
	{
		//This line needs tidying up - can use vector operations
		GSL_FUNCTION(vector,set)(WZ, i, GSL_FUNCTION(vector,get)(PnegP, i) * GSL_FUNCTION(vector,get)(XB, i) + 
		GSL_FUNCTION(vector,get)(pheno, i) - GSL_FUNCTION(vector,get)(p, i)) ;
	}
  	GSL_FUNCTION(vector,free)(XB); 
  	GSL_FUNCTION(vector,free)(p);
  	GSL_FUNCTION(vector,free)(PnegP);
  	//Need to do the final matrix operations
  	GSL_TYPE(vector) * tXWZ = GSL_FUNCTION(vector,alloc)(NPRED);
	BLAS_FUNCTION(gemv)(CblasTrans, 1.0, X, WZ, 0.0, tXWZ);
  	GSL_FUNCTION(vector,free)(WZ);
  	GSL_TYPE(matrix) * tXW = GSL_FUNCTION(matrix,alloc)(NPRED, NINDIV);
	// Can use CULA for this line
	BLAS_FUNCTION(gemm)(CblasTrans, CblasNoTrans, 1.0, X, W_in, 0.0, tXW);
  	GSL_TYPE(matrix) * tXWX = GSL_FUNCTION(matrix,alloc)(NPRED, NPRED);
  	BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, 1.0, tXW, X, 0.0, tXWX);
  	GSL_FUNCTION(matrix,free)(tXW);
  	GSL_FUNCTION(matrix,add)(tXWX, kI);
  	//solve tXWX
  	GSL_TYPE(matrix) * invtXWX_in = GSL_FUNCTION(matrix,calloc)(NPRED, NPRED);
  	MY_FUNCTION(solve)(tXWX, invtXWX_in);
  	BLAS_FUNCTION(gemv)(CblasNoTrans, 1.0, invtXWX_in, tXWZ, 0.0, beta);
  	GSL_FUNCTION(vector,free)(tXWZ);
  	if(DofF_flag)
    	{
      		GSL_FUNCTION(matrix,memcpy)(invtXWX_return, invtXWX_in);
      		GSL_FUNCTION(matrix,memcpy)(W_return, W_in);
    	}
	GSL_FUNCTION(matrix,free)(tXWX);
  	GSL_FUNCTION(matrix,free)(invtXWX_in);
  	GSL_FUNCTION(matrix,free)(W_in);
	return 0;
}

int getProb(GSL_TYPE(vector) * p, GSL_TYPE(vector) * XB)
{
	int i = 0;
	PREC tmp = 0.0;
	const size_t NINDIV = XB->size;
	// make expXB
	GSL_TYPE(vector) * expXB = GSL_FUNCTION(vector,alloc)(NINDIV);
	for (i = 0; i < NINDIV; ++i)
	{
		tmp = GSL_FUNCTION(vector,get)(XB, i);
		GSL_FUNCTION(vector,set)(expXB, i, MATHS_FUNCTION(exp)(tmp));
	}
	// make 1 + expXB (held in the vector expXB)
	GSL_FUNCTION(vector,memcpy)(p, expXB);
	GSL_FUNCTION(vector,add_constant)(expXB, 1.0);
	// make expXB/1+expXB
	GSL_FUNCTION(vector,div)(p, expXB);
	GSL_FUNCTION(vector,free)(expXB);
	return 0;
}


int my_gsl_solve(gsl_matrix * X,
	gsl_matrix * solvedX)
{
  /* NB only works for matrices of type double
     because the gsl_linalg_LU functions
     are only written for double precision */
  int nrow_x = X->size1;
  int ncol_x = X->size2;
  int nrow_solved_x = solvedX->size1;
  int ncol_solved_x = solvedX->size2;
  if(nrow_x != ncol_x || nrow_x != nrow_solved_x || ncol_x != ncol_solved_x)
    {
      error("ERROR: dimensions error in my_gsl_solve\n");
    } else {
    //solve tXWX
    int signum;
    gsl_permutation * perm = gsl_permutation_alloc(ncol_x);
    gsl_linalg_LU_decomp(X, perm, &signum);
    gsl_linalg_LU_invert(X, perm, solvedX);
    gsl_permutation_free(perm);
  }
  return 0;
}


int compute_XB_and_p(GSL_TYPE(matrix) * X,
	GSL_TYPE(vector) * B,
	GSL_TYPE(vector) * XB,
	GSL_TYPE(vector) * p)
{
	BLAS_FUNCTION(gemv)(CblasNoTrans, 1.0, X, B, 0.0, XB);
	getProb(p, XB);
	return 0;
}


int chooseHowManyK(GSL_TYPE(vector) * D)
{
  // int to return
  int howManyK = 0;
  // Vector to store D^2
  GSL_TYPE(vector) * D2 = GSL_FUNCTION(vector,alloc)(D->size);
  // Copy D2 into it
  GSL_FUNCTION(vector,memcpy)(D2, D);
  // Multiply by D
  GSL_FUNCTION(vector,mul)(D2, D);
  // Compute the total
  PREC tot = GSL_BLAS_FUNCTION(asum)(D2);
  // Vector view
  GSL_FUNCTION(vector,view) sub;
  int i = 0;
  sub = GSL_FUNCTION(vector,subvector)(D2, 0, i+1);
  // Sum of vector view
  PREC cumsum = GSL_BLAS_FUNCTION(asum)(&sub.vector);
  PREC tmp = cumsum / tot;
  while(tmp < 0.90)
    {
      i++;
      sub = GSL_FUNCTION(vector,subvector)(D2, 0, i+1);
      cumsum = GSL_BLAS_FUNCTION(asum)(&sub.vector);
      tmp = cumsum / tot;
    }
  GSL_FUNCTION(vector,free)(D2);
  howManyK=i+1;
  return howManyK;
}


/* return to original scale linear */
/* also puts variables back in their correct places */
int returnToOriginalScaleLinear(GSL_TYPE(vector) * betaOut, 
				GSL_TYPE(vector) * Bridge, 
				GSL_TYPE(vector) * means,
				GSL_TYPE(vector) * scales,
				PREC y_mean,
				int intercept_flag)
{
  /* Variable for iterating */
  int i = 0;
  /* vector to contain the fitted coefficients on the original scale */
  GSL_TYPE(vector) * scaledCoef = GSL_FUNCTION(vector,alloc)(Bridge->size);
  /* Copy the fitted coefficients to it */
  GSL_FUNCTION(vector,memcpy)(scaledCoef, Bridge);
  /* Divide */
  GSL_FUNCTION(vector,div)(scaledCoef, scales);
  /* Check whether the fitted coefficients contain any NAs */
  int contains_nans = 0;
  for(i = 0; i < scaledCoef->size; i++)
    {
      if(isnan(GSL_FUNCTION(vector, get)(scaledCoef, i)))
	{
	  contains_nans++;
	}
    }
  /* If there is an intercept flag */
  if(intercept_flag)
    {
      /* Adjust for other predictors */
      /* R code: */
      /* inter <- object$ym - scaledcoef %*% object$xm */
      /* scaledcoef corresponds to Bridge */
      /* object$xm corresponds to means */
      /* Double to hold product of Bridge and means */
      PREC factor_to_adjust_by = 0.0;
      /* compute the dot product */
      for(i = 0; i < scaledCoef->size; i++)
	{
	  factor_to_adjust_by = GSL_FUNCTION(vector, get)(scaledCoef, i) * GSL_FUNCTION(vector, get)(means, i);
	  /* Subtract */
	  y_mean = y_mean - factor_to_adjust_by;
	  factor_to_adjust_by = 0;
	}
      /* Put in the correct place in vector to be returned */
      GSL_FUNCTION(vector,set)(betaOut, 0, y_mean);
    }
  for(i = intercept_flag; i < betaOut->size; ++i)
    {
      GSL_FUNCTION(vector, set)(betaOut, i, GSL_FUNCTION(vector,get)(scaledCoef, i - intercept_flag));
    }
  GSL_FUNCTION(vector,free)(Bridge);
  GSL_FUNCTION(vector,free)(scaledCoef);
  return 0;
}


char setTrans(CBLAS_TRANSPOSE_t Trans)
{
  char trans = 0;
  if(Trans == CblasNoTrans)
    {
      trans = 'N';
    } else if(Trans == CblasTrans)
    {
      trans = 'T';
    }
  return trans;
}

/* Prepare for linear ridge regression */
int prepareForLinearRidge(GSL_TYPE(matrix) * X, 
			  GSL_TYPE(vector) * y, 
			  GSL_TYPE(matrix) * U, 
			  GSL_TYPE(matrix) * V, 
			  GSL_TYPE(vector) * D, 
			  GSL_TYPE(vector) * D2, 
			  GSL_TYPE(matrix) * Z,
			  GSL_TYPE(vector) * ahat)
{
  /* This function computes the canonical form non shrunken regression coefficeints ahat and the principal components in Z */
  /* Compute the single value decomposition of X */
  SVD_FUNCTION(X, U, V, D);
  /* Compute the principal components */
  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, 1.0, X, V, 0.0, Z);
  /* Make the vector of eigenvalues D2 */ 
  int i;
  for(i = 0; i < D->size; i++)
    {
      GSL_FUNCTION(vector,set)(D2, i, MATHS_FUNCTION(pow)(GSL_FUNCTION(vector,get)(D, i), 2));
    }
  /* Make a diagonal matrix of D2 */
  GSL_TYPE(matrix) * D2_diag = GSL_FUNCTION(matrix,calloc)(D2->size, D2->size);
  /* Set the diagonal elements of D2_diag */
  for(i = 0; i < D2->size; i++)
    {
      GSL_FUNCTION(matrix,set)(D2_diag, i, i, (1.0 / GSL_FUNCTION(vector,get)(D2, i)));
    } 
  /* Make the vector t(Z) %*% y */
  GSL_TYPE(vector) * tZy = GSL_FUNCTION(vector,calloc)(Z->size2);
  /* Fill it */
  BLAS_FUNCTION(gemv)(CblasTrans, 1.0, Z, y, 0.0, tZy);
  /* Compute ahat */
  BLAS_FUNCTION(gemv)(CblasNoTrans, 1.0, D2_diag, tZy, 0.0, ahat);
  GSL_FUNCTION(vector,free)(tZy);
  GSL_FUNCTION(matrix,free)(D2_diag);
  return 0;
}


int computeDofF(GSL_TYPE(vector) * D2,
		PREC Kr,
		PREC * DofF)
{
  // Computes tr(HH) as sum_{i = 1}^{D2->size} \left
  // PREC to store DofF
  // Numerator vector
  GSL_TYPE(vector) * numerator_vector = GSL_FUNCTION(vector,calloc)(D2->size);
  // Copy D2 to it
  GSL_FUNCTION(vector,memcpy)(numerator_vector, D2);
  // Multiply by D2
  GSL_FUNCTION(vector,mul)(numerator_vector, D2);
  // Denominator vector
  GSL_TYPE(vector) * denominator_vector = GSL_FUNCTION(vector,calloc)(D2->size);
  // Copy D2 to it
  GSL_FUNCTION(vector,memcpy)(denominator_vector, D2);
  // Add Kr to it
  GSL_FUNCTION(vector,add_constant)(denominator_vector, Kr);
  // Multiply it by itself
  GSL_FUNCTION(vector,mul)(denominator_vector, denominator_vector);
  /* Divide numerator vector by denominaotr vector */
  GSL_FUNCTION(vector,div)(numerator_vector, denominator_vector); 
  /* Sum it */
  * DofF = GSL_BLAS_FUNCTION(asum)(numerator_vector);
  GSL_FUNCTION(vector,free)(numerator_vector);
  GSL_FUNCTION(vector,free)(denominator_vector);
  return 0;
}

int invert_sum_of_matrices(PREC Ainv,
			   const GSL_TYPE(matrix) * B,
			   const GSL_TYPE(vector) * Dinv,
			   const GSL_TYPE(matrix) * tC,
			   GSL_TYPE(matrix) * out)
{
	/* Check dims */
	int p = B->size1;
	/* D is the diagonal of the square W and of length q */
	int q = Dinv->size;
	/* B is the matrix t(X) and is of dimensions p x q */
	int ncolB = B->size2;
	/* tC is the matrix X and is of dimensions q x p */
	int nrowtC = tC->size1;
	int ncoltC = tC->size2;
	/* out is the matrix inv(t(X) %*% W %*% X + kI) and is of dims p x p */
	int nrowout = out->size1;
	int ncolout = out->size1;
	if( q != ncolB || p != ncoltC || q != nrowtC || p != nrowout || p != ncolout )
	{
		error("ERROR: Wrong dimensions in invert_sum_of_matrices\n");
	} else {
	  // Counter for D
	  int i = 0;
	  // Make Ainv %*% B
	  GSL_TYPE(matrix) * AinvB = GSL_FUNCTION(matrix,calloc)(p, q);
	  GSL_FUNCTION(matrix,memcpy)(AinvB, B);
	  GSL_FUNCTION(matrix,scale)(AinvB,Ainv);
	  // Make tC %*% Ainv
	  GSL_TYPE(matrix) * tCAinv = GSL_FUNCTION(matrix,calloc)(q, p);
	  GSL_FUNCTION(matrix,memcpy)(tCAinv, tC);
	  GSL_FUNCTION(matrix,scale)(tCAinv, Ainv);
	  // Make tC %*% Ainv %*% B
	  // name is invDinvplustCAinvB
	  // because this is what the matrix stores later in the function
	  GSL_TYPE(matrix) * invDinvplustCAinvB = GSL_FUNCTION(matrix,calloc)(q, q);
	  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, 1.0, tCAinv, B, 0.0, invDinvplustCAinvB);
	  // Make D^(-1)
	  GSL_TYPE(matrix) * DinvplustCAinvB = GSL_FUNCTION(matrix,calloc)(q,q);
	  for(i = 0; i < q; i++)
	    {
	      GSL_FUNCTION(matrix,set)(DinvplustCAinvB, i, i, GSL_FUNCTION(vector,get)(Dinv, i));
	    }
	  // Make (D^(-1) + tC %*% Ainv %*% B)
	  GSL_FUNCTION(matrix,add)(DinvplustCAinvB, invDinvplustCAinvB);
	  // Solve (D^(-1) + tC %*% Ainv %*% B)
	  MY_FUNCTION(solve)(DinvplustCAinvB, invDinvplustCAinvB);
	  // Compute AinvB %*% invDinvplustCAinvB
	  GSL_TYPE(matrix) * first_two = GSL_FUNCTION(matrix,calloc)(p, q); // AinvB %*% invDinvplustCAinvB
	  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, 1.0, AinvB, invDinvplustCAinvB, 0.0, first_two);
	  // Multiply by tCAinv and make negative
	  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, -1.0, first_two, tCAinv, 0.0, out);
	  // Free matrices
	  GSL_FUNCTION(matrix,free)(AinvB);
	  GSL_FUNCTION(matrix,free)(tCAinv);
	  GSL_FUNCTION(matrix,free)(invDinvplustCAinvB);
	  GSL_FUNCTION(matrix,free)(DinvplustCAinvB);
	  GSL_FUNCTION(matrix,free)(first_two);
	  // Add Ainv to diagonal
	  GSL_FUNCTION(vector,view) diagOut = GSL_FUNCTION(matrix,diagonal)(out);
	  GSL_FUNCTION(vector,add_constant)(&diagOut.vector, Ainv);
	}
	return 0;
}

/* Compute unpenalized regression coefficients */
/* NB my_gsl_linear_fit only works on vectors of type double */
gsl_vector * my_gsl_linear_fit(gsl_matrix * X, 
		       gsl_vector * y, 
		       int NROW, 
		       int NCOL)
{
  // Allocate a vector for beta
  gsl_vector * beta = gsl_vector_calloc(NCOL);
  //allocate a workspace
  gsl_multifit_linear_workspace * workspace = gsl_multifit_linear_alloc(NROW, NCOL);
  //allocate some space for the matrix of covariances
  gsl_matrix * covar = gsl_matrix_calloc(NCOL, NCOL);
  //allocate a double for chisq
  double chisq = 0.0;
  gsl_multifit_linear(X, y, beta, covar, &chisq, workspace);
  //free the workspace after use!
  gsl_multifit_linear_free(workspace);
  /* Free the covariance matrix */
  GSL_FUNCTION(matrix,free)(covar);
  return beta;
}
#endif

typedef int make_iso_compilers_happy;

