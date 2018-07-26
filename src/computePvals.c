#include "computePvals.h"
#ifdef HAVE_GSL_HEADER
/* 
   computePvals.c
   All the computePvals functions we could need
   i.e. the approx ones and the permutation ones for 
   different sorts of models
*/
/* Compute Approx Ps */
// NB value of shrinkage parameter is already accounted for
// in the vector div (see the function computeLinearRidge)
// but we need it anyway for D2 + 2K

float my_ugaussian_function(float x)
{
	float res = 0.0;
	double x_double = (double)x;
	double res_double = gsl_cdf_ugaussian_P(x_double);
	res = (float)res_double;
	return res;
}

int computeApproxPsLinear(GSL_TYPE(vector) * B, 
			  GSL_TYPE(vector) * y,
			  GSL_TYPE(matrix) * U,
			  GSL_TYPE(vector) * D,
			  GSL_TYPE(vector) * D2,
			  GSL_TYPE(matrix) * V,
			  PREC k,
			  GSL_TYPE(vector) * pvals)
{
  int p = B->size;
  int n = y->size;
  int ddim = D->size;
  GSL_TYPE(vector) * div = GSL_FUNCTION(vector,calloc)(D2->size);
  GSL_FUNCTION(vector,memcpy)(div, D2);
  GSL_FUNCTION(vector,add_constant)(div, k);
  int i = 0;
  GSL_TYPE(vector) * diagvector = GSL_FUNCTION(vector,calloc)(ddim);
  GSL_TYPE(vector) * diagdenom = GSL_FUNCTION(vector,calloc)(ddim);
  PREC diagdenomval;
  /* Make sig2hat as as.numeric(crossprod(y - U %*% diag((D2)/(div)) %*% t(U) %*% y)) / (n - sum( (D2 * (D2 + 2*k)) / ((div)^2))) */
  /* Make the diagonal matrix diag(D2 / div) */
  GSL_FUNCTION(vector,memcpy)(diagvector, D2);
  GSL_FUNCTION(vector,div)(diagvector, div);
  GSL_TYPE(matrix) * diagmatrix = GSL_FUNCTION(matrix,calloc)(ddim, ddim);
  for(i = 0; i < ddim; i++)
    {
      GSL_FUNCTION(matrix,set)(diagmatrix, i, i, GSL_FUNCTION(vector,get)(diagvector, i));
    }
  /* Make U %*% diagmatrix */
  GSL_TYPE(matrix) * UD2div = GSL_FUNCTION(matrix,alloc)(U->size1, diagmatrix->size2);
  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, 1.0, U, diagmatrix, 0.0, UD2div);
  /* Make hat matrix */
  GSL_TYPE(matrix) * H = GSL_FUNCTION(matrix,alloc)(U->size1, U->size1);
  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasTrans, 1.0, UD2div, U, 0.0, H); 
  /* Make ypred */
  GSL_TYPE(vector) * ypred = GSL_FUNCTION(vector,alloc)(n);
  BLAS_FUNCTION(gemv)(CblasNoTrans, 1.0, H, y, 0.0, ypred);
  /* Make residuals */
  GSL_FUNCTION(vector,scale)(ypred, -1);
  GSL_FUNCTION(vector,add)(ypred, y);
  PREC numerator = 0.0;
  GSL_BLAS_FUNCTION(dot)(ypred, ypred, &numerator);
  /* Now make numerator */
  /* n - sum(D * (D2 + 2k) / div^2)*/
  /* Make div^2 */
  GSL_TYPE(vector) * div2 = GSL_FUNCTION(vector,alloc)(div->size);
  GSL_FUNCTION(vector,memcpy)(div2, div);
  GSL_FUNCTION(vector,mul)(div2, div);
  /* Make D2 * (D2 + 2k) */
  GSL_TYPE(vector) * DD22kdiv2 = GSL_FUNCTION(vector,alloc)(D2->size);
  GSL_FUNCTION(vector,memcpy)(DD22kdiv2, D2);
  GSL_FUNCTION(vector,add_constant)(DD22kdiv2, 2*k);
  GSL_FUNCTION(vector,mul)(DD22kdiv2, D2);
  /* Divide by div2 */
  GSL_FUNCTION(vector,div)(DD22kdiv2, div2);
  PREC denom = n - GSL_BLAS_FUNCTION(asum)(DD22kdiv2);
  PREC sig2hat;
  sig2hat = numerator / denom;
  /* Then make varmat */
  for(i = 0; i < ddim; i++)
    {
      diagdenomval = GSL_FUNCTION(vector,get)(div, i);
      diagdenomval = MATHS_FUNCTION(pow)(diagdenomval, 2);
      GSL_FUNCTION(vector,set)(diagdenom, i, diagdenomval);
    }
  GSL_FUNCTION(vector,memcpy)(diagvector, D2);
  GSL_FUNCTION(vector,div)(diagvector, diagdenom);
  for(i = 0; i < ddim; i++)
    {
      GSL_FUNCTION(matrix,set)(diagmatrix, i, i, GSL_FUNCTION(vector,get)(diagvector, i));
    }
  /* V %*% diag(D2 / (div^2)) %*% t(V)*/
  // Step 1 - V %*% diag(D2 / div^2)
  GSL_TYPE(matrix) * varMat = GSL_FUNCTION(matrix,alloc)(p, p);
  GSL_TYPE(matrix) * VD = GSL_FUNCTION(matrix,alloc)(V->size1, V->size2);
  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, 1.0, V, diagmatrix, 0.0, VD);
  /* VD %*% t(V) */
  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasTrans, 1.0, VD, V, 0.0, varMat);
  /* Make sig2hat as as.numeric(crossprod(y - U %*% diag((D2)/(div)) %*% t(U) %*% y)) / (n - sum( (D2 * (D2 + 2*k)) / ((div)^2))) */
  /* i.e. y - predy -> need to compute the hat matrix */
  /* Make the hat matrix */
  // BLAS_FUNCTION(gemv)();
  // PREC resid;
  // GSL_BLAS_FUNCTION(dot)(ypred, ypred, &resid);
  // printf("%f\n", resid);
  /* Make the crossproduct t(resid) %*% resid */
  GSL_FUNCTION(matrix,scale)(varMat, sig2hat);
  GSL_FUNCTION(vector,view) sdview = GSL_FUNCTION(matrix,diagonal)(varMat);
  GSL_TYPE(vector) * Zstat = GSL_FUNCTION(vector,alloc)(varMat->size1);
  GSL_FUNCTION(vector,memcpy)(Zstat,&sdview.vector);
  for(i = 0; i < Zstat->size; i++)
    {
      GSL_FUNCTION(vector,set)(Zstat, i, 1 / MATHS_FUNCTION(sqrt)(GSL_FUNCTION(vector,get)(Zstat, i)));
    }
  GSL_FUNCTION(vector,mul)(Zstat, B);
  for(i = 0; i < Zstat->size; i++)
    {
      GSL_FUNCTION(vector,set)(pvals, i,  2*(1 - UGAUSSIAN_FUNCTION(MATHS_FUNCTION(fabs)(GSL_FUNCTION(vector,get)(Zstat, i)))));
    }
  // Free the variables
  GSL_FUNCTION(matrix,free)(UD2div);
  GSL_FUNCTION(matrix,free)(H);
  GSL_FUNCTION(vector,free)(ypred);
  GSL_FUNCTION(vector,free)(div);
  GSL_FUNCTION(vector,free)(div2);
  GSL_FUNCTION(vector,free)(DD22kdiv2);
  GSL_FUNCTION(matrix,free)(VD);
  GSL_FUNCTION(vector,free)(diagdenom);
  GSL_FUNCTION(vector,free)(diagvector);
  GSL_FUNCTION(matrix,free)(diagmatrix);
  GSL_FUNCTION(matrix,free)(varMat);
  GSL_FUNCTION(vector,free)(Zstat);
  return 0;
}

int computeApproxPsLogistic(GSL_TYPE(vector) * B,
			    GSL_TYPE(matrix) * X,
			    GSL_TYPE(vector) * shrinkage,
			    int intercept_flag,
			    GSL_TYPE(vector) * approxPs)
{
  int NINDIV = X->size1;
  int NCOVAR = X->size2;
  int i;
  // Compute the standard deviations
  GSL_TYPE(vector) * XB = GSL_FUNCTION(vector,alloc)(NINDIV);
  BLAS_FUNCTION(gemv)(CblasNoTrans, 1.0, X, B, 0.0, XB);
  GSL_TYPE(vector) * p = GSL_FUNCTION(vector,alloc)(NINDIV);
  getProb(p, XB);
  GSL_TYPE(matrix) * W = GSL_FUNCTION(matrix,calloc)(NINDIV, NINDIV);
  for(i = 0; i < NINDIV; ++i)
    {
      GSL_FUNCTION(matrix,set)(W, i, i, GSL_FUNCTION(vector,get)(p, i)*(1 - GSL_FUNCTION(vector,get)(p, i)));
    }
  // V <- solve(t(X) %*% W %*% X + kI) %*% (t(X) %*% W %* X) %*% solve(t(X) %*% W %*% X + kI)
  // Varmat <- solve(t(X) %*% W %*% X
  GSL_TYPE(matrix) * Varmat = GSL_FUNCTION(matrix,alloc)(NCOVAR, NCOVAR);
  GSL_TYPE(matrix) * tmp = GSL_FUNCTION(matrix,alloc)(NCOVAR, NINDIV);
  BLAS_FUNCTION(gemm)(CblasTrans, CblasNoTrans, 1.0, X, W, 0.0, tmp);
  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, 1.0, tmp, X, 0.0, Varmat);
  GSL_TYPE(matrix) * VarmatRidge = GSL_FUNCTION(matrix,alloc)(NCOVAR, NCOVAR);
  GSL_FUNCTION(matrix,memcpy)(VarmatRidge, Varmat);
  GSL_TYPE(matrix) * kI = GSL_FUNCTION(matrix,calloc)(NCOVAR, NCOVAR);
  for(i = intercept_flag; i < NCOVAR; ++i)
    {
      GSL_FUNCTION(matrix,set)(kI, i, i, 2 * GSL_FUNCTION(vector,get)(shrinkage, i));
    }
  GSL_FUNCTION(matrix,add)(VarmatRidge, kI);
  GSL_TYPE(matrix) * invVarmatRidge = GSL_FUNCTION(matrix,alloc)(VarmatRidge->size1, VarmatRidge->size2);
	MY_FUNCTION(solve)(VarmatRidge,invVarmatRidge);
  GSL_TYPE(matrix) * V = GSL_FUNCTION(matrix,alloc)(NCOVAR, NCOVAR);
  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, 1.0, invVarmatRidge, Varmat, 0.0, V);
  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, 1.0, V, invVarmatRidge, 0.0, Varmat);
  GSL_TYPE(vector) * se = GSL_FUNCTION(vector,alloc)(NCOVAR);
  for(i = 0; i < NCOVAR; ++i)
    {
      GSL_FUNCTION(vector,set)(se, i, MATHS_FUNCTION(sqrt)(GSL_FUNCTION(matrix,get)(Varmat, i, i)));
    }
  GSL_TYPE(vector) * z = GSL_FUNCTION(vector,alloc)(NCOVAR);
  GSL_FUNCTION(vector,memcpy)(z, B);
  GSL_FUNCTION(vector,div)(z, se);
  PREC prob = 0.0;
  for(i = intercept_flag; i < NCOVAR; ++i)
    {
      prob = 2*(1 - UGAUSSIAN_FUNCTION(MATHS_FUNCTION(fabs)(GSL_FUNCTION(vector,get)(z, i))));
      GSL_FUNCTION(vector,set)(approxPs, i - intercept_flag, prob);
    }
  // Free memory
  GSL_FUNCTION(vector,free)(XB);
  GSL_FUNCTION(vector,free)(p);
  GSL_FUNCTION(matrix,free)(W);
  GSL_FUNCTION(matrix,free)(Varmat);
  GSL_FUNCTION(matrix,free)(tmp);
  GSL_FUNCTION(matrix,free)(VarmatRidge);
  GSL_FUNCTION(matrix,free)(kI);
  GSL_FUNCTION(matrix,free)(invVarmatRidge);
  GSL_FUNCTION(matrix,free)(V);
  GSL_FUNCTION(vector,free)(se);
  GSL_FUNCTION(vector,free)(z);
  return 0;
}

int computeApproxPsGeneralizedLinear(GSL_TYPE(vector) * B,
				     GSL_TYPE(matrix) * pred,
				     GSL_TYPE(vector) * y,
				     GSL_TYPE(vector) * shrinkage,
				     int intercept_flag,
				     GSL_TYPE(vector) * approxPs)
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
    // Make sig2hat:
       	// Make Yhat = X %*% B
    GSL_TYPE(vector) * XB = GSL_FUNCTION(vector,alloc)(NINDIV);
    BLAS_FUNCTION(gemv)(CblasNoTrans, 1.0, X, B, 0.0, XB);
    	// Make resid = y - Yhat
    GSL_TYPE(vector) * Yhat = GSL_FUNCTION(vector,alloc)(NINDIV);
    GSL_FUNCTION(vector,memcpy)(Yhat,y);
    GSL_FUNCTION(vector,sub)(Yhat,XB);
       	// Make sig2hat_numerator = crossprod(resid)
    PREC sig2hat_numerator = 0.0;
    GSL_BLAS_FUNCTION(dot)(Yhat,Yhat,&sig2hat_numerator);	
    	// Make hat matrix
    	// Make XX
    GSL_TYPE(matrix) * XX = GSL_FUNCTION(matrix,alloc)(NCOVAR,NCOVAR);
    BLAS_FUNCTION(gemm)(CblasTrans,CblasNoTrans,1.0,X,X,0.0,XX);
    	// Make XXk
    GSL_TYPE(matrix) * K = GSL_FUNCTION(matrix,calloc)(NCOVAR,NCOVAR);
    for(i = 0; i < NCOVAR; i++)
    {
    	GSL_FUNCTION(matrix,set)(K,i,i,GSL_FUNCTION(vector,get)(shrinkage_to_use,i));
    }
    GSL_TYPE(matrix) * XXk = GSL_FUNCTION(matrix,alloc)(NCOVAR,NCOVAR);
    GSL_FUNCTION(matrix,memcpy)(XXk,XX);
    GSL_FUNCTION(matrix,add)(XXk,K);
    GSL_FUNCTION(matrix,free)(K);
    	// make invXXk
    GSL_TYPE(matrix) * invXXk = GSL_FUNCTION(matrix,alloc)(NCOVAR,NCOVAR);
    MY_FUNCTION(solve)(XXk,invXXk);
    	// make XinvXXk
    GSL_TYPE(matrix) * XinvXXk = GSL_FUNCTION(matrix,alloc)(NINDIV,NCOVAR);
    BLAS_FUNCTION(gemm)(CblasNoTrans,CblasNoTrans,1.0,X,invXXk,0.0,XinvXXk);
    	// Make H
    GSL_TYPE(matrix) * H = GSL_FUNCTION(matrix,alloc)(NINDIV,NINDIV);
    BLAS_FUNCTION(gemm)(CblasNoTrans,CblasTrans,1.0,XinvXXk,X,0.0,H);
    	// Make 2H - HH'
    		// Make 2H (call it H2 because variable names cannot start with a number)
    GSL_TYPE(matrix) * H2 = GSL_FUNCTION(matrix,alloc)(NINDIV,NINDIV);
    GSL_FUNCTION(matrix,memcpy)(H2,H);
    GSL_FUNCTION(matrix,add)(H2,H);
    		// Make HH'
    GSL_TYPE(matrix) * HH = GSL_FUNCTION(matrix,alloc)(NINDIV,NINDIV);
    BLAS_FUNCTION(gemm)(CblasNoTrans,CblasTrans,1.0,H,H,0.0,HH);
    		// Make 2H - HH
    GSL_TYPE(matrix) * HHH = GSL_FUNCTION(matrix,alloc)(NINDIV,NINDIV);
    GSL_FUNCTION(matrix,memcpy)(HHH,H2);
    GSL_FUNCTION(matrix,sub)(HHH,HH);
    	// Make trace(2H - HH')
    GSL_FUNCTION(vector,view) Hdiag = GSL_FUNCTION(matrix,diagonal)(HHH);
	PREC trace2HHH = 0.0;
	sumVector(&Hdiag.vector,&trace2HHH);
    	// Make sig2hat_denominator = n - trace(2H - HH')
	PREC sig2hat_denominator = (PREC) NINDIV;
	sig2hat_denominator = sig2hat_denominator - trace2HHH;
    	// Make sig2hat = sig2hat_numerator / sig2hat_denominator
	PREC sig2hat = sig2hat_numerator / sig2hat_denominator;
    // Make varMat = sig2hat * invXXk %*% XX %*% invXXk:
    	// make invXXk %*% XX
    GSL_TYPE(matrix) * invXXkXX = GSL_FUNCTION(matrix,alloc)(NCOVAR,NCOVAR);
    BLAS_FUNCTION(gemm)(CblasNoTrans,CblasNoTrans,1.0,invXXk,XX,0.0,invXXkXX);
       	// make invXXk %*% XX %*% invXXk
    GSL_TYPE(matrix) * VarMat = GSL_FUNCTION(matrix,alloc)(NCOVAR,NCOVAR);
    BLAS_FUNCTION(gemm)(CblasNoTrans,CblasNoTrans,1.0,invXXkXX,invXXk,0.0,VarMat);
    	// make sig2hat * invXXk %*% XX %*% invXXk;
    GSL_FUNCTION(matrix,scale)(VarMat,sig2hat);
    // Extract square root of diagonal
    GSL_FUNCTION(vector,view) diag_varmat = GSL_FUNCTION(matrix,diagonal)(VarMat);
    // Make Zstat
    PREC Zstat = 0.0;
    for(i = 0; i < NCOVAR; i++)
    {
		Zstat = GSL_FUNCTION(vector,get)(B, i) / MATHS_FUNCTION(sqrt)(GSL_FUNCTION(vector,get)(&diag_varmat.vector,i));
		GSL_FUNCTION(vector,set)(approxPs,i,UGAUSSIAN_FUNCTION(Zstat));    	
    }
    // Compute pvalue
  	return 0;
}

/* Perm Ps - linear - all models */
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
		  char *model)
{
  /* Some counters */
  int i = 0;
  int j = 0;
  /* Number of individuals */
  int NINDIV = 0;
  /* Vector to hold tau */
  GSL_TYPE(vector) * tau_vector = GSL_FUNCTION(vector,calloc)(shrinkage->size);
  /* Get the number of individuals from the phenotype vector */
  if(strcmp(model, "linear") == 0) {
    NINDIV = pheno_linear->size;
  } else if(strcmp(model, "logistic") == 0) {
    NINDIV = pheno_logistic->size;
  }
  /* Get number of columns of pred */
  int NPRED = pred->size2;
  /* matrix to contain the permuted  betas */
  GSL_TYPE(matrix) * permbetamat = GSL_FUNCTION(matrix,calloc)(NPRED, NPERM);
  /* vector to contain that permuoted beta */
  GSL_TYPE(vector) *  permbeta = GSL_FUNCTION(vector,calloc)(NPRED);
  /* Create an instance of a random number generator */
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  /* Set the random seed for the random number generator */
  gsl_rng_set(r, SEED);
  /* Array of numbers to contain permuted Ys */
  /* Have to use an array for the permuting part with gsl_ran_shuffle */
  /* PREC (for linear) */
  PREC permYarray_linear[NINDIV];
  /* int (for logistic) */
  int permYarray_logistic[NINDIV];
  /* Create a vector view - PREC - for linear */
  GSL_FUNCTION(vector,view) permYview_linear;
  /* Create a vector view - int - for logistic */
  gsl_vector_int_view permYview_logistic;
  /* Fill the array */
  /* Linear case */
  if(strcmp(model, "linear") == 0)
    {
      /* For every individual */ 
      for(i = 0; i < NINDIV; i++)
  	{
	  /* Allocate an element of the phenotype vector */
  	  permYarray_linear[i] = GSL_FUNCTION(vector,get)(pheno_linear,i);
  	}
      /* logistic case */
    } else if (strcmp(model, "logistic")==0) {
    /* For every individual */
    for(i = 0; i < NINDIV; i++)
      {
	/* Allocate an element of the phenotype vector */
  	permYarray_logistic[i] = gsl_vector_int_get(pheno_logistic,i);
      }
    /* Prepare the penalty vector - shrinkage - logistic */
    PREC tau;
    tau = 2 * lambda;
    tau = 1 / tau;
    GSL_FUNCTION(vector,set_all)(tau_vector, tau);
    /* intercept */
    if(intercept_flag)
      {
  	GSL_FUNCTION(vector,set)(tau_vector, 0, 0.0);
      }
  }
  /* for every permutation iteration */
  for(i = 0; i < NPERM; ++i) 
    {
      /* linear case */ 
      if(strcmp(model, "linear") == 0) {
	/* permute the array permyarray_linear */
  	gsl_ran_shuffle(r, permYarray_linear, NINDIV, sizeof(PREC));
	/* Create a vector view of the permuted array */
  	permYview_linear = GSL_FUNCTION(vector,view_array)(permYarray_linear, NINDIV);
	/* logistic case */
      } else if (strcmp(model, "logistic") == 0) {
	/* permute the array permYarray_logistic */
  	gsl_ran_shuffle(r, permYarray_logistic, NINDIV, sizeof(int));
	/* Create a vector view of permYarray_logistic */
  	permYview_logistic = gsl_vector_int_view_array(permYarray_logistic, NINDIV);
      }
      if(strcmp(model, "linear") == 0)
  	{
  	  if(lambda != -1)
  	    {
  	      /* Allocate matrices for SVD of pred */
  	      GSL_TYPE(vector) * D = NULL;
  	      GSL_TYPE(vector) * D2 = NULL;
  	      GSL_TYPE(matrix) * U = NULL;
  	      GSL_TYPE(matrix) * V = NULL;
  	      GSL_TYPE(vector) * a = NULL;
  	      GSL_TYPE(matrix) * Z = NULL;
	      /* Get dims of pred */
	      int NROW_PRED = pred->size1;
	      int NCOL_PRED = pred->size2;
	      
	      V = GSL_FUNCTION(matrix,alloc)(NCOL_PRED, GSL_MIN(NROW_PRED, NCOL_PRED));
	      D = GSL_FUNCTION(vector,alloc)(GSL_MIN(NROW_PRED, NCOL_PRED));
	      D2 = GSL_FUNCTION(vector,alloc)(GSL_MIN(NROW_PRED, NCOL_PRED));
	      U = GSL_FUNCTION(matrix,calloc)(NINDIV, GSL_MIN(NROW_PRED, NCOL_PRED));
	      a = GSL_FUNCTION(vector,calloc)(GSL_MIN(NROW_PRED, NCOL_PRED));
	      Z = GSL_FUNCTION(matrix,calloc)(NINDIV, GSL_MIN(NROW_PRED, NCOL_PRED));
  
  	      /*
  		 Prepare for linear ridge regression
  		 prepares the matrix X and the vector y
  		 to be shrunk by amount k
  	      */

	      prepareForLinearRidge(pred, &permYview_linear.vector, U, V, D, D2, Z, a);

	      /* Compute linearRidge */

  	      computeLinearRidge(a, permbeta, D2, V, lambda);

	      /* Put the results in the appropriate place in permbeta */

	      /* Move results from perm_results to permbeta */

	      /* Free Z */
	      GSL_FUNCTION(matrix,free)(V);
	      GSL_FUNCTION(vector,free)(D);
	      GSL_FUNCTION(vector,free)(D2);
	      GSL_FUNCTION(matrix,free)(U);
	      GSL_FUNCTION(vector,free)(a);
	      GSL_FUNCTION(matrix,free)(Z);

  	    } else {
  	    computeLinearGeneralizedRidge(permbeta, pred, &permYview_linear.vector, shrinkage, intercept_flag);
  	  }
  	} else if (strcmp(model, "logistic") == 0)
  	{
  	  /* Make this part refer to coordinate descent logistic */
  	  coordinateDescentLogistic(permbeta,
  			      pred,
  			      &permYview_logistic.vector,
  			      tau_vector,
  			      intercept_flag,
  			      0,
  			      PREC_EPS);
  	}
      GSL_FUNCTION(matrix,set_col)(permbetamat, i, permbeta);
    }
  // compute the pval itself
  int tmp = 0;
  int start = 0; // start is the index of the predictors at which to start computing the permPvals.
  if(strcmp(model,"logistic") == 0)
    {
      start = intercept_flag;
    } else if (strcmp(model, "linear") == 0)
    {
      start = 0;
    }

  for(i = start; i < NPRED; ++i)
    {
      tmp = 0;
      for(j = 0; j < NPERM; ++j)
	{
	  if(MATHS_FUNCTION(fabs)(GSL_FUNCTION(matrix,get)(permbetamat, i, j)) > MATHS_FUNCTION(fabs)(GSL_FUNCTION(vector,get)(Bridge, i)))
	    {
	      tmp = tmp + 1;
	    } 
   	} 
      GSL_FUNCTION(vector,set)(permPs, i-start, (PREC) tmp/(PREC) NPERM);
    }
  /* Free the vector tau_vector */
  GSL_FUNCTION(vector,free)(tau_vector);
  /* Free the matrix permbetamat */
  GSL_FUNCTION(matrix,free)(permbetamat);
  /* Free the vector permbeta */
  GSL_FUNCTION(vector,free)(permbeta);
  /* Free the random number generator */
  gsl_rng_free(r);
  return 0;
}
#endif

typedef int make_iso_compilers_happy;

