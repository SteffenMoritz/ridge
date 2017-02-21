#include "linear.h"
#ifdef HAVE_GSL_HEADER

int linearMain(char * genofilename,
	       char * thinfilename,
	       char * phenofilename,
	       char * covarfilename,
	       char * betafilename,
	       char * lambdafilename,
	       char * lambdacovarfilename,
	       char * approxtestfilename,
	       char * permtestfilename,
	       PREC lambda,
	       PREC lambda_c,
	       unsigned long int seed,
	       int howManyK,
	       int individualK,
	       int intercept_flag,
	       int standardize_flag,
	       int standardize_c_flag,
	       int thinning_distance,
	       int NINDIV,
	       int NPRED,
	       int NCOVAR,
	       int NSNPS,
	       char ** SNPnames,
	       char ** COVARnames,
	       int predict_flag,
	       PREC convergence_threshold,
	       int verbose)
{

  /* Declare some variables */
  
  int i = 0; /* Needed for counting loops */

  GSL_TYPE(vector) * phenotypes = NULL; /* For storing genotypes, required unless we are predicting */

  int operationType = 0; /* operation type = 0 -> go straight to fitting using shorts, used for prediction and coordinate descent */  

  int automaticK = 0; /* SVD required */

  int singleK = 0; /* SVD required but K is specified as number of components */

  GSL_TYPE(vector) * betaOut = NULL;

  GSL_TYPE(matrix) * predictors = NULL; /* allocate matrix for predictors - this is never used if we are only reading in the variables as short */

  PREC y_mean = 0;

  if(!predict_flag)
    {
      if(verbose) 
	Rprintf("Reading phenotypes..."); 
      phenotypes = readLinearPhenotypes(phenofilename, NINDIV);
      /* Scale the phenotypes */
      scaley(phenotypes, &y_mean);
      if(verbose) 
	Rprintf("done\n"); 
    }

  GSL_TYPE(vector) * beta = NULL;

  operationType = checkOperationType(lambda,
				     lambda_c,
				     lambdafilename,
				     lambdacovarfilename,
				     approxtestfilename,
				     howManyK,
				     individualK,
				     &automaticK,
				     &singleK,
				     predict_flag);

  /* If computing K based on PCs is needed */ 
  if(operationType)
    {
      if(verbose) 
	Rprintf("Preparing predictors for automatic choice of shrinkage parameter...\n"); 
      
      int nThinnedSnps = 0;
      
      if(verbose) 
	Rprintf("Reading thinning file\n"); 
      
      gsl_vector_int * thin = readThinFile(thinfilename, 
					   SNPnames,
					   thinning_distance,
					   NINDIV,
					   NSNPS,
					   &nThinnedSnps,
					   verbose);
      
      if(verbose)
	{
	  Rprintf("Done\n");
	}
      
      GSL_TYPE(matrix) * Z = GSL_FUNCTION(matrix,calloc)(NINDIV, GSL_MIN(nThinnedSnps, NINDIV));
      GSL_TYPE(matrix) * thinnedGenotypes = GSL_FUNCTION(matrix, calloc)(NINDIV, nThinnedSnps);
      GSL_TYPE(vector) * D2 = GSL_FUNCTION(vector, calloc)(Z->size2);
      
      if(verbose) 
	Rprintf("Reading SNPs, Thinning and computing PCs\n"); 
      
      readSNPsThinAndComputePCs(genofilename,
				thin,
				Z,
				thinnedGenotypes,
				D2,
				&howManyK);
      /* Free the thin vector */
      gsl_vector_int_free(thin);
      GSL_FUNCTION(matrix, free)(thinnedGenotypes);
      
      if(verbose) 
	Rprintf("Done\n"); 
      
      GSL_FUNCTION(matrix,view) Zsub = GSL_FUNCTION(matrix, submatrix)(Z, 0, 0, Z->size1, howManyK);
      
      /* Compute regression coefficients a on all PCs */
      if(verbose) 
	Rprintf("Computing principal components regression coefficients...\n"); 
      GSL_TYPE(vector) * a = MY_FUNCTION(linear_fit)(&Zsub.matrix,
						     phenotypes,
						     Z->size1,
						     howManyK);
      /* coordinateDescentLinearFloat(Z, */
      /* 				   phenotypes, */
      /* 				   a, */
      /* 				   PREC_EPS); */
      if(verbose) 
	Rprintf("Done\n"); 
      
      if(!individualK)
      	{
      	  if(verbose) 
	    Rprintf("Automatically selecting shrinkage parameter\n"); 
      	} else if (individualK)
      	{
      	  if(verbose) 
	    Rprintf("Computing penalty based on fixed number of PCs: %d PCs\n", individualK); 
      	}
      
      
      PREC k[howManyK];
      PREC DofF_vector[howManyK];
      if(individualK)
	{
	  /* Compute linear kr with r = individualK  */
	  computeLinearKr(a,
			  Z,
	  		  phenotypes,
			  D2,
	  		  individualK,
	  		  &k[0],
			  &DofF_vector[0]);
	  lambda = k[0];
	} else {
	for(i = 0; i < howManyK; i++)
	  {
	    computeLinearKr(a,
			    Z,
	    		    phenotypes,
			    D2,
	    		    i+1,
	    		    &k[i],
			    &DofF_vector[i]);
	  }
	/* Free Z matrix */
	GSL_FUNCTION(matrix,free)(Z);
	/* Free D2 vector */
	GSL_FUNCTION(vector,free)(D2);
	if(verbose) 
 Rprintf("Choosing shrinkage parameter..."); 
	GSL_TYPE(vector) * chooseK = GSL_FUNCTION(vector,alloc)(howManyK);
	for(i = 0; i < howManyK; i++)
	  {
	    GSL_FUNCTION(vector,set)(chooseK, i, MATHS_FUNCTION(fabs)(DofF_vector[i] - (PREC)i - 1.0));
	  }
	int whichIndex = 0;
	{
	  whichIndex = GSL_FUNCTION(vector,min_index)(chooseK);
	}
	GSL_FUNCTION(vector,free)(chooseK);
	lambda = k[whichIndex];
	if(verbose) 
 Rprintf("done\n"); 
      } // ends else

      if(verbose) 
	Rprintf("Done\n"); 
      /* Free the vector a */
      GSL_FUNCTION(vector,free)(a);
    } // Ends if operation type

  /* Get short genotypes */
  if(verbose)
    {
      Rprintf("Reading short genotypes..."); 
    }
  gsl_matrix_int * genotypesShort = readShortGenotypes(genofilename,
  						       NINDIV,
  						       NSNPS);
  if(verbose)
    {
      Rprintf("done\n"); 
    }

  if(!predict_flag)
    {
      if(verbose) 
	{
	  Rprintf("Preparing genotypes..."); 
	}
      GSL_TYPE(vector) * means = GSL_FUNCTION(vector,calloc)(NSNPS);
      GSL_TYPE(vector) * scales = GSL_FUNCTION(vector,calloc)(NSNPS);
      beta = GSL_FUNCTION(vector,calloc)(NSNPS);
      getGenotypeInfo(genotypesShort,
		      standardize_flag,
		      1,
		      means,
		      scales,
		      SNPnames);
      if(verbose) 
	{
	  Rprintf("done\n"); 
	  Rprintf("Calling coordinate descent...\n"); 
	}
      coordinateDescentLinearGenotypes(genotypesShort,
				       phenotypes,
				       intercept_flag,
				       standardize_flag,
				       lambda,
				       means,
				       scales,
				       beta,
				       convergence_threshold);
      if(verbose) 
	{
	  Rprintf("Done\n"); 
	}
      
      if(approxtestfilename != NULL || permtestfilename != NULL)
	{
	  if(verbose) 
	    {
	      Rprintf("Computing p-values...\n"); 
	      /* Read in all the predictors and scale them to correlation form */
	      Rprintf("Preparing predictors (as double)...\n"); 
	    }
	  // allocate vector for all_means
	  GSL_TYPE(vector) * all_means = GSL_FUNCTION(vector,alloc)(NSNPS + NCOVAR);
	  // allocate vector for all_scales
	  GSL_TYPE(vector) * all_scales = GSL_FUNCTION(vector,alloc)(NSNPS + NCOVAR);
	  
	  predictors = preparePredictors(NINDIV,
					 NSNPS,
					 SNPnames,
					 NCOVAR,
					 COVARnames,
					 genofilename,
					 covarfilename,
					 intercept_flag,
					 standardize_flag,
					 standardize_c_flag,
					 all_means,
					 all_scales,
					 automaticK);
	  
	  GSL_TYPE(vector) * means = GSL_FUNCTION(vector,alloc)(NSNPS + NCOVAR);
	  GSL_TYPE(vector) * scales = GSL_FUNCTION(vector,alloc)(NSNPS + NCOVAR);
	  
	  PREC mean = 0.0;
	  PREC sd = 0.0;
	  
	  for(i = 0; i < means->size; i++)
	    {
	      mean = GSL_FUNCTION(vector,get)(all_means, i);
	      GSL_FUNCTION(vector,set)(means, i, mean);
	      sd = GSL_FUNCTION(vector,get)(all_scales, i);
	      GSL_FUNCTION(vector,set)(scales, i, sd);
	    }
      
	  GSL_FUNCTION(vector,free)(all_means);
	  GSL_FUNCTION(vector,free)(all_scales);
	  
	  /* Extract the matrix for fitting the genotypes */
	  GSL_FUNCTION(matrix,view) predictorsView = GSL_FUNCTION(matrix,submatrix)(predictors, 0, intercept_flag, NINDIV, predictors->size2 - intercept_flag);
	  
	  int NCOL = predictorsView.matrix.size2;
	  
	  GSL_TYPE(matrix) * V = GSL_FUNCTION(matrix,calloc)(NCOL, GSL_MIN(NCOL,NINDIV));
	  GSL_TYPE(vector) * D = GSL_FUNCTION(vector,calloc)(GSL_MIN(NCOL,NINDIV));
	  GSL_TYPE(vector) * D2 = GSL_FUNCTION(vector,calloc)(GSL_MIN(NCOL,NINDIV));
	  GSL_TYPE(matrix) * U = GSL_FUNCTION(matrix,calloc)(NINDIV, GSL_MIN(NCOL,NINDIV));
	  GSL_TYPE(vector) * a = GSL_FUNCTION(vector,calloc)(GSL_MIN(NCOL,NINDIV));
	  GSL_TYPE(matrix) * Z = GSL_FUNCTION(matrix,calloc)(NINDIV, GSL_MIN(NCOL,NINDIV));
	  
	  
	  prepareForLinearRidge(&predictorsView.matrix,
				phenotypes,
				U,
				V,
				D,
				D2,
				Z,
				a);
	  
	  GSL_TYPE(vector) * betaForTests = GSL_FUNCTION(vector,calloc)(predictorsView.matrix.size2);
	  
	  computeLinearRidge(a,
			     betaForTests,
			     D2,
			     V,
			     lambda);
	  /* Free the Z matrix */
	  GSL_FUNCTION(matrix,free)(Z);
	  
	  if(verbose)
	    {
	      Rprintf("done\n"); 
	    }
	  if(approxtestfilename != NULL)
	    {
	      if(verbose) 
		{
		  Rprintf("Computing approximate test p-values...\n"); 
		}
	      GSL_TYPE(vector) * approxPs = NULL;
	      if (singleK)
		{
		  approxPs = GSL_FUNCTION(vector,calloc)(beta->size);
		  computeApproxPsLinear(betaForTests,
					phenotypes,
					U,
					D,
					D2,
					V,
					lambda,
					approxPs);
		}

	      GSL_TYPE(vector) * approxPsOut = GSL_FUNCTION(vector, calloc)(NPRED);
	      GSL_TYPE(vector) * tmp_for_means = GSL_FUNCTION(vector, calloc)(approxPs->size);
	      GSL_TYPE(vector) * tmp_for_scales = GSL_FUNCTION(vector, calloc)(approxPs->size);
	      GSL_FUNCTION(vector, set_all)(tmp_for_scales, 1);
	      returnToOriginalScaleLinear(approxPsOut,
					  approxPs,
					  tmp_for_means,
					  tmp_for_scales,
					  0.0,
					  intercept_flag);
	      /* Free the temporary vector for means */
	      GSL_FUNCTION(vector, free)(tmp_for_means);
	      /* Free the temporary vector for scales*/
	      GSL_FUNCTION(vector, free)(tmp_for_scales);
	      writeOut(intercept_flag,
		       NSNPS,
		       NCOVAR,
		       SNPnames,
		       COVARnames,
		       approxtestfilename,
		       approxPsOut);
	      /* Free the allocated vector */
	      GSL_FUNCTION(vector,free)(approxPsOut);
	      if(verbose)
		{
		  Rprintf("done\n"); 
		}
	    } // Ends if approxtestfilename
	  
	  //  if permutation test is required
	  if(permtestfilename != NULL)
	    {
	      if(verbose) 
		{
		  Rprintf("Computing permutation test p-values..."); 
		}
	      GSL_TYPE(vector) * permPs = GSL_FUNCTION(vector, calloc)(betaForTests->size);
	      /* Temporary vector shrinkage */
	      GSL_TYPE(vector) * shrinkage = GSL_FUNCTION(vector, calloc)(predictors->size2);
	      int NCOL = predictors->size2 - intercept_flag;
	      GSL_FUNCTION(matrix,view) genotypes = GSL_FUNCTION(matrix,submatrix)(predictors, 0, intercept_flag, NINDIV, NCOL);
	      computePermPs(permPs,
			    &genotypes.matrix,
			    phenotypes,
			    NULL,
			    betaForTests,
			    lambda,
			    shrinkage,
			    1000,
			    1,
			    intercept_flag,
			    "linear");
	      /* free the vector shrinkage */
	      GSL_FUNCTION(vector,free)(shrinkage);
	      GSL_TYPE(vector) * permPsOut = GSL_FUNCTION(vector, calloc)(NPRED);
	      /* Temporary vector for means */
	      /* (Used for returnToOriginalScaleLinear) */
	      GSL_TYPE(vector) * tmp_for_means = GSL_FUNCTION(vector, calloc)(permPs->size);
	      /* Temporary vector for scales */
	      /* (Used for returnToOriginalScaleLinear) */
	      GSL_TYPE(vector) * tmp_for_scales = GSL_FUNCTION(vector, calloc)(permPs->size);
	      /* Set all elements of tmp_for_scales to 1 */
	      /* (So that p-values are not changed by call to returnToOriginalScaleLinear ) */
	      GSL_FUNCTION(vector, set_all)(tmp_for_scales, 1);
	      returnToOriginalScaleLinear(permPsOut,
					  permPs,
					  tmp_for_means,
					  tmp_for_scales,
					  0.0,
					  intercept_flag);
	      /* Free temporary vector for means */
	      GSL_FUNCTION(vector, free)(tmp_for_means);
	      /* Free temporary vector for scales */
	      GSL_FUNCTION(vector, free)(tmp_for_scales);
	      /* write out */
	      writeOut(intercept_flag,
		       NSNPS,
		       NCOVAR,
		       SNPnames,
		       COVARnames,
		       permtestfilename,
		       permPsOut);
	      GSL_FUNCTION(vector,free)(permPsOut);
	      if(verbose)
		{
		  Rprintf("done\n"); 
		}
	    } // Ends if permtestfilename
	  /* Free the variables */
	  GSL_FUNCTION(vector, free)(betaForTests);
	  GSL_FUNCTION(vector, free)(means);
	  GSL_FUNCTION(vector, free)(scales);
	  GSL_FUNCTION(matrix, free)(V);
	  GSL_FUNCTION(vector, free)(D);
	  GSL_FUNCTION(vector,free)(D2);
	  GSL_FUNCTION(matrix, free)(U);
	  GSL_FUNCTION(vector, free)(a);
	} // Ends if (approxtestfilename != NULL || permtestfilename != NULL)

      // However the beta were calculated, return them to their original scale here
      // Allocate vector for betaOut
      if(verbose)
	{ 
	  Rprintf("Returning to scale of original regression coefficients..."); 
	}
      betaOut = GSL_FUNCTION(vector,calloc)(NPRED);
      returnToOriginalScaleLinear(betaOut,
				  beta,
				  means,
				  scales,
				  y_mean,
				  intercept_flag);
      if(verbose) 
	{
	  Rprintf("done\n"); 
	}
      GSL_FUNCTION(vector,free)(means);
      GSL_FUNCTION(vector,free)(scales);
      safelyFreeMatrix(predictors);
      GSL_FUNCTION(vector,free)(phenotypes);
      
      if(verbose) 
	{
	  Rprintf("Writing out..."); 
	}
      writeOut(intercept_flag,
	       NSNPS,
	       NCOVAR,
	       SNPnames,
	       COVARnames,
	       betafilename,
	       betaOut);
      GSL_FUNCTION(vector,free)(betaOut);
      if(verbose) 
	{
	  Rprintf("done\n"); 
	}
      // Ends if !predict_flag 
    } else if (predict_flag) {
    if(verbose)
      {
	Rprintf("Predicting...\n");
      }
    // Read in genotypes as PREC
    GSL_TYPE(matrix) * genotypesDouble = readGenotypes(genofilename, NINDIV, NSNPS);
    // Read in betas as PREC
    // Flag to indicate whether there is an intercept 
    int intercept_flag_predict = 0;
    // Flag to hold the intercept coefficient
    PREC intercept_coefficient = 0;
    // Read the coefficeints
    GSL_TYPE(vector) * betas = readCoefficients(betafilename, &intercept_flag_predict, &intercept_coefficient);
    // Compute XB
    GSL_TYPE(vector) * XB = GSL_FUNCTION(vector, calloc)(NINDIV);
    // Set all the values to intercept_flag
    GSL_FUNCTION(vector, set_all)(XB, (PREC) intercept_flag_predict);
    // Compute XB (i.e. intercept + XB)
    BLAS_FUNCTION(gemv)(CblasNoTrans, 1.0, genotypesDouble, betas, intercept_coefficient, XB);
    // Free genotypesDouble
    GSL_FUNCTION(matrix,free)(genotypesDouble);
    // Free betas
    GSL_FUNCTION(vector, free)(betas);
    // Write to file:
    // Open a file pointer
    FILE * predfile = fopen(phenofilename, "w");
    // Write out the vector
    GSL_FUNCTION(vector, fprintf)(predfile, XB, "%f");
    // Close the file
    fclose(predfile);
    // Free the XB vector
    GSL_FUNCTION(vector,free)(XB);
    if(verbose)
      {
	Rprintf("done\n");
      }
  }
  /* Free the matrix of short genotyeps */
  gsl_matrix_int_free(genotypesShort);
  return 0;
}
#endif

typedef int make_iso_compilers_happy;

