#include "logistic.h"
#ifdef HAVE_GSL_HEADER

int logisticMain(char * genofilename,
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

  int i = 0; /* Needed for couting loops */

  PREC tau = 0; /* to store the elements for tau_vector */
  
  gsl_vector_int * phenotypes = NULL; /* For storing phenotypes, required unless we are predicting */

  int operationType = 0; /* operationType = 0 -> go straight to fitting using shorts, used for prediction and for coordinate descent
			    operationType = 1 -> need long genotypes, SVD. */
  int automaticK = 0; /* SVD required */

  int singleK = 0; /* SVD required but K is specified as number of components*/
  
  GSL_TYPE(matrix) * predictors; /* allocate matrix for predictors - this is never used if we are only reading in the variables as short */

  /*
    read phenotypes into int vector
    check that they are all int
    convert them to -1,1 (for coordinate descent)
  */
  
  if(!predict_flag)
    {
      if(verbose) 
	Rprintf("Reading phenotypes..."); 
      phenotypes = readLogisticPhenotypes(phenofilename, NINDIV);
      if(verbose) 
	Rprintf("done\n"); 
    }
  

  
  // Check what type of operation we need: Can we fit coordinateDescent using short genotypes
  // or do we need PRECs?
  

  // this function can go in commonFunctions.c / commonFunctions.h
  // this function also checks and sets automaticK and singleK
  

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
	{
	  Rprintf("\tReading thinning file..."); 
	}

      gsl_vector_int * thin = readThinFile(thinfilename,
					   SNPnames,
					   thinning_distance,
					   NINDIV,
					   NSNPS,
					   &nThinnedSnps,
					   verbose);

      if(verbose) 
	{
	  Rprintf("done\n"); 
	}

      GSL_TYPE(matrix) * Z = GSL_FUNCTION(matrix, calloc)(NINDIV, GSL_MIN(nThinnedSnps, NINDIV));
      GSL_TYPE(matrix) * thinnedGenotypes = GSL_FUNCTION(matrix, calloc)(NINDIV, nThinnedSnps);
      GSL_TYPE(vector) * D2 = GSL_FUNCTION(vector, calloc)(Z->size2);
    
      if(verbose)
	{
	  Rprintf("\tReading SNPs, Thinning and computing PCs..."); 
	}

      readSNPsThinAndComputePCs(genofilename,
				thin,
				Z,
				thinnedGenotypes,
				D2,
				&howManyK);
      /* Free the thin vector */
      gsl_vector_int_free(thin);
      /* Free the D2 vector */
      GSL_FUNCTION(vector,free)(D2);
      if(verbose)
	{
	  Rprintf("done\n"); 
	}

     
      if(automaticK)
  	{
	  if(!individualK)
	    {
	      if(verbose)
		{
		  Rprintf("\tAutomatically selecting shrinkage parameter\n"); 
		}
	    } else if (individualK)
	    {
	      if(verbose)
		{
		  Rprintf("\tComputing penalty based on fixed number of PCs: %d PCs\n", individualK); 
		}
	    }
	  
	  PREC doff_var = 0.0; // For storing doff
	  
	  PREC chosenK = 0.0; // For storing chosenK
	  
	  PREC k[howManyK];
	  PREC DofF_vector[howManyK];
	  PREC tmp = 0.0; // For storing sum(betaTmp^2)
	  
	  if(verbose) 
	    {
	      Rprintf("\tComputing penalty for each PC..."); 
	    }

	  int len = 0;
	  
	  for(i = 0; i < howManyK; i++)
	    {
	      if(individualK)
		{
		  len = individualK - 1; // -1 because GSL_TYPE vectors are 0-indexed
		} else {
		len = i;
	      }
	      GSL_TYPE(vector) * betaPC = GSL_FUNCTION(vector,calloc)(len+1);
	      GSL_TYPE(vector) * shrinkagePC = GSL_FUNCTION(vector,calloc)(len+1);
	      GSL_FUNCTION(matrix,view) Zsub = GSL_FUNCTION(matrix,submatrix)(Z, 0, 0, Z->size1, len+1);
	      
	      // Compute Logistic Ridge on first [i+1] PCs, no need to compute DofF
	      coordinateDescentLogistic(betaPC,
					&Zsub.matrix,
					phenotypes,
					shrinkagePC,
					0,
					1,
					PREC_EPS);
	      // Compute k
	      GSL_BLAS_FUNCTION(dot)(betaPC, betaPC, &tmp);
	      k[i] = ((PREC)len + 1.0) / tmp;
	      GSL_FUNCTION(vector,free)(betaPC);
	      GSL_FUNCTION(vector,free)(shrinkagePC);
	    } // Ends for i in 0..howManyK)
	  if(verbose) 
	    {
	      Rprintf("\ndone\n"); 
	    }
	  /* Free the Z matrix */
	  GSL_FUNCTION(matrix,free)(Z);

	  ////////////////////////////////////////////////////
	  // Compute d of f for each PC
	  if(verbose) 
	    {
	      Rprintf("computing DofF for each PC...\n"); 
	    }
	  // Compute Beta, DoF on full data set
	  GSL_TYPE(vector) * betaTmp2 = GSL_FUNCTION(vector,calloc)(thinnedGenotypes->size2);
	  GSL_TYPE(vector) * shrinkageTmp2 = GSL_FUNCTION(vector,calloc)(thinnedGenotypes->size2);
	  for(i = 0; i < howManyK; i++)
	    {
	      // Set shrinkage parameter
	      GSL_FUNCTION(vector,set_all)(shrinkageTmp2, 1 / (2 *k[i]));
	      GSL_FUNCTION(vector,set_zero)(betaTmp2);
	      // Initialise value of DofF to zero
	      doff_var = 0;
	      coordinateDescentLogistic(betaTmp2,
					thinnedGenotypes,
					phenotypes,
					shrinkageTmp2,
					0,
					0,
					PREC_EPS);
	      
	      doff_var = computeDofFLogistic(thinnedGenotypes,
					     betaTmp2,
					     k[i]);
	      
	      DofF_vector[i] = doff_var;
	    } // ends for i in 0..howManyK
	  /* Free thinnedGenotypes */
	  GSL_FUNCTION(matrix,free)(thinnedGenotypes);
	  GSL_FUNCTION(vector,free)(betaTmp2);
	  GSL_FUNCTION(vector,free)(shrinkageTmp2);
	  if(verbose)
	    {
	      Rprintf("done\n"); 
	    }
	  // Choose k
	  if(verbose) 
	    {
	      Rprintf("Choosing shrinkage parameter..."); 
	    }
	  
	  GSL_TYPE(vector) * chooseK = GSL_FUNCTION(vector,alloc)(howManyK);
	  for(i = 0; i < howManyK; i++)
	    {
	      GSL_FUNCTION(vector,set)(chooseK, i, MATHS_FUNCTION(fabs)(DofF_vector[i] - (PREC)i - 1.0));
	    }
	  int whichIndex;
	  whichIndex = GSL_FUNCTION(vector,min_index)(chooseK);
	  GSL_FUNCTION(vector,free)(chooseK);
	  chosenK = k[whichIndex];
	  if(verbose)
	    {
	      Rprintf("done\n"); 
	    }
	  if(verbose) 
	    {
	      Rprintf("Chosen shrinkage parameter is %f\n", chosenK); 
	    }
	  lambda = chosenK;
	  //	    GSL_FUNCTION(vector,free)(y_phen);
	  if(verbose) 
	    {
	      Rprintf("done\n"); 
	    }
	  /* GSL_FUNCTION(matrix, free)(predictors);  */
  	} // Matches if (automaticK)
    } // Matches if (operationType)
  
  // Now fitting: Genotypes as short needed
  if(verbose)
    {
      Rprintf("Reading genotypes as int..."); 
    }
  gsl_matrix_int * genotypesShort = readShortGenotypes(genofilename, NINDIV, NSNPS);
  if(verbose)
    {
      Rprintf("done\n"); 
    }
    /*   Prepare the predictors: */
  /*   If we are not going to standardize, then for each predictor get mean, scale */
  /*   Otherwise we have to scale the genotypes as well */
  /* *\/ */

  if(!predict_flag) {
    if(verbose)
      {
	Rprintf("Preparing genotypes..."); 
      }
    GSL_TYPE(vector) * means = GSL_FUNCTION(vector,calloc)(NSNPS);
    GSL_TYPE(vector) * scales = GSL_FUNCTION(vector,calloc)(NSNPS);
    getGenotypeInfo(genotypesShort,
  		    standardize_flag,
  		    1,
  		    means,
  		    scales,
		    SNPnames);
    if(verbose) 
      {
	Rprintf("done\n"); 
	/* Once we have k, genotypes stored as int, means, scales */
	/* then we can call coordinate descent on the genotype data */
	Rprintf("Calling coordinate descent..."); 
      }
    
    GSL_TYPE(vector) * beta = GSL_FUNCTION(vector,calloc)(intercept_flag + NSNPS);
    /* Prepare Shrinkage logistic */
    GSL_TYPE(vector) * tau_vector = GSL_FUNCTION(vector,calloc)(intercept_flag + NSNPS);
    
    tau = 2 * lambda;
    tau = 1 / tau;
    
    GSL_FUNCTION(vector,set_all)(tau_vector, tau);
    
    /* intercept */
    
    if(intercept_flag)
      {
	GSL_FUNCTION(vector,set)(tau_vector, 0, 0.0);
      }
    coordinateDescentLogisticGenotypes(genotypesShort,
				       phenotypes,
				       intercept_flag,
				       standardize_flag,
				       0,
				       tau_vector,
				       means,
				       scales,
				       beta,
				       convergence_threshold);
    
    if(verbose) 
      {
	Rprintf("done\n"); 
      }
    /* Free the tau_vector */
    GSL_FUNCTION(vector,free)(tau_vector);
      /* Check if p-values are required */
      /* If they are - then fill shrinkage vector */
      if((approxtestfilename != NULL) || (permtestfilename != NULL))
        {
    	/* We have to read in the predictors as PREC */
    	/* Read in all the predictors and scale them to correlation form */
    	if(verbose)
    	  {
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
	
  	/* Allocate shrinkage vector */
  	GSL_TYPE(vector) * shrinkage = GSL_FUNCTION(vector,calloc)(beta->size);
	
  	/* Fill shrinkage vector */
  	for(i = intercept_flag; i < beta->size; i++)
  	  {
  	    GSL_FUNCTION(vector,set)(shrinkage, i, lambda);
  	  }
	
  	/* If approx test - compute pvalues and write to file */
  	if (approxtestfilename != NULL)
  	  {
  	    if(verbose)
  	      {
  		Rprintf("Computing approximate test p-values...");
  	      }
  	    /* Allocate vector for approximate p-values */
  	    GSL_TYPE(vector) * approxPs = GSL_FUNCTION(vector,calloc)(beta->size - intercept_flag);
  	    /* Compute the approximate p-values */
  	    computeApproxPsLogistic(beta,
  				    predictors,
  				    shrinkage,
  				    intercept_flag,
  				    approxPs);
  	    /* Allocate vector for approximate p-values to write out */
  	    GSL_TYPE(vector) * approxPsOut = GSL_FUNCTION(vector, calloc)(NPRED);
  	    /* Temporary vector for means - leave all elements set to zero */
  	    GSL_TYPE(vector) * tmp_for_means = GSL_FUNCTION(vector, calloc)(approxPs->size);
  	    /* Temporary vector for scales */
  	    GSL_TYPE(vector) * tmp_for_scales = GSL_FUNCTION(vector, calloc)(approxPs->size);
  	    /* Set all of the temporary vector for scales to 1 */
  	    GSL_FUNCTION(vector, set_all)(tmp_for_scales, 1);
  	    /* Return to original scale - we use the Linear version because it puts things in the correct places */
  	    returnToOriginalScaleLinear(approxPsOut,
  					approxPs,
  					tmp_for_means,
  					tmp_for_scales,
  					0.0,
  					intercept_flag);
  	    /* Free the temporary vector for means */
  	    GSL_FUNCTION(vector, free)(tmp_for_means);
	    /* The vector for approxPs is already freed in returnToOriginalScaleLinear */
  	    /* Free the temporary vector for scales */
  	    GSL_FUNCTION(vector, free)(tmp_for_scales);
  	    /* Write out the p-values */
  	    writeOut(intercept_flag,
  		     NSNPS,
  		     NCOVAR,
  		     SNPnames,
  		     COVARnames,
  		     approxtestfilename,
  		     approxPsOut);
  	    /* Free the vector for approxPsOut */
  	    GSL_FUNCTION(vector, free)(approxPsOut);
  	    if(verbose)
  	      {
  		Rprintf("done\n");
  	      }
  	  } // Ends if (approxTestFileName != NULL)
  	/* If permutation test - compute pvalues and write to file  */
  	if (permtestfilename != NULL)
  	  {
  	    if(verbose)
  	      {
  		Rprintf("Computing permutation test p-values...\n");
  	      }
  	    /* Allocate vector for permutation p-values */
  	    GSL_TYPE(vector) * permPs = GSL_FUNCTION(vector,calloc)(beta->size - intercept_flag);
  	    /* Compute permutation p-values */
  	    computePermPs(permPs,
  			  predictors,
  			  NULL,
  			  phenotypes,
  			  beta,
  			  lambda,
  			  shrinkage,
  			  1000,
  			  1,
  			  intercept_flag,
  			  "logistic");
  	    /* Allocate vector for permutation p-values to write out */
  	    GSL_TYPE(vector) * permPsOut = GSL_FUNCTION(vector, calloc)(NPRED);
  	    /* Temporary vector for means - leave all elements set to zero */
  	    GSL_TYPE(vector) * tmp_for_means = GSL_FUNCTION(vector, calloc)(permPs->size);
  	    /* Temporary vector for scales */
  	    GSL_TYPE(vector) * tmp_for_scales = GSL_FUNCTION(vector, calloc)(permPs->size);
  	    /* Set all of the temporary vector for scales to 1 */
  	    GSL_FUNCTION(vector, set_all)(tmp_for_scales, 1);
  	    /* Return to original scale - we use the Linear version because it puts things in the correct places */
  	    returnToOriginalScaleLinear(permPsOut,
  					permPs,
  					tmp_for_means,
  					tmp_for_scales,
  					0.0,
  					intercept_flag);
  	    /* Free the temporary vector for means */
  	    GSL_FUNCTION(vector, free)(tmp_for_means);
  	    /* Free the temporary vector for scales */
  	    GSL_FUNCTION(vector, free)(tmp_for_scales);
  	    /* Write out the p-values */
  	    writeOut(intercept_flag,
  		     NSNPS,
  		     NCOVAR,
  		     SNPnames,
  		     COVARnames,
  		     permtestfilename,
  		     permPsOut);
  	    /* Free the vector for approxPsOut */
  	    GSL_FUNCTION(vector, free)(permPsOut);
  	  } // Ends if permtestfilename != null 
	/* Free the shrinkage vector */
  	GSL_FUNCTION(vector,free)(shrinkage);
	/* Free the predictors matrix */
	GSL_FUNCTION(matrix,free)(predictors);
	/* Free the means */
	GSL_FUNCTION(vector,free)(means);
	/* Free the scales */
	GSL_FUNCTION(vector,free)(scales);
	} // if((approxtestfilename != NULL) || (permtestfilename != NULL))

    /* Return the beta to the original scale */
    GSL_TYPE(vector) * betaOut = GSL_FUNCTION(vector,calloc)(intercept_flag + NSNPS);
    returnToOriginalScaleLogistic(betaOut,
  				  beta,
  				  means,
  				  scales,
  				  intercept_flag);
    /* And write them out */
    writeOut(intercept_flag,
    	     NSNPS,
    	     NCOVAR,
    	     SNPnames,
    	     COVARnames,
    	     betafilename,
    	     betaOut);
    /* free the phenotypes vector */
    gsl_vector_int_free(phenotypes);
    /* free the means */
    GSL_FUNCTION(vector,free)(means);
    /* free the scales */
    GSL_FUNCTION(vector,free)(scales);
    /* free beta */
    GSL_FUNCTION(vector,free)(beta);
    /* free betaOut */
    GSL_FUNCTION(vector,free)(betaOut);
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
    // Compute exp(xb):
    // Allocate vector
    GSL_TYPE(vector) * XB = GSL_FUNCTION(vector,alloc)(NINDIV);
    // Set all values to intercept_flag_predict
    GSL_FUNCTION(vector, set_all)(XB, (PREC) intercept_flag_predict);
    // Compute XB (i.e.e intercept + XB)
    BLAS_FUNCTION(gemv)(CblasNoTrans, 1.0, genotypesDouble, betas,intercept_coefficient, XB);
    GSL_TYPE(vector) * p = GSL_FUNCTION(vector, calloc)(NINDIV);
    getProb(p, XB);
    GSL_FUNCTION(matrix,free)(genotypesDouble);
    GSL_FUNCTION(vector, free)(betas);
    // Write to file
    FILE * predfile = fopen(phenofilename, "w");
    GSL_FUNCTION(vector, fprintf)(predfile, p, "%f");
    fclose(predfile);
    /* Free XB vector */
    GSL_FUNCTION(vector,free)(XB);
    /* Free p vector */
    GSL_FUNCTION(vector,free)(p);
    if(verbose)
      {
  	Rprintf("done\n");
      }
  } // Ends if (predict_flag)
  /* Free the short genotypes matrix */
  gsl_matrix_int_free(genotypesShort);
  return 0;
}
#endif

typedef int make_iso_compilers_happy;

