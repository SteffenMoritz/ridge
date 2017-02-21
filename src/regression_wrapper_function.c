/* Wrapper function to the regression function 
   to be called from within R 
   in the linearRidgeBig or logisticRidgeBig functions*/
#include "depends.h"

#ifdef HAVE_GSL_HEADER
#include "linear.h"
#include "logistic.h"
#include "ReadInData.h"

void regression_wrapper_function(char **g, 
				 char **p, 
				 char **b,
				 char **a,
				 char **perm,
				 char ** thin,
				 int * intercept,
				 double * l,
				 char **m,
				 int * predict,
				 int * v) {
  char * genofilename = *g;
  char * phenofilename = *p;
  char * betafilename = *b;

  double lambda = *l;
  char * model = *m;
  int predict_flag = *predict;

  int intercept_flag = *intercept;
  //  int standardize_flag = *standardize;
  int standardize_flag = 1;
  int standardize_c_flag = -1;

  int verbose = *v;

  /* Check a, perm and thin files */
  char * approxtestfilename = NULL;
  char * permtestfilename = NULL;
  char * thinfilename = NULL;
  if(strcmp(*a, "NULL") != 0)
    {
      approxtestfilename = *a;
    }
  if(strcmp(*perm, "NULL") != 0)
    {
      permtestfilename = *perm;
    }
  if(strcmp(*thin, "NULL") != 0)
    {
      thinfilename = *thin;
    }

  char * covarfilename = NULL;
  char * lambdafilename = NULL;
  char * lambdacovarfilename = NULL;

  PREC lambda_c = -1;
  PREC convergence_threshold = -1;
  unsigned long int seed = 0;
  int howManyK = 0;
  int individualK = 0;
  int thinning_distance = -1;

  /* Check for phenotypes file
     If it is supplied check it can be opened for reading */
  /* if(!predict_flag) */
  /*   { */
  /*     if(phenofilename == NULL) */
  /* 	{ */
  /* 	  error("You didn't supply a vector of phenotypes (dependent variables), please do so (use -p)\n"); */
  /* 	} else { */
  /* 	checkFileForReading(phenofilename); */
  /*     } */
  /*   } else if (predict_flag) { */
  /*   if(betafilename == NULL) */
  /*     { */
  /* 	error("You requested predict but didn't supply a file containing coefficients, please do so (use -b)\n"); */
  /*     } else { */
  /*     checkFileForReading(betafilename); */
  /*   } */
  /* } */
   
  /* if(howManyK != 0) */
  /*   { */
  /*     printf("howManyK set on command line: %d\n", howManyK); */
  /*   } */
  
  /* /\* Check for genotypes and/or covariates file */
  /*    If one or both are supplied check they can be opened for reading *\/ */
  /* if(genofilename == NULL && covarfilename == NULL) */
  /*   { */
  /*     printf("ERROR: You didn't supply a genotypes file or a covariates file, please do so (use -g and/or -c)\n"); */
  /*     exit(EXIT_FAILURE); */
  /*   } else { */
  /*   checkFileForReading(genofilename); */
  /*   checkFileForReading(covarfilename); */
  /* } */
  
  /* /\* Check beta file *\/ */
  /* if(betafilename == NULL && !predict_flag) */
  /*   { */
  /*     printf("WARNING: You did not supply a name for your beta file. Defaulting to beta.dat\n"); */
  /*     betafilename = strdup("beta.dat"); */
  /*   } */
  
  /* /\* Check prediction file *\/ */
  /* if(predictionfilename == NULL && predict_flag) */
  /*   { */
  /*     printf("WARNING: You did not supply a name for your prediction results file. Defaulting to prediction.dat\n"); */
  /*     predictionfilename = "prediction.dat"; */
  /*   } */
  
  /* // Check for a valid model */
  /* checkModel(model); */
  
  /* // handle the intercept */
  /* if(intercept_flag == -1) */
  /*   { */
  /*     printf("WARNING: No intercept specified. Use --intercept (the default) or --no-intercept\n"); */
  /*     intercept_flag = 1; */
  /*   } */
  
  /* // handle the standardization flag for genotypes */
  /* if(genofilename != NULL && standardize_flag == -1) */
  /*   { */
  /*     printf("WARNING: No standardization flag specified for genotypes. Use --standardize (the default) or --no-standardize\n"); */
  /*     standardize_flag = 1; */
  /*   } */
  
  /* // handle the standardization flag for covariates */
  /* if(covarfilename != NULL && standardize_c_flag == -1 && !predict_flag) */
  /*   { */
  /*     printf("WARNING: No standardization flag specified for covariates. Use --standardize-c (the default) or --no-standardize-c\n"); */
  /*     standardize_c_flag == 1; */
  /*   } */

  /* Get NINDIV: */
  /* Declare an integer variable to hold the number of individuals */
  int NINDIV = 0;
  /* If we are not predicting */
  if(!predict_flag)
    {
      /* File pointer to open the file containing phenotypes */
      FILE * phenofile;
      /* Open the fine */
      phenofile = fopen(phenofilename,"r");
      /* Get the number of individuals */
      NINDIV = getNROW(phenofile);
      /* Close the file */
      fclose(phenofile);
      /* If we are fitting and not predicting, we need to check for invariant predictors */
      checkForInvariantPredictors(genofilename,
				  NINDIV);
      /* Otherwise, if we are predicting */
    } else if (predict_flag) {
    /* File pointer to open the file containing genotypes */
    FILE * genofile;
    /* Open the file */
    genofile = fopen(genofilename, "r");
    /* Get the number of rows in the file */
    NINDIV = getNROW(genofile);
    /* Close th efile */
    fclose(genofile);
    /* Subtract 1 for the header row*/
    NINDIV = NINDIV - 1; // For the header row
  }
  
  /* Get number of SNPs, number of covariates, and names of each */
  
  int NSNPS = 0;
  char ** SNPnames = NULL;
  
  int NCOVAR = 0;
  char ** COVARnames = NULL;
  
  if(genofilename != NULL)
    {
      if(verbose)
  		{
  		 Rprintf("Getting SNP names...");
  		}
      SNPnames = getHeaderRow(genofilename, &NSNPS);
      if(verbose)
  		{
  			Rprintf("done\n");
  		}
    }


  /* if(covarfilename != NULL) */
  /*   { */
  /*     printf("Getting covariate names..."); */
  /*     COVARnames = getHeaderRow(covarfilename, &NCOVAR); */
  /*     printf("done\n"); */
  /*   } */
  
  /* /\* Check convergence threshold is set - if it was not set on the command line, set it to PREC_EPS *\/ */
  /* if(convergence_threshold == -1) */
  /*   { */
  /*     convergence_threshold = PREC_EPS; */
  /*   } else { */
  /*   printf("Using %f as threshold for convergence in coordinate descent\n", convergence_threshold); */
  /* } */

  /* total number of predictors */
  int NPRED = intercept_flag + NCOVAR + NSNPS;

  /* Call the function */
  /* If model is linear: */
  if(strcmp(model, "logistic") == 0)
    {
      /* We can use PREC_EPS as the convergence threshold */
      convergence_threshold = PREC_EPS;
      /* Call logisticMain */
      logisticMain(genofilename,
  		   thinfilename,
  		   phenofilename,
  		   covarfilename,
  		   betafilename,
  		   lambdafilename,
  		   lambdacovarfilename,
  		   approxtestfilename,
  		   permtestfilename,
  		   lambda,
  		   lambda_c,
  		   seed,
  		   howManyK,
  		   individualK,
  		   intercept_flag,
  		   standardize_flag,
  		   standardize_c_flag,
  		   thinning_distance,
  		   NINDIV,
  		   NPRED,
  		   NCOVAR,
  		   NSNPS,
  		   SNPnames,
  		   COVARnames,
  		   predict_flag,
  		   convergence_threshold,
  		   verbose);
  /* Else if model is logistic */
    } else if (strcmp(model, "linear") == 0)
    {
      /* We must set the convergence threshold to 0.000001 */
      convergence_threshold = 0.000001;
      /* Call linearMain */
      linearMain(genofilename,
  		 thinfilename,
  		 phenofilename,
  		 covarfilename,
  		 betafilename,
  		 lambdafilename,
  		 lambdacovarfilename,
  		 approxtestfilename,
  		 permtestfilename,
  		 lambda,
  		 lambda_c,
  		 seed,
  		 howManyK,
  		 individualK,
  		 intercept_flag,
  		 standardize_flag,
  		 standardize_c_flag,
  		 thinning_distance,
  		 NINDIV,
  		 NPRED,
  		 NCOVAR,
  		 NSNPS,
  		 SNPnames,
  		 COVARnames,
  		 predict_flag,
  		 convergence_threshold,
  		 verbose);
    }
  /* iterator for freeing SNP and COVAR names */
  int i = 0;
  /* Free SNP names */
  if(NSNPS > 0)
    {
      for(i = 0; i < NSNPS; i++)
	{
	  free(SNPnames[i]);
	}
      free(SNPnames);
    }
  /* Free covariate names */
  if(NCOVAR > 0)
    {
      for(i = 0; i < NCOVAR; i++)
	{
	  free(COVARnames[i]);
	}
      free(COVARnames);
    }
}

#endif

typedef int make_iso_compilers_happy;

