#include "logisticFunctions.h"

#ifdef HAVE_GSL_HEADER

/* Logistic Regression Functions */

gsl_vector_int * readLogisticPhenotypes(char * phenotypefilename, int NINDIV)
{
  gsl_vector_int * phen = gsl_vector_int_alloc(NINDIV);
  /* Open a file for reading */
  FILE * phenofile = fopen(phenotypefilename,"r");
  /* Turn off the error handler */
  gsl_set_error_handler_off();
  /* Scan the phenotypes into the int vector */
  int scancheck = 0;
  scancheck = gsl_vector_int_fscanf(phenofile, phen);
  if(scancheck)
    {
      if(scancheck == GSL_EFAILED)
	{
	  error("ERROR: phenotype file %s not formatted correctly\n", phenotypefilename);
	  
	} else  {
	error("failed, gsl_errno=%d\n", scancheck);
      }
    }
  /* Restore the error handler */
  gsl_set_error_handler(NULL);
  /* Close the file for reading */
  fclose(phenofile);
  /* 
     Check that the phenotypes are either zero or one
     If they are then change them to -1, 1 (for coordinate descent)
  */
  int i = 0;
  int tmp = 0;
  for(i = 0; i < NINDIV; i++)
    {
      tmp = gsl_vector_int_get(phen, i);
      if((tmp != 0) && (tmp != 1))
	{
	  error("ERROR: Phenotype value not permitted (must be 0 or 1)\n");
	} else {
	gsl_vector_int_set(phen, i, 2 * tmp - 1);
      }
    }
  return phen;
}

/* return to original scale */
int returnToOriginalScaleLogistic(GSL_TYPE(vector) * betaOut,
				  GSL_TYPE(vector) * Bridge, 
				  GSL_TYPE(vector) * means,
				  GSL_TYPE(vector) * scales, 
				  int intercept_flag)
{
  int i = 0;
  /* return to original scale */
  int NPRED = scales->size + intercept_flag;
  /* get vector view of beta excluding intercept */
  PREC beta = 0.0;
  // Use a vector view of Bridge 
  GSL_FUNCTION(vector,view) Bridge1 = GSL_FUNCTION(vector,subvector)(Bridge, intercept_flag, scales->size);
  // Divide it by scales
  GSL_FUNCTION(vector,div)(&Bridge1.vector, scales);
  // Fill betaOut
  for(i = intercept_flag; i < NPRED; i++)
    {
      beta = GSL_FUNCTION(vector,get)(Bridge, i);
      GSL_FUNCTION(vector,set)(betaOut, i, beta);
    }
  // Intercept
  if(intercept_flag)
    {
      beta = GSL_FUNCTION(vector,get)(Bridge, 0);
      GSL_FUNCTION(vector,mul)(&Bridge1.vector, means);
      PREC tmp = 0.0;
      for(i = 0; i < scales->size; i++)
      {
	tmp = tmp + GSL_FUNCTION(vector, get)(&Bridge1.vector, i);
      }
      beta = beta - tmp;
      GSL_FUNCTION(vector,set)(betaOut, 0, beta);
    }
  return 0;
}

#endif

typedef int make_iso_compilers_happy;

