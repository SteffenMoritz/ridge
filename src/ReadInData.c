#include "ReadInData.h"
#ifdef HAVE_GSL_HEADER

void writePhenotypes(int NINDIV)
{
  gsl_vector * pheno = gsl_vector_alloc(NINDIV);
  int i;
  int half = NINDIV/2;
  for(i = 0; i < half; ++i)
    {
      gsl_vector_set(pheno,i,0);
    }
  for(i = half; i < NINDIV; ++i)
    {
      gsl_vector_set (pheno,  i, 1);
    }
  FILE * phenofile = NULL;
  phenofile = fopen("phenotypes.dat","w");
  gsl_vector_fprintf(phenofile,pheno,"%1.0f");
  fclose(phenofile);
  gsl_vector_free(pheno);
}

int getNROW(FILE *fp)
{
  int NINDIV = 0;
  char chr = getc(fp);
  while (chr != EOF)
    {
      if (chr == '\n')
	{
	  NINDIV = NINDIV + 1;
	}
      chr = getc(fp);
    }
  return NINDIV;
}

char ** getHeaderRow(char * filename, int *N)
{
  int currlistpos = 0;
  FILE * file;
  int ch; 
  // First, get length of header row (N)
  file = fopen(filename, "r");
  ch = fgetc(file);
  int num = 0;
  while(ch != (int)'\n')
    {
      if(ch == (int)' ' || ch == (int)'\t')
	{
	  currlistpos++;
	}
      ch = fgetc(file); 
    }
  fseek(file, -2, SEEK_CUR);
  ch = fgetc(file);
  fclose(file);
  if(ch != (int)' ' && ch != (int)'\t')
    {
      num = 1;
      currlistpos++;
    }
  char ** headerlist = (char **)malloc(currlistpos * sizeof(char *));
  int currheadpos = -1;
  currlistpos = 0;
  file = fopen(filename, "r");
  ch = getc(file);
  while(ch != (int)'\n')
    {
      currheadpos = currheadpos + 1;
      if(ch == (int)' ' || ch == (int)'\t')
	{
	  headerlist[currlistpos] = (char *)realloc((char *)headerlist[currlistpos], (1 + currheadpos) * sizeof(char));
	  headerlist[currlistpos][currheadpos] = '\0';
	  currheadpos = -1;
	  currlistpos++;
	} else {
	if(currheadpos == 0)
	  {
	    headerlist[currlistpos] = (char *)malloc((1 + currheadpos) * sizeof(char));
	    headerlist[currlistpos][currheadpos] = (char)ch;
	  } else {
	  headerlist[currlistpos] = (char *)realloc((char *)headerlist[currlistpos], (1 + currheadpos) * sizeof(char));
	  headerlist[currlistpos][currheadpos] = (char)ch;
	}
      }
      ch = fgetc(file); 
    }
  // Terminate the final string
  currheadpos = currheadpos + 1;
  headerlist[currlistpos] = (char *)realloc((char *)headerlist[currlistpos], (1 + currheadpos) * sizeof(char));
  headerlist[currlistpos][currheadpos] = '\0';
  *N = currlistpos + num;
  // close the file
  fclose(file);
  return headerlist;
}

gsl_vector * concatenateTwoVectors(gsl_vector * result, gsl_vector * vec1, gsl_vector * vec2)
{
  //Check Dims
  int sizeRes = result->size;
  int sizeVec1 = vec1->size;
  int sizeVec2 = 0;
  if(vec2 != NULL)
    {
      sizeVec2 = vec2->size;
    }
  if(sizeRes != sizeVec1 + sizeVec2)
    {
		error("length of result vector must equal sum of length of input vectors"); 
    }
  int i;
  if(sizeVec1 > 0 )
    {
      for(i = 0; i < sizeVec1; ++i)
	{
	  gsl_vector_set(result, i, gsl_vector_get(vec1, i));
	}
    }
  if(sizeVec2 > 0)
    {
      for(i = 0; i < sizeVec2; ++i)
	{
	  gsl_vector_set(result, sizeVec1+i, gsl_vector_get(vec2, i)); 
	}	
    }
  return result;
}

void printMatrix(GSL_TYPE(matrix) * mat)
{
  int i, j;
  int SIZE1 = mat->size1;
  int SIZE2 = mat->size2;
  for(i = 0; i < SIZE1; ++i)
    {
      for(j = 0; j < SIZE2; ++j)
	{
	  Rprintf("%f ", GSL_FUNCTION(matrix,get)(mat, i, j)); 
	}
      Rprintf("\n"); 
    }
  return;
}

void printMatrixTen(GSL_TYPE(matrix) * mat)
{
  int i, j;
  int SIZE1 = 10;
  int SIZE2 = 10;
  int SIZE1_MAT = mat->size1;
  int SIZE2_MAT = mat->size2;
  if(SIZE1_MAT < 10)
    {
      SIZE1 = SIZE1_MAT;
    }
  if(SIZE2_MAT < 10)
    {
      SIZE2 = SIZE2_MAT;
    }
  Rprintf("\n"); 
  for(i = 0; i < SIZE1; ++i)
    {
      for(j = 0; j < SIZE2; ++j)
	{
	  Rprintf("%f ", GSL_FUNCTION(matrix,get)(mat, i, j)); 
	}
      Rprintf("\n"); 
    }
  return;
}

void printIntMatrix(gsl_matrix_int * mat)
{
  int i, j;
  int SIZE1 = mat->size1;
  int SIZE2 = mat->size2;
  for(i = 0; i < SIZE1; ++i)
    {
      for(j = 0; j < SIZE2; ++j)
	{
	  Rprintf("%d ", gsl_matrix_int_get(mat, i, j)); 
	}
      Rprintf("\n"); 
    }
  return;
}

void printIntMatrixTen(gsl_matrix_int * mat)
{
  int i, j;
  int SIZE1 = 10;
  int SIZE2 = 10;
  int SIZE1_MAT = mat->size1;
  int SIZE2_MAT = mat->size2;
  if(SIZE1_MAT < 10)
    {
      SIZE1 = SIZE1_MAT;
    }
  if(SIZE2_MAT < 10)
    {
      SIZE2 = SIZE2_MAT;
    }
  Rprintf("\n"); 
  for(i = 0; i < SIZE1; ++i)
    {
      for(j = 0; j < SIZE2; ++j)
	{
	  Rprintf("%d ", gsl_matrix_int_get(mat, i, j)); 
	}
      Rprintf("\n"); 
    }
  return;
}

void printVector(GSL_TYPE(vector) * Vec)
{
  int i;
  int SIZE = Vec->size;
  Rprintf("\n"); 
  for(i = 0; i < SIZE; ++i)
    {
      Rprintf("%f ", GSL_FUNCTION(vector,get)(Vec, i)); 
      Rprintf("\n"); 
    }
  return;
}

void printIntVector(gsl_vector_int * Vec)
{
  int i;
  int SIZE = Vec->size;
  Rprintf("\n"); 
  for(i = 0; i < SIZE; ++i)
    {
      Rprintf("%d ", gsl_vector_int_get(Vec, i)); 
      Rprintf("\n"); 
    }
  return;
}

void printIntVectorTen(gsl_vector_int * Vec)
{
  int i;
  int SIZE = GSL_MIN(Vec->size,10);
  Rprintf("\n"); 
  for(i = 0; i < SIZE; ++i)
    {
      Rprintf("%d ", gsl_vector_int_get(Vec, i)); 
      Rprintf("\n"); 
    }
  return;
}

void printVectorTen(GSL_TYPE(vector) * Vec)
{
  int i;
  int SIZE = 10;
  int SIZE_VEC = Vec->size;
  if(SIZE_VEC < SIZE)
    {
      SIZE = SIZE_VEC;
    }
  Rprintf("\n"); 
  for(i = 0; i < SIZE; ++i)
    {
      Rprintf("%f ", GSL_FUNCTION(vector,get)(Vec, i)); 
      Rprintf("\n"); 
    }
  return;
}

int sumVector(GSL_TYPE(vector) * vector, PREC * sum)
{
  int length = vector->size;
  int i;
  for(i=0; i < length; ++i)
    {
      * sum = * sum + GSL_FUNCTION(vector,get)(vector, i);
    }
  return 0;
} 

int sumVectorDouble(gsl_vector * vector, double * sum)
{
  int length = vector->size;
  int i;
  for(i=0; i < length; ++i)
    {
      * sum = * sum + gsl_vector_get(vector, i);
    }
  return 0;
} 

/* 
   prepare shrinkage vector of length intercept_flag + NSNPS + NCOVAR 
   if ((lambda != -1 && !NCOVAR) || (lambda_c != -1 && !NSNPS) || (lambda == lambda_c)) then shrinkage is null
   else it is allocated as follows:
   if intercept_flag == 1 then shrinkage[0] = 0
   then
   shrinkage[intercept_flag:NSNPS] are either lambda or contents of lambdafile
   shrinkage[intercept_flag+NSNPS:intercept_flag+NSNPS+NCOVAR] are either lambda_c or contents of lambdacovarfile
*/

/* Not used in R package */

/* int prepareShrinkage(char * model, */
/* 		     PREC * lambda, */
/* 		     PREC * lambda_c,  */
/* 		     char * lambdafilename,  */
/* 		     char * lambdacovarfilename,  */
/* 		     int NSNPS,  */
/* 		     int NCOVAR,  */
/* 		     int intercept_flag,  */
/* 		     GSL_TYPE(vector) * shrinkage,  */
/* 		     int * automaticK, */
/* 		     int * singleK) */
/* { */
/*   if(verbose)  */
/*     Rprintf("Model is %s\n", model);  */
/*   // Check dims of arguemnts */
/*   int shrinkageSize = shrinkage->size; */
/*   int desiredShrinkageSize = intercept_flag + NSNPS + NCOVAR; */
/*   if(shrinkageSize != desiredShrinkageSize) */
/*     { */
/*       error("ERROR: Error in prepareShrinkage - arguments have incorrect dimensinos\n"); */
/*     } */
/*   if(verbose)  */
/*     Rprintf("Checking shrinkage parameters...");  */
/*   if(!NSNPS) */
/*     { */
/*       if(*lambda != -1 && lambdafilename != NULL) */
/* 	{ */
/* 	  error("\nERROR: There are no genotypes, but you specified a value of lambda and a lambdafile name.\n"); */
/* 	} else if (*lambda != -1 && lambdafilename == NULL) */
/* 	{ */
/* 	  error("\nERROR: There are no genotypes, but you specified a value of lambda.\n"); */
/* 	} else if (*lambda == -1 && lambdafilename != NULL) */
/* 	{ */
/* 	  error("\nERROR: There are no genotypes, but you specified a lambdafile name.\n"); */
/* 	} */
/*     } // ends if !NSNPS loop */
/*   // NCOVAR */
/*   if(!NCOVAR) */
/*     { */
/*       if(*lambda_c != -1 && lambdacovarfilename != NULL) */
/* 	{ */
/* 	  error("\nERROR: There are no non-genetic covariates, but you specified a value of lambda_c and a lambdacovarfile name.\n"); */
/* 	} else if (*lambda_c != -1 && lambdacovarfilename == NULL) */
/* 	{ */
/* 	  error("\nERROR: There are no non-genetic covariates, but you specified a value of lambda_c.\n"); */
/* 	} else if (*lambda_c == -1 && lambdacovarfilename != NULL) */
/* 	{ */
/* 	  error("\nERROR: There are no non-genetic covariates, but you specified a lambdacovarfilefile name.\n"); */
/* 	} */
/*     } // ends if !NCOVAR loo */
/*   GSL_TYPE(vector) * lambdavec = NULL; */
/*   GSL_TYPE(vector) * lambdacovarvec = NULL; */
/*   if(NSNPS) */
/*     { */
/*       if(*lambda != -1 && lambdafilename != NULL) */
/* 	{ */
/* 	  error("\nERROR: You supplied both lambda and lambdafile, please supply one only.\n"); */
/* 	} else if(*lambda == -1 && lambdafilename == NULL) */
/* 	{ */
/* 	  warning("There are genetic covariates but you didn't supply a shrinkage paramter (-l) or file of shrinkage paramters (-f).\n"); */
/* 	  warning("The shrinkage parameter will be selected automatically. For unpenalised regression supply -l 0\n");  */
/* 	  * automaticK = 1; */
/* 	  //  */
/* 	} else { */
/* 	lambdavec = GSL_FUNCTION(vector,calloc)(NSNPS); */
/* 	if (*lambda != -1 && lambdafilename == NULL) */
/* 	  { */
/* 	    GSL_FUNCTION(vector,set_all)(lambdavec, *lambda); */
/* 	  } else if (*lambda == -1 && lambdafilename != NULL) */
/* 	  { */
/* 	    FILE * lambdafile; */
/* 	    lambdafile = fopen(lambdafilename, "r"); */
/* 	    GSL_FUNCTION(vector,fscanf)(lambdafile, lambdavec); */
/* 	    fclose(lambdafile); */
/* 	  } */
/*       } */
/*     } // ends if(NSNPS) loop */
/*   if(NCOVAR) */
/*     { */
/*       if(*lambda_c != -1 && lambdacovarfilename != NULL) */
/* 	{ */
/* 	  error("\nERROR: You supplied both lambda_c and lambdacovarfile, please supply one only.\n"); */
/* 	} else if(*lambda_c == -1 && lambdacovarfilename == NULL) */
/* 	{ */
/* 	  if(NSNPS) */
/* 	    { */
/* 	      warning(" There are non-genetic covariates but you did not supply a shrinkage parameter (-d) or file of shrinkage parameters (-v). The non-genetic covariates will not be penalised\n");  */
/* 	      *lambda_c = 0; */
/* 	    } else { */
/* 	    error("\nERROR: There are non-genetic covariates but you didn't supply a shrinkage paramter (-d) or file of shrinkage paramters (-v). For unpenalised regression supply -d 0\n"); */
/* 	  } */
/* 	} else { */
/* 	lambdacovarvec = GSL_FUNCTION(vector,calloc)(NCOVAR); */
/* 	if (*lambda_c != -1 && lambdacovarfilename == NULL) */
/* 	  { */
/* 	    GSL_FUNCTION(vector,set_all)(lambdacovarvec, *lambda_c); */
/* 	  } else if (*lambda_c == -1 && lambdacovarfilename != NULL) */
/* 	  { */
/* 	    FILE * lambdacovarfile; */
/* 	    lambdacovarfile = fopen(lambdacovarfilename, "r"); */
/* 	    GSL_FUNCTION(vector,fscanf)(lambdacovarfile, lambdacovarvec); */
/* 	    fclose(lambdacovarfile); */
/* 	  } */
/*       } */
/*     } // ends if(NCOVAR) loop */
/*   // copy the lambdavector and the lambdacovarvector into the shrinkage vector */
/*   if(intercept_flag) */
/*     { */
/*       GSL_FUNCTION(vector,set)(shrinkage, 0, 0); */
/*     } */
/*   int i = 0; */
/*   if(! * automaticK) */
/*     { */
/*       if(NSNPS) */
/* 	{ */
/* 	  for(i = intercept_flag; i < intercept_flag + NSNPS; i++) */
/* 	    { */
/* 	      GSL_FUNCTION(vector,set)(shrinkage, i, GSL_FUNCTION(vector,get)(lambdavec, i - intercept_flag)); */
/* 	    } */
/* 	} */
/*       if(NCOVAR) */
/* 	{ */
/* 	  for(i = intercept_flag + NSNPS; i < intercept_flag + NSNPS + NCOVAR; i++) */
/* 	    { */
/* 	      GSL_FUNCTION(vector,set)(shrinkage, i, GSL_FUNCTION(vector,get)(lambdacovarvec, i - intercept_flag -  */
/* 	      NSNPS)); */
/* 	    } */
/* 	} */
/*     } */
/*   // set the singleK flag */
/*   if(NSNPS && *lambda != -1 && lambdafilename == NULL && !NCOVAR) */
/*     { */
/*       *singleK = 1; */
/*     } */
/*   if(NCOVAR && *lambda_c != -1 && lambdacovarfilename == NULL && !NSNPS) */
/*     { */
/*       *singleK = 1; */
/*     } */
/*   if(*automaticK && *lambda_c == -1 && lambdacovarfilename == NULL) */
/*     { */
/*       *singleK = 1; */
/*     } */
/*   if(verbose)  */
/*     { */
/*       Rprintf("done\n");  */
/*     } */
/*   if(lambdavec != NULL) */
/*     { */
/*       GSL_FUNCTION(vector,free)(lambdavec); */
/*     } */
/*   if(lambdacovarvec != NULL) */
/*     { */
/*       GSL_FUNCTION(vector,free)(lambdacovarvec); */
/*     } */
/*   return 0; */
/* } // ends function */

/* Check model */
int checkModel(char * model)
{
  if( model == NULL || ( (strcmp(model, "linear") != 0) && (strcmp(model, "logistic") != 0) ))
    {
      error("ERROR: please specify a valid model linear or logistic\n");
    }
  return 0;
}

/* Prepare matrix of predictors */
GSL_TYPE(matrix) *  preparePredictors(int NINDIV, 
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
				int automaticK)
{
  int i = 0;
  int j = 0;
  PREC tmp;
  GSL_TYPE(matrix) * genotypes = NULL;
  GSL_TYPE(matrix) * covariates = NULL;

  genotypes = getDataWithoutHeaderRow(genofilename, NINDIV, NSNPS);

  // Check the genotypes file for valid genotypes
 for(i = 0; i < NINDIV; i++)
   {
     for(j = 0; j < NSNPS; j++)
 	{
 	  tmp = GSL_FUNCTION(matrix,get)(genotypes, i, j);
 	  if(tmp != 0.0 && tmp != 1.0 && tmp !=2.0)
 	    {
 	      error("ERROR: Invalid genotype. Valid genotypes are 0, 1, 2\n");
 	    }
 	}
   }

  covariates = getDataWithoutHeaderRow(covarfilename, NINDIV, NCOVAR);
  // 1 - Fill the predictors matrix

  // 1.1 - include the genotypes
  GSL_TYPE(vector) * means_g = NULL;
  GSL_TYPE(vector) * scales_g = NULL;
  if(NSNPS)
    {
      means_g = GSL_FUNCTION(vector,alloc)(NSNPS);
      scales_g = GSL_FUNCTION(vector,alloc)(NSNPS);
      checkForInvAndStandardize(genotypes,
				intercept_flag,
				intercept_flag + NSNPS,
				standardize_flag,
				1,
				means_g,
				scales_g,
				SNPnames);
    }


  // 1.2 include the covariates 
  GSL_TYPE(vector) * means_c = NULL;
  GSL_TYPE(vector) * scales_c = NULL;
  if(NCOVAR)
    {
      means_c = GSL_FUNCTION(vector,alloc)(NCOVAR);
      scales_c = GSL_FUNCTION(vector,alloc)(NCOVAR);
      checkForInvAndStandardize(covariates,
				intercept_flag + NSNPS,
				intercept_flag + NSNPS + NCOVAR,
				standardize_c_flag,
				0,
				means_c,
				scales_c,
				COVARnames);
    } // Ends if NCOVAR

  GSL_TYPE(matrix) * predictors = GSL_FUNCTION(matrix,alloc)(NINDIV, intercept_flag + NSNPS + NCOVAR);

  // 2.1 - include the intercept
	if(intercept_flag)
	{
		for(i = 0; i < NINDIV; ++i)
		{
	  		GSL_FUNCTION(matrix,set)(predictors, i, 0, 1.0);
		}
   }

  // For genotypes
  GSL_FUNCTION(vector,view) column;
  if(NSNPS)
    {
      for(i=0;i<NSNPS; i++)
  	{
	  column = GSL_FUNCTION(matrix,column)(genotypes, i);
	  GSL_FUNCTION(matrix,set_col)(predictors, intercept_flag + i, &column.vector);
	  GSL_FUNCTION(vector,set)(means, i, GSL_FUNCTION(vector,get)(means_g, i));
	  GSL_FUNCTION(vector,set)(scales, i, GSL_FUNCTION(vector,get)(scales_g, i));
  	} // ends for i=0;i<NSNPS;i++
      GSL_FUNCTION(matrix,free)(genotypes);
      GSL_FUNCTION(vector,free)(means_g);
      GSL_FUNCTION(vector,free)(scales_g);
    } // ends if(NSNPS)

  // 2.3 - include the covariants in predictors, means and scales
  if(NCOVAR)
    {
      for(i=0;i<NCOVAR; i++)
	{
	  // Is this Covariate variant? 
	  column = GSL_FUNCTION(matrix,column)(covariates, i);
	  GSL_FUNCTION(matrix,set_col)(predictors, intercept_flag + i, &column.vector); 
	  GSL_FUNCTION(vector,set)(means, i, GSL_FUNCTION(vector,get)(means_c, i));
	  GSL_FUNCTION(vector,set)(scales, i, GSL_FUNCTION(vector,get)(scales_c, i));
	} // ends for i=0;i<NSNPS;i++
      GSL_FUNCTION(matrix,free)(covariates);
      GSL_FUNCTION(vector,free)(means_c);
      GSL_FUNCTION(vector,free)(scales_c);
    } // ends if(NCOVAR)
  return predictors;
}

/* Function not used in R code */
/* int prepareGenotypes(gsl_matrix * genotypes, int NINDIV, int NSNPS, char * genofilename) */
/* { */
/*   int i, j; */
/*   double tmp; */
/*   if(genofilename != NULL) */
/*     { */
/*       FILE * genofile; */
/*       genofile = fopen(genofilename,"r"); */
/*       gsl_matrix_fscanf(genofile, genotypes); */
/*       fclose(genofile); */
/*       // Check the genotypes file for valid genotypes */
/*       // Commented out for error testing purposes */
/*       //      for(i = 0; i < NINDIV; i++) */
/*       //	{ */
/*       //	  for(j = 0; j < NSNPS; j++) */
/*       //	    { */
/*       //	      tmp = gsl_matrix_get(genotypes, i, j); */
/* 	      //	      if(tmp != 0.0 && tmp != 1.0 && tmp !=2.0) */
/* 	      //		{ */
/* 	      //		  error("ERROR: Invalid genotype. Valid genotypes are 0, 1, 2\n"); */
/* 	      //		} */
/*       //	    } */
/*       //	} */
/*     } */
/*   return 0; */
/* } */

int scaley(GSL_TYPE(vector) * y, PREC * y_mean)
{
  PREC mean = GSL_STATS_FUNCTION(mean)(y->data, y->stride, y->size);
  * y_mean = mean;
  mean = mean * -1;
  GSL_FUNCTION(vector,add_constant)(y, mean);
  return 0;
}

/* scaleX not used in R package */
/* int scaleX(gsl_matrix * X, int n2, int NSNPS) */
/* { */
/*   int j = 0; */
/*   double sqrtn1 = sqrt(n2 - 1); */
/*   double mean, sd; */
/*   gsl_vector * genocol = gsl_vector_alloc(n2); */
/*   for(j = 0; j < NSNPS; j++) */
/*     { */
/*       // get the column */
/*       gsl_matrix_get_col(genocol, X, j); */
/*       // get the mean */
/*       mean = gsl_stats_mean(genocol->data, genocol->stride, genocol->size); */
/*       mean = mean * -1; */
/*       // get the SD */
/*       sd = gsl_stats_sd(genocol->data, genocol->stride, genocol->size); */
/*       sd = sd * sqrtn1; */
/*       sd = 1 / sd; */
/*       // subdract the mean */
/*       gsl_vector_add_constant(genocol, mean); */
/*       // divide by the standard deviation times sqrtn1 */
/*       gsl_vector_scale(genocol, sd); */
/*       // return the column to the matrix X */
/*       gsl_matrix_set_col(X, j, genocol); */
/*     } */
/*   gsl_vector_free(genocol); */
/* return 0;*/
/* } */


/* Append char to string */
int appendToString(char * ptr, int *currentpos, char charToAppend)
{
  (*currentpos)++;
  ptr = (char *)realloc((char *)ptr, (1 + *currentpos) * sizeof(char));
  ptr[*currentpos] = charToAppend;
  return 0;
}

/* Check file exists and can be opened for reading */
void checkFileForReading(char * filename)
{
  if(filename != NULL && access(filename, R_OK) != 0)
    {
      error("Cannot open file %s for reading\n", filename);
    }
}

GSL_TYPE(matrix) * getDataWithoutHeaderRow(char * filename, int NROW, int NCOL)
{
  GSL_TYPE(matrix) * mat = NULL;
  if(filename != NULL)
    {
      FILE * file;
      char ch;
      file = fopen(filename, "r");
      ch = fgetc(file);
      while(ch != (int)'\n')
	{
	  ch = fgetc(file);
	}
      mat = GSL_FUNCTION(matrix,calloc)(NROW, NCOL);
      GSL_FUNCTION(matrix,fscanf)(file, mat);
      fclose(file);
    }
  return mat;
}

/* Check for invariant SNPs and covariates, standardize */
int checkForInvAndStandardize(GSL_TYPE(matrix) * mat, 
			      int START, 
			      int END, 
			      int standardize_flag,
			      int corr_form_flag,
			      GSL_TYPE(vector) * means, 
			      GSL_TYPE(vector) * scales, 
			      char ** names)
{
  int i;
  int LEN = END - START;
  int NINDIV = mat->size1;
  PREC mean, sd;
  if(!standardize_flag)
    {
      GSL_FUNCTION(vector,set_all)(means, 0.0);
      GSL_FUNCTION(vector,set_all)(scales, 1.0);
    }
  for(i = 0; i < LEN; i++)
    {
      // Get the column
      GSL_FUNCTION(vector,view) column = GSL_FUNCTION(matrix,column)(mat, i);
      // Compute the mean
      mean = GSL_STATS_FUNCTION(mean)(column.vector.data, column.vector.stride, column.vector.size);
      // Store the mean
      if(standardize_flag)
	GSL_FUNCTION(vector,set)(means, i, mean);
      // Compute the sd
      sd = GSL_STATS_FUNCTION(sd)(column.vector.data, column.vector.stride, column.vector.size);
      if(sd == 0)
	{
	  error(" %s is invariant in the sample, please remove it from the data\n", names[i]); 
	}
      if(standardize_flag)
	{
	  GSL_FUNCTION(vector,set)(scales, i, sd);
	  // Subtract the mean
	  GSL_FUNCTION(vector,add_constant)(&column.vector, -1 * mean);
	  // Divide by the standard deviation
	  GSL_FUNCTION(vector,scale)(&column.vector, 1/sd);
	  // Return to the matrix genotypes
	} // Ends if(standardize_flag)
    } // Ends for i = 0; i < NCOVAR; i++)
  // Get genotypes in correlation form
  if(corr_form_flag)
    {
    PREC sqrtn1 = MATHS_FUNCTION(sqrt)((PREC)NINDIV - 1.0);
    GSL_FUNCTION(matrix,scale)(mat, 1 / sqrtn1);
    GSL_FUNCTION(vector,scale)(scales, sqrtn1);
    }
  return 0;
}

int sumIntVec(gsl_vector_int * vec)
{
  int i;
  int count = 0;
  int length = vec->size;
  for(i = 0; i < length; i++)
    {
      count = count + gsl_vector_int_get(vec, i);
    }
  return count;
}

int writeOut(int intercept_flag,
	     int NSNPS,
	     int NCOVAR,
	     char ** SNPnames,
	     char ** COVARnames,
	     char * betafilename,
	     GSL_TYPE(vector) * betaOut)
{
  int i = 0;
  char * name = NULL;
  PREC beta = 0;
  FILE * betafile;
  betafile = fopen(betafilename,"w");
  if(intercept_flag)
    {
      name = "Intercept";
      beta = GSL_FUNCTION(vector,get)(betaOut, 0);
      printBeta(name, beta, betafile);
    }
  if(NSNPS)
    {
      for(i = intercept_flag; i < intercept_flag+NSNPS; i++)
      	{
	  
	  beta = GSL_FUNCTION(vector,get)(betaOut, i);
	  printBeta(SNPnames[i - intercept_flag], beta, betafile);
      	}
    }
  if(NCOVAR)
    {
      for(i = intercept_flag + NSNPS; i < intercept_flag + NSNPS + NCOVAR; i++)
      	{
	  name = COVARnames[i - NSNPS - intercept_flag];
	  beta = GSL_FUNCTION(vector,get)(betaOut, i);
	  printBeta(name, beta, betafile);
      	}
    }
  fclose(betafile);
  return 0;
}

int printBeta(char * name, 
	      PREC beta, 
	      FILE * file)
{
  fprintf(file, "%s\t", name);
  if(isnan(beta))
    {
       fprintf(file, "NA\n");
    } else {
    fprintf(file, "%f\n", beta);
  }
  return 0;
}

int safelyFreeVector(GSL_TYPE(vector) * vec)
{
  if (vec != NULL)
    {
      GSL_FUNCTION(vector,free)(vec);
    }
  return 0;
}

int safelyFreeMatrix(GSL_TYPE(matrix) * mat)
{
  if(mat != NULL)
    {
      GSL_FUNCTION(matrix,free)(mat);
    }
  return 0;
}

/* Print Opening Blurb */
/* Function not uesd in R package*/
/* void printOpening(void) */
/* { */
/*  Rprintf("\n\n");  */
/*  Rprintf("\t##############\n");  */
/*  Rprintf("\t# Regression #\n");  */
/*  Rprintf("\t##############\n");  */
/*  Rprintf("\n");  */
/*  Rprintf("\tWritten by Erika Cule 2012\n\n");  */
/* } */
#endif

typedef int make_iso_compilers_happy;

