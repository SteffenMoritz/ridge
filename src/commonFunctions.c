#include "commonFunctions.h"
#ifdef HAVE_GSL_HEADER

/* Read short genotypes */
/* Includes error checking to be sure that genotypes are formatted as short */
gsl_matrix_int * readShortGenotypes(char * filename,
	int NROW,
	int NCOL)
{
  /* Allocate a pointer to hold the matrix */
  gsl_matrix_int * mat = NULL;
  /* If the filename is supplied */
  if(filename != NULL)
    {
      /* A file pointer */
      FILE * file;
      /* A character */
      char ch;
      /* Open the file */
      file = fopen(filename, "r");
      /* Read a character */
      ch = fgetc(file);
      /* This part to skip the header row */
      while(ch != (int)'\n')
	{
	  ch = fgetc(file);
	}
      /* Allocate a matrix to contain the contents of the file */
      mat = gsl_matrix_int_alloc(NROW, NCOL);
      /* Turn off the error handler */
      gsl_set_error_handler_off();
      /* integer to store the error handler number */
      int scancheck = 0;
      /* Scan the genotypes into the int matrix */
      scancheck = gsl_matrix_int_fscanf(file, mat);;
      /* Check the returned value */
      if(scancheck)
	{
	  if(scancheck == GSL_EFAILED)
	    {
	      error("ERROR: phenotype file %s not formatted correctly\n", filename);
	      
	    } else  {
	    error("failed, gsl_errno=%d\n", scancheck);
	  }
	}
      /* Restore the error handler */
      gsl_set_error_handler(NULL);
      /* Close the file */
      fclose(file);
    }
  /* Check the genotypes */
  checkGenotypes(mat);
  /* Return the matrix */
  return mat;
}

GSL_TYPE(matrix) * readGenotypes(char * filename,
	int NROW,
	int NCOL)
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
      mat = GSL_FUNCTION(matrix,alloc)(NROW, NCOL);
      GSL_FUNCTION(matrix, fscanf)(file, mat);
      fclose(file);
    }
  return mat;
}

GSL_TYPE(vector) * readCoefficients(char * filename,
				    int * intercept_flag,
				    PREC * intercept_coefficient)
{
  // A file pointer
  FILE * file;
  // Open the file
  file = fopen(filename, "r");
  // Count number of rows
  int NROWS = getNROW(file);
  // Close the file
  fclose(file);
  // Read the first row of filename into a string
  // Max length of line
  int MY_LINE_MAX = 256;
  // Character array to hold the row
  char row_of_file[MY_LINE_MAX];
  // Open the file
  file = fopen(filename, "r");
  // Get the line
  if(fgets(row_of_file, MY_LINE_MAX, file) == NULL)
  	error("error reading from %s\n", filename);
  // Close the file
  fclose(file);
  // Split the line on whitespace 
  // Declare the delimiter
  char * delim = " \t\n";
  // Pointer to string to hold the first word
  char * first_word = NULL;
  // Split the line
  first_word = strtok(row_of_file, delim);
  // If the first word is "Intercept"
  if(strcmp(first_word, "Intercept") == 0)
    {
      // Set intercept_flag to 1
      *intercept_flag = 1;
      // Store the intercept coefficient: 
      // Read it as string; convert it to PREC
      * intercept_coefficient = (PREC) atof(strtok(NULL, delim));
    }
  // Allocate an array 
  PREC coefficients_array[NROWS - *intercept_flag];
  // Open the file
  file = fopen(filename, "r");
  // If intercept_flag
  if(*intercept_flag)
    {
      // Skip the first line
	if(fgets(row_of_file, MY_LINE_MAX, file) == NULL)
  		error("error reading from %s\n", filename);
    }
  // A counter to keep track of lines
  int count = 0;
  // Then for every line
  while(fgets(row_of_file, MY_LINE_MAX, file))
    {
      // Don't need the first word
      first_word = strtok(row_of_file, delim);
      // Get the coefficient
      coefficients_array[count] = (PREC) atof(strtok(NULL, delim));
      // Increase the count
      count = count + 1;
    }
  // Close file
  fclose(file);
  // Create a vector to hold the coefficients
  GSL_TYPE(vector) * coefs = GSL_FUNCTION(vector, alloc)(NROWS - *intercept_flag );
  // Copy the array to the coefs vector:
  // // Create a vector view of the array
  GSL_FUNCTION(vector,view) coefs_view = GSL_FUNCTION(vector, view_array)(coefficients_array, NROWS - *intercept_flag);
  // Copy the vector view of the array to coefs 
  GSL_FUNCTION(vector,memcpy)(coefs, &coefs_view.vector);
  return coefs;
}

int getGenotypeInfo(gsl_matrix_int * genotypes,
	int standardize_flag,
	int corr_form_flag,
	GSL_TYPE(vector) * means,
	GSL_TYPE(vector) * scales,
	char ** names)
{
	int i = 0;
	int NINDIV = genotypes->size1;
	int NSNPS = genotypes->size2;
	PREC mean, sd;
	if(!standardize_flag)
	{
		GSL_FUNCTION(vector,set_all)(means, 0.0);
		GSL_FUNCTION(vector,set_all)(scales, 1.0);
    }
	for(i = 0; i < NSNPS; i++)
    {
      // Get the column
      gsl_vector_int_view column = gsl_matrix_int_column(genotypes, i);
      // Compute the mean
      mean = gsl_stats_int_mean(column.vector.data, column.vector.stride, column.vector.size);
      // Store the mean
      if(standardize_flag)
		GSL_FUNCTION(vector,set)(means, i, mean);
      // Compute the sd
      sd = gsl_stats_int_sd(column.vector.data, column.vector.stride, column.vector.size);
      if(sd == 0)
		{
		  error("%s (SNP %d) is invariant in the sample and will have an NA coefficient\n Please remove it from your input data", names[i], i); 
		}
      if(standardize_flag)
	{
	  GSL_FUNCTION(vector,set)(scales, i, sd);
	} // Ends if(standardize_flag)
    } // Ends for i = 0; i < NCOVAR; i++)
	// Get genotypes in correlation form
	if(corr_form_flag && standardize_flag)
	  {
	    PREC sqrtn1 = MATHS_FUNCTION(sqrt)((PREC)NINDIV - 1.0);
	    GSL_FUNCTION(vector,scale)(scales, sqrtn1);
	  }
	return 0;
}

int convert_int_vector(const gsl_vector_int * src, GSL_TYPE(vector) * dest)
{
	int src_size = src->size;
	int dest_size = dest->size;
	if(src_size != dest_size)
	{
		error("ERROR: Mismatched lengths in convert_int_vector_to_float\n");
	}
	int i = 0;
	for(i = 0; i < src_size; i++)
	{
	  GSL_FUNCTION(vector,set)(dest, i, (PREC) gsl_vector_int_get(src, i));
	}
	return 0;
}

int checkOperationType(PREC lambda,
		       PREC lambda_c,
		       char * lambdafilename,
		       char * lambdacovarfilename,
		       char * approxfilename,
		       int howManyK,
		       int individualK,
		       int * automaticK,
		       int * singleK,
		       int predict_flag)
{
  // Function to decide whether to read coefficients in from file
  // as int or PREC
  // int -> 0
  // PREC -> 1
  // int to store the operationType
  int operationType = 0;
  // This function also sets automaticK and singleK
  // If predict_flag, then operationType is 0
  if(!predict_flag)  
    {
    if(lambda == -1 && lambda_c == -1 && lambdafilename == NULL && lambdacovarfilename == NULL)
      {
	*automaticK = 1;
      }
    if(*automaticK == 1 || individualK || ( lambda != -1 && lambda_c == -1 && lambdafilename == NULL && lambdacovarfilename == NULL))
      {
	*singleK = 1;
      }
    if(*automaticK || individualK)
      {
	operationType = 1;
      }
  }
  return operationType;
}

int checkGenotypes(gsl_matrix_int * mat)
{
  /* Number of rows of mat */
  int NROW = mat->size1;
  /* Number of columns of mat */
  int NCOL = mat->size2;
  /* Counters for iterating */
  int i = 0;
  int j = 0;
  /* int tmp to hold each genotype */
  int tmp = 0;
  for(i = 0; i < NROW; i++)
    {
      for(j = 0; j < NCOL; j++)
	{
	  tmp = gsl_matrix_int_get(mat, i, j);
	  if(tmp != 0 && tmp !=1 && tmp !=2)
	    error("Genotypes must be coded as 0, 1, 2\n");
	}
    }
  return 0;
}

int  checkForInvariantPredictors(char * genofilename, 
				 int NINDIV)
{
  /* Get number of predictors and predictor names */

  int NPRED = 0;

  char ** PREDNAMES = NULL;

  PREDNAMES = getHeaderRow(genofilename, &NPRED);

  /* Read in the entire matrix */

  gsl_matrix_int * genotypes = gsl_matrix_int_calloc(NINDIV, NPRED);

  char ch;

  FILE * genofile = fopen(genofilename, "r");

  /* Skip the header row */
  ch = fgetc(genofile);
  while(ch != (int)'\n')
    {
      ch = fgetc(genofile);
    }
  /* Read in data, checking that the data are formatted as integers */

  /* Turn off the error handler */

  gsl_set_error_handler_off();

  /* integer to store the error handler number */

  int scancheck = 0;

  /* Scan the genotypes into the int matrix */

  scancheck = gsl_matrix_int_fscanf(genofile, genotypes);;

  /* Check the returned value */
  if(scancheck)
    {
      if(scancheck == GSL_EFAILED)
    {
      error("ERROR: phenotype file %s not formatted correctly\n", genofilename);
      
    } else  {
    error("failed, gsl_errno=%d\n", scancheck);
      }
    }
  /* Restore the error handler */
  gsl_set_error_handler(NULL);

  fclose(genofile);

  checkGenotypes(genotypes);

  int i = 0;
  PREC sd = 0;

  /* For each predictor */
  for(i = 0; i < NPRED; i++)
    /* Read the predictor into a vector of PREC */
    {
      /* Temporary int vector for reading column of genotypes as int */
      gsl_vector_int_view predictorsColFloat = gsl_matrix_int_column(genotypes, i);
      /* Temporary vector for holding column of genotypes as PREC */
      GSL_TYPE(vector) * predictorsCol = GSL_FUNCTION(vector, calloc)(NINDIV);
      /* Convert int to PREC */
      convert_int_vector(&predictorsColFloat.vector, predictorsCol);
      /* Compute the standard deviation */
      sd = GSL_STATS_FUNCTION(sd)(predictorsCol->data, predictorsCol->stride, predictorsCol->size);
      /* If the standard deviation is 0, then print an error message */
      if(sd == 0)
  	{
	  error("predictor %d in %s (%s) is invariant in your input data\nPlease remove invariant predictors and re-run\n", 
		i + 1, genofilename, PREDNAMES[i]); 
  	}
      /* Free the temporary vector for the column of predictors*/
      GSL_FUNCTION(vector, free)(predictorsCol);
    }
  /* Free the PREDNAMES */
  /* They are not needed anymore */
  if(NPRED > 0)
    {
      for(i = 0; i < NPRED; i++)
	{
	  free(PREDNAMES[i]);
	}
      free(PREDNAMES);
    }
  /* Free the genotypes matrix */
  gsl_matrix_int_free(genotypes);
  return 0;
}

#endif

typedef int make_iso_compilers_happy;

