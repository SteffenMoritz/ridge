#include "thin.h"
#ifdef HAVE_GSL_HEADER

gsl_vector_int * readThinFile(char * thinfilename,
			      char ** SNPNAMES, 
			      int thinning_distance,
			      int NINDIV,
			      int NSNPS,
			      int * nThinnedSnps,
			      int verbose)
{

  /* Allocate vector of int to indicate which SNPs to get for the thinned SNPs */
  gsl_vector_int * thin = gsl_vector_int_calloc(NSNPS);

  int count = 0;

  /* if thinfilename has not been provided */
  if(thinfilename == NULL)
    {
      if(thinning_distance == -1)
	{
	  thinning_distance = GSL_MAX(1, NSNPS / NINDIV);
	}
      /* temporary int to store next pos */
      int next_pos = 0;
      for(count = 0; count < NSNPS; count++)
	{
	  if(count == next_pos)
	    {
	      gsl_vector_int_set(thin, count, 1);
	      next_pos = count + thinning_distance;
	    } else {
	    // Do nothing
	  }
	}
      /* else if thinfilename has been provided */
    } else {
    /* Check thinning distance */
    if(thinning_distance == -1)
      {
	if(verbose) 
	  {
	    Rprintf("Thinning SNPs using default distance of 100000 bp\n"); 
	  }
	thinning_distance = 100000;
      } else {
      if(verbose)
	{
	  Rprintf("Thinning SNPs using distance of %d bp\n", thinning_distance); 
	}
    }
    
    /* Allocate vector of int for chromosomes */
    gsl_vector_int * chromosomes = gsl_vector_int_alloc(NSNPS);
    /* Allocate vector of int for bp positions */
    gsl_vector_int * positions = gsl_vector_int_alloc(NSNPS);
    
    /* Length of line */
    int maxcharssnpname = 256; // Maximum number of characters of a SNP name
    int maxcharschrom = 2; // Maximum number of characters for a chromosome is 2 (i.e. for chromosomes 10..22)
    int maxcharspos = 11; // Maximum number of characters for chromosome position is 11 (i.e. 10,000,000,000 - nb human csome 1 has ~ 250,000,000 bp)
    char line[maxcharssnpname + maxcharschrom + maxcharspos + 3]; // The +3 are for the spaces and end of line character
    char * tmp;
    /* Make a file pointer for thinfile */
    FILE * thinfile = fopen(thinfilename, "r");
    count = -1;
    /* Read in thinfile */
    if(thinfile != NULL)
      {
	while ( fgets ( line, sizeof line, thinfile ) != NULL ) /* read a line */
	  {
	    count++;
	    /* Extract SNP name */
	    tmp = strtok(line, " ");
	    if(strcmp(SNPNAMES[count],tmp) != 0)
	      {
		error("SNPnames in genotype file and thinfile do not match (%s vs %s)\n", SNPNAMES[count], tmp);
	      }
	    /* Extract chromosome */
	    tmp = strtok(NULL, " ");
	    gsl_vector_int_set(chromosomes, count, atoi(tmp));
	    /* Extract SNP position */
	    tmp = strtok(NULL, "\n");
	    gsl_vector_int_set(positions, count, atoi(tmp));
	  }
      } else {
      error("could not open %s for reading\n", thinfilename);
    }
    /* Close thinfile */
    fclose(thinfile);
    
    /* Choose which SNPs to thin by: follow method in plinkcomp */
    int current_chr = -1;
    int start_of_interval = 0;
    int end_of_interval = 0;
    int first_snp = 0; // boolean
    int this_snp_chr = 0;
    
    for(count = 0; count < NSNPS; count++)
      {
	this_snp_chr = gsl_vector_int_get(chromosomes, count);
	if(current_chr != this_snp_chr)
	  {
	    current_chr = this_snp_chr;
	    start_of_interval = 0;
	    end_of_interval = start_of_interval + thinning_distance;
	    first_snp = 1; // boolean
	  }
	if(current_chr == 0)
	  {
	    // Do nothing
	  }
	else
	  {
	    if(first_snp)
	      {
		
		gsl_vector_int_set(thin, count, 1);
		start_of_interval = gsl_vector_int_get(positions, count);
		end_of_interval = start_of_interval + thinning_distance;
		first_snp = 0;
	      }
	    if((gsl_vector_int_get(positions, count) >= end_of_interval))
	      {
		gsl_vector_int_set(thin, count, 1);
		start_of_interval = gsl_vector_int_get(positions, count);
		end_of_interval = start_of_interval + thinning_distance;
	      } else {
	      // Do nothing
	    }
	    
	  } // ends else 
      } // ends for loop
  } // ends else (i.e. if thinfilename != NULL)
  * nThinnedSnps = sumIntVec(thin);
  return thin;
}


int readSNPsThinAndComputePCs(char * genofilename,
			      gsl_vector_int * thin,
			      GSL_TYPE(matrix) * Z,
			      GSL_TYPE(matrix) * thinnedGenotypes,
			      GSL_TYPE(vector) * D2,
			      int * howManyK)
{

  /* Allocate memory for the thinned genotypes */
  int NSNPS = thin->size;
  int nThinnedSnps = sumIntVec(thin);
  int NINDIV = Z->size1;
  int i = 0;
  int j = 0;
  int count = 0;

  /* Allocate a bit integer matrix to read all of the SNPs in */
  gsl_matrix_int * genotypes = gsl_matrix_int_calloc(NINDIV, NSNPS);

  GSL_TYPE(matrix) * U = GSL_FUNCTION(matrix, calloc)(Z->size1, Z->size2);
  GSL_TYPE(matrix) * V = GSL_FUNCTION(matrix, calloc)(nThinnedSnps, Z->size2);
  GSL_TYPE(vector) * D = GSL_FUNCTION(vector, calloc)(Z->size2);

  /* open genofile for reading */
  FILE * genofile = NULL;
  genofile = fopen(genofilename, "r");
  /* skip the header row */
  char ch;
  ch = fgetc(genofile);
  while(ch != (int)'\n')
    {
      ch = fgetc(genofile);
    }
  /* Read in genotypes as int */ 
  gsl_matrix_int_fscanf(genofile, genotypes);
  fclose(genofile);
  
  /* Temporary variables for counting */
  int tmp = 0;
  for(i = 0; i < NSNPS; i++)
   /* Read the predictor into a vector of PREC */
    {
      tmp = gsl_vector_int_get(thin, i);
      if(tmp == 1)
	{
	  /* vector view of genotypes */
	  gsl_vector_int_view predictorsCol = gsl_matrix_int_column(genotypes, i);
	  /* vector view of thinned genotypes */
	  GSL_FUNCTION(vector,view) thinnedGenotypesCol = GSL_FUNCTION(matrix,column)(thinnedGenotypes, count);
	  /* Convert int to PREC and put in thinnedGenotypes matrix */
	  convert_int_vector(&predictorsCol.vector, &thinnedGenotypesCol.vector);
	  count++;
	}
    }

  /* Free the genotypes matrix */
  gsl_matrix_int_free(genotypes);

  /* scale the genotypes that have been read in */
  /* some temporary variables for the mean and sd */
  PREC mean = 0;
  PREC sd = 0;
  PREC sqrtn1 = MATHS_FUNCTION(sqrt)((PREC)NINDIV - 1.0);
  /* vector view for the column of genotypes matrix */
  for(j = 0; j < nThinnedSnps; j++)
    {
      /* Get a vector view of the column of the genotypes */
      GSL_FUNCTION(vector,view) genotype_vector = GSL_FUNCTION(matrix, column)(thinnedGenotypes, j);
      /* Compute the mean */
      mean = GSL_STATS_FUNCTION(mean)(genotype_vector.vector.data, genotype_vector.vector.stride, genotype_vector.vector.size);
      /* Compute the standard deviation */
      sd = GSL_STATS_FUNCTION(sd)(genotype_vector.vector.data, genotype_vector.vector.stride, genotype_vector.vector.size);
      /* Subtract the mean */
      GSL_FUNCTION(vector, add_constant)(&genotype_vector.vector, -1*mean);
      /* Divide by the standard deviation corrected to correlation form */
      GSL_FUNCTION(vector, scale)(&genotype_vector.vector, 1 / (sd * sqrtn1));
    }
  /* SVD them */ 
  SVD_FUNCTION(thinnedGenotypes, U, V, D);
  /* Compute D2 */
  GSL_FUNCTION(vector,memcpy)(D2, D);
  GSL_FUNCTION(vector,mul)(D2, D);

  /* Compute Z */
  BLAS_FUNCTION(gemm)(CblasNoTrans, CblasNoTrans, 1.0, thinnedGenotypes, V, 0.0, Z);
  if(*howManyK == 0)
    {
      *howManyK = chooseHowManyK(D);
    }

  /* Free U, V and D */
  GSL_FUNCTION(matrix, free)(U);
  GSL_FUNCTION(matrix, free)(V);
  GSL_FUNCTION(vector, free)(D);

  return 0;
}

#endif

typedef int make_iso_compilers_happy;

