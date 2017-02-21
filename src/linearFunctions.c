#include "linearFunctions.h"
#ifdef HAVE_GSL_HEADER

GSL_TYPE(vector) * readLinearPhenotypes(char * phenotypefilename,
					int NINDIV)
{
  GSL_TYPE(vector) * phen = GSL_FUNCTION(vector,calloc)(NINDIV);
  /* Open a file for reading */
  FILE * phenofile = fopen(phenotypefilename,"r");
  /* Scan the phenotypes into the int vector */
  GSL_FUNCTION(vector,fscanf)(phenofile, phen);
  /* Close the file for reading */
  fclose(phenofile);
  return phen;
}

/*  */

/* compute Kr based on a, Z, y and r */

int computeLinearKr(GSL_TYPE(vector) * a,
		    GSL_TYPE(matrix) * Z,
		    GSL_TYPE(vector) * y,
		    GSL_TYPE(vector) * D2,
		    int r,
		    PREC * kr,
		    PREC * DofF)
{
  // n is the number of individuals
  int n = y->size;
  // A vector view of a
  GSL_FUNCTION(vector,view) a_view = GSL_FUNCTION(vector,subvector)(a, 0, r);
  // The denominator of Kr
  PREC denom = 0.0;
  // Compute it
  GSL_BLAS_FUNCTION(dot)(&a_view.vector, &a_view.vector, &denom);
  // A matrix view of W
  GSL_FUNCTION(matrix,view) Zview = GSL_FUNCTION(matrix,submatrix)(Z, 0, 0, n, r);
  // Make Wview * a_view
  // A vector for the residuals
  GSL_TYPE(vector) * resid = GSL_FUNCTION(vector,calloc)(n);
  // Make the fitted ys and put them in the resid vector
  BLAS_FUNCTION(gemv)(CblasNoTrans, 1.0, &Zview.matrix, &a_view.vector, 0.0, resid);
  // Make the residual vector 
  GSL_FUNCTION(vector,scale)(resid, -1);
  GSL_FUNCTION(vector,add)(resid, y);
  // Make the crossproduct
  PREC crossprod = 0.0;
  GSL_BLAS_FUNCTION(dot)(resid, resid, &crossprod);
  /* Divide it by n - r */
  crossprod = crossprod / ((PREC)n - (PREC) r);
  // times it by i
  PREC numerator = ((PREC) r) * crossprod;
  *kr = numerator / denom;
  GSL_FUNCTION(vector,free)(resid);
  /* Compute DofF */
  computeDofF(D2, *kr, DofF);
  return 0;
}

#endif

typedef int make_iso_compilers_happy;
