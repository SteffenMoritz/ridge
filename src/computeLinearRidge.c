#include "computeLinearRidge.h"
#ifdef HAVE_GSL_HEADER


/* SVD of any matrix */
int svdAnyMat(gsl_matrix * X, 
	      gsl_matrix * U, 
	      gsl_matrix * V, 
	      gsl_vector * D)
{
  // SVD part
  gsl_vector * work;
  int n = X->size1;
  int p = X->size2;
  if(p > n)
    {
      work = gsl_vector_alloc(n);
      gsl_matrix * tmpV = gsl_matrix_alloc(n, n);
      gsl_matrix * tmpU = gsl_matrix_alloc(p, n);
      // Transpose the matrix X
      gsl_matrix_transpose_memcpy(tmpU, X);
      // TmpU contians t(X)
      // There is linear dependence in X
      // Need to replace the NaNs with zeros in the matrix
      // Or something 
      gsl_linalg_SV_decomp(tmpU, tmpV, D, work);
      gsl_vector_free(work);
      // Swap back
      gsl_matrix * tmp1 = gsl_matrix_alloc(tmpU->size1, tmpU->size2);
      gsl_matrix * tmp2 = gsl_matrix_alloc(tmpV->size1, tmpV->size2);
      gsl_matrix_memcpy(tmp1, tmpU);
      gsl_matrix_memcpy(tmp2, tmpV);
      gsl_matrix_free(tmpU);
      gsl_matrix_free(tmpV);
      gsl_matrix_memcpy(V, tmp1);
      gsl_matrix_memcpy(U, tmp2);
      gsl_matrix_free(tmp1);
      gsl_matrix_free(tmp2);
    } else {
    work = gsl_vector_alloc(p);
    gsl_matrix_memcpy(U, X);
    gsl_linalg_SV_decomp(U, V, D, work);
    gsl_vector_free(work);
  }
  return 0;
}

int prepareLambdas(gsl_vector * y, 
		   gsl_matrix * U, 
		   gsl_vector * D2, 
		   gsl_vector * lambdaVeckHKB, 
		   char * skhkbfilename, 
		   char * sklwfilename, 
		   gsl_vector * lambdaVeckLW, 
		   int randomized, 
		   int s)
{
  double kHKB;
  double kLW;
  double crossprod;
  double numerator;
  double denominatorkHKB;
  double denominatorkLW;
  int lengthLambdaVec = lambdaVeckHKB->size;
  gsl_matrix_view Uview; // a matrix view
  int n = y->size;
  int i, j;
  gsl_vector * resid = gsl_vector_alloc(n);
  gsl_matrix * H = gsl_matrix_alloc(n, n);
  for(i = 0; i < lengthLambdaVec; i++)
    {
      gsl_matrix * diag = gsl_matrix_calloc((i+1), (i+1));
      Uview = gsl_matrix_submatrix(U, 0, 0, n, (i + 1));
      // Make the hat matrix
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &Uview.matrix, &Uview.matrix, 0.0, H);
      // make the fitted ys - put in the resid vector
      gsl_blas_dgemv(CblasNoTrans, 1.0, H, y, 0.0, resid);
      // make the denominaotor for kLW
      if(sklwfilename != NULL)
	{
	  gsl_blas_ddot(y, resid, &denominatorkLW);
	}
      // Make the residual vector 
      gsl_vector_scale(resid, -1);
      gsl_vector_add(resid, y);
      // make the crossproduct
      gsl_blas_ddot(resid, resid, &crossprod);
      // times it by i
      numerator = crossprod * ((float) i + 1.0);
      // this gives the numerator
      // Make the denominator for kHKB
      // Make the diagonal matrix
      for(j = 0; j < diag->size1; j++)
	{
	  gsl_matrix_set(diag, j, j, 1.0 / gsl_vector_get(D2, j));
	}
      // 
      // Make the matrix U diag D2
      gsl_matrix * UD2 = gsl_matrix_alloc(n, (i + 1));
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &Uview.matrix, diag, 0.0, UD2);
      // Make the matrix U diag D2 U' - put it into H
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, UD2, &Uview.matrix, 0.0, H);
      // Make the matrix U diag D2 U' y - put it into resid
      gsl_blas_dgemv(CblasNoTrans, 1.0, H, y, 0.0, resid);
      // Make the dot product
      gsl_blas_ddot(y, resid, &denominatorkHKB);
      // put in the matrix
      if(skhkbfilename != NULL)
	{
	  gsl_blas_ddot(y, resid, &denominatorkHKB);
	  denominatorkHKB = ((float) n - (float) i - 1.0) * denominatorkHKB;
	  kHKB = numerator / denominatorkHKB;
	  gsl_vector_set(lambdaVeckHKB, i, kHKB);
	}
      if(sklwfilename != NULL)
	{
	  denominatorkLW = ((float) n - (float) i - 1.0) * denominatorkLW;
	  kLW = numerator / denominatorkLW;
	  gsl_vector_set(lambdaVeckLW, i, kLW);
	}
      gsl_matrix_free(UD2);
      gsl_matrix_free(diag);
    }
  if(randomized)
    {
      gsl_rng * rndm = gsl_rng_alloc(gsl_rng_mt19937);
      double weight;
      gsl_rng_set(rndm, s);
      for(i=0; i<lambdaVeckHKB->size; i++)
	{
	  weight = gsl_ran_flat(rndm, 0.2, 1.0);
	  gsl_vector_set(lambdaVeckHKB, i, weight * gsl_vector_get(lambdaVeckHKB, i));
	  weight = gsl_ran_flat(rndm, 0.2, 1.0);
	  gsl_vector_set(lambdaVeckLW, i, weight * gsl_vector_get(lambdaVeckLW, i));
	}
      gsl_rng_free(rndm);
    }
  gsl_vector_free(resid);
  gsl_matrix_free(H);
  return 0;
}

/* /\* and one function for which parts are different *\/  */
/* gsl_vector * computeLinearRidge(gsl_matrix * VoneOverS, gsl_vector * d2, gsl_vector * ahat, gsl_vector * Bridge, double shrinkage) */
/* { */
/*   int i; */
/*   gsl_vector * aridge = gsl_vector_alloc(ahat->size); */
/*   gsl_vector_memcpy(aridge, ahat);  */
/*   gsl_vector * d2tmp = gsl_vector_alloc(d2->size); */
/*   gsl_vector_memcpy(d2tmp, d2); */
/*   gsl_vector_mul(aridge, d2); */
/*   gsl_vector_add_constant(d2tmp, shrinkage); */
/*   gsl_vector_div(aridge, d2tmp); */
/*   gsl_blas_dgemv(CblasNoTrans, 1.0, VoneOverS, aridge, 0.0, Bridge); */
/*   gsl_vector_free(aridge); */
/*   gsl_vector_free(d2tmp); */
/*   return Bridge; */
/* } */

#endif

typedef int make_iso_compilers_happy;

