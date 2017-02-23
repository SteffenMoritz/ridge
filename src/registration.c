#include <stdlib.h> // for NULL
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void regression_wrapper_function(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"regression_wrapper_function", (DL_FUNC) &regression_wrapper_function, 11},
    {NULL, NULL, 0}
};

void R_init_ridge(DllInfo *dll)
{
  
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  
    R_useDynamicSymbols(dll, FALSE);
    
    R_forceSymbols(dll, TRUE);
    
    
}

