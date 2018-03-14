#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _spatialcluster_rcpp_get_component_vector(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_spatialcluster_rcpp_get_component_vector", (DL_FUNC) &_spatialcluster_rcpp_get_component_vector, 1},
    {NULL, NULL, 0}
};

void R_init_spatialcluster(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
