#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _spatialcluster_rcpp_alk(SEXP, SEXP, SEXP);
extern SEXP _spatialcluster_rcpp_clk(SEXP, SEXP, SEXP, SEXP);
extern SEXP _spatialcluster_rcpp_cut_tree(SEXP, SEXP, SEXP, SEXP);
extern SEXP _spatialcluster_rcpp_full_initial(SEXP, SEXP);
extern SEXP _spatialcluster_rcpp_full_merge(SEXP, SEXP, SEXP);
extern SEXP _spatialcluster_rcpp_mst(SEXP);
extern SEXP _spatialcluster_rcpp_slk(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_spatialcluster_rcpp_alk",          (DL_FUNC) &_spatialcluster_rcpp_alk,          3},
    {"_spatialcluster_rcpp_clk",          (DL_FUNC) &_spatialcluster_rcpp_clk,          4},
    {"_spatialcluster_rcpp_cut_tree",     (DL_FUNC) &_spatialcluster_rcpp_cut_tree,     4},
    {"_spatialcluster_rcpp_full_initial", (DL_FUNC) &_spatialcluster_rcpp_full_initial, 2},
    {"_spatialcluster_rcpp_full_merge",   (DL_FUNC) &_spatialcluster_rcpp_full_merge,   3},
    {"_spatialcluster_rcpp_mst",          (DL_FUNC) &_spatialcluster_rcpp_mst,          1},
    {"_spatialcluster_rcpp_slk",          (DL_FUNC) &_spatialcluster_rcpp_slk,          4},
    {NULL, NULL, 0}
};

void R_init_spatialcluster(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
