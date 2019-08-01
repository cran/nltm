#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern void profileLik(double *beta, double *x1, double *x2, int *status, int *dd, 
		       int *rr, double *s0, char **survModel, int *cure, double *tol, 
		       int *nvar1, int *nvar2, int *ntime, int *nobs, int *npred, 
		       int *verbose, double *plik);

extern void informationMatrix(double *beta, double *x1, double *x2, int *status, 
		       int *dd, int *rr, double *s0, char **survModel, 
		       int *cure, int *nvar1, int *nvar2, int *ntime, 
		       int *nobs, int *npred, int *verbose, double *imat);

static R_NativePrimitiveArgType profileLikT[] =
  { REALSXP, REALSXP, REALSXP, INTSXP, INTSXP,
    INTSXP, REALSXP, STRSXP, INTSXP, REALSXP,
    INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
    INTSXP, REALSXP};

static R_NativePrimitiveArgType informationMatrixT[] =
  { REALSXP, REALSXP, REALSXP, INTSXP, 
    INTSXP, INTSXP, REALSXP, STRSXP, 
    INTSXP, INTSXP, INTSXP, INTSXP,
    INTSXP, INTSXP, INTSXP, REALSXP};

static const R_CMethodDef CEntries[] = {
  {"profileLik",        (DL_FUNC) &profileLik,        17, profileLikT},
  {"informationMatrix", (DL_FUNC) &informationMatrix, 16, informationMatrixT},
  {NULL, NULL, 0, NULL}
};

void R_init_nltm(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
