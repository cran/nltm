#include <R.h>

// array=(row1, row2, ...)
double **dmat(double *array, int nrow, int ncol)
{
  int i;
  double **pointer;
  
  pointer=(double **)R_alloc(nrow+1, sizeof(double *));
  for(i=0; i<nrow; ++i)
    pointer[i]=(double *)R_alloc(ncol+1,sizeof(double));
  for (i=0; i<nrow; i++) {
    pointer[i] = array;
    array += ncol;
  }
  return(pointer);
}
