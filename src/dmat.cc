#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;


// array=(row1, row2, ...)
double **dmat(double *array, int nrow, int ncol)
{
  register int i;
  register double **pointer;
  
  if((pointer=(double **)malloc(sizeof(double *)*nrow))==NULL) {
    cerr<<"ERROR: malloc failed."<<endl;
    exit(1);
  }
  for(i=0; i<nrow; ++i)
    if((pointer[i]=(double *)malloc(sizeof(double)*ncol))==NULL) {
      cerr<<"ERROR: malloc failed."<<endl;
      exit(1);
    }
  //  pointer = (double **) malloc(nrow, sizeof(double *));
  for (i=0; i<nrow; i++) {
    pointer[i] = array;
    array += ncol;
  }
  return(pointer);
}
