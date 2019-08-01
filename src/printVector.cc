#include <vector>
#include "Rostream.h"

using namespace std;

#define VF 999 // vectorial function

void printDV(vector<double> &a)
{
  int i;

  for(i=0; i<int(a.size()); i++)
    Rcout<<a[i]<<" ";

  Rcout<<endl;
}

void printDM(vector<vector<double> > &a)
{
  int i;

  for(i=0; i<int(a.size()); i++)
    printDV(a[i]);
}


void printIV(vector<int> &a)
{
  int i;

  for(i=0; i<int(a.size()); i++)
    Rcout<<a[i]<<" ";
  Rcout<<endl;
}

void printIM(vector<vector<int> > &a)
{
  int i;

  for(i=0; i<int(a.size()); i++)
    printIV(a[i]);
}




void printDVector(double *a, int n)
{
  int i;

  for(i=0; i<n; i++)
    Rcout<<a[i]<<" ";

  Rcout<<endl;
}

void printDMatrix(double **a, int nrow, int ncol)
{
  int i;

  for(i=0; i<nrow; i++)
    printDVector(a[i], ncol);
}


void printIVector(int *a, int n)
{
  int i;

  for(i=0; i<n; i++)
    Rcout<<a[i]<<" ";
  Rcout<<endl;
}

void printIMatrix(int **a, int nrow, int ncol)
{
  int i;

  for(i=0; i<nrow; i++)
    printIVector(a[i], ncol);
}




// Print in R format
void printDMRformat(vector<vector<double> > &a)
{
  int i, j;

  Rcout<<"matrix(c(";
  for(i=0; i<int(a.size())-1; i++)
    for(j=0; j<int(a[0].size()); j++)
      Rcout<<a[i][j]<<", ";
  for(j=0; j<int(a[0].size())-1; j++)
    Rcout<<a[a.size()-1][j]<<", ";
  Rcout<<a[a.size()-1][a[0].size()-1]<<"), nrow="<<a.size()<<", ncol="
      <<a[0].size()<<", byrow=TRUE)"<<endl;
}


void printDMatrixRformat(double **a, int nrow, int ncol)
{
  int i, j;

  Rcout<<"matrix(c(";
  for(i=0; i<nrow-1; i++)
    for(j=0; j<ncol; j++)
      Rcout<<a[i][j]<<", ";
  for(j=0; j<ncol-1; j++)
    Rcout<<a[nrow-1][j]<<", ";
  Rcout<<a[nrow-1][ncol-1]<<"), nrow="<<nrow<<", ncol="<<ncol<<", byrow=TRUE)"
      <<endl;

}

void printModelFunction(string fName, vector<double> &pred, 
			double s, int cc, double resU, vector<double> resVF)
{
  int i;

  Rcout<<fName<<" ";
  for(i=0; i<int(pred.size()); i++)
    Rcout<<pred[i]<<" ";
  Rcout<<s<<" "<<cc<<" ";
  if(resU!=VF)
    Rcout<<resU<<" ";
  else{
    for(i=0; i<int(resVF.size()); i++)
      Rcout<<resVF[i]<<" ";
  }
  Rcout<<endl;
}
