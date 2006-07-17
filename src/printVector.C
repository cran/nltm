#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

#define VF 999 // vectorial function

void printDV(ofstream *ofs, vector<double> &a)
{
  int i;

  for(i=0; i<int(a.size()); i++)
    *ofs<<a[i]<<" ";

  *ofs<<endl;
}

void printDM(ofstream *ofs, vector<vector<double> > &a)
{
  int i;

  for(i=0; i<int(a.size()); i++)
    printDV(ofs, a[i]);
}


void printIV(ofstream *ofs, vector<int> &a)
{
  int i;

  for(i=0; i<int(a.size()); i++)
    *ofs<<a[i]<<" ";
  *ofs<<endl;
}

void printIM(ofstream *ofs, vector<vector<int> > &a)
{
  int i;

  for(i=0; i<int(a.size()); i++)
    printIV(ofs, a[i]);
}




void printDVector(ofstream *ofs, double *a, int n)
{
  int i;

  for(i=0; i<n; i++)
    *ofs<<a[i]<<" ";

  *ofs<<endl;
}

void printDMatrix(ofstream *ofs, double **a, int nrow, int ncol)
{
  int i;

  for(i=0; i<nrow; i++)
    printDVector(ofs, a[i], ncol);
}


void printIVector(ofstream *ofs, int *a, int n)
{
  int i;

  for(i=0; i<n; i++)
    *ofs<<a[i]<<" ";
  *ofs<<endl;
}

void printIMatrix(ofstream *ofs, int **a, int nrow, int ncol)
{
  int i;

  for(i=0; i<nrow; i++)
    printIVector(ofs, a[i], ncol);
}




// Print in R format
void printDMRformat(ofstream *ofs, vector<vector<double> > &a)
{
  int i, j;

  *ofs<<"matrix(c(";
  for(i=0; i<int(a.size())-1; i++)
    for(j=0; j<int(a[0].size()); j++)
      *ofs<<a[i][j]<<", ";
  for(j=0; j<int(a[0].size())-1; j++)
    *ofs<<a[a.size()-1][j]<<", ";
  *ofs<<a[a.size()-1][a[0].size()-1]<<"), nrow="<<a.size()<<", ncol="
      <<a[0].size()<<", byrow=TRUE)"<<endl;
}


void printDMatrixRformat(ofstream *ofs, double **a, int nrow, int ncol)
{
  int i, j;

  *ofs<<"matrix(c(";
  for(i=0; i<nrow-1; i++)
    for(j=0; j<ncol; j++)
      *ofs<<a[i][j]<<", ";
  for(j=0; j<ncol-1; j++)
    *ofs<<a[nrow-1][j]<<", ";
  *ofs<<a[nrow-1][ncol-1]<<"), nrow="<<nrow<<", ncol="<<ncol<<", byrow=TRUE)"
      <<endl;

}

void printModelFunction(ofstream *ofs, string fName, vector<double> &pred, 
			double s, int cc, double resU, vector<double> resVF)
{
  int i;

  *ofs<<fName<<" ";
  for(i=0; i<int(pred.size()); i++)
    *ofs<<pred[i]<<" ";
  *ofs<<s<<" "<<cc<<" ";
  if(resU!=VF)
    *ofs<<resU<<" ";
  else{
    for(i=0; i<int(resVF.size()); i++)
      *ofs<<resVF[i]<<" ";
  }
  *ofs<<endl;
}
