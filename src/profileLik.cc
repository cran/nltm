#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <R.h>
#include <stdlib.h>


using namespace std;

// dd: number of deaths
// rr: number of censored or dead

#define CENSOR 0
#define FAILURE 1
#define ITMAX 100000 // Maximum number of iterations for calculation of hh
#define VERYBIG 1e30

double Theton(vector<double> &pred, double x, int cc, int model);
double ThetonCure(vector<double> &pred, double s, int cc, int model);
double vtheta(vector<double> &pred, double x, int cc, int model);
double vthetaCure(vector<double> &pred, double s, int cc, int model);
void printDV(ofstream *ofs, vector<double> &a);
void printDM(ofstream *ofs, vector<vector<double> > &a);
void printDVector(ofstream *ofs, double *a, int n);
void printDMatrix(ofstream *ofs, double **a, int nrow, int ncol);
void printIVector(ofstream *ofs, int *a, int n);
double **dmat(double *array, int nrow, int ncol);
int nmodel(string model);

ofstream ofsDebug;

// pred: nn x npred
// beta=(beta_theta, beta_eta, beta_cure)
void predictor(double **xx1, double **xx2, int nvar1, int nvar2, double *beta, 
	       int cure, vector<vector<double> > &pred)
{
  int i, j, icure=nvar1+nvar2;
 
  for(i=0; i<int(pred.size()); i++){
    pred[i][0]=0;
    for(j=0; j<nvar1; j++)
      pred[i][0]+=xx1[i][j]*beta[j];
    pred[i][0]=exp((cure ? pred[i][0]+beta[icure] : pred[i][0]));
  }

  if(pred[0].size()>1){
    for(i=0; i<int(pred.size()); i++){
      pred[i][1]=0;
      for(j=0; j<nvar2; j++)
	pred[i][1]+=xx2[i][j]*beta[j+nvar1];
      pred[i][1]=exp(pred[i][1]);
    }
  }
}


// s0=exp(-h)
void survivalJump(int *status, int *dd, int *rr, vector<vector<double> > &pred,
		  int model, int cure, vector<double> &s0, int verbose)
{
  int i=0, j, k, nn=pred.size(), nt=s0.size();
  double ss=1, sum=0;
  vector<double> ThetonVal(nn);

  for(j=0; j<nt-cure; j++){
    ss*=s0[j];
    for(k=0; k<rr[j]; k++){
      ThetonVal[i]=Theton(pred[i], ss, status[i], model);
      i++;
    }
  }
  if(cure){
    for(k=0; k<rr[nt-1]; k++){
      ThetonVal[i]=ThetonCure(pred[i], ss, status[i], model);
      i++;
    }
  }
  
  i=nn-1;
  if(cure){
    s0[nt-1]=0;
    for(k=0; k<rr[nt-1]; k++)
      sum+=ThetonVal[i--];
  }
  for(j=nt-1-cure; j>=0; j--){
    for(k=0; k<rr[j]; k++)
      sum+=ThetonVal[i--];
    s0[j]=exp(-(double)dd[j]/sum);
  }
}


void checkSelfConsistency(int *status, int *dd, int *rr, 
			  vector<vector<double> > &pred, int model, int cure,
			  vector<double> &s0)
{
  int i=0, j, k, nn=pred.size(), nt=s0.size();
  double ss=1, sum=0, max=-1;
  vector<double> ThetonVal(nn);

  for(j=0; j<nt-cure; j++){
    ss*=s0[j];
    for(k=0; k<rr[j]; k++){
      ThetonVal[i]=Theton(pred[i], ss, status[i], model);
      i++;
    }
  }
  if(cure){
    for(k=0; k<rr[nt-1]; k++){
      ThetonVal[i]=ThetonCure(pred[i], ss, status[i], model);
      i++;
    }
  }

  i=nn-1;
  if(cure){
    s0[nt-1]=0;
    for(k=0; k<rr[nt-1]; k++)
      sum+=ThetonVal[i--];
  }
  for(j=nt-1-cure; j>=0; j--){
    for(k=0; k<rr[j]; k++)
      sum+=ThetonVal[i--];
    s0[j]=fabs(exp(-(double)dd[j]/sum)-s0[j]);
    max=(s0[j]>max ? s0[j] : max);
  }
}


// Find hh satisfying the self-consistency equation, return exp(-hh) 
// Function FitFfixedP in [SP] 
void fitSurvival(int *status, int *dd, int *rr, vector<vector<double> > &pred, 
		 int model, int cure, double tol, double *s0, int nt, 
		 int verbose)
{
  int i, j;
  double sum=2*tol;
  vector<double> s0Aux;

  s0Aux.resize(nt);
  for(j=0; j<nt; j++)
    s0Aux[j]=s0[j];

  for(i=0; i<ITMAX & sum>tol; i++){ // Some warning should come if i==ITMAX
    sum=0;
    survivalJump(status, dd, rr, pred, model, cure, s0Aux, verbose);

    for(j=0; j<nt; j++)
      sum+=fabs(s0Aux[j]-s0[j]);

    for(j=0; j<nt; j++)
      s0[j]=s0Aux[j];
  }
}


// log likelihood
double likelihood(int *status, int *dd, int *rr, int model, int cure, 
		  double *s0, vector<vector<double> > &pred, int nt)
{
  int i=0, j, k;
  double ss=1, lik=0, vt;

  for(j=0; j<nt-cure; j++){
    ss*=s0[j];
    for(k=0; k<rr[j]; k++){
      vt=vtheta(pred[i], ss, status[i], model);
      lik+=log(vt);
      i++;
    }
    lik+=dd[j]*log(-log(s0[j]));
  } 
  
  if(cure){
    for(k=0; k<rr[nt-1]; k++){
      vt=vthetaCure(pred[i], ss, status[i], model);
      lik+=log(vt);
      i++;
    }
  }
  return(lik);
}


// randon number generator that works both in Windows and Unix
double jrnd(void) 
{
  double rv = ((double) rand())/RAND_MAX;
  return(rv);
}


// x1: covariates matrix for long term predictor
// x2: covariates matrix for short term predictor
// nvar: number of variables once categorical variables have been dichotomized
extern "C" {
// profile likelihood
void profileLik(double *beta, double *x1, double *x2, int *status, int *dd, 
		int *rr, double *s0, char **survModel, int *cure, double *tol, 
		int *nvar1, int *nvar2, int *ntime, int *nobs, int *npred, 
		int *verbose, double *plik) 
{
  int i, nt, nn, model, nbeta;
  double **xx1, **xx2;
  vector<vector<double> > pred;

  nbeta=(*nvar1)+(*nvar2)+(*cure);
  nt=*ntime;
  nn=*nobs;
  model=nmodel(*survModel);

  if(*nvar1>0)
    xx1=dmat(x1, nn, *nvar1);
  if(*npred>1 && *nvar2>0)
    xx2=dmat(x2, nn, *nvar2);

  if(verbose){
    ofsDebug<<"nn: "<<nn<<" nvar1: "<<*nvar1<<" nvar2: "<<*nvar2<<endl;
    ofsDebug<<"beta "<<nbeta<<endl;
    printDVector(&ofsDebug, beta, nbeta);
    
    if(*nvar1>0){
      ofsDebug<<"xx1"<<endl;
      printDMatrix(&ofsDebug, xx1, nn, *nvar1);
    }
    if(*npred>1 && *nvar2>0){
      ofsDebug<<"xx2"<<endl;
      printDMatrix(&ofsDebug, xx2, nn, *nvar2);
    }
    ofsDebug<<"dd"<<endl;
    printIVector(&ofsDebug, dd, nt);
    ofsDebug<<"rr"<<endl;
    printIVector(&ofsDebug, rr, nt);
    ofsDebug<<"status"<<endl;
    printIVector(&ofsDebug, status, nt);
    
    ofsDebug<<"s0"<<endl;
    printDVector(&ofsDebug, s0, nt);
  }
    
  pred.resize(nn);
  for(i=0; i<int(pred.size()); i++)
    pred[i].resize(*npred);
  predictor(xx1, xx2, *nvar1, *nvar2, beta, *cure, pred);
  
  fitSurvival(status, dd, rr, pred, model, *cure, *tol, s0, nt, *verbose);

  if(verbose){
    ofsDebug<<"s0 ";
    printDVector(&ofsDebug, s0, nt);
  }
  
  *plik=likelihood(status, dd, rr, model, *cure, s0, pred, nt);
    
  if(verbose)
    ofsDebug<<"plik: "<<*plik<<endl;

// Note: drand48 doens't work in Windows  
//   *plik=(isinf(*plik) || isnan(*plik) ? -VERYBIG*(1+drand48()*0.1) : *plik);
//  *plik=(isinf(*plik) || isnan(*plik) ? -VERYBIG*(1+jrnd()*0.1) : *plik);
  *plik = !R_FINITE(*plik) ? -VERYBIG*(1+jrnd()*0.1) : *plik;

  if(verbose)
    ofsDebug<<"plik: "<<*plik<<endl;
}
}


extern "C" {
void openDebug(char **fileDebug)
{
  ofsDebug.open(*fileDebug);
  if(!ofsDebug) {
    cerr<<"ERROR: couldn't create debug file."<<endl;
    exit(1);
  }
  ofsDebug<<"start"<<endl;
}
}


extern "C" {
void closeDebug()
{
  ofsDebug.close();
}
}
