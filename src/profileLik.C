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
void predictor(double **xx, int nvar, double *beta, int cure,
	       vector<vector<double> > &pred)
{
  int i, j, icure=pred[0].size()*nvar;
 
  //  ofsDebug<<"begin pred "<<nvar<<" "<<cure<<endl;

  for(i=0; i<int(pred.size()); i++){
    pred[i][0]=0;
    for(j=0; j<nvar; j++)
      pred[i][0]+=xx[i][j]*beta[j];
    pred[i][0]=exp((cure ? pred[i][0]+beta[icure] : pred[i][0]));
  }

  if(pred[0].size()>1){
    for(i=0; i<int(pred.size()); i++){
      pred[i][1]=0;
      for(j=0; j<nvar; j++)
	pred[i][1]+=xx[i][j]*beta[j+nvar];
      pred[i][1]=exp(pred[i][1]);
    }
  }
}



// s0=exp(-h)
void survivalJump(int *status, int *dd, int *rr, vector<vector<double> > &pred,
		  int model, int cure, vector<double> &s0)
{
  int i=0, j, k, nn=pred.size(), nt=s0.size();
  double ss=1, sum=0;
  vector<double> ThetonVal(nn);

#if 0
  ofsDebug<<"survivalJump"<<endl;
  ofsDebug<<"cure: "<<cure<<" nn: "<<nn<<endl;
  ofsDebug<<"s0: ";
  printDV(&ofsDebug, s0);
  ofsDebug<<"dd"<<endl;
  printIVector(&ofsDebug, dd, nt);
  ofsDebug<<"rr"<<endl;
  printIVector(&ofsDebug, rr, nt);
  ofsDebug<<"status"<<endl;
  printIVector(&ofsDebug, status, nn);
#endif
  //  int l;
  for(j=0; j<nt-cure; j++){
    ss*=s0[j];
    for(k=0; k<rr[j]; k++){
      ThetonVal[i]=Theton(pred[i], ss, status[i], model);
//       ofsDebug<<i<<" "<<status[i]<<" "<<ss<<" ";
//       for(l=0; l<int(pred[i].size()); l++)
// 	ofsDebug<<pred[i][l]<<" ";
//       ofsDebug<<ThetonVal[i]<<endl;
      i++;
    }
  }
  if(cure){
    for(k=0; k<rr[nt-1]; k++){
      ThetonVal[i]=ThetonCure(pred[i], ss, status[i], model);
//       ofsDebug<<i<<" "<<status[i]<<" "<<ss<<" ";
//       for(l=0; l<int(pred[i].size()); l++)
// 	ofsDebug<<pred[i][l]<<" ";
//       ofsDebug<<ThetonVal[i]<<" cure"<<endl;
      i++;
    }
  }

//   ofsDebug<<"survivalJump: ThetonVal"<<endl;
//   printDV(&ofsDebug, ThetonVal);

  
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
//   ofsDebug<<"survivalJump: s0"<<endl;
//   printDV(&ofsDebug, s0);
//   ofsDebug<<"end survivalJump"<<endl;
}


void checkSelfConsistency(int *status, int *dd, int *rr, 
			  vector<vector<double> > &pred, int model, int cure,
			  vector<double> &s0)
{
  int i=0, j, k, nn=pred.size(), nt=s0.size();
  double ss=1, sum=0, max=-1;
  vector<double> ThetonVal(nn);

#if 0
  ofsDebug<<"checkSelfConsistency"<<endl;
  ofsDebug<<"cure: "<<cure<<" nn: "<<nn<<endl;
  ofsDebug<<"s0: ";
  printDV(&ofsDebug, s0);
  ofsDebug<<"dd"<<endl;
  printIVector(&ofsDebug, dd, nt);
  ofsDebug<<"rr"<<endl;
  printIVector(&ofsDebug, rr, nt);
  ofsDebug<<"status"<<endl;
  printIVector(&ofsDebug, status, nn);
#endif
  //  int l;
  for(j=0; j<nt-cure; j++){
    ss*=s0[j];
    for(k=0; k<rr[j]; k++){
      ThetonVal[i]=Theton(pred[i], ss, status[i], model);
//       ofsDebug<<i<<" "<<status[i]<<" "<<ss<<" ";
//       for(l=0; l<int(pred[i].size()); l++)
// 	ofsDebug<<pred[i][l]<<" ";
//       ofsDebug<<ThetonVal[i]<<endl;
      i++;
    }
  }
  if(cure){
    for(k=0; k<rr[nt-1]; k++){
      ThetonVal[i]=ThetonCure(pred[i], ss, status[i], model);
//       ofsDebug<<i<<" "<<status[i]<<" "<<ss<<" ";
//       for(l=0; l<int(pred[i].size()); l++)
// 	ofsDebug<<pred[i][l]<<" ";
//       ofsDebug<<ThetonVal[i]<<" cure"<<endl;
      i++;
    }
  }

//   ofsDebug<<"checkSelfConsistency: ThetonVal"<<endl;
//   printDV(&ofsDebug, ThetonVal);

  
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

#if 0
  ofsDebug<<"checkSelfConsistency: s0"<<endl;
  printDV(&ofsDebug, s0);
  ofsDebug<<"Max diff: "<<max<<endl;
  ofsDebug<<"end checkSelfConsistency"<<endl;
#endif
}

// Find hh satisfying the self-consistency equation, return exp(-hh) 
// Function FitFfixedP in [SP] 
void fitSurvival(int *status, int *dd, int *rr, vector<vector<double> > &pred, 
		 int model, int cure, double tol, double *s0, int nt)
{
  int i, j;
  double sum=2*tol;
  vector<double> s0Aux;

  s0Aux.resize(nt);
  for(j=0; j<nt; j++)
    s0Aux[j]=s0[j];

#if 0
  ofsDebug<<"fitSurvival"<<endl;
  ofsDebug<<"s0Aux: ";
  printDV(&ofsDebug, s0Aux);
  ofsDebug<<"dd"<<endl;
  printIVector(&ofsDebug, dd, nt);
  ofsDebug<<"rr"<<endl;
  printIVector(&ofsDebug, rr, nt);
  ofsDebug<<"status"<<endl;
  printIVector(&ofsDebug, status, nt);
#endif

  for(i=0; i<ITMAX & sum>tol; i++){ // Some warning should come if i==ITMAX
    sum=0;
    survivalJump(status, dd, rr, pred, model, cure, s0Aux);
//     ofsDebug<<"fitSurvival: s0Aux"<<endl;
//     printDV(&ofsDebug, s0Aux);

    for(j=0; j<nt; j++)
      sum+=fabs(s0Aux[j]-s0[j]);

    //    ofsDebug<<"i: "<<i<<" diff: "<<sum<<endl;

    for(j=0; j<nt; j++)
      s0[j]=s0Aux[j];
  }
#if 0
  ofsDebug<<"fitSurvival iter: "<<i<<" tol: "<<tol<<endl;
  checkSelfConsistency(status, dd, rr, pred, model, cure, s0Aux);
#endif
}



// likelihood
double likelihood(int *status, int *dd, int *rr, int model, int cure, 
		  double *s0, vector<vector<double> > &pred, int nt)
{
  int i=0, j, k;
  double ss=1, lik=0, vt;

  //  ofsDebug<<"likelihood: vt"<<endl;
  for(j=0; j<nt-cure; j++){
    ss*=s0[j];
    for(k=0; k<rr[j]; k++){
      vt=vtheta(pred[i], ss, status[i], model);
      //      ofsDebug<<j<<" "<<k<<" "<<vt<<endl;
      lik+=log(vt);
      i++;
    }
    lik+=dd[j]*log(-log(s0[j]));
  } 
  
  if(cure){
    for(k=0; k<rr[nt-1]; k++){
      vt=vthetaCure(pred[i], ss, status[i], model);
      //      ofsDebug<<j<<" "<<k<<" "<<vt<<endl;
      lik+=log(vt);
      i++;
    }
  }


  //  cout<<"lik: "<<lik<<endl;

  return(lik);
}


// x: covariates matrix
// nvar: number of variables once categorical variables have been dichotomized
extern "C" {
// profile likelihood
void profileLik(double *beta, double *x, int *status, int *dd, int *rr,
		double *s0, char **survModel, int *cure, double *tol, 
		int *nvar1, int *ntime, int *nobs, int *npred, int *verbose,
		double *plik) 
{
  int i, nvar, nt, nn, model, nbeta;//, isinfLik;
  double **xx;
  vector<vector<double> > pred;

  //  cout<<model<<endl;
  nvar=*nvar1;
  nbeta=nvar*(*npred)+(*cure);
  nt=*ntime;
  nn=*nobs;
  model=nmodel(*survModel);

  xx=dmat(x, nn, nvar);

  if(verbose){
    ofsDebug<<"beta"<<endl;
    printDVector(&ofsDebug, beta, nbeta);
    
//     ofsDebug<<"x"<<endl;
//     printDVector(&ofsDebug, x, nn*nvar);
    ofsDebug<<"xx"<<endl;
    printDMatrix(&ofsDebug, xx, nn, nvar);
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
  predictor(xx, nvar, beta, *cure, pred);
  
  if(verbose){
    ofsDebug<<"pred"<<endl;
    printDM(&ofsDebug, pred);
    ofsDebug<<"s0"<<endl;
    printDVector(&ofsDebug, s0, nt);
  }

  fitSurvival(status, dd, rr, pred, model, *cure, *tol, s0, nt);

  if(verbose){
    ofsDebug<<"s0 ";
    printDVector(&ofsDebug, s0, nt);
  }
  
  *plik=likelihood(status, dd, rr, model, *cure, s0, pred, nt);
    
  if(verbose)
    ofsDebug<<"plik: "<<*plik<<endl;
  
  //  isinfLik=isinf(*plik);
  *plik=(isinf(*plik) || isnan(*plik) ? -VERYBIG*(1+drand48()*0.1) : *plik);

  if(verbose)
    ofsDebug<<"plik: "<<*plik<<endl;
}
}

#if 0
// use the Nelson-Aalen etimator (Klein, Moeschberger, page 86) to
// initialize the hazard jumps
// Note: exp(-cumsum(hh)) is the KM in [SP]
void initSurvival(int *dd, int *rr, int *nt, double *s0)
{
  int j;
  int sum=0;

  for(j=*nt-1; j>=0; j--){
    sum+=rr[j];
    s0[j]=exp(-dd[j]/(double)sum);
  }
}
#endif

extern "C" {
void openDebug(char **survModel)
{
  string fileDebug;
    
  fileDebug=*survModel;
  fileDebug="debugLikelihood"+fileDebug;
  ofsDebug.open(fileDebug.c_str());
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
