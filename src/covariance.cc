#include <string>
#include <vector>
#include <math.h>
#include <R.h>
#include "Rostream.h"

using namespace std;

#define TINY 1e-10
#define SPECIALN 0
#define SPECIALY 1
#define VF 999 // vectorial function

double vtheta(vector<double> &pred, double s, int cc, int model);
void vtheta_pred(vector<double> &pred, double x, int cc, int model,
		 vector<double> &der1);
void vtheta_2pred(vector<double> &pred, double s, int cc, int model,
		  vector<double> &der1);
double vthetaCure(vector<double> &pred, double s, int cc, int model);
void vthetaCure_pred(vector<double> &pred, double s, int cc, int model,
		     vector<double> &der1);
void vthetaCure_2pred(vector<double> &pred, double s, int cc, int model,
		      vector<double> &der2);
void Theton_pred(vector<double> &pred, double s, int cc, int model, 
		 vector<double> &der1);
double Theton_h(vector<double> &pred, double s, int cc, int model);
void ThetonCure_pred(vector<double> &pred, double s, int cc, int model,
		     vector<double> &der1);
double ThetonCure_h(vector<double> &pred, double s, int cc, int model);
void printDV(vector<double> &a); 
void printDM(vector<vector<double> > &a); 
void printDVector(double *a, int n);
void printDMatrix(double **a, int nrow, int ncol);
void printIVector(int *a, int n);
void printDMRformat(vector<vector<double> > &a);
void printDMatrixRformat(double **a, int nrow, int ncol);
void printModelFunction(string fName, vector<double> &pred, 
			double s, int cc, double resU, vector<double> resVF);
void predictor(double **xx1, double **xx2, int nvar1, int nvar2, double *beta, 
	       int cure, vector<vector<double> > &pred);
void fitSurvival(int *status, int *dd, int *rr, vector<vector<double> > &pred, 
		 int model, int cure, double tol, double *s0, int nt);
double **dmat(double *array, int nrow, int ncol);
int nmodel(string model);


void der1vthetabeta(double *x1i, double *x2i, int nvar1, int nvar2, 
		    vector<double> &predi, int statusi, double ss, int model, 
		    int cure, int special, vector<double> &der1, int verbose)
{  
  int k, npred;
  vector<double> d1;

  npred=predi.size();
  d1.resize(npred);

  switch(special){
  case SPECIALN:
    vtheta_pred(predi, ss, statusi, model, d1);
    if(verbose)
      printModelFunction("vtheta_pred", predi, ss, statusi, VF, d1);
    break;
  case SPECIALY:
    vthetaCure_pred(predi, ss, statusi, model, d1);
    if(verbose)
      printModelFunction("vthetaCure_pred", predi, ss, statusi, VF, 
			 d1);
    break;
  default:
    Rcerr<<"der1vthetabeta: incorrect special value "<<special<<endl;
  }

  for(k=0; k<nvar1; k++){
    der1[k]=d1[0]*x1i[k]*predi[0];
  }
  for(k=0; k<nvar2; k++){
    der1[k+nvar1]=d1[1]*x2i[k]*predi[1];
  }

  if(cure){
    der1[nvar1+nvar2]=d1[0]*predi[0];
  }
}


// Second partial derivative of the log likelihood with respect to beta
// nbeta: nvar1+nvar2+cure
// pred: nn x npred
// der2: nbeta x nbeta
// The vector of second derivatives of vtheta with respect to theta
// has to be organized like this (vtheta_00, vtheta_11, vtheta_01)
void der2likBeta(double **xx1, double **xx2, vector<vector<double> > &pred, 
		 int *rr, int *status, vector<double> &ss, int model, int cure,
		 int nvar1, int nvar2, double **der2, int verbose)
{  
  int i, j, k, l1, l2, npred, nbeta, nt;
  double vt, vt2, aux1, aux2, aux3;
  vector<double> d1, d2;

  nt=ss.size();
  npred=pred[0].size();
  nbeta=nvar1+nvar2+cure;
  d1.resize(nbeta);
  d2.resize(npred*(npred+1)/2);

  for(i=0; i<nbeta; i++)
    for(j=0; j<nbeta; j++)
      der2[i][j]=0;

  i=0;
  for(k=0; k<nt-cure; k++){
    for(j=0; j<rr[k]; j++){
      vt=vtheta(pred[i], ss[k], status[i], model);
      if(verbose){
	Rcout<<i<<" ";
	printModelFunction("vtheta", pred[i], ss[k], status[i], 
			   vt, d1);
      }
      vt=(fabs(vt)<TINY ? (vt<0 ? -TINY : TINY) : vt);
      vt2=vt*vt;
      if(nvar2>0)
	der1vthetabeta(xx1[i], xx2[i], nvar1, nvar2, pred[i], status[i], ss[k],
		       model, cure, SPECIALN, d1, verbose);
      else
	der1vthetabeta(xx1[i], 0, nvar1, nvar2, pred[i], status[i], ss[k], 
		       model, cure, SPECIALN, d1, verbose);
      vtheta_2pred(pred[i], ss[k], status[i], model, d2);
      if(verbose)
	printModelFunction("vtheta_2pred", pred[i], ss[k], 
			   status[i], VF, d2);

      // derivative with respect to beta_k1, beta_k2 both contributing
      // to predictor theta, if cure model then also derivatives with
      // respect to cure parameter
      aux1=d2[0]*pred[i][0]*pred[i][0];
      for(l1=0; l1<nvar1; l1++){
	aux2=d1[l1]/vt2;
	aux3=(aux1*xx1[i][l1]+d1[l1])/vt;
	for(l2=l1; l2<nvar1; l2++)
	  der2[l1][l2]+=-aux2*d1[l2]+aux3*xx1[i][l2];
	if(cure)
	  der2[l1][nbeta-1]+=-aux2*d1[nbeta-1]+aux3;
      }

      if(cure)
	der2[nbeta-1][nbeta-1]+=-d1[nbeta-1]/vt2*d1[nbeta-1]+
	  (aux1+d1[nbeta-1])/vt;

      // derivative with respect to beta_k1, beta_k2 both contributing
      // to predictor eta
      if(npred>1){
	aux1=d2[1]*pred[i][1]*pred[i][1];
	for(l1=0; l1<nvar2; l1++){
	  aux2=d1[l1+nvar1]/vt2;
	  aux3=(aux1*xx2[i][l1]+d1[l1+nvar1])/vt;
	  for(l2=l1; l2<nvar2; l2++)
	    der2[l1+nvar1][l2+nvar1]+=-aux2*d1[l2+nvar1]+aux3*xx2[i][l2];
	}

	// off-diagonal submatrix, i.e. derivative with respect to beta_k1,
	// beta_k2 corresponding to different predictors
	aux1=d2[2]*pred[i][0]*pred[i][1]/vt;
	for(l1=0; l1<nvar1; l1++){
	  aux2=d1[l1]/vt2;
	  aux3=aux1*xx1[i][l1];
	  for(l2=0; l2<nvar2; l2++)
	    der2[l1][l2+nvar1]+=-aux2*d1[l2+nvar1]+aux3*xx2[i][l2];
	}
	if(cure)
	  for(l2=0; l2<nvar2; l2++)
	    der2[l2+nvar1][nbeta-1]+=-d1[nbeta-1]/vt2*d1[l2+nvar1]+
	      aux1*xx2[i][l2];
      }
      i++;
    }
  }

//   Rcout<<"##################### cure terms #####################"<<endl;
  if(cure){
    for(j=0; j<rr[nt-1]; j++){
      vt=vthetaCure(pred[i], ss[nt-2], status[i], model);
      if(verbose)
	printModelFunction("vthetaCure", pred[i], ss[nt-2], 
			   status[i], vt, d1);
      vt=(fabs(vt)<TINY ? (vt<0 ? -TINY : TINY) : vt);
      vt2=vt*vt;
      if(nvar2>0)
	der1vthetabeta(xx1[i], xx2[i], nvar1, nvar2, pred[i], status[i], 
		       ss[nt-2], model, cure, SPECIALY, d1, verbose);
      else
	der1vthetabeta(xx1[i], 0, nvar1, nvar2, pred[i], status[i], ss[nt-2], 
		       model, cure, SPECIALY, d1, verbose);
      vthetaCure_2pred(pred[i], ss[nt-2], status[i], model, d2);
      if(verbose)
	printModelFunction("vthetaCure_2pred", pred[i], ss[nt-2], 
			   status[i], VF, d2);

      // derivative with respect to beta_k1, beta_k2 both contributing
      // to predictor theta, includes derivatives with respect to cure
      // parameter
      aux1=d2[0]*pred[i][0]*pred[i][0];
      for(l1=0; l1<nvar1; l1++){
	aux2=d1[l1]/vt2;
	aux3=(aux1*xx1[i][l1]+d1[l1])/vt;
	for(l2=l1; l2<nvar1; l2++)
	  der2[l1][l2]+=-aux2*d1[l2]+aux3*xx1[i][l2];
	der2[l1][nbeta-1]+=-aux2*d1[nbeta-1]+aux3;
      }
      der2[nbeta-1][nbeta-1]+=-d1[nbeta-1]/vt2*d1[nbeta-1]+
	(aux1+d1[nbeta-1])/vt;
           
      // derivative with respect to beta_k1, beta_k2 both contributing
      // to predictor eta
      if(npred>1){
	aux1=d2[1]*pred[i][1]*pred[i][1];
	for(l1=0; l1<nvar2; l1++){
	  aux2=d1[l1+nvar1]/vt2;
	  aux3=(aux1*xx2[i][l1]+d1[l1+nvar1])/vt;
	  for(l2=l1; l2<nvar2; l2++)
	    der2[l1+nvar1][l2+nvar1]+=-aux2*d1[l2+nvar1]+aux3*xx2[i][l2];
	}

	// off-diagonal submatrix, i.e. derivative with respect to beta_k1,
	// beta_k2 corresponding to different predictors
	aux1=d2[2]*pred[i][0]*pred[i][1]/vt;
	for(l1=0; l1<nvar1; l1++){
	  aux2=d1[l1]/vt2;
	  aux3=aux1*xx1[i][l1];
	  for(l2=0; l2<nvar2; l2++)
	    der2[l1][l2+nvar1]+=-aux2*d1[l2+nvar1]+aux3*xx2[i][l2];
	}
	for(l2=0; l2<nvar2; l2++)
	  der2[l2+nvar1][nbeta-1]+=-d1[nbeta-1]/vt2*d1[l2+nvar1]+
	    aux1*xx2[i][l2];
      }
      i++;
    }
  }

  for(l1=0; l1<nbeta; l1++)
    for(l2=l1+1; l2<nbeta; l2++)
      der2[l2][l1]=der2[l1][l2];
}


// Second partial derivative of the log likelihood with respect to
// beta and the hazard jumps h
// nbeta: nvar1+nvar2+cure
// pred: nn x npred
// der2: nh x nbeta
// nh: nt-cure
// Note: In cure models h_nt doesn't show up
void der2likBetah(double **xx1, double **xx2, int nvar1, int nvar2, 
		  vector<vector<double> > &pred, int *rr, int *status, 
		  vector<double> &ss, int model, int cure,
		  vector<vector<double> > &der2, int verbose)
{
  int i, j, k, l, nt, nn, npred, nbeta;
  vector<double> d1;

  nt=ss.size();
  nn=pred.size();
  npred=pred[0].size();
  nbeta=der2[0].size();

  d1.resize(npred);

  for(j=0; j<nbeta; j++)
    der2[nt-cure-1][j]=0;

  i=nn-1;
  if(cure){
    for(j=0; j<rr[nt-1]; j++){
      ThetonCure_pred(pred[i], ss[nt-2], status[i], model, d1);
      if(verbose)
	printModelFunction("ThetonCure_pred", pred[i], ss[nt-2], 
			   status[i], VF, d1);
      for(l=0; l<nvar1; l++)
	  der2[nt-2][l]-=d1[0]*pred[i][0]*xx1[i][l];
      for(l=0; l<nvar2; l++)
	der2[nt-2][l+nvar1]-=d1[1]*pred[i][1]*xx2[i][l];
      der2[nt-2][nbeta-1]-=d1[0]*pred[i][0];
      i--;
    }     
  }else{
    for(j=0; j<rr[nt-1]; j++){
      Theton_pred(pred[i], ss[nt-1], status[i], model, d1);
      if(verbose)
	printModelFunction("Theton_pred", pred[i], ss[nt-1], 
			   status[i], VF, d1);
      for(l=0; l<nvar1; l++)
	der2[nt-1][l]-=d1[0]*pred[i][0]*xx1[i][l];
      for(l=0; l<nvar2; l++)
	der2[nt-1][l+nvar1]-=d1[1]*pred[i][1]*xx2[i][l];
      i--;
    }
  }

  for(k=nt-2; k>=0; k--){
    // The case k=nt-1 has to be done separately because of the loop
    // below and the cure model
    if(k<nt-2 || !cure)
      for(j=0; j<nbeta; j++)
	der2[k][j]=der2[k+1][j];
    for(j=0; j<rr[k]; j++){
      Theton_pred(pred[i], ss[k], status[i], model, d1);
      if(verbose)
	printModelFunction("Theton_pred", pred[i], ss[k], status[i],
			   VF, d1);
      for(l=0; l<nvar1; l++)
	der2[k][l]-=d1[0]*pred[i][0]*xx1[i][l];
      for(l=0; l<nvar2; l++)
	der2[k][l+nvar1]-=d1[1]*pred[i][1]*xx2[i][l];
      if(cure)
	der2[k][nbeta-1]-=d1[0]*pred[i][0];
      i--;
    }    
  }
}


// Terms used to calculate the second derivative of the log-likelihood
// with respect to h
// a_k=sum_{j:t_j>=t_k} der_h Theton(S_j)
// diag_k=D_k/h_k
void der1ThetonhDiag(vector<vector<double> > &pred, int *rr, int *dd, 
		     int *status, double *s0, vector<double> &ss, int model, 
		     int cure, vector<double> &aa, vector<double> &diag, 
		     int verbose)
{
  int i, j, k, nn, nt;
  double hh, aux;

  nt=ss.size();
  nn=pred.size();
  aa[nt-cure-1]=0;

  i=nn-1;
  if(cure){
    for(j=0; j<rr[nt-1]; j++){
      aux=ThetonCure_h(pred[i], ss[nt-2], status[i], model);
      aa[nt-2]+=aux;
      if (verbose)
	printModelFunction("ThetonCure_h", pred[i], ss[nt-2], 
			   status[i], aux, aa); 
      // aa is passed just to put something
      i--;
    }     
  }else{
    for(j=0; j<rr[nt-1]; j++){
      aux=Theton_h(pred[i], ss[nt-1], status[i], model);
      aa[nt-1]+=aux;
      if (verbose)
	printModelFunction("Theton_h", pred[i], ss[nt-1], status[i],
			   aux, aa);
      i--;
    }     
  }

  for(k=nt-2; k>=0; k--){
    // The case k=nt-1 has to be done separately because of the line
    // below
    if(k<nt-2 || !cure)
      aa[k]=aa[k+1];
    for(j=0; j<rr[k]; j++){
      aux=Theton_h(pred[i], ss[k], status[i], model);
      aa[k]+=aux;
      if (verbose)
	printModelFunction("Theton_h", pred[i], ss[k], status[i],
			   aux, aa);
      i--;
    }    
  }

  for(k=0; k<nt-cure; k++){
    hh=-log(s0[k]);
    diag[k]=(hh<TINY ? double(dd[k])/TINY : double(dd[k])/(hh*hh)); 
  }
}


// Solve the linear system in order to find the derivatives of hazard
// jumps h(beta) with respect to beta.
// The notation refers to notes A1 and N1
// The linear system is of the form (D+R)x=bb where D is a diagonal
// matrix with elements diag, R_kj=-sum{i>=max{k,j}} aa_i 
// fiVec and solveLinearSystem are involved

// nh=nt-cure
// Calculates x_nh,...,x_1 given x_*=sum_{k=1}^nh x_i; i.e
// diag_k x_k=
// bb_k+sum{i=k}^nn aa_i x_* -sum_{j=k+1}^nh sum_{i=k}^{j-1} aa_i x_j
void fiVec(double xSum, vector<double> &aa, vector<double> &bb, 
	   vector<double> &diag, vector<double> &x)
{
  int j, k, nh=x.size();
  double sum=0, aux;

  for(k=nh-1; k>=0; k--){
    sum=sum+aa[k];
    x[k]=bb[k]+sum*xSum;
    aux=aa[k];
    for(j=k+1; j<nh; j++){
      x[k]=x[k]-aux*x[j];
      aux=aux+aa[j];
    }
    x[k]=(fabs(diag[k])<TINY ? x[k]/TINY : x[k]/diag[k]);
  }
}


// Solves a linear system of equations (D+R)x=bb, where D and R are as
// described above.
// Let fi_*(y)=sum(fiVec(y)), y in R, fi_* is a linear function of y.
// Steps:
// Solve fi_*(y)=y. To do this let g(y)=fi_*(y)-y, since fi_* is linear
// g(y)=ay+b. Find a and b by calculating g(0) and g(1). Now, the solution
// of fi_*(y)=y is y=-b/a
// Find x=fiVec(-b/a) 
void solveLinearSystem(vector<double> &aa, vector<double> &bb, 
		       vector<double> &diag, vector<double> &x)
{
  int j;
  double fi0=0, fi1=0;

  fiVec(0.0, aa, bb, diag, x);
  for(j=0; j<int(x.size()); j++)
    fi0+=x[j];

  fiVec(1.0, aa, bb, diag, x);
  for(j=0; j<int(x.size()); j++)
    fi1+=x[j];
  
  if(fabs(fi0+1-fi1)>TINY)
    fiVec(fi0/(fi0+1-fi1), aa, bb, diag, x);
  else
    Rcerr<<"solveLinearSystem: fi0+1-fi1=0"<<endl; 

// This shouldn't happen, if it happens, the solution would be to
// evaluate fiVec at other points, instead of 0 and 1
}


// Finds the derivative of h(beta) with respect to beta
// The second derivative of likelihood with respect to h and beta,
// is equal to -derivative of Theton with respect to beta, d1Tb
// d1Tb: nt-cure x nbeta
// diag, aa: nt-cure
// der1: nt-cure x nbeta
void der1Hbeta(vector<double> &diag, vector<double> &aa, 
	       vector<vector<double> > &d2likbh, vector<vector<double> > &der1)
{
  int j, k, nh, nbeta;
  vector<double> bb, x, rr;

  nh=diag.size();
  nbeta=d2likbh[0].size();

  // columns of first derivative of Theton with respect to beta
  bb.resize(nh);
  x.resize(nh);
  rr.resize(nh);
  
  rr[nh-1]=-aa[nh-1];
  for(k=nh-2; k>=0; k--)
    rr[k]=-aa[k]+aa[k+1];
  
  for(j=0; j<nbeta; j++){
    for(k=0; k<nh; k++)
      bb[k]=d2likbh[k][j];
    solveLinearSystem(rr, bb, diag, x);
    for(k=0; k<nh; k++)
      der1[k][j]=x[k];
  }
}


// Second partial derivative of likelihood times derivative of dHb with respect
// to beta plus the transpose of this product
// res: nbeta x nbeta
void term23(vector<vector<double> > &d1hb, vector<vector<double> > &d2likbh,
	    vector<vector<double> > &res)
{
  int i, j, k, nh, nbeta;

  nh=d1hb.size();
  nbeta=d1hb[0].size();
  
  for(i=0; i<nbeta; i++)
    for(j=0; j<nbeta; j++){
      res[i][j]=0;
      for(k=0; k<nh; k++)
	res[i][j]+=d2likbh[k][i]*d1hb[k][j];
    }
  
  for(i=0; i<nbeta; i++)
    for(j=0; j<=i; j++){
      res[i][j]+=res[j][i];
      res[j][i]=res[i][j];
    }
}


// First derivative of h(beta) with respect to beta times second
// derivative of likelihood with respet to h times first derivative of
// h(beta) with respect to beta
void term4(vector<double> &aa, vector<double> &diag, 
	   vector<vector<double> > &d1hb, vector<vector<double> > &res)
{
  int i, j, k, nh, nbeta;
  double sum1, sum2;

  nh=d1hb.size();
  nbeta=d1hb[0].size();

  for(i=0; i<nbeta; i++)
    for(j=0; j<nbeta; j++){
      res[i][j]=0;
      sum1=0.0;
      sum2=0.0;
      for(k=1; k<nh; k++)
	sum2+=d1hb[k][i]*aa[k];
      for(k=0; k<nh-1; k++){
	sum1+=d1hb[k][i];
	res[i][j]-=(sum1*aa[k]+d1hb[k][i]*diag[k]+sum2)*d1hb[k][j];
	sum2-=d1hb[k+1][i]*aa[k+1];
      }
      // last term has to be calculated separately because of sum2
      sum1+=d1hb[nh-1][i];
      res[i][j]-=(sum1*aa[nh-1]+d1hb[nh-1][i]*diag[nh-1]+sum2)*d1hb[nh-1][j];
    }
}


void der2likh(vector<double> &aa, vector<double> &diag, 
	      vector<vector<double > > &d2lh) 
{
  int i, j, nh;

  nh=aa.size();

  for(i=0; i<nh; i++)
    d2lh[i][i]=-aa[i]-diag[i];

  for(i=1; i<nh; i++)
    for(j=0; j<i; j++){
      d2lh[i][j]=-aa[i];
      d2lh[j][i]=d2lh[i][j];
    }
}

double checkIs0(vector<vector<double> > &d2lh, vector<vector<double> > &d1hb, 
		vector<vector<double> > &d2likbh)
{
  int i, j, k;
  double max=-1, sum;

  for(i=0; i<int(d2lh.size()); i++)
    for(j=0; j<int(d1hb[0].size()); j++){
      sum=0;
      for(k=0; k<int(d2lh[0].size()); k++)
	sum+=d2lh[i][k]*d1hb[k][j];
      sum+=d2likbh[i][j];
      sum=fabs(sum);
      max=(max>sum ? max : sum);
    }
  return(max);
}

extern "C"{
void informationMatrix(double *beta, double *x1, double *x2, int *status, 
		       int *dd, int *rr, double *s0, char **survModel, 
		       int *cure, int *nvar1, int *nvar2, int *ntime, 
		       int *nobs, int *npred, int *verbose, double *imat)
{
  int i, j, nt, nn, nbeta, nh, model;
  double **xx1, **xx2, **infMat;
  vector<double> ss, aa, diag;
  vector<vector<double> > auxMat, d2likbh, d1hb, d2lh;
 
  nt=*ntime;
  nh=nt-*cure;
  nn=*nobs;
  model=nmodel(*survModel);
  nbeta=(*nvar1)+(*nvar2)+(*cure);

  if(*verbose){
    Rcout<<"information matrix"<<endl;
    Rcout<<"nn: "<<nn<<" nvar1: "<<*nvar1<<" nvar2: "<<*nvar2<<endl;
    Rcout<<"beta "<<nbeta<<endl;
    printDVector(beta, nbeta);
  }
  xx1=dmat(x1, nn, *nvar1);
  xx2=dmat(x2, nn, *nvar2);
  infMat=dmat(imat, nbeta, nbeta);

  vector<vector<double> > pred(*nobs, std::vector<double>(*npred, 0.0));
  predictor(xx1, xx2, *nvar1, *nvar2, beta, *cure, pred);

  ss.resize(nt);
  ss[0]=s0[0];
  for(j=1; j<nt; j++)
    ss[j]=ss[j-1]*s0[j];
  if(*verbose){
    Rcout<<"pred"<<endl;
    printDM(pred);
    Rcout<<"s0"<<endl;
    printDVector(s0, nt);
  }

   // Partial second derivative of likelihood with respect to beta
  der2likBeta(xx1, xx2, pred, rr, status, ss, model, *cure, *nvar1, *nvar2, 
	      infMat, *verbose);
  if(*verbose){
    Rcout<<"der2likbeta start"<<endl;
    Rcout<<"d2lbeta <- ";
    printDMatrixRformat(infMat, nbeta, nbeta);
  }

  // Partial second derivative of likelihood with respect to beta and h
  d2likbh.resize(nh);
  for(i=0; i<nh; i++)
    d2likbh[i].resize(nbeta);
  der2likBetah(xx1, xx2, *nvar1, *nvar2, pred, rr, status, ss, model, *cure, 
	       d2likbh, *verbose);
  if(*verbose){
    Rcout<<"d2lbetah <- ";
    printDMRformat(d2likbh);
  }

  // Elements for construction of partial second derivative of
  // likelihood with respect to h, it is also used for calculation of
  // derivative of h(beta) with respet to beta
  aa.resize(nh);
  diag.resize(nh);
  der1ThetonhDiag(pred, rr, dd, status, s0, ss, model, *cure, aa, diag, 
		  *verbose);
  
  if(*verbose){
    Rcout<<"aa"<<endl;
    printDV(aa);
    Rcout<<"diag"<<endl;
    printDV(diag);
  }

  // First derivative of h(beta) with respect to beta
  d1hb.resize(nh);
  for(i=0; i<nh; i++)
    d1hb[i].resize(nbeta);
  der1Hbeta(diag, aa, d2likbh, d1hb);
  if(*verbose){
    Rcout<<"der1Hbeta <- ";
    printDMRformat(d1hb);
  }

  // Second partial derivative of likelihood times derivative of dHb
  // with respect to beta plus the transpose of this product
  auxMat.resize(nbeta);
  for(i=0; i<nbeta; i++)
    auxMat[i].resize(nbeta);
  term23(d1hb, d2likbh, auxMat);
  if(*verbose){
    Rcout<<"term23 <- ";
   printDMRformat(auxMat);
  }

  for(i=0; i<nbeta; i++)
    for(j=0; j<nbeta; j++)
      infMat[i][j]+=auxMat[i][j];

  // First derivative of h(beta) with respect to beta times second
  // derivative of likelihood with respet to h times first derivative
  // of h(beta) with respect to beta
  term4(aa, diag, d1hb, auxMat);
  if(*verbose){
    Rcout<<"term4 <- ";
    printDMRformat(auxMat);
  }
  for(i=0; i<nbeta; i++)
    for(j=0; j<nbeta; j++){
      infMat[i][j]+=auxMat[i][j];
      infMat[i][j]*=(-1);
    }

  if(*verbose){
    Rcout<<"infMat <- ";
    printDMatrixRformat(infMat, nbeta, nbeta);
  }
}
} // extern "C"



