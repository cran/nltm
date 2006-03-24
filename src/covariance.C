#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <R.h>

using namespace std;

#define TINY 1e-10
#define SPECIALN 0
#define SPECIALY 1

extern ofstream ofsDebug;

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
void printDV(ofstream *ofs, vector<double> &a); 
void printDM(ofstream *ofs, vector<vector<double> > &a); 
void printDVector(ofstream *ofs, double *a, int n);
void printDMatrix(ofstream *ofs, double **a, int nrow, int ncol);
void printIVector(ofstream *ofs, int *a, int n);
void printDMRformat(ofstream *ofs, vector<vector<double> > &a);
void printDMatrixRformat(ofstream *ofs, double **a, int nrow, int ncol);
void predictor(double **xx, int nvar, double *beta, int cure, 
 	       vector<vector<double> > &pred);
void fitSurvival(int *status, int *dd, int *rr, vector<vector<double> > &pred, 
		 int model, int cure, double tol, double *s0, int nt);
double **dmat(double *array, int nrow, int ncol);
int nmodel(string model);


void der1vthetabeta(double *xi, vector<double> &predi, int statusi, 
		    double ss, int model, int cure, int special, 
		    vector<double> &der1)
{  
  int j, k, npred, nvar;
  vector<double> d1;

  npred=predi.size();
  nvar=(der1.size()-cure)/npred;
  d1.resize(npred);

  //  ofsDebug<<special<<" "<<ss<<" "<<statusi<<" "<<nvar<<" "<<npred<<" ";
  switch(special){
  case SPECIALN:
    vtheta_pred(predi, ss, statusi, model, d1);
    break;
  case SPECIALY:
    vthetaCure_pred(predi, ss, statusi, model, d1);
    break;
  default:
    cerr<<"der1vthetabeta: incorrect special value "<<special<<endl;
  }

  for(j=0; j<npred; j++){
    //    ofsDebug<<d1[j]<<" "<<predi[j]<<" ";
    for(k=0; k<nvar; k++){
      der1[k+j*nvar]=d1[j]*xi[k]*predi[j];
      //      ofsDebug<<xi[k]<<" "<<der1[k+j*nvar]<<" ";
    }
  }
  if(cure){
    der1[nvar*npred]=d1[0]*predi[0];
    //    ofsDebug<<"cure "<<der1[nvar*npred]<<" ";
  }

  //  ofsDebug<<endl;
}


// Second partial derivative of the log likelihood with respect to beta
// nbeta: nvar*npred
// pred: nn x npred
// der2: nbeta x nbeta
// The vector of second derivatives of vtheta with respect to theta
// has to be organized like this (vtheta_00, vtheta_11, vtheta_01)
void der2likBeta(double **xx, vector<vector<double> > &pred, int *rr, 
		 int *status, vector<double> &ss, int model, int cure, 
		 int nvar, double **der2)
{  
  int i, j, k, l1, l2, npred, nn, nbeta, nt;
  double vt, vt2, aux1, aux2, aux3;
  vector<double> d1, d2;

  nt=ss.size();
  nn=pred.size();
  npred=pred[0].size();
  nbeta=nvar*npred+cure;
  d1.resize(nbeta);
  d2.resize(npred*(npred+1)/2);

  for(i=0; i<nbeta; i++)
    for(j=0; j<nbeta; j++)
      der2[i][j]=0;

//   ofsDebug<<"der2likBeta "<<model<<" "<<cure<<" "<<nvar<<" "<<nbeta<<" "
// 	  <<nt<<" "<<nn<<" "<<npred<<" "<<endl;
//  int i1;
  i=0;
  for(k=0; k<nt-cure; k++){
    for(j=0; j<rr[k]; j++){
      vt=vtheta(pred[i], ss[k], status[i], model);
//       for(i1=0; i1<npred; i1++)
// 	ofsDebug<<pred[i][i1]<<" ";
//       ofsDebug<<status[i]<<" "<<ss[k]<<" "<<vt<<endl;
      vt=(fabs(vt)<TINY ? (vt<0 ? -TINY : TINY) : vt);
      vt2=vt*vt;
//       ofsDebug<<"der1vthetabeta"<<endl;
      der1vthetabeta(xx[i], pred[i], status[i], ss[k], model, cure, SPECIALN, 
		     d1);
//       printDV(&ofsDebug, d1);
      vtheta_2pred(pred[i], ss[k], status[i], model, d2);
//       ofsDebug<<"vtheta_2pred "<<pred[i][0]<<" "<<ss[k]<<" "<<status[i]<<" "
// 	      <<endl;
//       printDV(&ofsDebug, d2);


      // derivative with respect to beta_k1, beta_k2 both contributing
      // to predictor theta, if cure model then also derivatives with
      // respect to cure parameter
//       ofsDebug<<"enter beta_theta, beta_theta loop "<<nvar<<endl;
      aux1=d2[0]*pred[i][0]*pred[i][0];
      for(l1=0; l1<nvar; l1++){
	aux2=d1[l1]/vt2;
	aux3=(aux1*xx[i][l1]+d1[l1])/vt;
	for(l2=l1; l2<nvar; l2++){
	  der2[l1][l2]+=-aux2*d1[l2]+aux3*xx[i][l2];
// 	  ofsDebug<<l1<<" "<<l2<<" "<<status[i]<<" "<<ss[k]<<" "<<xx[i][l1]
// 		  <<" "<<xx[i][l2]<<" "<<d1[l1]<<" "<<d1[l2]<<" "<<pred[i][0]
// 		  <<" "<<-aux2*d1[l2]+aux3*xx[i][l2]<<" "<<der2[l1][l2]<<endl;
	}
	if(cure){
	  der2[l1][nbeta-1]+=-aux2*d1[nbeta-1]+aux3;
// 	  ofsDebug<<l1<<" "<<nbeta-1<<" cure "<<" "<<xx[i][l1]<<" "<<d1[l1]
// 		  <<" "<<d1[nbeta-1]<<" "<<pred[i][0]<<" "
// 		  <<-aux2*d1[nbeta-1]+aux3<<" "<<der2[l1][nbeta-1]<<endl;
	  
	}
      }

      if(cure){
	der2[nbeta-1][nbeta-1]+=-d1[nbeta-1]/vt2*d1[nbeta-1]+
	  (aux1+d1[nbeta-1])/vt;
// 	ofsDebug<<nbeta-1<<" "<<nbeta-1<<" cure "<<d1[nbeta-1]<<" "<<pred[i][0]
// 		<<" "<<-d1[nbeta-1]/vt2*d1[nbeta-1]+(aux1+d1[nbeta-1])/vt<<" "
// 		<<der2[nbeta-1][nbeta-1]<<endl;
      }
//       ofsDebug<<"exit beta_theta, beta_theta loop "<<nvar<<endl;

      // derivative with respect to beta_k1, beta_k2 both contributing
      // to predictor eta
      if(npred>1){
// 	ofsDebug<<"enter beta_eta, beta_eta loop "<<nvar<<endl;
	aux1=d2[1]*pred[i][1]*pred[i][1];
	for(l1=0; l1<nvar; l1++){
	  aux2=d1[l1+nvar]/vt2;
	  aux3=(aux1*xx[i][l1]+d1[l1+nvar])/vt;
	  for(l2=l1; l2<nvar; l2++){
	    der2[l1+nvar][l2+nvar]+=-aux2*d1[l2+nvar]+aux3*xx[i][l2];
// 	    ofsDebug<<l1+nvar<<" "<<l2+nvar<<" "<<status[i]<<" "<<ss[k]
// 		    <<" "<<xx[i][l1]<<" "<<xx[i][l2]<<" "<<d1[l1+nvar]<<" "
// 		    <<d1[l2+nvar]<<" "<<pred[i][1]<<" "
// 		    <<-aux2*d1[l2+nvar]+aux3*xx[i][l2]<<" "
// 		    <<der2[l1+nvar][l2+nvar]<<endl;
	  }
	}
// 	ofsDebug<<"exit beta_eta, beta_eta loop "<<nvar<<endl;

	// off-diagonal submatrix, i.e. derivative with respect to beta_k1,
	// beta_k2 corresponding to different predictors
// 	ofsDebug<<"enter beta_theta, beta_eta loop "<<nvar<<endl;
	aux1=d2[2]*pred[i][0]*pred[i][1]/vt;
	for(l1=0; l1<nvar; l1++){
	  aux2=d1[l1]/vt2;
	  aux3=aux1*xx[i][l1];
	  for(l2=0; l2<nvar; l2++){
	    der2[l1][l2+nvar]+=-aux2*d1[l2+nvar]+aux3*xx[i][l2];
// 	    ofsDebug<<l1<<" "<<l2+nvar<<" "<<status[i]<<" "<<ss[k]<<" "
// 		    <<xx[i][l1]<<" "<<xx[i][l2]<<" "<<d1[l1]<<" "<<d1[l2+nvar]
// 		    <<" "<<pred[i][0]<<" "<<pred[i][1]<<" "
// 		    <<-aux2*d1[l2+nvar]+aux3*xx[i][l2]<<" "
// 		    <<der2[l1][l2+nvar]<<endl;
	  }
	  if(cure){
	    der2[l1+nvar][nbeta-1]+=-d1[nbeta-1]/vt2*d1[l1+nvar]+
	      aux1*xx[i][l1];
// 	    ofsDebug<<l1+nvar<<" "<<nbeta-1<<" cure "
// 		    <<xx[i][l1]<<" "<<d1[nbeta-1]<<" "<<d1[l1+nvar]<<" "
// 		    <<" "<<pred[i][0]<<" "<<pred[i][1]
// 		    <<-d1[nbeta-1]/vt2*d1[l1+nvar]+aux1*xx[i][l2]<<" "
// 		    <<der2[l1+nvar][nbeta-1]<<endl;
	  }
	}
// 	ofsDebug<<"exit beta_theta, beta_eta loop "<<nvar<<endl;
      }
      i++;
    }
  }

//   ofsDebug<<"##################### cure terms #####################"<<endl;
  if(cure){
    for(j=0; j<rr[nt-1]; j++){
      vt=vthetaCure(pred[i], ss[nt-2], status[i], model);
//       for(i1=0; i1<npred; i1++)
//  	ofsDebug<<pred[i][i1]<<" ";
//       ofsDebug<<status[i]<<" "<<ss[k]<<" "<<vt<<endl;
      vt=(fabs(vt)<TINY ? (vt<0 ? -TINY : TINY) : vt);
      vt2=vt*vt;
      der1vthetabeta(xx[i], pred[i], status[i], ss[nt-2], model, cure, 
		     SPECIALY, d1);
//       ofsDebug<<"der1vthetabeta"<<endl;
//       printDV(&ofsDebug, d1);
      vthetaCure_2pred(pred[i], ss[nt-2], status[i], model, d2);
//       ofsDebug<<"vtheta_2pred"<<endl;
//       printDV(&ofsDebug, d2);


      // derivative with respect to beta_k1, beta_k2 both contributing
      // to predictor theta, includes derivatives with respect to cure
      // parameter
//       ofsDebug<<"enter beta_theta, beta_theta loop "<<nvar<<endl;
      aux1=d2[0]*pred[i][0]*pred[i][0];
      for(l1=0; l1<nvar; l1++){
	aux2=d1[l1]/vt2;
	aux3=(aux1*xx[i][l1]+d1[l1])/vt;
	for(l2=l1; l2<nvar; l2++){
	  der2[l1][l2]+=-aux2*d1[l2]+aux3*xx[i][l2];
// 	  ofsDebug<<l1<<" "<<l2<<" "<<status[i]<<" "<<ss[k]<<" "<<xx[i][l1]
// 		  <<" "<<xx[i][l2]<<" "<<d1[l1]<<" "<<d1[l2]<<" "<<pred[i][0]
// 		  <<" "<<-aux2*d1[l2]+aux3*xx[i][l2]<<" "<<der2[l1][l2]<<endl;
	}
	der2[l1][nbeta-1]+=-aux2*d1[nbeta-1]+aux3;
// 	ofsDebug<<l1<<" "<<nbeta-1<<" cure "<<" "<<xx[i][l1]<<" "<<d1[l1]
// 		<<" "<<d1[nbeta-1]<<" "<<pred[i][0]<<" "
// 		<<-aux2*d1[nbeta-1]+aux3<<" "
// 		<<der2[l1][nbeta-1]<<endl;
      }
      der2[nbeta-1][nbeta-1]+=-d1[nbeta-1]/vt2*d1[nbeta-1]+
	(aux1+d1[nbeta-1])/vt;
//       ofsDebug<<nbeta-1<<" "<<nbeta-1<<" cure "<<d1[nbeta-1]<<" "<<pred[i][0]
// 	      <<" "<<-d1[nbeta-1]/vt2*d1[nbeta-1]+(aux1+d1[nbeta-1])/vt<<" "
// 	      <<der2[nbeta-1][nbeta-1]<<endl;
//       ofsDebug<<"exit beta_theta, beta_theta loop "<<nvar<<endl;
           
      // derivative with respect to beta_k1, beta_k2 both contributing
      // to predictor eta
      if(npred>1){
// 	ofsDebug<<"enter beta_eta, beta_eta loop "<<nvar<<endl;
	aux1=d2[1]*pred[i][1]*pred[i][1];
	for(l1=0; l1<nvar; l1++){
	  aux2=d1[l1+nvar]/vt2;
	  aux3=(aux1*xx[i][l1]+d1[l1+nvar])/vt;
	  for(l2=l1; l2<nvar; l2++){
	    der2[l1+nvar][l2+nvar]+=-aux2*d1[l2+nvar]+aux3*xx[i][l2];
// 	    ofsDebug<<l1+nvar<<" "<<l2+nvar<<" "<<status[i]<<" "<<ss[k]
// 		    <<" "<<xx[i][l1]<<" "<<xx[i][l2]<<" "<<d1[l1+nvar]<<" "
// 		    <<d1[l2+nvar]<<" "<<pred[i][1]<<" "
// 		    <<-aux2*d1[l2+nvar]+aux3*xx[i][l2]<<" "
// 		    <<der2[l1+nvar][l2+nvar]<<endl;
	  }
	}
// 	ofsDebug<<"exit beta_eta, beta_eta loop "<<nvar<<endl;

	// off-diagonal submatrix, i.e. derivative with respect to beta_k1,
	// beta_k2 corresponding to different predictors
	aux1=d2[2]*pred[i][0]*pred[i][1]/vt;
// 	ofsDebug<<"enter beta_theta, beta_eta loop "<<nvar<<endl;
	for(l1=0; l1<nvar; l1++){
	  aux2=d1[l1]/vt2;
	  aux3=aux1*xx[i][l1];
	  for(l2=0; l2<nvar; l2++){
	    der2[l1][l2+nvar]+=-aux2*d1[l2+nvar]+aux3*xx[i][l2];
// 	    ofsDebug<<l1<<" "<<l2+nvar<<" "<<status[i]<<" "<<ss[k]<<" "
// 		    <<xx[i][l1]<<" "<<xx[i][l2]<<" "<<d1[l1]<<" "<<d1[l2+nvar]
// 		    <<" "<<pred[i][0]<<" "<<pred[i][1]<<" "
// 		    <<-aux2*d1[l2+nvar]+aux3*xx[i][l2]<<" "
// 		    <<der2[l1][l2+nvar]<<endl;
	  }
	  der2[l1+nvar][nbeta-1]+=-d1[nbeta-1]/vt2*d1[l1+nvar]+aux1*xx[i][l1];
// 	  ofsDebug<<l1+nvar<<" "<<nbeta-1<<" cure "
// 		  <<xx[i][l1]<<" "<<d1[nbeta-1]<<" "<<d1[l1+nvar]<<" "
// 		  <<" "<<pred[i][0]<<" "<<pred[i][1]
// 		  <<-d1[nbeta-1]/vt2*d1[l1+nvar]+aux1*xx[i][l1]<<" "
// 		  <<der2[l1+nvar][nbeta-1]<<endl;
	}
// 	ofsDebug<<"exit beta_theta, beta_eta loop "<<nvar<<endl;
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
// nbeta: nvar*npred
// pred: nn x npred
// der2: nh x nbeta
// nh: nt-cure
// Note: In cure models h_nt doesn't show up
void der2likBetah(double **xx, vector<vector<double> > &pred, int *rr, 
		  int *status, vector<double> &ss, int model, int cure,
		  vector<vector<double> > &der2)
{
  int i, j, k, l, n1, nt, nn, npred, nbeta, nvar;
  vector<double> d1;

  nt=ss.size();
  nn=pred.size();
  npred=pred[0].size();
  nbeta=der2[0].size();
  nvar=(nbeta-cure)/npred;

  d1.resize(npred);

  for(j=0; j<nbeta; j++)
    der2[nt-cure-1][j]=0;

  //  ofsDebug<<"Theton_pred"<<endl;
  i=nn-1;
  if(cure){
    for(j=0; j<rr[nt-1]; j++){
      ThetonCure_pred(pred[i], ss[nt-2], status[i], model, d1);
      //    printDV(&ofsDebug, d1);
      for(n1=0; n1<npred; n1++)
	for(l=0; l<nvar; l++)
	  der2[nt-2][l+n1*nvar]-=d1[n1]*pred[i][n1]*xx[i][l];
      der2[nt-2][nbeta-1]-=d1[0]*pred[i][0];
      i--;
    }     
  }else{
    for(j=0; j<rr[nt-1]; j++){
      Theton_pred(pred[i], ss[nt-1], status[i], model, d1);
      //    printDV(&ofsDebug, d1);
      for(n1=0; n1<npred; n1++)
	for(l=0; l<nvar; l++)
	  der2[nt-1][l+n1*nvar]-=d1[n1]*pred[i][n1]*xx[i][l];
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
      //      printDV(&ofsDebug, d1);
      for(n1=0; n1<npred; n1++)
	for(l=0; l<nvar; l++)
	  der2[k][l+n1*nvar]-=d1[n1]*pred[i][n1]*xx[i][l];
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
void der1ThetonhDiag(double **xx, vector<vector<double> > &pred, int *rr, 
		     int *dd, int *status, double *s0, vector<double> &ss, 
		     int model, int cure, vector<double> &aa, 
		     vector<double> &diag)
{
  int i, j, k, nn, nt;
  double hh;

  nt=ss.size();
  nn=pred.size();
  aa[nt-cure-1]=0;

  //  ofsDebug<<"enter der1ThetonhDiag "<<nt<<" "<<nn<<endl;

  i=nn-1;
  if(cure){
    for(j=0; j<rr[nt-1]; j++){
      aa[nt-2]+=ThetonCure_h(pred[i], ss[nt-2], status[i], model);
      i--;
    }     
  }else{
    for(j=0; j<rr[nt-1]; j++){
      aa[nt-1]+=Theton_h(pred[i], ss[nt-1], status[i], model);
      i--;
    }     
  }

  for(k=nt-2; k>=0; k--){
    // The case k=nt-1 has to be done separately because of the line
    // below
    if(k<nt-2 || !cure)
      aa[k]=aa[k+1];
    for(j=0; j<rr[k]; j++){
      aa[k]+=Theton_h(pred[i], ss[k], status[i], model);
      i--;
    }    
  }

  for(k=0; k<nt-cure; k++){
    hh=-log(s0[k]);
    diag[k]=(hh<TINY ? double(dd[k])/TINY : double(dd[k])/(hh*hh)); 
  }
#if 0
  ofsDebug<<"der1ThetonhDiag ss"<<endl;
  printDV(&ofsDebug, ss);

  ofsDebug<<"der1ThetonhDiag s0"<<endl;
  printDVector(&ofsDebug, s0, nt);

  ofsDebug<<"der1ThetonhDiag dd"<<endl;
  printIVector(&ofsDebug, dd, nt);

  ofsDebug<<"der1ThetonhDiag diag"<<endl;
  printDV(&ofsDebug, diag);

  ofsDebug<<"exit der1ThetonhDiag "<<endl;
#endif
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
    cerr<<"solveLinearSystem: fi0+1-fi1=0"<<endl; 

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

  //  ofsDebug<<"term4"<<endl;
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
	//	ofsDebug<<k<<" "<<sum1<<" "<<sum2<<" "<<res[i][j]<<endl;
	sum2-=d1hb[k+1][i]*aa[k+1];
      }
      // last term has to be calculated separately because of sum2
      sum1+=d1hb[nh-1][i];
      res[i][j]-=(sum1*aa[nh-1]+d1hb[nh-1][i]*diag[nh-1]+sum2)*d1hb[nh-1][j];
      //      ofsDebug<<nh-1<<" "<<sum1<<" "<<sum2<<" "<<res[i][j]<<endl;
    }
}


void der2likh(vector<double> &aa, vector<double> &diag, 
	      vector<vector<double > > &d2lh) 
{
  int i, j, nh;

  nh=aa.size();

  for(i=0; i<nh; i++)
    d2lh[i][i]=-aa[i]-diag[i];

  //  ofsDebug<<"der2likh"<<endl;
  for(i=1; i<nh; i++)
    for(j=0; j<i; j++){
      d2lh[i][j]=-aa[i];
      d2lh[j][i]=d2lh[i][j];
//       ofsDebug<<i<<" "<<j<<" "<<d2lh[i][j]<<" "<<d2lh[j][i]<<" "<<aa[i]<<endl;
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
void informationMatrix(double *beta, double *x, int *status, int *dd, 
		       int *rr, double *s0, char **survModel, int *cure,
		       int *nvar1, int *ntime, int *nobs, int *npred, 
		       int *verbose, double *imat)
{
  int i, j, nvar, nt, nn, nbeta, nh, model;
  double **xx, **infMat;
  vector<double> ss, aa, diag;
  vector<vector<double> > pred, auxMat, d2likbh, d1hb, d2lh;
 
  nvar=*nvar1;
  nt=*ntime;
  nh=nt-*cure;
  nn=*nobs;
  model=nmodel(*survModel);
  nbeta=nvar*(*npred)+(*cure);

#if 0
  fileDebug=*survModel;
  fileDebug="result/cov"+fileDebug+(beta[0]==0 ? "0" : "MLE");
  ofsDebug.open(fileDebug.c_str());
  if(!ofsDebug) {
    cerr<<"ERROR: couldn't create debug file."<<endl;
    exit(1);
  }
  ofsDebug<<"information matrix"<<endl;
  ofsDebug<<"beta "<<nvar<<" "<<nt<<" "<<nn<<" "<<*npred<<" "<<nbeta<<" "
	  <<*cure<<endl;
  printDVector(&ofsDebug, beta, nbeta);
#endif

  xx=dmat(x, nn, nvar);
  infMat=dmat(imat, nbeta, nbeta);

  pred.resize(nn);
  for(i=0; i<int(pred.size()); i++)
    pred[i].resize(*npred);
//   ofsDebug<<"pred"<<endl;
  predictor(xx, nvar, beta, *cure, pred);
//   printDM(&ofsDebug, pred);

  ss.resize(nt);
  ss[0]=s0[0];
  for(j=1; j<nt; j++)
    ss[j]=ss[j-1]*s0[j];
#if 0
  ofsDebug<<"s0"<<endl;
  printDVector(&ofsDebug, s0, nt);
  ofsDebug<<"ss"<<endl;
  printDV(&ofsDebug, ss);
#endif
   // Partial second derivative of likelihood with respect to beta
//   ofsDebug<<"der2likbeta start"<<endl;
  der2likBeta(xx, pred, rr, status, ss, model, *cure, nvar, infMat);
//   ofsDebug<<"d2lbeta <- ";
//   printDMatrixRformat(&ofsDebug, infMat, nbeta, nbeta);

  // Partial second derivative of likelihood with respect to beta and h
  d2likbh.resize(nh);
  for(i=0; i<nh; i++)
    d2likbh[i].resize(nbeta);
  der2likBetah(xx, pred, rr, status, ss, model, *cure, d2likbh);
//   ofsDebug<<"d2lbetah <- ";
//   printDMRformat(&ofsDebug, d2likbh);

  // Elements for construction of partial second derivative of
  // likelihood with respect to h, it is also used for calculation of
  // derivative of h(beta) with respet to beta
  aa.resize(nh);
  diag.resize(nh);
  der1ThetonhDiag(xx, pred, rr, dd, status, s0, ss, model, *cure, aa, diag);
//   ofsDebug<<"aa"<<endl;
//   printDV(&ofsDebug, aa);
//   ofsDebug<<"diag"<<endl;
//   printDV(&ofsDebug, diag);



  // First derivative of h(beta) with respect to beta
  d1hb.resize(nh);
  for(i=0; i<nh; i++)
    d1hb[i].resize(nbeta);
//   ofsDebug<<"enter der1Hbeta"<<endl;
  der1Hbeta(diag, aa, d2likbh, d1hb);
//   ofsDebug<<"der1Hbeta <- ";
//   printDMRformat(&ofsDebug, d1hb);

  /////////////////////////////////////
  // DELETE
#if 0
  d2lh.resize(nh);
  for(i=0; i<nh; i++)
    d2lh[i].resize(nh);
  der2likh(aa, diag, d2lh);
  ofsDebug<<"d2lh <- ";
  printDMRformat(&ofsDebug, d2lh);

  // Check d2lh*d1hb+=0
  ofsDebug<<"max(abs(d2lh*d1hb+d2likbh))"<<checkIs0(d2lh, d1hb, d2likbh)<<endl;
#endif
  /////////////////////////////////////


  // Second partial derivative of likelihood times derivative of dHb
  // with respect to beta plus the transpose of this product
  auxMat.resize(nbeta);
  for(i=0; i<nbeta; i++)
    auxMat[i].resize(nbeta);
  term23(d1hb, d2likbh, auxMat);
//   ofsDebug<<"term23 <- ";
//   printDMRformat(&ofsDebug, auxMat);

  for(i=0; i<nbeta; i++)
    for(j=0; j<nbeta; j++)
      infMat[i][j]+=auxMat[i][j];

  // First derivative of h(beta) with respect to beta times second
  // derivative of likelihood with respet to h times first derivative
  // of h(beta) with respect to beta
  term4(aa, diag, d1hb, auxMat);
//   ofsDebug<<"term4 <- ";
//   printDMRformat(&ofsDebug, auxMat);

  for(i=0; i<nbeta; i++)
    for(j=0; j<nbeta; j++){
      infMat[i][j]+=auxMat[i][j];
      infMat[i][j]*=(-1);
    }

//   ofsDebug<<"infMat <- ";
//   printDMatrixRformat(&ofsDebug, infMat, nbeta, nbeta);
  //  ofsDebug.close();
}
} // extern "C"



