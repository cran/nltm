#include <string>
#include <math.h>
#include <vector>
#include "Rostream.h"

using namespace std;

#define CENSOR 0
#define FAILURE 1
#define PH 0
#define PHC 1
#define PO 2
#define PHPHC 3
#define PHPOC 4
#define GFM 5
#define PHPO 6
#define ERROR_CODE -1
#define BIG 100

// s: cumulative survival function

// pred: predictor. 
//  single predictor: pred=exp(x*beta)
//  single predictor cure: pred=exp(x*beta1+betac), beta=(beta1, betac)
//  two predictors: pred=(exp(x*beta1), exp(x*beta2)), beta=(beta1, beta2)
//  two predictors cure: pred=(exp(x*beta1+betac), exp(x*beta2)), 
//  beta=(beta1, beta2, betac)

// gammaD1 didn't seem to be used except when defining vtheta for
// failure observations. vtheta for failure is gammaD1*s so it is
// better to define it, similarly for the derivatives

// gammaD2 doesn't seem to be in use - was not checked

#define gammaPH(pred, s)(s<=0 ? 0 : (s>=1 ? 1 : pow(s, pred)))
#define vthetafPH(pred,s)(s<=0 ? 0 : (s>=1 ? pred : pred*pow(s, pred)))
#define gammaD2PH(pred, s)(s<=0 ? 0 : (s>=1 ? pred*(pred-1) : pred*(pred-1)*pow(s, pred-2)))

#define gammaPHC(pred, s)(s<=0 ? exp(-pred) : (s>=1 ? 1 : exp(-pred*(1-s))))
#define vthetafPHC(pred, s)(s<=0 ? 0 : (s>=1 ? pred : pred*exp(-pred*(1-s))*s))
#define gammaD2PHC(pred, s)(s<=0 ? pred*pred*exp(-pred) : (s>=1 ? pred*pred : pred*pred*exp(-pred*(1-s))))
#define gammaPHC_pred(pred, s)(s<=0 ? -exp(-pred) : (s>=1 ? 0 : -exp(-pred*(1-s))*(1-s)))
#define vthetafPHC_pred(pred, s)(s<=0 ? 0 : (s>=1 ? 1 : (1-pred*(1-s))*exp(-pred*(1-s))*s))
#define gammaPHC_2pred(pred, s)(s<=0 ? exp(-pred) : (s>=1 ? 0 : (1-s)*(1-s)*exp(-pred*(1-s))))
#define vthetafPHC_2pred(pred, s)(s<=0 || s>=1 ? 0 : (1-s)*(-2+pred*(1-s))*exp(-pred*(1-s))*s)

#define gammaPO(pred, s)(s<=0 ? 0 : (s>=1 ? 1 : pred/(pred-log(s))))

///////////////////////////////
//         PH model          //
///////////////////////////////
// First derivative of gammaPH with respect to pred
double gammaPH_pred(double pred, double s)
{
  double res;

  if(s<=0 || s>=1)
    return(0);
  else{
    res=pow(s,pred)*log(s);
    return((isnan(res) ? 0 : res));
  }
}

double vthetafPH_pred(double pred, double s)
{
  double res;

  if(s<=0)
    return(0);
  else{
    if(s>=1)
      return(1);
    else{
      res=pow(s,pred)*(1+pred*log(s));
      return((isnan(res) ? 0 : res));
    }
  }
}

// Second derivative of gammaPH with respect to pred
double gammaPH_2pred(double pred, double s)
{
  double aux;

  if(s<=0 || s>=1)
    return(0);
  else{
    aux=log(s);
    aux=pow(s, pred)*aux*aux;
    return((isnan(aux) ? 0 : aux));
  }
}

// s times second derivative of gammaD1PH with respect to pred
double vthetafPH_2pred(double pred, double s)
{
  double aux;

  if(s<=0 || s>=1)
    return(0);
  else{
    aux=log(s);
    aux=pow(s, pred)*aux*(2+pred*aux);
    return((isnan(aux) ? 0 : aux));
  }
}

double ThetonPH(double pred, double s, int cc)
{
  switch(cc){
  case CENSOR:
    return(pred);
    break;
  case FAILURE:
    return(pred);
    break;
  default:
    Rcerr<<"ThetonPH: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

double ThetonPH_pred(double pred, double s, int cc)
{
  switch(cc){
  case CENSOR:
    return(1);
    break;
  case FAILURE:
    return(1);
    break;
  default:
    Rcerr<<"ThetonPH_pred: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

double ThetonPH_h(double pred, double s, int cc)
{
  switch(cc){
  case CENSOR:
    return(0);
    break;
  case FAILURE:
    return(0);
    break;
  default:
    Rcerr<<"ThetonPH_h: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}


///////////////////////////////
//         PHC model         //
///////////////////////////////
double ThetonPHC(double pred, double s, int cc)
{
  switch(cc){
  case CENSOR:
    return((s<=0 ? 0 : (s>=1 ? pred : pred*s)));
    break;
  case FAILURE:
    return((s<=0 ? 1 : (s>=1 ? 1+pred : 1+pred*s)));
    break;
  default:
    Rcerr<<"ThetonPHC: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

double ThetonPHC_pred(double pred, double s, int cc)
{
  switch(cc){
  case CENSOR:
    return((s<=0 ? 0 : (s>=1 ? 1 : s)));
    break;
  case FAILURE:
    return((s<=0 ? 0 : (s>=1 ? 1 : s)));
    break;
  default:
    Rcerr<<"ThetonPHC: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

double ThetonPHC_h(double pred, double s, int cc)
{
  switch(cc){
  case CENSOR:
    return((s<=0 ? 0 : (s>=1 ? -pred : -pred*s)));
    break;
  case FAILURE:
    return((s<=0 ? 0 : (s>=1 ? -pred : -pred*s)));
    break;
  default:
    Rcerr<<"ThetonPHC: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

double ThetonCurePHC(double pred, double s)
{
  double aux;

  if(s<=0)
    return(1);
  else{
    if(s>=1)
      return(pred/(1-exp(-pred)));
    else{
      aux=s*pred/(1-exp(-pred*s));
      return((isnan(aux) ? 1 : aux));
    }
  }
}

double ThetonCurePHC_pred(double pred, double s)
{
  double aux;

  if(s<=0)
    return(0);
  else{
    if(s>=1){
      aux=exp(-pred);
      return((1-aux*(1+pred))/(1-aux)/(1-aux));
    }else{
      aux=exp(-pred*s);
      aux=s*(1-aux*(1+s*pred))/(1-aux)/(1-aux);
      return((isnan(aux) ? 0 : aux));
    }
  }
}

double ThetonCurePHC_h(double pred, double s)
{
  double aux;

  if(s<=0)
    return(0);
  else{
    if(s>=1){
      aux=exp(-pred);
      return(-pred*(1-aux*(1+pred))/(1-aux)/(1-aux));
    }else{
      aux=exp(-pred*s);
      aux=-pred*s*(1-aux*(1+s*pred))/(1-aux)/(1-aux);
      return((isnan(aux) ? 0 : aux));
    }
  }
}


///////////////////////////////
//          PO model         //
///////////////////////////////
// s times first derivative of gammaPO with respect to s
double vthetafPO(double pred, double s)
{
  double aux;

  if(s<=0)
    return(0);
  else{
    if(s>=1)
      return(1/pred);
    else{
      aux=pred-log(s);
      return(pred/aux/aux);      
    }
  }
}

// Second derivative of gammaPO with respect to s
double gammaD2PO(double pred, double s)
{
  double aux;

  if(s<=0)
    return(0);
  else{
    if(s>=1)
      return((pred-2)/(pred*pred));
    else{
      aux=pred-log(s);
      return(-pred*(aux-2)/(s*s*aux*aux*aux));
    }
  }
}

// First derivative of gammaPO with respect to pred
double gammaPO_pred(double pred, double s)
{
  double aux;

  if(s<=0 || s>=1)
    return(0);
  else{
    aux=-log(s);	
    return((aux>BIG ? 1/(pred*pred/aux+2*pred+aux) : 
	    aux/(pred+aux)/(pred+aux)));
  }
}   

// s times first derivative of gammaD1PO with respect to pred
double vthetafPO_pred(double pred, double s)
{
  double aux;

  if(s<=0)
    return(0);
  else{
    if(s>=1)
      return(-1/(pred*pred));
    else{
      aux=pred-log(s);
      aux=(aux-2*pred)/aux/aux/aux;
      return((isnan(aux) ? 0 : aux));
    }
  }
}

// Second derivative of gammaPO with respect to pred
double gammaPO_2pred(double pred, double s)
{
  double aux;

  if(s<=0 || s>=1)
    return(0);
  else{
    aux=pred-log(s);
    aux=-2*(aux-pred)/aux/aux/aux;
    return((isnan(aux) ? 0 : aux));
  }
}

// s times second derivative of gammaD1PO with respect to pred
double vthetafPO_2pred(double pred, double s)
{
  double aux1, aux2;

  if(s<=0)
    return(0);
  else{
    if(s>=1)
      return(2/(pred*pred*pred));
    else{
      aux1=log(s);
      aux2=pred-aux1;
      aux2=2*(pred+2*aux1)/aux2/aux2/aux2/aux2;
      return((isnan(aux2) ? 0 : aux2));
    }
  }
}

double ThetonPO(double pred, double s, int cc)
{
  switch(cc){
  case CENSOR:
    return((s<=0 ? 0 : (s>=1 ? 1/pred : 1/(pred-log(s)))));
    break;
  case FAILURE:
    return((s<=0 ? 0 : (s>=1 ? 2/pred : 2/(pred-log(s)))));
    break;
  default:
    Rcerr<<"ThetonPO: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

// ThetonPO_h equals ThetonPO_pred
double ThetonPO_pred(double pred, double s, int cc)
{
  double aux;
  
  switch(cc){
  case CENSOR:
    if(s<=0)
      return(0);
    else{
      if(s>=1)
	return(-1/pred/pred);
      else{
	aux=1/(pred-log(s));
	return(-aux*aux);
      }
    }
    break;
  case FAILURE:
    if(s<=0)
      return(0);
    else{
      if(s>=1)
	return(-2/pred/pred);
      else{
	aux=1/(pred-log(s));
	return(-2*aux*aux);
      }
    }
    break;
  default:
    Rcerr<<"ThetonPO_pred and TheronPO_h: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}


///////////////////////////////
//        PHPHC model        //    
///////////////////////////////
double gammaPHPHC(vector<double> &pred, double s)
{
  return((s<=0 ? exp(-pred[0]) : (s>=1 ? 1 : 
				  exp(-pred[0]*(1-pow(s, pred[1]))))));
}

// s times gammaD1PHPHC
double vthetafPHPHC(vector<double> &pred, double s)
{
  double aux;
  
  if(s<=0)
    return(0);
  else{
    if(s>=1)
      return(pred[0]*pred[1]);
    else{
      aux=pow(s, pred[1]);
      return(exp(-pred[0]*(1-aux))*pred[0]*pred[1]*aux);
    }
  }
}

double gammaD2PHPHC(vector<double> &pred, double s)
{
  double aux1, aux2, prod;

  if(s<=0)
    return(0);
  else{
    prod=pred[0]*pred[1];
    if(s>=1)
      return(prod*(prod+pred[1]-1));
    else{
      aux1=pow(s, pred[1]-2);
      aux2=aux1*s;
      return(exp(-pred[0]*(1-aux2*s))*prod*(prod*aux2*aux2+(pred[1]-1)*aux1));
    }
  }
}

// First derivative of gammaPHPHC with respect to pred
void gammaPHPHC_pred(vector<double> &pred, double s, vector<double> &der1)
{
  double aux1, aux2;

  if(s<=0){
    der1[0]=-exp(-pred[0]);
    der1[1]=0;
  }else{
    if(s>=1){
      der1[0]=0;
      der1[1]=0;
    }else{
      aux1=pow(s, pred[1]); 
      aux2=exp(-pred[0]*(1-aux1));
      der1[0]=-aux2*(1-aux1);
      der1[1]=aux2*pred[0]*aux1*log(s);
      der1[1]=(isnan(der1[1]) ? 0 : der1[1]);
    }
  }
}

// s times first derivative of gammaD1PHPHC with respect to pred
void vthetafPHPHC_pred(vector<double> &pred, double s, vector<double> &der1)
{
  double aux1, aux2, aux3;

  if(s<=0){
    der1[0]=0;
    der1[1]=0;
  }else{
    if(s>=1){
      der1[0]=pred[1];
      der1[1]=pred[0];
    }else{
      aux1=pow(s, pred[1]);
      aux2=exp(-pred[0]*(1-aux1))*aux1; 
      aux3=log(s);
      der1[0]=aux2*pred[1]*(-pred[0]*(1-aux1)+1);
      der1[1]=aux2*pred[0]*(pred[0]*pred[1]*aux1*aux3+1+pred[1]*aux3);
    }
  }
}

// Second derivative of gammaPHPHC with respect to pred
// der2[0]: gamma_{theta,theta}
// der2[1]: gamma_{eta,eta}
// der2[2]: gamma_{theta,eta}
void gammaPHPHC_2pred(vector<double> &pred, double s, vector<double> &der2)
{
  double aux1, aux2, aux3;

  if(s<=0){
    der2[0]=exp(-pred[0]);
    der2[1]=0;
    der2[2]=0;
  }else{
    if(s>=1){
      der2[0]=0;
      der2[1]=0;
      der2[2]=0;
    }else{
      aux1=pow(s, pred[1]);
      aux2=exp(-pred[0]*(1-aux1));
      aux3=log(s);
   
      der2[0]=aux2*(1-aux1)*(1-aux1);
      der2[1]=aux2*pred[0]*aux3*aux3*aux1*(pred[0]*aux1+1);
      der2[1]=(isnan(der2[1]) ? 0 : der2[1]);
      der2[2]=aux2*aux1*aux3*(-pred[0]*(1-aux1)+1);
      der2[2]=(isnan(der2[2]) ? 0 : der2[2]);
    }
  }
}

// s times second derivative of gammaD1PHPHC with respect to pred
void vthetafPHPHC_2pred(vector<double> &pred, double s, vector<double> &der2)
{
  double aux1, aux2, aux3, aux4;

  if(s<=0){
    der2[0]=0;
    der2[1]=0;
    der2[2]=0;
  }else{
    if(s>=1){
      der2[0]=0;
      der2[1]=0;
      der2[2]=1;
    }else{
      aux1=pow(s, pred[1]);
      aux2=aux1*exp(-pred[0]*(1-aux1));
      aux3=log(s);
      aux4=pred[0]*pred[1]*aux3;
   
      der2[0]=-aux2*pred[1]*(1-aux1)*(2-pred[0]*(1-aux1));
      der2[0]=(isnan(der2[0]) ? 0 : der2[0]);
      der2[1]=aux2*pred[0]*aux3*(pred[0]*aux1*aux1*aux4+2*pred[0]*aux1+
				 3*aux1*aux4+pred[1]*aux3+2);
      der2[1]=(isnan(der2[1]) ? 0 : der2[1]);
      der2[2]=aux2*(1+3*aux1*aux4+pred[1]*aux3-pred[0]*aux1*(1-aux1)*aux4-
		    pred[0]*(1-aux1)-aux4);
      der2[2]=(isnan(der2[2]) ? 0 : der2[2]);
    }
  }
}

double ThetonPHPHC(vector<double> &pred, double s, int cc)
{
  switch(cc){
  case CENSOR:
    return((s<=0 ? 0 : (s>=1 ? pred[0]*pred[1] : 
			pred[0]*pred[1]*pow(s, pred[1])))); 
    break;
  case FAILURE:
    return((s<=0 ? pred[1] : (s>=1 ? pred[1]*(pred[0]+1) :
			      pred[1]*(pred[0]*pow(s, pred[1])+1))));
    break;
  default:
    Rcerr<<"ThetonPHPHC: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

void ThetonPHPHC_pred(vector<double> &pred, double s, int cc, 
		      vector<double> &der1)
{
  double aux1;

  switch(cc){
  case CENSOR:
    if(s<=0){
      der1[0]=0;
      der1[1]=0;
    }else{
      if(s>=1){
	der1[0]=pred[1];
	der1[1]=pred[0];
      }else{
	aux1=pow(s,pred[1]);
	der1[0]=pred[1]*aux1;
	der1[1]=pred[0]*aux1*(1+pred[1]*log(s));
	der1[1]=(isnan(der1[1]) ? 0 : der1[1]);
      }
    }
    break;
  case FAILURE:
    if(s<=0){
      der1[0]=0;
      der1[1]=1;
    }else{
      if(s>=1){
	der1[0]=pred[1];
	der1[1]=pred[0]+1;
      }else{
	aux1=pow(s,pred[1]);
	der1[0]=pred[1]*aux1;
	der1[1]=pred[0]*aux1*(1+pred[1]*log(s))+1;
	der1[1]=(isnan(der1[1]) ? 1 : der1[1]);
      }
    }
    break;
  default:
    Rcerr<<"ThetonPHPHC_pred: Observation not censored or failure"<<endl;
  }
}

double ThetonPHPHC_h(vector<double> &pred, double s, int cc)
{
  switch(cc){
  case CENSOR:
    return((s<=0 ? 0 : (s>=1 ? -pred[0]*pred[1]*pred[1] : 
			-pred[0]*pred[1]*pred[1]*pow(s, pred[1]))));
    break;
  case FAILURE:
    return((s<=0 ? 0 : (s>=1 ? -pred[0]*pred[1]*pred[1] : 
			-pred[0]*pred[1]*pred[1]*pow(s, pred[1]))));
    break;
  default:
    Rcerr<<"ThetonPHPHC_h: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

double ThetonCurePHPHC(vector<double> &pred, double s)
{
  double aux;

  if(s<=0)
    return(pred[1]);
  else{
    if(s>=1)
      return(pred[0]*pred[1]/(1-exp(-pred[0])));
    else{
      aux=pow(s, pred[1])*pred[0];
      aux=pred[1]*aux/(1-exp(-aux));
      return((isnan(aux) ? pred[1] : aux));
    }
  }
}

void ThetonCurePHPHC_pred(vector<double> &pred, double s, vector<double> &der1)
{
  double aux1, aux2;

  if(s<=0){
    der1[0]=0;
    der1[1]=1;
  }else{
    if(s>=1){
      aux1=exp(-pred[0]);
      der1[0]=pred[1]*(1-aux1*(1+pred[0]))/(1-aux1)/(1-aux1);
      der1[1]=pred[0]/(1-aux1);
    }else{
      aux1=pow(s, pred[1]);
      aux2=exp(-aux1*pred[0]);
      der1[0]=pred[1]*aux1*(1-aux2*(1+aux1*pred[0]))/(1-aux2)/(1-aux2);
      der1[0]=(isnan(der1[0]) ? 0 : der1[0]);
      aux1*=pred[0];
      der1[1]=aux1/(1-aux2)*(1+pred[1]*log(s)*(1-aux2*(1+aux1))/(1-aux2));
      der1[1]=(isnan(der1[1]) ? 1 : der1[1]);
    }
  }
}

double ThetonCurePHPHC_h(vector<double> &pred, double s)
{
  double aux1, aux2;

  if(s<=0)
    return(0);
  else{
    if(s>=1){
      aux1=exp(-pred[0]);
      return(-pred[0]*pred[1]*pred[1]*(1-aux1*(1+pred[0]))/(1-aux1)/(1-aux1));
    }else{
      aux1=pow(s, pred[1])*pred[0];
      aux2=exp(-aux1);
      aux1=-pred[1]*pred[1]*aux1*(1-aux2*(1+aux1))/(1-aux2)/(1-aux2);
      return((isnan(aux1) ? 0 : aux1));
    }
  }
}


///////////////////////////////
//        PHPOC model        //    
///////////////////////////////
double gammaPHPOC(vector<double> &pred, double s)
{
  return((s<=0 ? exp(-pred[0]) : (s>=1 ? 1 : 
				  exp(-pred[0]*(1-s)/(1-(1-pred[1])*s)))));
}

// s times gammaD1PHPOC
double vthetafPHPOC(vector<double> &pred, double s)
{
  double aux;

  if(s<=0)
    return(0);
  else{
    if(s>=1)
      return(pred[0]/pred[1]);
    else{
      aux=1-(1-pred[1])*s;
      return(exp(-pred[0]*(1-s)/aux)*pred[0]*pred[1]/(aux*aux)*s);
    }
  }
}

double gammaD2PHPOC(vector<double> &pred, double s)
{
  double aux, prod=pred[0]*pred[1];

  if(s<=0)
    return(exp(-pred[0])*prod*prod+2*prod*(1-pred[1]));
  else{
    if(s>=1)
      return(prod*prod+2*prod*pred[1]*(1-pred[1]));
    else{
      aux=1-(1-pred[1])*s;
      return(exp(-pred[0]*(1-s)/aux)*prod/pow(aux, 4)*
	     (prod+2*(1-pred[1])-2*(1-pred[1])*(1-pred[1])*s));
    }
  }
}

// First derivative of gammaPHPOC with respect to pred
void gammaPHPOC_pred(vector<double> &pred, double s, vector<double> &der1)
{
  double aux;

  if(s<=0){
    der1[0]=-exp(-pred[0]);
    der1[1]=0;
  }else{
    if(s>=1){
      der1[0]=0;
      der1[1]=0;
    }else{
      aux=1-(1-pred[1])*s;
      der1[0]=-exp(-pred[0]*(1-s)/aux)*(1-s)/aux;
      der1[1]=-der1[0]*pred[0]*s/aux;
    }
  }
}

// s times first derivative of gammaD1PHPOC with respect to pred
void vthetafPHPOC_pred(vector<double> &pred, double s, vector<double> &der1)
{
  double aux1, aux2;

  if(s<=0){
    der1[0]=0;
    der1[1]=0;
  }else{
    if(s>=1){
      der1[0]=1/pred[1];
      der1[1]=-pred[0]/pred[1]/pred[1];
    }else{
      aux1=1-(1-pred[1])*s;
      aux2=exp(-pred[0]*(1-s)/aux1)*s;
      der1[0]=aux2*pred[1]/(aux1*aux1)*(1-pred[0]*(1-s)/aux1);
      der1[1]=aux2*pred[0]/(aux1*aux1*aux1)*
	(pred[0]*pred[1]*s*(1-s)/aux1+1-(1+pred[1])*s);
    }
  }
}

// Second derivative of gammaPHPOC with respect to pred
// der2[0]: gamma_{theta,theta}
// der2[1]: gamma_{eta,eta}
// der2[2]: gamma_{theta,eta}
void gammaPHPOC_2pred(vector<double> &pred, double s, vector<double> &der2)
{
  double aux1, aux2;

  if(s<=0){
    der2[0]=exp(-pred[0]);
    der2[1]=0;
    der2[2]=0;
  }else{
    if(s>=1){
      der2[0]=0;
      der2[1]=0;
      der2[2]=0;
    }else{
      aux1=1-(1-pred[1])*s;
      aux2=exp(-pred[0]*(1-s)/aux1)*(1-s)/aux1;
      der2[0]=aux2*(1-s)/aux1;
      der2[1]=aux2*pred[0]*s/aux1*s/aux1*(pred[0]*(1-s)/aux1-2);
      der2[2]=aux2*s/aux1*(1-pred[0]*(1-s)/aux1);
    }
  }
}

// s times second derivative of gammaD1PHPOC with respect to pred
void vthetafPHPOC_2pred(vector<double> &pred, double s, vector<double> &der2)
{
  double aux1, aux2, aux3, aux4, aux5;

  if(s<=0){
    der2[0]=0;
    der2[1]=0;
    der2[2]=0;
  }else{
    if(s>=1){
      der2[0]=0;
      aux1=pred[1]*pred[1];
      der2[1]=2*pred[0]/aux1/pred[1];
      der2[2]=-1/aux1;
    }else{
      aux1=1-(1-pred[1])*s;
      aux2=aux1*aux1;
      aux3=aux1*aux2;
      aux5=pred[0]*(1-s); 
      aux4=exp(-aux5/aux1)*s;
      der2[0]=aux4*pred[1]*(1-s)/aux3/aux1*(aux5-2+2*(1-pred[1])*s);
      der2[1]=aux4*pred[0]/aux3*
	((aux5*s-3*s*aux1)/aux2*(pred[1]*s*aux5/aux1+1-(1+pred[1])*s)+
	 aux5*s*(1-s)/aux2-s);
      der2[2]=aux4/aux3*(-pred[1]*aux5*aux5*s/aux1/aux1-
			 aux5/aux1*(1-s*(1+3*pred[1]))+1-(1+pred[1])*s);
    }
  }
}

double ThetonPHPOC(vector<double> &pred, double s, int cc)
{
  double aux;

  switch(cc){
  case CENSOR:
    if(s<=0)
      return(0);
    else{
      if(s>=1)
	return(pred[0]/pred[1]);
      else{
	aux=1-(1-pred[1])*s;
	return(pred[0]*pred[1]*s/aux/aux);
      }
    }
    break;
  case FAILURE:
    if(s<=0)
      return(1);
    else{
      if(s>=1)
	return((2+pred[0]-pred[1])/pred[1]);
      else{
	aux=1-(1-pred[1])*s;
	return((pred[0]*pred[1]*s/aux+(1+(1-pred[1])*s))/aux);
      }
    }
    break;
  default:
    Rcerr<<"ThetonPHPOC: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

void ThetonPHPOC_pred(vector<double> &pred, double s, int cc, 
		      vector<double> &der1)
{
  double aux1, aux2;

  switch(cc){
  case CENSOR:
    if(s<=0){
      der1[0]=0;
      der1[1]=0;
    }else{
      if(s>=1){
	der1[0]=1/pred[1];
	der1[1]=-pred[0]/pred[1]/pred[1];
      }else{
	aux1=1-(1-pred[1])*s;
	aux2=aux1*aux1;
	der1[0]=pred[1]*s/aux2;
	der1[1]=pred[0]*s*(1-s-pred[1]*s)/aux2/aux1;
	der1[1]=(isnan(der1[1]) ? 0 : der1[1]);
      }
    }
    break;
  case FAILURE:
    if(s<=0){
      der1[0]=0;
      der1[1]=0;
    }else{
      if(s>=1){
	der1[0]=1/pred[1];
	der1[1]=-(2+pred[0])/pred[1]/pred[1];
      }else{
	aux1=1-(1-pred[1])*s;
	aux2=aux1*aux1;
	der1[0]=pred[1]*s/aux2;
	der1[1]=pred[0]*s*(1-s-pred[1]*s)/aux2/aux1-2*s/aux2;
	der1[1]=(isnan(der1[1]) ? 0 : der1[1]);
      }
    }
    break;
  default:
    Rcerr<<"ThetonPHPOC_pred: Observation not censored or failure"<<endl;
  }
}

double ThetonPHPOC_h(vector<double> &pred, double s, int cc)
{
  double aux1;

  switch(cc){
  case CENSOR:
    if(s<=0)
      return(0);
    else{
      if(s>=1)
	return(-pred[0]*(2-pred[1])/pred[1]/pred[1]);
      else{
	aux1=1-(1-pred[1])*s;
	return(-pred[0]/aux1*pred[1]/aux1*s*(1+(1-pred[1])*s)/aux1);
      }
    }
    break;
  case FAILURE:
    if(s<=0)
      return(0);
    else{
      if(s>=1)
	return((-2*pred[0]+pred[0]*pred[1]-2+2*pred[1])/pred[1]/pred[1]);
      else{
	aux1=1-(1-pred[1])*s;
	return(s*(-pred[0]*pred[1]*(1+(1-pred[1])*s)-2*(1-pred[1])*aux1)/aux1/
	       aux1/aux1);
      }
    }
    break;
  default:
    Rcerr<<"ThetonPHPOC_h: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

double ThetonCurePHPOC(vector<double> &pred, double s)
{
  double aux1, aux2, aux3;
  
  if(s<=0)
    return(1);
  else{
    if(s>=1)
      return(pred[0]/pred[1]/(1-exp(-pred[0])));
    else{
      aux1=pred[0]*pred[1];
      aux2=aux1*s;
      aux3=1-(1-pred[1])*s;
      aux3=aux2/aux3/aux3/(1-exp(-aux2/aux3));
      return((isnan(aux3) ? 1 : aux3));
    }
  }
}


void ThetonCurePHPOC_pred(vector<double> &pred, double s, vector<double> &der1)
{
  double aux1, aux2, aux3, aux4;
  
  if(s<=0){
    der1[0]=0;
    der1[1]=0;
  }else{
    if(s>=1){
      aux1=exp(-pred[0]);
      der1[0]=(1-aux1*(1+pred[0]))/pred[1]/(1-aux1)/(1-aux1);
      der1[1]=-pred[0]/pred[1]/pred[1]/(1-exp(-pred[0]));
    }else{
      aux1=pred[0]*pred[1]*s;
      aux2=1-(1-pred[1])*s;
      aux3=exp(-aux1/aux2);
      aux4=aux2*aux2*aux2;
      der1[0]=pred[1]*s*(aux2-aux3*(aux2+aux1))/aux4/(1-aux3)/(1-aux3);
      der1[0]=(isnan(der1[0]) ? 0 : der1[0]);
      der1[1]=pred[0]*s*((1-aux3)*(1-(1+pred[1])*s)-aux3*aux1*(1-s)/aux2)/
	aux4/(1-aux3)/(1-aux3);
      der1[1]=(isnan(der1[1]) ? 0 : der1[1]);
    }
  }
}


double ThetonCurePHPOC_h(vector<double> &pred, double s)
{
  double aux1, aux2, aux3;
  
  if(s<=0)
    return(0);
  else{
    if(s>=1){
      aux1=exp(-pred[0]);
      aux2=(1-aux1);
      return(-pred[0]*((2-pred[1])*aux2-pred[0]*aux1)/pred[1]/pred[1]/aux2
	     /aux2);
    }else{
      aux1=pred[0]*pred[1]*s;
      aux2=1-(1-pred[1])*s;
      aux3=exp(-aux1/aux2);
      aux3=-aux1*((1+(1-pred[1])*s)*(1-aux3)-aux1*aux3/aux2)/aux2/aux2/aux2
	/(1-aux3)/(1-aux3);
      return((isnan(aux3) ? 0 : aux3));
    }
  }
}


///////////////////////////////
//          GF model         //
///////////////////////////////
double gammaGF(vector<double> &pred, double s)
{
  return((s<=0 ? 0 : (s>=1 ? 1 : pow(pred[0]/(pred[0]-log(s)), pred[1]))));
}

// s times first derivative of gammaGF with respect to s
double vthetafGF(vector<double> &pred, double s)
{
  double aux;

  if(s<=0)
    return(0);
  else{
    if(s>=1)
      return(pred[1]/pred[0]);
    else{
      aux=pred[0]-log(s);
      return(pred[1]*pow(pred[0]/aux, pred[1])/aux);
    }
  }	
}

// Second derivative of gammaGF with respect to s
double gammaD2GF(vector<double> &pred, double s)
{
  double aux;

  if(s<=0)
    return(0);
  else{
    if(s>=1)
      return(pred[1]/pred[0]*((pred[1]+1)/pred[0]-1));
    else{
      aux=pred[0]-log(s);
      return(pred[1]*pow(pred[0]/aux, pred[1])/(aux*s*s)*((pred[1]+1)/aux-1));
    }
  }
}

// First derivative of gammaGF with respect to pred
void gammaGF_pred(vector<double> &pred, double s, vector<double> &der1)
{
  double aux1, aux2, aux3;

  if(s<=0 || s>=1){
    der1[0]=0;
    der1[1]=0;
  }else{
    aux1=log(s); 
    aux2=pred[0]/(pred[0]-aux1);
    aux3=pow(aux2, pred[1]-1);
    der1[0]=-aux1*pred[1]*aux3/(pred[0]-aux1)/(pred[0]-aux1);
    der1[0]=(isnan(der1[0]) ? 0 : der1[0]);
    der1[1]=aux3*aux2*log(aux2);
    der1[1]=(isnan(der1[1]) ? 0 : der1[1]);
  } 
}

// s times first derivative of gammaD1GF with respect to pred
void vthetafGF_pred(vector<double> &pred, double s, vector<double> &der1)
{
  double aux1, aux2, aux3;

  if(s<=0){
    der1[0]=0;
    der1[1]=0;
  }else{
    if(s>=1){
      der1[0]=-pred[1]/(pred[0]*pred[0]);
      der1[1]=1/pred[0];
    }else{
      aux1=log(s);
      aux2=pred[0]-aux1;
      aux3=pow(pred[0]/aux2, pred[1]-1);
      der1[0]=-pred[1]*aux3/aux2/aux2/aux2*(pred[0]+pred[1]*aux1);
      der1[0]=(isnan(der1[0]) ? 0 : der1[0]);
      aux1=pred[0]/aux2;
      der1[1]=aux3*aux1/aux2*(1+pred[1]*log(aux1));
      der1[1]=(isnan(der1[1]) ? 0 : der1[1]);
    }
  }
}

// Second derivative of gammaGF with respect to pred
// der2[0]: gamma_{theta,theta}
// der2[1]: gamma_{eta,eta}
// der2[2]: gamma_{theta,eta}
void gammaGF_2pred(vector<double> &pred, double s, vector<double> &der2)
{
  double aux1, aux2, aux3, aux4;

  if(s<=0 || s>=1){
    der2[0]=0;
    der2[1]=0;
    der2[2]=0;
  }else{
    aux1=log(s);
    aux2=pred[0]-aux1;
    aux3=pred[0]/aux2;
    aux4=pow(aux3, pred[1]-2);
    der2[0]=aux1*pred[1]*aux4/pow(aux2, 4)*(2*pred[0]+(pred[1]-1)*aux1);
    der2[0]=(isnan(der2[0]) ? 0 : der2[0]);
    aux4*=aux3;
    der2[2]=-aux1;
    aux1=log(aux3);
    der2[2]*=aux4/(aux2*aux2)*(1+pred[1]*aux1);
    der2[2]=(isnan(der2[2]) ? 0 : der2[2]);
    der2[1]=aux4*aux3*aux1*aux1;
    der2[1]=(isnan(der2[1]) ? 0 : der2[1]);
  }
}

// s times second derivative of gammaD1GF with respect to pred 
void vthetafGF_2pred(vector<double> &pred, double s, vector<double> &der2)
{
  double aux1, aux2, aux3, aux4, aux5;

  if(s<=0){
    der2[0]=0;
    der2[1]=0;
    der2[2]=0;
  }else{
    if(s>=1){
      aux1=pred[0]*pred[0];
      der2[0]=2*pred[1]/(aux1*pred[0]);
      der2[1]=0;
      der2[2]=-1/aux1;
    }else{
      aux1=log(s); 
      aux2=log(pred[0]-aux1);
      aux3=log(pred[0]);
      aux4=pred[0]/(pred[0]-aux1);
      aux5=pow(aux4, pred[1]-2);
      der2[0]=pred[1]*aux5/pow(pred[0]-aux1, 5)*
	(2*pred[0]*pred[0]+4*pred[0]*pred[1]*aux1+
	 pred[1]*(pred[1]-1)*aux1*aux1);
      der2[0]=(isnan(der2[0]) ? 0 : der2[0]);
      aux5*=aux4;
      der2[2]=-aux5/(pred[0]-aux1)/(pred[0]-aux1)/(pred[0]-aux1)*
	(pred[0]*(1+pred[1]*aux3)+pred[1]*aux1*(2+pred[1]*aux3)-
	 pred[1]*(pred[0]+pred[1]*aux1)*aux2);
      der2[2]=(isnan(der2[2]) ? 0 : der2[2]);
      aux5*=aux4;
      der2[1]=aux5/(pred[0]-aux1)*
	(2*aux3+pred[1]*aux3*aux3-2*(pred[1]*aux3+1)*aux2+pred[1]*aux2*aux2);
      der2[1]=(isnan(der2[1]) ? 0 : der2[1]);

    }
  }
}

double ThetonGF(vector<double> &pred, double s, int cc)
{
  switch(cc){
  case CENSOR:
    return((s<=0 ? 0 : (s>=1 ? pred[1]/pred[0] : pred[1]/(pred[0]-log(s)))));
    break;
  case FAILURE:
    return((s<=0 ? 0 : (s>=1 ? (pred[1]+1)/pred[0] : 
			(pred[1]+1)/(pred[0]-log(s)))));
    break;
  default:
    Rcerr<<"ThetonGF: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

void ThetonGF_pred(vector<double> &pred, double s, int cc, 
		   vector<double> &der1)
{
  double aux;

  switch(cc){
  case CENSOR:
    if(s<=0){
      der1[0]=0;
      der1[1]=0;
    }else{
      if(s>=1){
	der1[0]=-pred[1]/(pred[0]*pred[0]);
	der1[1]=1/pred[0];
      }else{
	aux=pred[0]-log(s);
	der1[0]=-pred[1]/(aux*aux);
	der1[1]=1/aux;
      }
    }
    break;
  case FAILURE:
    if(s<=0){
      der1[0]=0;
      der1[1]=0;
    }else{
      if(s>=1){
	der1[0]=-(pred[1]+1)/(pred[0]*pred[0]);
	der1[1]=1/pred[0];
      }else{
	aux=pred[0]-log(s);
	der1[0]=-(pred[1]+1)/(aux*aux);
	der1[1]=1/aux;
      }
    }
    break;
  default:
    Rcerr<<"ThetonGF_pred: Observation not censored or failure"<<endl;
  }
}

// ThetonGF_h equals ThetonGF_pred0
double ThetonGF_h(vector<double> &pred, double s, int cc)
{
  double aux;

  switch(cc){
  case CENSOR:
    if(s<=0)
      return(0);
    else{
      if(s>=1)
	return(-pred[1]/(pred[0]*pred[0]));
      else{
	aux=pred[0]-log(s);
	return(-pred[1]/(aux*aux));
      }
    }
    break;
  case FAILURE:
    if(s<=0)
      return(0);
    else{
      if(s>=1)
	return(-(pred[1]+1)/(pred[0]*pred[0]));
      else{
	aux=pred[0]-log(s);
	return(-(pred[1]+1)/(aux*aux));
      }
    }
    break;
  default:
    Rcerr<<"ThetonGF_h: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}


///////////////////////////////
//        PHPO model         //    
///////////////////////////////
double gammaPHPO(vector<double> &pred, double s)
{
  return((s<=0 ? 0 : (s>=1 ? 1 : pred[0]/(pow(s, -pred[1])-(1-pred[0])))));
}

// s times first derivative of gammaPHPO with respect to s
double vthetafPHPO(vector<double> &pred, double s)
{
  double aux1, aux2;

  if(s<=0)
    return(0);
  else{
    if(s>=1)
      return(pred[1]/pred[0]);
    else{
      aux1=pow(s,pred[1]);
      aux2=1-(1-pred[0])*aux1;
      return(pred[0]*pred[1]*aux1/aux2/aux2);
    }
  }
}

// Second derivative of gammaPHPO with respect to s
double gammaD2PHPO(vector<double> &pred, double s)
{
  double aux1, aux2;

  if(s<=0)
    return(0);
  else{
    if(s>=1)
      return(-pred[0]*(1-1/pred[1]));
    else{
      aux1=pow(s, pred[1]-2);
      aux2=aux1*s*s*(1-pred[0]);
      return(pred[0]*pred[1]*aux1/pow(1-aux2, 3)*
	     (pred[1]-1+(pred[1]+1)*aux2));
    }
  }
}

// First derivative of gammaPHPO with respect to pred
void gammaPHPO_pred(vector<double> &pred, double s, vector<double> &der1)
{
  double aux1, aux2;

  if(s<=0 || s>=1){
    der1[0]=0;
    der1[1]=0;
  }else{
    aux1=pow(s, pred[1]);
    aux2=1-(1-pred[0])*aux1;
    aux2=aux1/(aux2*aux2);
    der1[0]=aux2*(1-aux1);
    der1[1]=pred[0]*aux2*log(s);
    der1[1]=(isnan(der1[1]) ? 0 : der1[1]);
  }
}

// s times first derivative of gammaD1PHPO with respect to pred
void vthetafPHPO_pred(vector<double> &pred, double s, vector<double> &der1)
{
  double aux1, aux2, aux3;

  if(s<=0){
    der1[0]=0;
    der1[1]=0;
  }else{
    if(s>=1){
      der1[0]=-pred[1]/(pred[0]*pred[0]);
      der1[1]=1/pred[0];
    }else{
      aux1=pow(s, pred[1]); 
      aux2=1-(1-pred[0])*aux1; 
      aux2=aux1/aux2/aux2/aux2;
      aux3=pred[1]*log(s);
      der1[0]=pred[1]*aux2*(1-(1+pred[0])*aux1);
      der1[1]=pred[0]*aux2*(1-(1-pred[0])*aux1*(1-aux3)+aux3);
    }
  }
}

// Second derivative of gammaPHPO with respect to pred
// der2[0]: gamma_{theta,theta}
// der2[1]: gamma_{eta,eta}
// der2[2]: gamma_{theta,eta}
void gammaPHPO_2pred(vector<double> &pred, double s, vector<double> &der2)
{
  double aux1, aux2, aux3;

  if(s<=0 || s>=1){
    der2[0]=0;
    der2[1]=0;
    der2[2]=0;
  }else{
    aux1=pow(s, pred[1]); 
    aux2=1-(1-pred[0])*aux1; 
    aux3=log(s);
    aux2=aux1/(aux2*aux2*aux2);
    der2[0]=-2*aux2*aux1*(1-aux1);
    aux2*=aux3;
    der2[1]=pred[0]*aux2*aux3*(1+(1-pred[0])*aux1);
    der2[1]=(isnan(der2[1]) ? 0 : der2[1]);
    der2[2]=aux2*(1-(1+pred[0])*aux1);
    der2[2]=(isnan(der2[2]) ? 0 : der2[2]);
  }
}

// s times second derivative of gammaD1PHPO with respect to pred
void vthetafPHPO_2pred(vector<double> &pred, double s, vector<double> &der2)
{
  double aux1, aux2, aux3, aux4;

  if(s<=0){
    der2[0]=0;
    der2[1]=0;
    der2[2]=0;
  }else{
    if(s>=1){
      der2[0]=2*pred[1]/(pred[0]*pred[0]*pred[0]);
      der2[1]=0;
      der2[2]=-1/(pred[0]*pred[0]);
    }else{
      aux1=pow(s, pred[1]); 
      aux2=1-(1-pred[0])*aux1; 
      aux2=aux1/aux2/aux2/aux2/aux2;
      aux3=log(s);
      der2[0]=-2*pred[1]*aux2*aux1*(2-(2+pred[0])*aux1);
      der2[1]=pred[0]*aux2*aux3;
      aux3*=pred[1];
      aux4=(1-pred[0])*aux1;
      der2[1]*=(2+aux3*(1+4*aux4+aux4*aux4)-2*aux4*aux4);
      der2[1]=(isnan(der2[1]) ? 0 : der2[1]);
      der2[2]=aux2*(1-2*aux1+(1-pred[0]*pred[0])*aux1*aux1+
		    aux3*(1-4*pred[0]*aux1-(1-pred[0]*pred[0])*aux1*aux1));
      der2[2]=(isnan(der2[2]) ? 0 : der2[2]);
    }
  }
}

double ThetonPHPO(vector<double> &pred, double s, int cc)
{
  double aux;

  switch(cc){
  case CENSOR:
    return((s<=0 ? pred[1] : (s>=1 ? pred[1]/pred[0] :
			      pred[1]/(1-(1-pred[0])*pow(s, pred[1])))));
    break;
  case FAILURE:
    if(s<=0)
      return(pred[1]);
    else{
      if(s>=1)
	return(pred[1]*(2/pred[0]-1));
      else{
	aux=(1-pred[0])*pow(s, pred[1]);
	return(pred[1]*(1+aux)/(1-aux));
      }
    }
    break;
  default:
    Rcerr<<"ThetonPHPO: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

void ThetonPHPO_pred(vector<double> &pred, double s, int cc, 
		     vector<double> &der1)
{
  double aux1, aux2, aux3;

  switch(cc){
  case CENSOR:
    if(s<=0){
      der1[0]=0;
      der1[1]=1;
    }else{
      if(s>=1){
	der1[0]=-pred[1]/(pred[0]*pred[0]);
	der1[1]=1/pred[0];
      }else{
	aux1=pow(s,pred[1]);
	aux2=1-(1-pred[0])*aux1;
	aux3=(1-pred[0])*aux1;
	aux2*=aux2;

	der1[0]=-pred[1]*aux1/aux2;
	der1[1]=(1-(1-pred[0])*aux1*(1-pred[1]*log(s)))/aux2;
      }
    }
    break;
  case FAILURE:
    if(s<=0){
      der1[0]=0;
      der1[1]=1;
    }else{
      if(s>=1){
	der1[0]=-2*pred[1]/(pred[0]*pred[0]);
	der1[1]=2/pred[0]-1;
      }else{
	aux1=pow(s,pred[1]);
	aux2=1-(1-pred[0])*aux1;
	aux3=(1-pred[0])*aux1;
	aux2*=aux2;

	der1[0]=-2*pred[1]*aux1/aux2;
	der1[1]=(1-aux3*aux3+2*pred[1]*aux3*log(s))/aux2;
      }
    }
    break;
  default:
    Rcerr<<"ThetonPHPO_pred: Observation not censored or failure"<<endl;
  }
}

double ThetonPHPO_h(vector<double> &pred, double s, int cc)
{
  double aux1, aux2;

  switch(cc){
  case CENSOR:
    if(s<=0)
      return(0);
    else{
      if(s>=1)
	return(-pred[1]/pred[0]*pred[1]/pred[0]*(1-pred[0]));
      else{
	aux1=(1-pred[0])*pow(s, pred[1]);
	aux2=pred[1]/(1-aux1);
	return(-aux2*aux2*aux1);
      }
    }
    break;
  case FAILURE:
    if(s<=0)
      return(0);
    else{
      if(s>=1)
	return(-2*pred[1]/pred[0]*pred[1]/pred[0]*(1-pred[0]));
      else{
	aux1=(1-pred[0])*pow(s, pred[1]);
	aux2=pred[1]/(1-aux1);
	return(-2*aux2*aux2*aux1);
      }
    }
    break;
  default:
    Rcerr<<"ThetonPHPO_h: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}


///////////////////////////////
//     General functions     //
///////////////////////////////
double gamma(vector<double> &pred, double s, int model)
{
  switch(model){
  case PH:
    return(gammaPH(pred[0], s));
    break;
  case PHC:
    return(gammaPHC(pred[0], s));
    break;
  case PO:
    return(gammaPO(pred[0], s));
    break;
  case PHPHC:
    return(gammaPHPHC(pred, s));
    break;
  case PHPOC:
    return(gammaPHPOC(pred, s));
    break;
  case GFM:
    return(gammaGF(pred, s));
    break;
  case PHPO:
    return(gammaPHPO(pred, s));
    break;
  default:
    Rcerr<<"gamma: Not one of the supported models"<<endl;
  }
  return(ERROR_CODE);
}

double vthetaf(vector<double> &pred, double s, int model)
{
  switch(model){
  case PH:
    return(vthetafPH(pred[0], s));
    break;
  case PHC:
    return(vthetafPHC(pred[0], s));
    break;
  case PO:
    return(vthetafPO(pred[0], s));
    break;
  case PHPHC:
    return(vthetafPHPHC(pred, s));
    break;
  case PHPOC:
    return(vthetafPHPOC(pred, s));
    break;
  case GFM:
    return(vthetafGF(pred, s));
    break;
  case PHPO:
    return(vthetafPHPO(pred, s));
    break;
  default:
    Rcerr<<"vthetaf: Not one of the supported models"<<endl;
  }
  return(ERROR_CODE);
}

// Doesn't seem to be in use - was not checked
double gammaD2(vector<double> &pred, double s, int model)
{
  switch(model){
  case PH:
    return(gammaD2PH(pred[0], s));
    break;
  case PHC:
    return(gammaD2PHC(pred[0], s));
    break;
  case PO:
    return(gammaD2PO(pred[0], s));
    break;
  case PHPHC:
    return(gammaD2PHPHC(pred, s));
    break;
  case PHPOC:
    return(gammaD2PHPOC(pred, s));
    break;
  case GFM:
    return(gammaD2GF(pred, s));
    break;
  case PHPO:
    return(gammaD2PHPO(pred, s));
    break;
  default:
    Rcerr<<"Not one of the supported models"<<endl;
  }
  return(ERROR_CODE);
}

void gamma_pred(vector<double> &pred, double s, int model, 
		vector<double> &der1)
{
  switch(model){
  case PH:
    der1[0]=gammaPH_pred(pred[0], s);
    break;
  case PHC:
    der1[0]=gammaPHC_pred(pred[0], s);
    break;
  case PO:
    der1[0]=gammaPO_pred(pred[0], s);
    break;
  case PHPHC:
    gammaPHPHC_pred(pred, s, der1);
  break;
  case PHPOC:
    gammaPHPOC_pred(pred, s, der1);
  break;
  case GFM:
    gammaGF_pred(pred, s, der1);
    break;
  case PHPO:
    gammaPHPO_pred(pred, s, der1);
    break;
  default:
    Rcerr<<"gamma_pred: Not one of the supported models"<<endl;
  }
}

void vthetaf_pred(vector<double> &pred, double s, int model, 
		  vector<double> &der1)
{
  switch(model){
  case PH:
    der1[0]=vthetafPH_pred(pred[0], s);
    break;
  case PHC:
    der1[0]=vthetafPHC_pred(pred[0], s);
    break;
  case PO:
    der1[0]=vthetafPO_pred(pred[0], s);
    break;
  case PHPHC:
    vthetafPHPHC_pred(pred, s, der1);
    break;
  case PHPOC:
    vthetafPHPOC_pred(pred, s, der1);
    break;
  case GFM:
    vthetafGF_pred(pred, s, der1);
    break;
  case PHPO:
    vthetafPHPO_pred(pred, s, der1);
    break;
  default:
    Rcerr<<"vthetaf_pred: Not one of the supported models"<<endl;
  }
}

void gamma_2pred(vector<double> &pred, double s, int model, 
		 vector<double> &der2)
{
  switch(model){
  case PH:
    der2[0]=gammaPH_2pred(pred[0], s);
    break;
  case PHC:
    der2[0]=gammaPHC_2pred(pred[0], s);
    break;
  case PO:
    der2[0]=gammaPO_2pred(pred[0], s);
    break;
  case GFM:
    gammaGF_2pred(pred, s, der2);
    break;
  case PHPHC:
    gammaPHPHC_2pred(pred, s, der2);
    break;
  case PHPOC:
    gammaPHPOC_2pred(pred, s, der2);
    break;
  case PHPO:
    gammaPHPO_2pred(pred, s, der2);
    break;
  default:
    Rcerr<<"gamma_2pred: Not one of the supported models"<<endl;
  }
}

void vthetaf_2pred(vector<double> &pred, double s, int model, 
		   vector<double> &der2)
{
  switch(model){
  case PH:
    der2[0]=vthetafPH_2pred(pred[0], s);
    break;
  case PHC:
    der2[0]=vthetafPHC_2pred(pred[0], s);
    break;
  case PO:
    der2[0]=vthetafPO_2pred(pred[0], s);
    break;
  case PHPHC:
    vthetafPHPHC_2pred(pred, s, der2);
    break;
  case PHPOC:
    vthetafPHPOC_2pred(pred, s, der2);
    break;
  case GFM:
    vthetafGF_2pred(pred, s, der2);
    break;
  case PHPO:
    vthetafPHPO_2pred(pred, s, der2);
    break;
  default:
    Rcerr<<"dammaD1_2pred: Not one of the supported models"<<endl;
  }
}

double vtheta(vector<double> &pred, double s, int cc, int model)
{
  switch(cc){
  case CENSOR:
    return(gamma(pred, s, model));
    break;
  case FAILURE:
    return(vthetaf(pred, s, model));
  default:
    Rcerr<<"vtheta: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

// First derivative of vtheta with respect to pred
void vtheta_pred(vector<double> &pred, double s, int cc, int model,
		 vector<double> &der1)
{
  switch(cc){
  case CENSOR:
    gamma_pred(pred, s, model, der1);
    break;
  case FAILURE:
    vthetaf_pred(pred, s, model, der1);
    break;
  default:
    Rcerr<<"vtheta_pred: Observation not censored or failure"<<endl;
  }
}

// Second derivative of vtheta with respect to pred
void vtheta_2pred(vector<double> &pred, double s, int cc, int model,
		  vector<double> &der2)
{
  switch(cc){
  case CENSOR:
    gamma_2pred(pred, s, model, der2);
    break;
  case FAILURE:
    vthetaf_2pred(pred, s, model, der2);
    break;
  default:
    Rcerr<<"vtheta_2pred: Observation not censored or failure"<<endl;
  }
}

double vthetaCure(vector<double> &pred, double s, int cc, int model)
{
  switch(cc){
  case CENSOR:
    return(gamma(pred, 0, model));
    break;
  case FAILURE:
    return(gamma(pred, s, model)-gamma(pred, 0, model));
    break;
  default:
    Rcerr<<"vthetaCure: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

void vthetaCure_pred(vector<double> &pred, double s, int cc, int model,
		     vector<double> &der1)
{
  int i;
  vector<double> d1;

  switch(cc){
  case CENSOR:
    gamma_pred(pred, 0, model, der1);
    break;
  case FAILURE:
    d1.resize(der1.size());
    gamma_pred(pred, s, model, der1);
    gamma_pred(pred, 0, model, d1);
    for(i=0; i<int(der1.size()); i++)
      der1[i]-=d1[i];
    break;
  default:
    Rcerr<<"vthetaCure_pred: Observation not censored or failure"<<endl;
  }
}

void vthetaCure_2pred(vector<double> &pred, double s, int cc, int model,
		      vector<double> &der2)
{
  int i;
  vector<double> d2;

  switch(cc){
  case CENSOR:
    gamma_2pred(pred, 0, model, der2);
    break;
  case FAILURE:
    d2.resize(der2.size());
    gamma_2pred(pred, s, model, der2);
    gamma_2pred(pred, 0, model, d2);
    for(i=0; i<int(der2.size()); i++)
      der2[i]-=d2[i];
    break;
  default:
    Rcerr<<"vthetaCure_2pred: Observation not censored or failure"<<endl;
  }
}

double Theton(vector<double> &pred, double s, int cc, int model)
{
  switch(model){
  case PH:
    return(ThetonPH(pred[0], s, cc));
    break;
  case PHC:
    return(ThetonPHC(pred[0], s, cc));
    break;
  case PO:
    return(ThetonPO(pred[0], s, cc));
    break;
  case PHPHC:
    return(ThetonPHPHC(pred, s, cc));
  break;
  case PHPOC:
    return(ThetonPHPOC(pred, s, cc));
  break;
  case GFM:
    return(ThetonGF(pred, s, cc));
    break;
  case PHPO:
    return(ThetonPHPO(pred, s, cc));
  break;
  default:
    Rcerr<<"Theton: Not one of the supported models"<<endl;
  }
  return(ERROR_CODE);
}

void Theton_pred(vector<double> &pred, double s, int cc, int model, 
		 vector<double> &der1)
{
  switch(model){
  case PH:
    der1[0]=ThetonPH_pred(pred[0], s, cc);
    break;
  case PHC:
    der1[0]=ThetonPHC_pred(pred[0], s, cc);
    break;
  case PO:
    der1[0]=ThetonPO_pred(pred[0], s, cc);
    break;
  case PHPHC:
    ThetonPHPHC_pred(pred, s, cc, der1);
    break;
  case PHPOC:
    ThetonPHPOC_pred(pred, s, cc, der1);
    break;
  case GFM:
    ThetonGF_pred(pred, s, cc, der1);
    break;
  case PHPO:
    ThetonPHPO_pred(pred, s, cc, der1);
    break;
  default:
    Rcerr<<"Theton_pred: Not one of the supported models"<<endl;
  }
}

double Theton_h(vector<double> &pred, double s, int cc, int model)
{
  switch(model){
  case PH:
    return(ThetonPH_h(pred[0], s, cc));
    break;
  case PHC:
    return(ThetonPHC_h(pred[0], s, cc));
    break;
  case PO:
    // ThetonPO_pred and Theton_h are the same
    return(ThetonPO_pred(pred[0], s, cc));
    break;
  case PHPHC:
    return(ThetonPHPHC_h(pred, s, cc));
    break;
  case PHPOC:
    return(ThetonPHPOC_h(pred, s, cc));
    break;
  case GFM:
    return(ThetonGF_h(pred, s, cc));
    break;
  case PHPO:
    return(ThetonPHPO_h(pred, s, cc));
    break;
  default:
    Rcerr<<"Theton_h: Not one of the supported models"<<endl;
  }
  return(ERROR_CODE);
}

double ThetonCure(vector<double> &pred, double s, int cc, int model)
{
  switch(cc){
  case CENSOR:
    return(0);
    break;
  case FAILURE:
    switch(model){
    case PHC:
      return(ThetonCurePHC(pred[0], s));
      break;
    case PHPHC:
      return(ThetonCurePHPHC(pred, s));
      break;
    case PHPOC:
      return(ThetonCurePHPOC(pred, s));
      break;
    default:
      Rcerr<<"ThetonCure: Not one of the supported models or not a cure model"
	  <<endl;
    }
    break;
  default:
    Rcerr<<"ThetonCure: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

void ThetonCure_pred(vector<double> &pred, double s, int cc, int model,
		     vector<double> &der1)
{
  int i;

  switch(cc){
  case CENSOR:
    for(i=0; i<int(der1.size()); i++)
      der1[i]=0;
    break;
  case FAILURE:
    switch(model){
    case PHC:
      der1[0]=ThetonCurePHC_pred(pred[0], s);
      break;
    case PHPHC:
      ThetonCurePHPHC_pred(pred, s, der1);
      break;
    case PHPOC:
      ThetonCurePHPOC_pred(pred, s, der1);
      break;
    default:
      Rcerr<<"ThetonCure_pred: Not one of the supported models"
	  <<"or not a cure model"<<endl;
    }
    break;
  default:
    Rcerr<<"ThetonCure_pred: Observation not censored or failure"<<endl;
  }
}

double ThetonCure_h(vector<double> &pred, double s, int cc, int model)
{
  switch(cc){
  case CENSOR:
    return(0);
    break;
  case FAILURE:
    switch(model){
    case PHC:
      return(ThetonCurePHC_h(pred[0], s));
      break;
    case PHPHC:
      return(ThetonCurePHPHC_h(pred, s));
      break;
    case PHPOC:
      return(ThetonCurePHPOC_h(pred, s));
      break;
    default:
      Rcerr<<"ThetonCure_h: Not one of the supported models or not a cure model"
	  <<endl;
    }
    break;
  default:
    Rcerr<<"ThetonCure_h: Observation not censored or failure"<<endl;
  }
  return(ERROR_CODE);
}

int nmodel(string model)
{
  if(model=="PH") return(PH);
  else{ 
    if(model=="PHC") return(PHC);
    else{
      if(model=="PO") return(PO);
      else{
	if(model=="PHPHC") return(PHPHC);
	else{
	  if(model=="PHPOC") return(PHPOC);
	  else{
	    if(model=="GFM") return(GFM);
	    else{
	      if(model=="PHPO") return(PHPO);
	      else
		Rcerr<<"nmodel: Not one of the supported models"<<endl;
	    }
	  }
	}
      }
    }
  }
  return(ERROR_CODE);
}

