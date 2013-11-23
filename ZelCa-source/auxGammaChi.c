/*
    Copyright (c) 2013       Svetlin Tassev
                             Princeton University
 
    This file is part of Zeldovich Calculator (ZelCa).

    ZelCa is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ZelCa is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ZelCa.  If not, see <http://www.gnu.org/licenses/>.
*/





//////////////////////////////////////////////////
///
///
/// This file is a mess for now, mixing GSL and CUBA integrators. 
///
/// All calculations in this code are for linear theory, so nothing new or exciting.
/// No effort was made to make those calculations fast. 
/// So, precomputing linear theory, or Gamma and Chi will be slow!
/// Once precomputed, the interpolation of Gamma and Chi in ChiGamma() is fast.
///
/// The only thing from this file one needs to know to do the ZA calculations is that:
/// 1. This file contains a function ChiGamma(Q,&chi,&gamma) which for a given Q
///    returns chi(Q) and gamma(Q) as defined in the paper. 
///    The code tries to be fast about it and not smart.
/// 2. Returns a 3x3 diagonal matrix PsiDiag with diag on its diagonal,
///    where diag=( sigma,sigma,sigma) for real space; or
///    diag=( sigma,sigma,sigma*RSfactor) for redshift space.
///    sigma is in fact \sigma^2 in the paper
/// 
///
//////////////////////////////////////////////////















#include <iostream>
#include <Eigen/LU>
using namespace std;
using namespace Eigen;


#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <cuba.h>
#include <stdio.h>

#define NTABLE 20000
#define EPS_SIGMA 1.e-4
extern double ns;
extern int CALCULATE_CHI_GAMMA;
double minK,maxK;
double minKP,maxKP;
double minQ,maxQ;
double GAMMA9999,CHI9999;
double sigma;
extern double  GROWTH;

extern string FILE_IN_PK;
extern string FILE_IN_OR_OUT_CHIGAMMAXI;
extern string FILE_OUT_LINEAR_THEORY;


gsl_interp_accel *accP,*accChi,*accGamma ;
gsl_spline *splineP,*splineChi,*splineGamma;
double   CHI(double Q);
double   GAMMA(double Q);
double   SIGMA(void);
double powerK(double K);

void   ChiGamma(double Q,double*chi,double*gamma);
void   initChiGamma(void);
void   finalizeSplines(int sw);
int LINEAR_RS;
double   XIL(double Q);
double   XILrs4(double Q);
double   XILrs2(double Q);
gsl_integration_workspace * w;
double x[NTABLE];
double y[NTABLE];
double z[NTABLE];
double xy[NTABLE];
double xz[NTABLE];
double   ChiT(double Q);
double   GammaT(double Q);

void   calcChiGamma(void);
void   calcLinearRS(void);

Matrix3d PsiDiag;
extern int REDSHIFT_SPACE;
extern double RSfactor;

double  INTd02(double Q);
double  INTd1(double Q);
double  INTd1k(double Q);
double  INTd2(double Q);


void initSplines(string filename)
{
    
    double x[NTABLE];
    double y[NTABLE];
    int i;
    
    FILE *input= fopen(filename.c_str(), "r");

    i=0;
    while(!feof(input))
     {    fscanf(input,"%lf%lf", x+i, y+i);
          x[i]=log(x[i]);
          y[i]=log(y[i]); 
          i++;
     }
    i--; 
     
    fclose(input);
    
    accP = gsl_interp_accel_alloc ();
    splineP = gsl_spline_alloc (gsl_interp_cspline, i);
    
    gsl_spline_init(splineP, x, y, i);
    minK=x[0];
    maxK=x[i-1];
    minKP=y[0];
    maxKP=y[i-1];
    
}

void initChiGamma(void)
{
    
    int i=0,n=NTABLE;

    double xi[n];
    
    if (CALCULATE_CHI_GAMMA){
        initSplines(FILE_IN_PK);
        w = gsl_integration_workspace_alloc (5000000);
        FILE *file;
        double minlQ,maxlQ;
        
        if (LINEAR_RS) {
            file= fopen(FILE_OUT_LINEAR_THEORY.c_str(), "w");
            minlQ=log10(1.);
            maxlQ=log10(250.);
        }
        else {
            file= fopen(FILE_IN_OR_OUT_CHIGAMMAXI.c_str(), "w");
            minlQ=-3;
            maxlQ=4;
        }

        
        maxQ=pow(10.,maxlQ);
        minQ=pow(10.,minlQ);
    
        n=((int)((maxlQ-minlQ)*200.0+0.5));

        sigma=SIGMA();
        CHI9999= CHI(maxQ*0.9999);
        GAMMA9999= GAMMA(maxQ*0.9999);
        fprintf(file,"%i %.8e %.8e\n",n,CHI9999,GAMMA9999);
        fprintf(file,"%.8le %.8le %.8le\n",minQ,maxQ,sigma);
        double q;
        for (i=0;i<n;i++){
            q=pow(10., ((double) i)/((double)n-1.0)*(maxlQ-minlQ)+minlQ);
            x[i]=q;
            y[i]=CHI(q);
            z[i]=GAMMA(q);
            if (LINEAR_RS)fprintf(file,"%.8le %.8le %.8le %.8le %.8le %.8le %.8le %.8le %.8le %.8le\n",x[i],y[i],z[i],XIL(q),XILrs2(q),XILrs4(q),INTd1k(q),INTd1(q),INTd2(q),INTd02(q));
            else fprintf(file,"%.8le %.8le %.8le %.8le\n",x[i],y[i],z[i],XIL(q));
        }
        fclose(file);
        gsl_integration_workspace_free (w);
        finalizeSplines(1);
    }
    else
    {
        FILE *input= fopen(FILE_IN_OR_OUT_CHIGAMMAXI.c_str(), "r");
        fscanf(input,"%i%lf%lf", &n,&CHI9999,&GAMMA9999);
        fscanf(input,"%lf%lf%lf", &minQ,&maxQ,&sigma);
        
        while(!feof(input))
        {    fscanf(input,"%lf%lf%lf%lf", x+i, y+i,z+i,xi+i);
            i++;
        }
     
        fclose(input);
    }
    
    sigma*=GROWTH*GROWTH;

    
    accChi = gsl_interp_accel_alloc ();
    splineChi = gsl_spline_alloc (gsl_interp_cspline, n);
    gsl_spline_init(splineChi, x, y, n);

    accGamma = gsl_interp_accel_alloc ();
    splineGamma = gsl_spline_alloc (gsl_interp_cspline, n);
    gsl_spline_init(splineGamma, x, z, n);

    x[0]=0;
    y[0]=sigma-EPS_SIGMA;
    z[0]=0;
    double q;
    for (i=1;i<=2000;i++){
        q=((double)i)*0.001;
        x[i]=q;
        y[i]=ChiT(q);
        z[i]=GammaT(q);
    }
    for (i=1;i<=2000;i++){
        q=((double)i)*0.01+2.0;
        x[i+2000]=q;
        y[i+2000]=ChiT(q);
        z[i+2000]=GammaT(q);
    }
    for (i=1;i<=2000;i++){
        q=((double)i)*0.1+22.0;
        x[i+4000]=q;
        y[i+4000]=ChiT(q);
        z[i+4000]=GammaT(q);
    }
    for (i=1;i<=2000;i++){
        q=((double)i)*1.0+222.0;
        x[i+6000]=q;
        y[i+6000]=ChiT(q);
        z[i+6000]=GammaT(q);
    }
    for (i=1;i<=2000;i++){
        q=((double)i)*10.0+2222.0;
        x[i+8000]=q;
        y[i+8000]=ChiT(q);
        z[i+8000]=GammaT(q);
    }
    
    for (i=1;i<=10000;i++){
        xy[i-1]=(y[i]-y[i-1])/(x[i]-x[i-1]);
        xz[i-1]=(z[i]-z[i-1])/(x[i]-x[i-1]);
    }
    xy[10000]=0;
    xz[10000]=0;
    
    Vector3d diag( sigma,sigma,sigma);
    if (REDSHIFT_SPACE)diag(2)*=RSfactor;
    PsiDiag=diag.asDiagonal();
    
    finalizeSplines(0);
}

//This is called millions of time. MUST BE FAST! So, simple linear interpolation by hand. NO GSL!
void   ChiGamma(double Q,double*chi,double*gamma){
       if (Q<=0.001) {
           chi[0]=sigma-EPS_SIGMA;
           gamma[0]=0;
           return;
       }
       if (Q>=22222.0){
            chi[0]=0;
            gamma[0]=0;
            return;
        }
       int i,i0;
       double sl,q0;
       sl=10.0;
       q0=2222.0;
       i0=8000;
       if(Q<=2222.0){
           sl=1.0;
           q0=222.0;
           i0=6000;
       }
       if(Q<=222.0){
           sl=0.1;
           q0=22.0;
           i0=4000;
       }
       if(Q<=22.0){
           sl=0.01;
           q0=2.0;
           i0=2000;
       }
       if(Q<=2.0){
           sl=0.001;
           q0=0.0;
           i0=0;
       }
       i=((int) ((Q-q0)/sl+0.5))+i0;
       chi[0]= xy[i]*(Q-x[i])+y[i];
       gamma[0]=xz[i]*(Q-x[i])+z[i];
       return ;
}





void finalizeSplines(int sw)
{
    if (sw==1){
        gsl_spline_free(splineP);
        gsl_interp_accel_free(accP);
    }
    else{
        gsl_spline_free(splineChi);
        gsl_interp_accel_free(accChi);
        
        gsl_spline_free(splineGamma);
        gsl_interp_accel_free(accGamma);
    }
}







double   ChiT(double Q){
       if (Q<minQ) return sigma-EPS_SIGMA;
       if (Q>maxQ) return GROWTH*GROWTH*CHI9999*pow(Q/(maxQ*0.9999),-1.0-ns);
       double res= GROWTH*GROWTH*gsl_spline_eval (splineChi, Q, accChi); 
       
       if (res>sigma-GammaT(Q)-EPS_SIGMA) return sigma-GammaT(Q)-EPS_SIGMA;
       return res;
}


double   GammaT(double Q){
       if (Q<minQ) return 0;
       if (Q>maxQ) return GROWTH*GROWTH*GAMMA9999*pow(Q/(maxQ*0.9999),-1.0-ns);
       return GROWTH*GROWTH*gsl_spline_eval (splineGamma, Q, accGamma); 
       
}


double powerK(double lK)
{
    if (lK<=minK) return  exp((lK-minK)*ns)*exp(minKP);
    if (lK>=maxK) return  0.0;
    return exp(gsl_spline_eval (splineP, lK, accP)); 

}



double min(double a, double b){
    if (a>b) return b;
    return a;
}

double max(double a, double b){
    if (a>b) return a;
    return b;
}

double j0(double x){
    if(x > 1E-6) return sin(x)/x;
    return 1.0-x*x/6.0;
}


double j02(double x){
    if(x > 1E-6) return 3.0*(sin(x)-cos(x)*x)/(x*x*x);
    return 1.0-x*x/10.0;
}


double j2(double x){
    if(x > 1E-6) return 3.0*(sin(x)-cos(x)*x)/(x*x*x)  - sin(x)/x;
    return x*x/15.0*(1.0-x*x/14.0);
}


double chiInt (double lK, void * params) {
       double K=exp(lK);
       double q = *(double *) params;
       
       return K*powerK(lK)*( j02(q*K) );
     }



double   CHI(double Q)
     {
         
       double result, error;
       double alpha=Q;

       double relerr=1.e-5;

       gsl_function F;
       F.params = &alpha;
       
       double b=max(exp(maxK),1.0/Q*100.);
       double a=exp(minK);
       const size_t limit=50000;       
       

       
       relerr=1.e-5;
       if (Q>0.01)relerr=1.e-7;
       F.function = &chiInt;
       gsl_integration_qag (&F, log(a),log(b), 0, relerr, limit,6,
                               w, &result, &error); 
                               
       fprintf(stderr,"CHI   at Q=%.5f: \t%.8f +- %.8f\n",Q,result,error);
       return result/(6.0*M_PI*M_PI);
     }


double sigmaInt (double lK, void * ) {
       double K=exp(lK);
       return K*powerK(lK);
     }
     
double   SIGMA(void)
     {
       double result, error;
       double alpha;     
       gsl_function F;
       F.function = &sigmaInt;
       F.params = &alpha;
       gsl_integration_qag (&F, min(minK,log(1.0/10000.)), max(maxK,log(10000.)), 0, 1e-7, 15000,6,
                             w, &result, &error); 
     
       fprintf(stderr,"SIGMA: \t%.8f +- %.8f\n",result,error);
       return result/(6.0*M_PI*M_PI);
     }





double xiLInt (double lK, void * params) {
       double K=exp(lK);
       double q = *(double *) params;
       
       return K*powerK(lK)*( j0(q*K) )*K*K;
     }



double   XIL(double Q)
     {
       if (Q >250.) return 0;  
       double result, error;
       double alpha=Q;

       double relerr=1.e-5;

       gsl_function F;
       F.params = &alpha;
       
       double b=exp(maxK);
       double a=exp(minK);
       const size_t limit=50000;       
       

       
       relerr=1.e-4;
       if (Q>0.01)relerr=1.e-5;
       if (Q>50)relerr=1.e-3;
       F.function = &xiLInt;
       gsl_integration_qag (&F, log(a),log(b), 1.e-5, relerr, limit,6,
                               w, &result, &error); 
                               
       fprintf(stderr,"Xi_linear   at X=%.5f: \t%.8f +- %.8f\n",Q,result,error);
       return result/(2.0*M_PI*M_PI);
     }






double xiLrs2Int (double lK, void * params) {
       double K=exp(lK);
       double q = *(double *) params;
       
       return K*powerK(lK)*( j2(q*K) )*K*K;
     }

double   XILrs2(double Q)
     {
       if (Q >250.) return 0;  
       double result, error;
       double alpha=Q;

       double relerr=2.e-3;

       gsl_function F;
       F.params = &alpha;
       
       double b=exp(maxK);
       double a=exp(minK);
       const size_t limit=500000;       
       

       

       F.function = &xiLrs2Int;
       gsl_integration_qag (&F, log(a),log(b), 1.e-5, relerr, limit,6,
                               w, &result, &error); 
                               
       fprintf(stderr,"Xi_L,Qaud   at X=%.5f: \t%.8f +- %.8f\n",Q,result,error);
       return -result/(2.0*M_PI*M_PI);
     }






double j4(double x){
    if(x > 1E-6) return (5.0*x*(-21.0 + 2.0*x*x)*cos(x) + (105.0 - 45.0*x*x + x*x*x*x)*sin(x))/pow(x,5);
    return pow(x,4)/945.0*(1.0-x*x/22.);
}



double xiLrs4Int (double lK, void * params) {
       double K=exp(lK);
       double q = *(double *) params;
       
       return K*powerK(lK)*( j4(q*K) )*K*K;
     }

double   XILrs4(double Q)
     {
       if (Q >250.) return 0;  
       double result, error;
       double alpha=Q;

       double relerr=2.e-3;

       gsl_function F;
       F.params = &alpha;
       
       double b=exp(maxK);
       double a=exp(minK);
       const size_t limit=500000;       
       

       

       F.function = &xiLrs4Int;
       gsl_integration_qag (&F, log(a),log(b), 1.e-5, relerr, limit,6,
                               w, &result, &error); 
                               
       fprintf(stderr,"Xi_L,Hexa   at X=%.5f: \t%.8f +- %.8f\n",Q,result,error);
       return result/(2.0*M_PI*M_PI);
     }







/*********************************************************************/


static int Integrand(const int *, const double x[], const int *, double f[], void *userdata){
        double *d =(double*)  userdata;
        double q =d[0];
        double a=d[1];

        
        double K;
        K=x[0]*(1.0/q-a)+a;
        double xx=K*q;

        
        f[0]= powerK(log(xx/q))*j2(xx);
        f[0]+= powerK(log(1.0/q/xx))*j2(1.0/xx)/xx/xx;
        f[0]*=(1.0/q-a);

         return 0;
     }

double gammaInt (double K, void * params) {
       double q = *(double *) params;
       double x=K*q;
       
       double res= powerK(log(x/q))*j2(x);
       res+= powerK(log(1.0/q/x))*j2(1.0/x)/x/x;
       return res;
     }


#define NDIM 1
#define NCOMP 1
//#define USERDATA NULL
#define EPSREL 1e-4
#define EPSABS 0.0
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 5000000

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL

#define NNEW 1000
#define FLATNESS 25.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

double  GAMMA(double Q)
{
  int verbose, fail; //comp, nregions,
  long long int neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
 double integ[1],err[1];
       double alpha=Q;
       


       double alphaN[2];
       
       double a=min(exp(minK),1.0/Q/1000.);
       
       alphaN[0]=Q;
       alphaN[1]=a;
       void * USERDATA = (void*) &alphaN;



       if (Q>0.1){

			double relerr=1.e-5;
            if (Q>1.0) relerr=1.e-7;
            gsl_function F;
			F.function = &gammaInt;
            F.params = &alpha;
			gsl_integration_qag (&F, a,1.0/Q, 0, relerr, 50000,6,
									w, integ, err); 
            integral[0]=integ[0];
            error[0]=err[0];
            fprintf(stderr,"GAMMA at Q=%.5f: \t%.8f +- %.8f\n",Q,integral[0],error[0]);
       }
       else {
             const char *env = getenv("CUBAVERBOSE");
             verbose = 0;
             if( env ) verbose = atoi(env);             
             llVegas(NDIM, NCOMP, Integrand, USERDATA,
                 EPSREL, EPSABS, verbose, SEED,
                 MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                 GRIDNO, STATEFILE,
                 &neval, &fail, integral, error, prob);
  
            //Suave(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose | LAST, SEED,
                //MINEVAL, MAXEVAL, NNEW, FLATNESS,
                //&nregions, &neval, &fail, integral, error, prob);
            
            //Divonne(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose, SEED,
                //MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
                //BORDER, MAXCHISQ, MINDEVIATION,
                //NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
                //&nregions, &neval, &fail, integral, error, prob);

            //Cuhre(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose | LAST,
                //MINEVAL, MAXEVAL, KEY,
                //&nregions, &neval, &fail, integral, error, prob);
            fprintf(stderr,"GAMMA at Q=%.5f:\t%.8f +- %.8f\tp = %.3f\n",Q,integral[0],error[0],prob[0]);
       }
       
       


  return integral[0]/(6.0*M_PI*M_PI);
}


void calcChiGamma(void)
{
	LINEAR_RS=0;
    CALCULATE_CHI_GAMMA=1;    
    initChiGamma();


}

void calcLinearRS(void)
{
	LINEAR_RS=1;
    CALCULATE_CHI_GAMMA=1;    
    initChiGamma();


}



/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////

//int j2*P*k^2 dk





static int IntegrandD2(const int *, const double x[], const int *, double f[], void *userdata){
        double *d =(double*)  userdata;
        double q =d[0];
        double a=d[1];

        
        double K;
        K=x[0]*(1.0/q-a)+a;
        double xx=K*q;

        
        f[0]= powerK(log(xx/q))*j2(xx)*(xx/q)*(xx/q);
        f[0]+= powerK(log(1.0/q/xx))*j2(1.0/xx)/xx/xx*(1.0/q/xx)*(1.0/q/xx);
        f[0]*=(1.0/q-a);

         return 0;
     }

double gslIntD2 (double K, void * params) {
       double q = *(double *) params;
       double x=K*q;
       
       double res= powerK(log(x/q))*j2(x)*(x/q)*(x/q);
       res+= powerK(log(1.0/q/x))*j2(1.0/x)/x/x*(1.0/q/x)*(1.0/q/x);
       return res;
     }



double  INTd2(double Q)
{
  int verbose, fail; //comp, nregions,
  long long int neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
 double integ[1],err[1];
       double alpha=Q;
       


       double alphaN[2];
       
       double a=min(exp(minK),1.0/Q/100000.);
       
       alphaN[0]=Q;
       alphaN[1]=a;
       void * USERDATA = (void*) &alphaN;



       if (Q>0.1){

			double relerr=1.e-5;
            //if (Q>1.0) relerr=1.e-7;
            gsl_function F;
			F.function = &gslIntD2;
            F.params = &alpha;
			gsl_integration_qag (&F, a,1.0/Q, 0, relerr, 50000,6,
									w, integ, err); 
            integral[0]=integ[0];
            error[0]=err[0];
            fprintf(stderr,"INTd2 at Q=%.5f: \t%.8f +- %.8f\n",Q,integral[0],error[0]);
       }
       else {
             const char *env = getenv("CUBAVERBOSE");
             verbose = 0;
             if( env ) verbose = atoi(env);             
             llVegas(NDIM, NCOMP, IntegrandD2, USERDATA,
                 EPSREL, EPSABS, verbose, SEED,
                 MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                 GRIDNO, STATEFILE,
                 &neval, &fail, integral, error, prob);
  
            //Suave(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose | LAST, SEED,
                //MINEVAL, MAXEVAL, NNEW, FLATNESS,
                //&nregions, &neval, &fail, integral, error, prob);
            
            //Divonne(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose, SEED,
                //MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
                //BORDER, MAXCHISQ, MINDEVIATION,
                //NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
                //&nregions, &neval, &fail, integral, error, prob);

            //Cuhre(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose | LAST,
                //MINEVAL, MAXEVAL, KEY,
                //&nregions, &neval, &fail, integral, error, prob);
            fprintf(stderr,"INTd2 at Q=%.5f:\t%.8f +- %.8f\tp = %.3f\n",Q,integral[0],error[0],prob[0]);
       }
       
       


  return integral[0]/(6.0*M_PI*M_PI);
}



/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////

//int j02*P*k^2 dk





static int IntegrandD02(const int *, const double x[], const int *, double f[], void *userdata){
        double *d =(double*)  userdata;
        double q =d[0];
        double a=d[1];

        
        double K;
        K=x[0]*(1.0/q-a)+a;
        double xx=K*q;

        
        f[0]= powerK(log(xx/q))*j02(xx)*(xx/q)*(xx/q);
        f[0]+= powerK(log(1.0/q/xx))*j02(1.0/xx)/xx/xx*(1.0/q/xx)*(1.0/q/xx);
        f[0]*=(1.0/q-a);

         return 0;
     }

double gslIntD02 (double K, void * params) {
       double q = *(double *) params;
       double x=K*q;
       
       double res= powerK(log(x/q))*j02(x)*(x/q)*(x/q);
       res+= powerK(log(1.0/q/x))*j02(1.0/x)/x/x*(1.0/q/x)*(1.0/q/x);
       return res;
     }



double  INTd02(double Q)
{
  int verbose, fail; //comp, nregions,
  long long int neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
 double integ[1],err[1];
       double alpha=Q;
       


       double alphaN[2];
       
       double a=min(exp(minK),1.0/Q/100000.);
       
       alphaN[0]=Q;
       alphaN[1]=a;
       void * USERDATA = (void*) &alphaN;



       if (Q>0.1){

			double relerr=1.e-5;
            //if (Q>1.0) relerr=1.e-7;
            gsl_function F;
			F.function = &gslIntD02;
            F.params = &alpha;
			gsl_integration_qag (&F, a,1.0/Q, 0, relerr, 50000,6,
									w, integ, err); 
            integral[0]=integ[0];
            error[0]=err[0];
            fprintf(stderr,"INTd02 at Q=%.5f: \t%.8f +- %.8f\n",Q,integral[0],error[0]);
       }
       else {
             const char *env = getenv("CUBAVERBOSE");
             verbose = 0;
             if( env ) verbose = atoi(env);             
             llVegas(NDIM, NCOMP, IntegrandD02, USERDATA,
                 EPSREL, EPSABS, verbose, SEED,
                 MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                 GRIDNO, STATEFILE,
                 &neval, &fail, integral, error, prob);
  
            //Suave(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose | LAST, SEED,
                //MINEVAL, MAXEVAL, NNEW, FLATNESS,
                //&nregions, &neval, &fail, integral, error, prob);
            
            //Divonne(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose, SEED,
                //MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
                //BORDER, MAXCHISQ, MINDEVIATION,
                //NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
                //&nregions, &neval, &fail, integral, error, prob);

            //Cuhre(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose | LAST,
                //MINEVAL, MAXEVAL, KEY,
                //&nregions, &neval, &fail, integral, error, prob);
            fprintf(stderr,"INTd02 at Q=%.5f:\t%.8f +- %.8f\tp = %.3f\n",Q,integral[0],error[0],prob[0]);
       }
       
       


  return integral[0]/(6.0*M_PI*M_PI);
}





/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////

//int j1*P*k^3 dk

double j1(double x){
    if(x > 1E-6) return -(cos(x)/x) + sin(x)/(x*x);
    return x/3. - x*x*x/30.;
}



static  int IntegrandD1(const  int *, const double x[], const  int *, double f[], void *userdata){
        double *d =(double*)  userdata;
        double q =d[0];
        double a=d[1];

        
        double K;
        K=x[0]*(1.0/q-a)+a;
        double xx=K*q;

        
        f[0]= powerK(log(xx/q))*j1(xx)*(xx/q)*(xx/q)*(xx/q);
        f[0]+= powerK(log(1.0/q/xx))*j1(1.0/xx)/xx/xx*(1.0/q/xx)*(1.0/q/xx)*(1.0/q/xx);
        f[0]*=(1.0/q-a);

         return 0;
     }

double gslIntD1 (double K, void * params) {
       double q = *(double *) params;
       double x=K*q;
       
       double res= powerK(log(x/q))*j1(x)*(x/q)*(x/q)*(x/q);
       res+= powerK(log(1.0/q/x))*j1(1.0/x)/x/x*(1.0/q/x)*(1.0/q/x)*(1.0/q/x);
       return res;
     }



double  INTd1(double Q)
{
  int verbose,   fail; //comp, nregions,
  long long int neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
 double integ[1],err[1];
       double alpha=Q;
       


       double alphaN[2];
       
       double a=min(exp(minK),1.0/Q/1000000.);
       
       alphaN[0]=Q;
       alphaN[1]=a;
       void * USERDATA = (void*) &alphaN;



       if (Q>0.1){

			double relerr=1.e-5;
            //if (Q>1.0) relerr=1.e-7;
            gsl_function F;
			F.function = &gslIntD1;
            F.params = &alpha;
			gsl_integration_qag (&F, a,1.0/Q, 0, relerr, 50000,6,
									w, integ, err); 
            integral[0]=integ[0];
            error[0]=err[0];
            fprintf(stderr,"INTd1 at Q=%.5f: \t%.8f +- %.8f\n",Q,integral[0],error[0]);
       }
       else {
             const char *env = getenv("CUBAVERBOSE");
             verbose = 0;
             if( env ) verbose = atoi(env);             
             llVegas(NDIM, NCOMP, IntegrandD1, USERDATA,
                 EPSREL, EPSABS, verbose, SEED,
                 MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                 GRIDNO, STATEFILE,
                 &neval, &fail, integral, error, prob);
  
            //Suave(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose | LAST, SEED,
                //MINEVAL, MAXEVAL, NNEW, FLATNESS,
                //&nregions, &neval, &fail, integral, error, prob);
            
            //Divonne(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose, SEED,
                //MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
                //BORDER, MAXCHISQ, MINDEVIATION,
                //NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
                //&nregions, &neval, &fail, integral, error, prob);

            //Cuhre(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose | LAST,
                //MINEVAL, MAXEVAL, KEY,
                //&nregions, &neval, &fail, integral, error, prob);
            fprintf(stderr,"INTd1 at Q=%.5f:\t%.8f +- %.8f\tp = %.3f\n",Q,integral[0],error[0],prob[0]);
       }
       
       


  return integral[0]/(6.0*M_PI*M_PI);
}



/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////
/////////////////////////////////////

//int j1*P*k dk

double j1k(double x){
    if(x > 1E-6) return -(cos(x)/x) + sin(x)/(x*x);
    return x/3. - x*x*x/30.;
}



static  int IntegrandD1k(const  int *, const double x[], const  int *, double f[], void *userdata){
        double *d =(double*)  userdata;
        double q =d[0];
        double a=d[1];

        
        double K;
        K=x[0]*(1.0/q-a)+a;
        double xx=K*q;

        
        f[0]= powerK(log(xx/q))*j1(xx)*(xx/q);
        f[0]+= powerK(log(1.0/q/xx))*j1(1.0/xx)/xx/xx*(1.0/q/xx);
        f[0]*=(1.0/q-a);

         return 0;
     }

double gslIntD1k (double K, void * params) {
       double q = *(double *) params;
       double x=K*q;
       
       double res= powerK(log(x/q))*j1(x)*(x/q);
       res+= powerK(log(1.0/q/x))*j1(1.0/x)/x/x*(1.0/q/x);
       return res;
     }



double  INTd1k(double Q)
{
  int verbose,   fail; //comp, nregions,
  long long int neval;
  double integral[NCOMP], error[NCOMP], prob[NCOMP];
 double integ[1],err[1];
       double alpha=Q;
       


       double alphaN[2];
       
       double a=min(exp(minK),1.0/Q/100000.);
       
       alphaN[0]=Q;
       alphaN[1]=a;
       void * USERDATA = (void*) &alphaN;



       if (Q>0.1){

			double relerr=1.e-5;
            //if (Q>1.0) relerr=1.e-7;
            gsl_function F;
			F.function = &gslIntD1k;
            F.params = &alpha;
			gsl_integration_qag (&F, a,1.0/Q, 0, relerr, 50000,6,
									w, integ, err); 
            integral[0]=integ[0];
            error[0]=err[0];
            fprintf(stderr,"INTd1k at Q=%.5f: \t%.8f +- %.8f\n",Q,integral[0],error[0]);
       }
       else {
             const char *env = getenv("CUBAVERBOSE");
             verbose = 0;
             if( env ) verbose = atoi(env);             
             llVegas(NDIM, NCOMP, IntegrandD1k, USERDATA,
                 EPSREL, EPSABS, verbose, SEED,
                 MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
                 GRIDNO, STATEFILE,
                 &neval, &fail, integral, error, prob);
  
            //Suave(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose | LAST, SEED,
                //MINEVAL, MAXEVAL, NNEW, FLATNESS,
                //&nregions, &neval, &fail, integral, error, prob);
            
            //Divonne(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose, SEED,
                //MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
                //BORDER, MAXCHISQ, MINDEVIATION,
                //NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
                //&nregions, &neval, &fail, integral, error, prob);

            //Cuhre(NDIM, NCOMP, Integrand, USERDATA,
                //EPSREL, EPSABS, verbose | LAST,
                //MINEVAL, MAXEVAL, KEY,
                //&nregions, &neval, &fail, integral, error, prob);
            fprintf(stderr,"INTd1k at Q=%.5f:\t%.8f +- %.8f\tp = %.3f\n",Q,integral[0],error[0],prob[0]);
       }
       
       


  return integral[0]/(6.0*M_PI*M_PI);
}
