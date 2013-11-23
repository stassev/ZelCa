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
/// Computes the real space matter 2-pt function in the Zel'dovich approximation.
/// 
///
//////////////////////////////////////////////////

#include <iostream>
#include <Eigen/LU>
using namespace std;
using namespace Eigen;
#include <stdio.h>
#include <cuba.h>
extern int REDSHIFT_SPACE;
extern string  OUTPUT_FILE;
extern double REDSHIFT;
extern double GROWTH;
extern double RATE_OF_GROWTH;
extern double sigma;


void initChiGamma(void);
void createNij33(double *q,Matrix3d *Nij);
void printMatrix(double *A,int N);
void XI(double *X,double *res,double *err,double *probability);
 

static double det0;
static double qmin,eps;
static double sigmasqrt;
static double tanhcoef;

void calcXi(void)
{
    REDSHIFT_SPACE=0;
    initChiGamma();
    
    double x,res,err,prob;
    int i;
    FILE *file= fopen(OUTPUT_FILE.c_str(), "w");
    
    double maxX=150;
    double minX=5;
    
    int n=((int)((maxX-minX)+0.5));

    fprintf(file,"%i\n",n);
    fprintf(file,"%.8e %.8e %.8e %.8e %.8e\n",minX,maxX,REDSHIFT,GROWTH,RATE_OF_GROWTH);
    for (i=0;i<=n;i++){
        x= ((double) i)/((double)n)*(maxX-minX)+minX;
        XI(&x,&res,&err,&prob);
        fprintf(file,"%.8e %.8e %.8e %.8e\n",x,res,err,prob);
    }
    fclose(file);
    
}




static int Integrand2pt(const int *, const double q[], const int *, double f[], void *userdata){
        double *x = (double *) userdata;
        double Q[3];
        double det,res;
        double sum,sum0,mu,Qabs;
        Matrix3d m;
        Vector3d luV;
        Matrix3d mInv,mInv0;

        double vv[3];
        double ex;
        double ex0;
        bool invertible;
        double A,B;
        Vector3d diag(1,1,1);
        double q0,jacQ,jacMu;
        
        //point vector X along \hat x
        double X[]={x[0],0,0};
        
        //The Cuba library integrates over the unit hypercube. So, we need to rescale q[] to Q[] and include the jacobians.
        
        //Set the magnitude of Q to sample around |X|=x[0] with a width of sigmasqrt. Do that with a atanh.
        q0=q[0]*(1.0-2.0*qmin)+qmin;
        Qabs=x[0] + 2.0*sigmasqrt*atanh(q0 + (q0-1.0)*tanhcoef);
        jacQ=((2.0*sigmasqrt)/(1.0-q0*q0 - (q0-1.0)*(q0-1.0)*tanhcoef)); // = d Qabs / d q0
        
        //phi integral is trivial, so it was already done. 
        //In principle one can then avoid using 3x3 matrices and 3-vectors and simplify.
        //But the code will no longer be as transparent.
        jacMu=2.0; // d \mu/ d q[1] =2
        mu=q[1]*2.0-1.0;
        Q[0]=Qabs*mu;
        Q[1]=Qabs*sqrt(1.0-mu*mu);
        Q[2]=0.0; 
        
        vv[0]=Q[0]-X[0];
        vv[1]=Q[1]-X[1];
        vv[2]=Q[2]-X[2];
        
        Map<Vector3d> v(vv);
        createNij33(Q,&m);
        
        m.computeInverseAndDetWithCheck(mInv,det,invertible);
        
        res=0;
        if (invertible){
            //2pt:
            sum=-0.5*v.transpose()*mInv*v;
            
            //disconnected piece:
            mInv0=0.5/sigma*diag.asDiagonal();
            sum0=-0.5*v.transpose()*mInv0*v;
            //end 
            
            sum-=sum0;
            
            if (sum<-300)ex=0.0;
            else ex=exp(sum);
                    
            if (sum0<-300)ex0=0.0;
            else ex0=exp(sum0);
            
            
            det=sqrt(det);
            A=(ex*det0-det)*ex0 *  Qabs*Qabs  *  jacMu * jacQ;
            B=det0*det;
            res=0;
            if ((fabs(B)*eps<fabs(A)) && (fabs(B) > eps) )
                res  =  A/B;
            
            if (isnan(res) || fabs(res)>1.e6) {
                printf("Warning! Possible loss of precision or NaN! %4g %4g %4g %4g %4g %4g %4g %4g\n",A,B,res,Qabs,ex0,ex,sum0,sum);
                printMatrix(&mInv(0),3);//This may produce weird output if running on multiple threads.
            }
    
        }
        else{
            printf("Matrix non-invertible! *********************************************************************\n");
            printf("%4e %4e %4e %4e %4e %4e\n",q[0],q[1],q[2],q[3],q[4],q[5]);
        }


        f[0]=res;

        return 0;
     }



//*********************************************************************/
//*********************************************************************/
//*********************************************************************/
//**********************SOME STANDARD CUBA DEFINITIONS*****************/
//*********************************************************************/
//*********************************************************************/
//*********************************************************************/


#define NDIM 2
#define NCOMP 1
#define EPSREL 1e-3
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




void  XI(double *X,double *res,double *err,double *probability)
{
       int verbose, nregions,  fail;
       long long int neval;
       double integral[NCOMP], error[NCOMP], prob[NCOMP];
       double x=X[0];
       void * USERDATA = (void*) &x;
       sigmasqrt=sqrt(2.*sigma);
       tanhcoef=tanh(x/(2.*sigmasqrt));
       det0=sqrt(8.*sigma*sigma*sigma);
       eps=1.e-13;
       qmin=1.e-10;
        
       const char *env = getenv("CUBAVERBOSE");
       verbose = 0;
       if( env ) verbose = atoi(env);             
       //Vegas(NDIM, NCOMP, Integrand2pt, USERDATA,
       //    EPSREL, EPSABS, verbose, SEED,
       //    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
       //    GRIDNO, STATEFILE,
       //    &neval, &fail, integral, error, prob);
  
       //Suave(NDIM, NCOMP, Integrand2pt, USERDATA,
           //EPSREL, EPSABS, verbose | LAST, SEED,
           //MINEVAL, MAXEVAL, NNEW, FLATNESS,
           //&nregions, &neval, &fail, integral, error, prob);
       
       //Divonne(NDIM, NCOMP, Integrand2pt, USERDATA,
           //EPSREL, EPSABS, verbose, SEED,
           //MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
           //BORDER, MAXCHISQ, MINDEVIATION,
           //NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
           //&nregions, &neval, &fail, integral, error, prob);

       llCuhre(NDIM, NCOMP, Integrand2pt, USERDATA,
           EPSREL, EPSABS, verbose | LAST,
           MINEVAL, MAXEVAL, KEY,NULL,
           &nregions, &neval, &fail, integral, error, prob);

       integral[0]/=sqrt(2.0*M_PI);
       error[0]/= sqrt(2.0*M_PI);
       
       fprintf(stderr,"XI at X=%.5e:\t%.8e +- %.3e\tp = %.3f\n",X[0],integral[0],error[0],prob[0]);

       res[0]=integral[0];
       err[0]= error[0];
       probability[0]=prob[0];
}






