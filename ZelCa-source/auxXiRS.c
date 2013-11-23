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
/// Computes the redshift space matter 2-pt function in the Zel'dovich approximation.
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
extern string OUTPUT_FILE;
extern double REDSHIFT; 
extern double GROWTH;
extern double RATE_OF_GROWTH;
extern double sigma;
extern double RSfactor;
extern int MULTIPOLE;

extern double RSfactorOver;

void initChiGamma(void);
void createNij33(double *q,Matrix3d *Nij);
void printMatrix(double *A,int N);
void XIRS(double *X,double *res,double *err,double *probability);
void calcXiRS(void);
Vector3d diagRS;

static double det0;
static double eps;
static double qmin;
static double sigmasqrt;
static double tanhcoef;



void calcXiRS(void)
{
    REDSHIFT_SPACE=1;
    initChiGamma();
    
    double x,res,err,prob;
    int i;
    FILE *file= fopen(OUTPUT_FILE.c_str(), "w");
    
    double minX=5.;
    double maxX=150.;

    int n=((int)((maxX-minX)+0.5));

    fprintf(file,"%i\n",n);
    fprintf(file,"%.8e %.8e %.8e %.8e %.8e\n",minX,maxX,REDSHIFT,GROWTH,RATE_OF_GROWTH);
    fprintf(file,"%i %i\n",REDSHIFT_SPACE,MULTIPOLE);
    for (i=0;i<=n;i++){
        x= ((double) i)/((double)n)*(maxX-minX)+minX;
        XIRS(&x,&res,&err,&prob);
        fprintf(file,"%.8e %.8e %.8e %.8e\n",x,res,err,prob);
    }
    fclose(file);
    
}




static int Integrand2ptRS(const int *, const double q[], const int *, double f[], void *userdata){
        double *x = (double *) userdata;
        double Xabs=x[0];
        double X[3];
        double Q[3];
        double det,res;

        double sum,sum0,muX,Qabs,muQ,phiQ;
        double vv[3];

        Matrix3d m;
        Vector3d luV;
        Matrix3d mInv,mInv0;
        double ex;
        double ex0;
        bool invertible;
        double A,B;
        double jacMuQ,jacMuX,jacPhi,jacQ,q0;
        double Qy,Qz,sqQ,sqX;
        double JACOBIAN;
        Map<Vector3d> v(vv);

        muX=q[1]*2.0-1.0;
        jacMuX=2.0; // d muX/d q[1]
        sqX=sqrt(1.0-muX*muX);
        
        X[0]=0.;
        X[1]=Xabs*sqX;
        X[2]=Xabs*muX;
        
        //The Cuba library integrates over the unit hypercube. So, we need to rescale q[] to Q[] and include the jacobians.
        
        //Set the magnitude of Q to sample around |X|=Xabs with a width of sigmasqrt. Do that with a atanh.
        q0=q[0]*(1.0-2.0*qmin)+qmin;
        Qabs=Xabs + 2.0*sigmasqrt*atanh(q0 + (q0-1.0)*tanhcoef);
        jacQ=((2.0*sigmasqrt)/(1.0-q0*q0 - (q0-1.0)*(q0-1.0)*tanhcoef)); // = d Qabs / d q0
        
        phiQ=(q[3]*2.0-1.0)*M_PI; 
        jacPhi=2*M_PI;  // = d phiQ/ d q[3]
        muQ=1.0 - 2.0*q[2]*q[2]*q[2]; 
        jacMuQ=6.0*q[2]*q[2]; // = d \mu/ d q[2]
        
        sqQ=sqrt(1.0-muQ*muQ);  
        
        Q[0]=Qabs*sqQ*sin(phiQ);
        Qy=Qabs*sqQ*cos(phiQ);
        Qz=Qabs*muQ;
        
        //Now rotate vector {Q[0],Qy,Qz} such that z axis points along X. This speeds up things ...
        //Don't remember if I tested this with a shift of Q->Q+X. Might be faster as for 3pt.
        Q[1]=Qy*muX+sqX*Qz;
        Q[2]=Qz*muX-sqX*Qy; //note that jacobian of rot = 1
        
        vv[0]=Q[0]-X[0];
        vv[1]=Q[1]-X[1];
        vv[2]=Q[2]-X[2];
        
        createNij33(Q,&m);

        m.computeInverseAndDetWithCheck(mInv,det,invertible);
        det=sqrt(det);
        B=det0*det;
        
        
        
        res=0;
        if (invertible && (fabs(B) > eps)){
        
                //2pt:
                sum=-0.5*v.transpose()*mInv*v;
                
                //disconnected piece
                mInv0=0.5/sigma*diagRS.asDiagonal();
                sum0=-0.5*v.transpose()*mInv0*v;
                
                sum-=sum0;
                
                if (sum<-300)ex=0.0;
                else ex=exp(sum);
                if (sum0<-300)ex0=0.0;
                else ex0=exp(sum0);
                
                JACOBIAN=jacMuQ*jacPhi*jacQ*jacMuX ; 
                A=(ex*det0-det)*ex0 *  Qabs*Qabs  *  JACOBIAN;
                res  =  A/B;
                
                if (isnan(res) || fabs(res)>1.e6) {
                    printf("Warning! Possible loss of precision or NaN! %4g %4g %4g %4g %4g %4g %4g %4g\n",A,B,res,Qabs,ex0,ex,sum0,sum);
                    printMatrix(&mInv(0),3);//This may produce weird output if running on multiple threads.
                }
        }
        else{
            printf("Matrix non-invertible or too small denominator! *********************************************************************\n");
            printf("%4e %4e %4e %4e %4e %4e\n",q[0],q[1],q[2],q[3],q[4],q[5]);
        }

        f[0]=res;

        if (MULTIPOLE==0) f[0]*=0.5;
        if (MULTIPOLE==2) f[0]*=5.0/2.0*(-1.0 + 3.0*muX*muX)/2.;
        if (MULTIPOLE==4) f[0]*=9.0/2.0*((3.0 - 30.0*muX*muX + 35.0*muX*muX*muX*muX)/8.);


        return 0;
     }

//*********************************************************************/
//*********************************************************************/
//*********************************************************************/
//**********************SOME STANDARD CUBA DEFINITIONS*****************/
//*********************************************************************/
//*********************************************************************/
//*********************************************************************/


#define NDIM 4
#define NCOMP 1
//#define USERDATA NULL
#define EPSREL 9e-4
#define EPSABS 2e-5
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 1000000000

#define NSTART 50000
#define NINCREASE 5000
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL

#define NNEW 10000
#define FLATNESS 125.

//#define KEY1 512
//#define KEY1 1024
#define KEY1 4096
//#define KEY1 512
//#define KEY1 41
#define KEY2 1
//#define KEY3 256
#define KEY3 512
//#define KEY3 1024
//#define KEY3 1


#define MAXPASS 50
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0




void  XIRS(double *X,double *res,double *err,double *probability)
{
       int verbose, nregions,  fail;
       long long int neval;
       double integral[NCOMP], error[NCOMP], prob[NCOMP];
       double XGIVEN[NGIVEN*LDXGIVEN];
       double x=X[0];
       void * USERDATA = (void*) &x;
       sigmasqrt=sqrt(2.*sigma);
       tanhcoef=tanh(x/(2.*sigmasqrt));
       diagRS << 1,1,1.0/RSfactor;
       det0=sqrt(8.*sigma*sigma*sigma*RSfactor); 
       eps=1.e-13;
       qmin=1.e-10;


       const char *env = getenv("CUBAVERBOSE");
       verbose = 0;
       if( env ) verbose = atoi(env);     
    
       //llVegas(NDIM, NCOMP, Integrand2ptRS, USERDATA,
       //   EPSREL, EPSABS, verbose, SEED,
       //   MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
       //   GRIDNO, STATEFILE,
       //   &neval, &fail, integral, error, prob);
  
       //llSuave(NDIM, NCOMP, Integrand2ptRS, USERDATA,
           //EPSREL, EPSABS, verbose | LAST, SEED,
           //MINEVAL, MAXEVAL, NNEW, FLATNESS,
           //&nregions, &neval, &fail, integral, error, prob);
       
       llDivonne(NDIM, NCOMP, Integrand2ptRS, USERDATA,
           EPSREL, EPSABS, verbose, SEED,
           MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
           BORDER, MAXCHISQ, MINDEVIATION,
           NGIVEN, LDXGIVEN, XGIVEN, NEXTRA, NULL,NULL,
           &nregions, &neval, &fail, integral, error, prob);

       //llCuhre(NDIM, NCOMP, Integrand2ptRS, USERDATA,
           //EPSREL, EPSABS, verbose | LAST,
           //MINEVAL, MAXEVAL, KEY,
           //&nregions, &neval, &fail, integral, error, prob);

       integral[0]/=pow(2.0*M_PI,1.5);
       error[0]/= pow(2.0*M_PI,1.5);
               
       fprintf(stderr,"XI_RS l=%i at X=%.5e:\t%.8e +- %.3e\tp = %.3f\n",MULTIPOLE,X[0],integral[0],error[0],prob[0]);

       res[0]=integral[0];
       err[0]= error[0];
       probability[0]=prob[0];
}



