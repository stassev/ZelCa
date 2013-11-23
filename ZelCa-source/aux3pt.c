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
/// Computes the real space matter 3-pt function in the Zel'dovich approximation.
/// 
///
//////////////////////////////////////////////////


#include <iostream>
#include <Eigen/LU>
using namespace std;
using namespace Eigen;
#include <stdio.h>
#include <cuba.h>
#include <cmath>
#include <Eigen/Cholesky> 
extern int REDSHIFT_SPACE;
extern string OUTPUT_FILE;
extern double REDSHIFT; 
extern double sigma;
void initChiGamma(void);
void printMatrix(double *A,int N);

typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6,1> Vector6d;
void createNij66(Matrix6d *Nij,bool *c,Matrix3d *psi1,Matrix3d *psi2,Matrix3d *psi12);
void createNij66init(double *q1,double *q2,Matrix3d *psi1,Matrix3d *psi2,Matrix3d *psi12);
void  ThreePt(double X[3],double *res,double *err,double *probability);
extern Matrix3d PsiDiag;

static Matrix6d mInit;

void calc3pt(void){
    double r1;
    double x[3],res,err,prob;
    REDSHIFT_SPACE=0;
    initChiGamma();
    
    //input triangle is set by the two sides x[0] and x[1] [Mpc/h]. The cos() of the angle between them is set in x[2]
    
    
    FILE *file= fopen(OUTPUT_FILE.c_str(), "w");
    
    x[2]=0.5; // angle=60deg
    
    int i;
    fprintf(file,"%.8e\n",REDSHIFT);
    for (i=0;i<=33;i++){
        r1=((double )i)/33.0*130.0+20.0;
        
        x[0]=r1;// set equilateral
        x[1]=r1;// set equilateral
        ThreePt(x,&res,&err,&prob);
        fprintf(file,"%.8e %.8e %.8e %.8e %.8e %.8e\n",x[0],x[1],x[2],res,err,prob);
        fflush(stdout);
        fflush(file);
    }
    fclose(file);

}










static  int Integrand3pt(const  int *, const double q[], const  int *, double f[], void *userdata){
        double *x = (double *) userdata;
        double X1abs=x[0];
        double X2abs=x[1];
        double mu12=x[2];
        

        double res;
        double Q1[3],Q2[3];
        double X1[3],X2[3];
        double sq12,sq;
        bool c[9];
        
//Declare 3pt:
        Matrix6d m012=mInit; // Initalize with mInit.
        FullPivLU<Matrix6d> lu012;
        
        


        
//Declare disconnected pieces:
//      2pt times a 1pt (3 symmetric):
        Matrix6d m01=mInit; // Initalize with mInit.
        FullPivLU<Matrix6d> lu01;

        Matrix6d m02=mInit; // Initalize with mInit.
        FullPivLU<Matrix6d> lu02;
        
        Matrix6d m12=mInit; // Initalize with mInit.
        FullPivLU<Matrix6d> lu12;
        
//      1pt times a 1pt times a 1pt:
        Matrix6d m0x1x2=mInit; // Initalize with mInit.
        FullPivLU<Matrix6d> lu0x1x2;
//end disconnected pieces

        double sum,Q1abs,muQ1,phiQ1,Q2abs,muQ2,phiQ2;
        double vv[6];
        double ex;

        Map<Vector6d> v(vv);
        Vector6d luV;
        
        //The Cuba library integrates over the unit hypercube. So, we need to rescale q[] to Q[] and include the jacobians.
        
        double jacMu;
        double sigmasqrt=3.0*sqrt(sigma);
        double q1=q[0];
        double q2=q[3];
        Q1abs= 2.0*sigmasqrt*atanh(q1);
        Q2abs= 2.0*sigmasqrt*atanh(q2);
        

        
//X1
        X1[0]=0;
        X1[1]=0;
        X1[2]=X1abs;
//X1
        sq12=sqrt(1.0-mu12*mu12);
        X2[0]=0;
        X2[1]=X2abs*sq12;
        X2[2]=X2abs*mu12;
//Q1
        muQ1=2.0*q[1]-1.0; 
        jacMu=2.0; 
        phiQ1=(q[2]*2.0-1.0)*M_PI;
        sq=sqrt(1.0-muQ1*muQ1);
        Q1[0]=Q1abs*sq*sin(phiQ1);   
        Q1[1]=Q1abs*sq*cos(phiQ1);   
        Q1[2]=Q1abs*muQ1;            

//Q2
        muQ2= 2.0*q[4]-1.0; 
        jacMu*=2.0;
        phiQ2=2.0*q[5]*M_PI;
        sq=sqrt(1.0-muQ2*muQ2);       
        Q2[0]=Q2abs*sq*sin(phiQ2); 
        Q2[1]=Q2abs*sq*cos(phiQ2); 
        Q2[2]=Q2abs*muQ2;          
    
        
//This contains v=X-Q, but we shift the Q's next. So, we end up with
        vv[0]=Q1[0];
        vv[1]=Q1[1];
        vv[2]=Q1[2];
        
        vv[3]=Q2[0];
        vv[4]=Q2[1];
        vv[5]=Q2[2];
        
//Center Q's around X's
        Q1[0]=X1[0]-Q1[0];
        Q1[1]=X1[1]-Q1[1];
        Q1[2]=X1[2]-Q1[2];
             
        Q2[0]=X2[0]-Q2[0];
        Q2[1]=X2[1]-Q2[1];
        Q2[2]=X2[2]-Q2[2];
        
////////////////////////////
////////////////////////////Initialize the psi's
////////////////////////////

        Matrix3d psi12;
        Matrix3d psi1 ;
        Matrix3d psi2 ;
        createNij66init(Q1,Q2,&psi1,&psi2,&psi12);

////////////////////////////
////////////////////////////
////////////////////////////
        //The c[] array tells createNij66 what are the connects and what are the disconnected pieces. 
        //Only 3 elements of c[] are used.
        //Do 3pt:
        c[0+3*1]=1;  c[0+3*2]=1;  
                     c[1+3*2]=1;                  
        createNij66(&m012,c,&psi1,&psi2,&psi12); 
        lu012.compute(m012);

        res=0;
        if (lu012.isInvertible()){
            ///////////////////////
            ///////////////////////
            //Do disconnected pieces:
                            //01
                            c[0+3*1]=1;  c[0+3*2]=0;
                                        c[1+3*2]=0; 
                            createNij66(&m01,c,&psi1,&psi2,&psi12); 
                            lu01.compute(m01);
                            //02      
                            c[0+3*1]=0;  c[0+3*2]=1;
                                        c[1+3*2]=0; 
                            createNij66(&m02,c,&psi1,&psi2,&psi12); 
                            lu02.compute(m02);
                            //12
                            c[0+3*1]=0;  c[0+3*2]=0;
                                        c[1+3*2]=1; 
                            createNij66(&m12,c,&psi1,&psi2,&psi12); 
                            lu12.compute(m12);
                        
                    //
                            //0x1x2
                            c[0+3*1]=0;  c[0+3*2]=0;
                                        c[1+3*2]=0;                
                            createNij66(&m0x1x2,c,&psi1,&psi2,&psi12); 
                            lu0x1x2.compute(m0x1x2);
            
                            
                        
            ///////////////////////
            ///////////////////////
            //3pt:                
                            ex=0;
                            
                            luV=lu012.solve(v);
                            sum = -0.5*v.transpose()*luV;
                            if (sum>-1500) ex+=exp(sum)/sqrt(lu012.determinant());
            //disconnected:
                            //01,2
                            luV=lu01.solve(v);
                            sum = -0.5*v.transpose()*luV;
                            if (sum>-1500) ex-=exp(sum)/sqrt(lu01.determinant());
                            //02,1
                            luV=lu02.solve(v);
                            sum = -0.5*v.transpose()*luV;
                            if (sum>-1500) ex-=exp(sum)/sqrt(lu02.determinant());
                            //0,12
                            luV=lu12.solve(v);
                            sum = -0.5*v.transpose()*luV;
                            if (sum>-1500) ex-=exp(sum)/sqrt(lu12.determinant());
                            //0,1,2
                            luV=lu0x1x2.solve(v);
                            sum = -0.5*v.transpose()*luV;
                            if (sum>-1500) ex+=2.0*exp(sum)/sqrt(lu0x1x2.determinant());
              //done              
                            
                            
                            double JACOBIAN=jacMu/(2.0*M_PI); //jacMu*(2pi)^2/(2pi)^(3)
                            JACOBIAN*=((2.0*sigmasqrt)/(1.0-q1*q1));
                            JACOBIAN*=((2.0*sigmasqrt)/(1.0-q2*q2));
                    
                            res=ex *  Q1abs*Q1abs*Q2abs*Q2abs  *  JACOBIAN;  
            
                            if (isnan(fabs(res)) || fabs(res)>1.e10) {
                                printf("Warning: Possible loss of precision or NaN! ");
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

#define NDIM 6
#define NCOMP 1
//#define USERDATA NULL
#define EPSREL 1.1e-2//1.e-7
#define EPSABS 1e-7//0.1//0.0
#define LAST 0
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 3000000000

#define NSTART 800000
#define NINCREASE 0
#define NBATCH 100000
#define GRIDNO 0
#define STATEFILE NULL

#define NNEW 10000
#define FLATNESS 4000.

//#define KEY1 512
//#define KEY1 1024
//#define KEY1 4096 //4096
#define KEY1 8096
//#define KEY1 512
//#define KEY1 41
#define KEY2 1
//#define KEY3 256
#define KEY3 3
//#define KEY3 1024
//#define KEY3 1


#define MAXPASS 50
#define BORDER 1.e-10
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 4096




void  ThreePt(double X[3],double *res,double *err,double *probability)
{
       int verbose, nregions, fail;
       long long int neval;
       double integral[NCOMP], error[NCOMP], prob[NCOMP];
       double x[]={X[0],X[1],X[2]};
       
       void * USERDATA = (void*) &x;
    
//////////////////
//////////////////
       mInit.block(0,0,3,3) = 2*PsiDiag;
       mInit.block(3,3,3,3) = 2*PsiDiag;
       mInit.block(0,3,3,3) =   PsiDiag;
       mInit.block(3,0,3,3) =   PsiDiag;
//////////////////
//////////////////

       const char *env = getenv("CUBAVERBOSE");
       verbose = 0;
       if( env ) verbose = atoi(env);     
    
       //llVegas(NDIM, NCOMP, Integrand3pt, USERDATA,
           //EPSREL, EPSABS, verbose, SEED,
           //MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
           //GRIDNO, STATEFILE,
           //&neval, &fail, integral, error, prob);
  
       //llSuave(NDIM, NCOMP, Integrand3pt, USERDATA,
           //EPSREL, EPSABS, verbose | LAST, SEED,
           //MINEVAL, MAXEVAL, NNEW, FLATNESS,NULL,
           //&nregions, &neval, &fail, integral, error, prob);
       
       llDivonne(NDIM, NCOMP, Integrand3pt, USERDATA,
           EPSREL, EPSABS,  verbose , SEED,
           MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
           BORDER, MAXCHISQ, MINDEVIATION,
           NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,NULL,
           &nregions, &neval, &fail, integral, error, prob);

       //llCuhre(NDIM, NCOMP, Integrand3pt, USERDATA,
           //EPSREL, EPSABS, verbose | LAST,
           //MINEVAL, MAXEVAL, KEY,NULL,
           //&nregions, &neval, &fail, integral, error, prob);

       
               
       fprintf(stderr,"ZETA at X=%.5e, %.5e, %.5e:\t%.8e +- %.3e\tp = %.3f\n",X[0],X[1],X[2],integral[0],error[0],prob[0]);

       res[0]=integral[0];
       err[0]= error[0];
       probability[0]=prob[0];
}





