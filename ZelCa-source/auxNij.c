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


#include <iostream>
#include <Eigen/LU>
using namespace std;
using namespace Eigen;
#include <stdio.h>
#include <cuba.h>

extern double RSfactorOver;
extern double RSfactor;
extern int REDSHIFT_SPACE;
extern int CALCULATE_CHI_GAMMA;
extern double sigma;
void   ChiGamma(double Q,double*chi,double*gamma);
void createPsi(double *q,Matrix3d *psi);

extern Matrix3d PsiDiag;
typedef Matrix<double, 6, 6> Matrix6d;

void initChiGamma(void);
void printMatrix(double *A,int N);



/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/




void createNij66init(double *q1,double *q2,Matrix3d *psi1,Matrix3d *psi2,Matrix3d *psi12){


    double q12[3];

    q12[0]=q2[0]-q1[0];
    q12[1]=q2[1]-q1[1];
    q12[2]=q2[2]-q1[2];

    createPsi(q12, psi12);
    createPsi(q1,  psi1);
    createPsi(q2,  psi2);

     
} 



void createNij66(Matrix6d *Nij,bool *c,Matrix3d *psi1,Matrix3d *psi2,Matrix3d *psi12){
//11
    if (c[0+3*1]){
        Nij[0].block(0,0,3,3)-= 2*psi1[0];
        Nij[0].block(0,3,3,3)-=   psi1[0];
        Nij[0].block(3,0,3,3)-=   psi1[0];
    }
//22
    if (c[0+3*2]){
        Nij[0].block(3,3,3,3)-= 2*psi2[0];
        Nij[0].block(0,3,3,3)-=   psi2[0];
        Nij[0].block(3,0,3,3)-=   psi2[0];
    }
//12
    if (c[1+3*2]) {
        Nij[0].block(0,3,3,3)+= psi12[0];
        Nij[0].block(3,0,3,3)+= psi12[0];
    }
} 
 
 


void createNij33(double *q,Matrix3d *Nij){

    Matrix3d psi;
    
    createPsi(q,&psi);
    Nij[0]  = PsiDiag;
    Nij[0] -=  psi;
    Nij[0] *=  2;
}



void createPsi(double *q, Matrix3d *psi){
    double Q2=q[0]*q[0]+q[1]*q[1]+q[2]*q[2];
    double Q=sqrt(Q2);
    double chi,gamma;
    
    ChiGamma(Q,&chi,&gamma);
    Vector3d diagChi( chi,chi,chi);


    psi[0]=diagChi.asDiagonal();
    if (Q >1.e-8) {
        Map<Vector3d> v(q);
        psi[0]-=3.0*gamma/Q2*v*v.transpose();
    }
    if (REDSHIFT_SPACE) {
        Vector3d diagRS( 1,1,RSfactorOver);
        psi[0]=diagRS.asDiagonal()*psi[0]*diagRS.asDiagonal();
    }
    
}




/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
