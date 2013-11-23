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

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;
#include <string>

int REDSHIFT_SPACE;
int CALCULATE_CHI_GAMMA;
int MULTIPOLE;
double REDSHIFT;

double  Om; 
double  ns;
double  GROWTH;
double  RATE_OF_GROWTH;
double  RSfactor;
extern void  XI(double *X,double *res,double *err,double *probability);
extern void  XIRS(double *X,double *res,double *err,double *probability);


string OUTPUT_FILE;
string FILE_IN_PK;
string FILE_IN_OR_OUT_CHIGAMMAXI;
string FILE_OUT_LINEAR_THEORY;

double RSfactorOver;
void doChecks(void);
void doChecksRS(void);
void init(void);
double  growthD(double a);
double growthRate(double a);
void   initChiGamma(void);
void   calcChiGamma(void);
void   calcLinearRS(void);
void   calc3pt(void);
extern void calcXi(void);
extern void calcXiRS(void);
void  ThreePt(double X[3],double *res,double *err,double *probability);
extern void calcXiRS2d(void);
