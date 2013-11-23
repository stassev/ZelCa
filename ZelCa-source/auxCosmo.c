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


#include <math.h>
#include <gsl/gsl_sf_hyperg.h> 
#include <gsl/gsl_errno.h>

extern double Om;
double growthRate(double a);
double growthD(double a);
extern double GROWTH;

double growthDtemp(double a){
    // Decided to use the analytic expression for LCDM. More transparent if I change this to numerical integration?
    double x=-Om/(Om-1.0)/(a*a*a);
    double hyperP=0,hyperM=0;
    
    if (fabs(x-1.0) < 1.e-3) {
        hyperP=0.859596768064608 - 0.1016599912520404*(-1.0 + x) + 0.025791094277821357*pow(-1.0 + x,2) - 0.008194025861121475*pow(-1.0 + x,3) + 0.0029076305993447644*pow(-1.0 + x,4) - 0.0011025426387159761*pow(-1.0 + x,5) + 0.00043707304964624546*pow(-1.0 + x,6) - 0.0001788889964687831*pow(-1.0 + x,7);
        hyperM=1.1765206505266006 + 0.15846194123099624*(-1.0 + x) - 0.014200487494738975*pow(-1.0 + x,2) + 0.002801728034399257*pow(-1.0 + x,3) - 0.0007268267888593511*pow(-1.0 + x,4) + 0.00021801569226706922*pow(-1.0 + x,5) - 0.00007163321597397065*pow(-1.0 + x,6) +    0.000025063737576245116*pow(-1.0 + x,7);
    }
    else {
        if (x < 1.0) {
            hyperP=gsl_sf_hyperg_2F1(1.0/2.0,2.0/3.0,5.0/3.0,-x);
            hyperM=gsl_sf_hyperg_2F1(-1.0/2.0,2.0/3.0,5.0/3.0,-x);
        }
        x=1.0/x;
        if ((x < 1.0) && (x>1.0/30)) {
            
            hyperP=gsl_sf_hyperg_2F1(-1.0/6.0,0.5,5.0/6.0,-x);
            hyperP*=4*sqrt(x);
            hyperP+=-3.4494794123063873799*pow(x,2.0/3.0);
        
            hyperM=gsl_sf_hyperg_2F1(-7.0/6.0,-0.5,-1.0/6.0,-x);
            hyperM*=4.0/7.0/sqrt(x);
            hyperM+=pow(x,2.0/3.0)*(-1.4783483195598803057); //-(Gamma[-7/6]*Gamma[5/3])/(2*sqrt[Pi])
        }
        if (x<=1.0/30.0){
            hyperP=3.9999999999999996*sqrt(x) - 3.4494794123063865*pow(x,0.6666666666666666) + 0.3999999999999999*pow(x,1.5) -    0.13636363636363635*pow(x,2.5) + 0.07352941176470587*pow(x,3.5) - 0.04755434782608695*pow(x,4.5) +    0.033943965517241374*pow(x,5.5) - 0.02578125*pow(x,6.5) + 0.020436356707317072*pow(x,7.5) -    0.01671324384973404*pow(x,8.5) + 0.013997779702240564*pow(x,9.5) - 0.011945562847590041*pow(x,10.5) + 0.01035003662109375*pow(x,11.5) - 0.009080577904069926*pow(x,12.5);
            hyperM=0.5714285714285715/sqrt(x) + 2.000000000000001*sqrt(x) - 1.4783483195598794*pow(x,0.66666666666666666) +    0.10000000000000002*pow(x,1.5) - 0.022727272727272735*pow(x,2.5) + 0.009191176470588237*pow(x,3.5) -    0.004755434782608697*pow(x,4.5) + 0.002828663793103449*pow(x,5.5) - 0.0018415178571428578*pow(x,6.5) +    0.0012772722942073172*pow(x,7.5) - 0.0009285135472074472*pow(x,8.5) + 0.0006998889851120285*pow(x,9.5) -    0.0005429801294359111*pow(x,10.5) + 0.0004312515258789064*pow(x,11.5) - 0.00034925299631038194*pow(x,12.5);
        }
    }
    
    
    if (a>0.2) return sqrt(1.0 + (-1.0 + pow(a,-3))*Om)*(3.4494794123063873799*pow(-1.0 + 1.0/Om,0.666666666666666666666666666) + (hyperP*(4*pow(a,3)*(-1.0 + Om) - Om) - 7.0*pow(a,3)*hyperM*(-1.0 + Om))/(pow(a,5)*(-1.0+ Om) - pow(a,2)*Om));
    return (a*pow(1 - Om,1.5)*(1291467969*pow(a,12)*pow(-1 + Om,4) + 1956769650*pow(a,9)*pow(-1 + Om,3)*Om + 8000000000*pow(a,3)*(-1 + Om)*pow(Om,3) + 37490640625*pow(Om,4)))/(1.5625e10*pow(Om,5));
    
}






double decayD(double a){ // D_{-}, the decaying mode
    return sqrt(Om/(a*a*a)+1.0-Om);
}



double DprimeQ(double a,double nGrowth){ // returns Q*d(D_{+}^nGrowth*D_{-}^nDecay)/da, where Q=Qfactor(a)
        float nDecay=0.0;// not interested in decay modes in this code.
        float Nn=6.0*pow(1.0 - Om,1.5)/growthDtemp(1.0);
        return (pow(decayD(a),-1.0 + nDecay)*pow(growthD(a),-1.0 + nGrowth)*(nGrowth*Nn- (3.0*(nDecay + nGrowth)*Om*growthD(a))/(2.*a)));
        
    }
    

double Qfactor(double a){ // Q\equiv a^3 H(a)/H0.
    return sqrt(Om/(a*a*a)+1.0-Om)*a*a*a;
}


double growthD(double a){ // growth factor for LCDM
    return growthDtemp(a)/growthDtemp(1.0);
}

double growthRate(double a){
    return DprimeQ(a,1)/Qfactor(a)*a/growthD(a);
}

