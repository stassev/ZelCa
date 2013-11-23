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

#include "stuff.h"
#include <iostream>
using namespace std;
#include <string>

int main()
{

   
    cout << "\nEnter Omega_M (CDM+baryons) for z=0 (e.g. 0.274): ";
    cin >> Om;
    
    cout << "\nEnter scalar spectral index, n_s (e.g. 0.95): ";
    cin >> ns;
    
    cout << "\nEnter file to input matter power spectrum (provided example: camb.dat): ";
    cin >> FILE_IN_PK;
    
    int CalculateChiGamma;
    cout << "\nDo you need to precompute chi and gamma (it will be slow)?";
    cout << "\nOnly needed for different cosmologies, not for different";
    cout << "\nredshift. (yes=1; no=0): ";
    cin >> CalculateChiGamma;
    if (CalculateChiGamma==1){
        cout << "\nEnter file to output chi, gamma: ";
        cin >> FILE_IN_OR_OUT_CHIGAMMAXI;
        calcChiGamma(); //~4min
        
        // It seems I have to exit after each calculation. Otherwise weird things happen.
        // The problem must be either in calling initChiGamma() multiple times; or more  
        // likely because the Cuba forks of the different calculations seem to cross-talk
        // (which is probably due to some shared variable, but who knows...).
        // This has nothing to do with the ZA pieces of the code, which is great!
        return 0;
    }
    else   cout << "\nEnter file to input chi, gamma (provided example: chigammaxi.dat): ";
    cin >> FILE_IN_OR_OUT_CHIGAMMAXI;
    
    

    
    
    int CalculateLinearTheory;
    cout << "\nDo you want to calculate the linear theory quantities (it will be slow)? (yes=1; no=0): ";
    cin >> CalculateLinearTheory;
    if (CalculateLinearTheory==1){
        cout << "\nEnter file to output linear theory results: ";
        cin >> FILE_OUT_LINEAR_THEORY;
        
        calcLinearRS(); //~4min
        
        return 0;
    }

    int space;

    cout << "\nDo you want to do calculation for ZA 2pt in real or redshift space? (real=0; redshift space (multipoles)=1; redshift space (2-dim)=2; skip this=-1): ";
    cin >> space;
    if (space==0){
        cout << "\nEnter file to output real space 2-pt function: ";
        cin >> OUTPUT_FILE;
        cout << "\nEnter redshift (e.g. 0.55): ";
        cin >> REDSHIFT;
        cout << "\n";
        init();
        
        calcXi();            //~10sec
        
        return 0;
    }
    if (space==1) {
        cout << "\nWhich multipole do you want? (0, 2 or 4): ";
        cin >> MULTIPOLE;
        cout << "\nEnter file to output redshift space 2-pt function: ";
        cin >> OUTPUT_FILE;
        cout << "\nEnter redshift (e.g. 0.55): ";
        cin >> REDSHIFT;
        cout << "\n";
        init();
        
        calcXiRS();      //l=0: ~2min; l=2: ~4min; l=4: ~10min
        
        return 0;
    }
    if (space==2){
        cout << "\nEnter file to output 2d redshift space 2-pt function: ";
        cin >> OUTPUT_FILE;
        cout << "\nEnter redshift (e.g. 0.35): ";
        cin >> REDSHIFT;
        cout << "\n";
        init();
        
        calcXiRS2d();            //~several minutes
        
        return 0;
    }
    
    int CalculateZetaZA;
    cout << "\nDo you want to do calculation for ZA 3pt in real space? (Yes=1; No=0): ";
    cin >> CalculateZetaZA;
    if (CalculateZetaZA==1){
        cout << "\nEnter file to output real space 3-pt function: ";
        cin >> OUTPUT_FILE;
        cout << "\nEnter redshift (e.g. 0.55): ";
        cin >> REDSHIFT;
        cout << "\n";
        init();
        
        calc3pt();           // ~2-20min/point
        
        return 0;
    }
    
    
    return 0;
}






void init(void){
    GROWTH=growthD(1.0/(1.0+REDSHIFT)) ;
    RATE_OF_GROWTH=growthRate(1.0/(1.0+REDSHIFT)) ;
    RSfactor=(1.0+RATE_OF_GROWTH)*(1.0+RATE_OF_GROWTH);
    RSfactorOver=1.0+RATE_OF_GROWTH;

    CALCULATE_CHI_GAMMA=0;   // this will be overwritten if calcChiGamma() is called.

}


void printMatrix(double *A,int N){
     int i,j;
     printf("------------\n");
     for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            printf("%9.4g ",A[j+N*i]);
        }
        printf("\n");
     }
     printf("------------\n");
}


