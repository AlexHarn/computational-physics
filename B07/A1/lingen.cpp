#include <eigen3/Eigen/Core>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include "lingen.h"

using namespace std;
//-------------------------------------------------------------------------------
//                             CLASS LinConGen
//                        PRIVATE IMPLEMENTATION SECTION
//-------------------------------------------------------------------------------
void LinConGen::calccongruentgen(int64_t r0, int64_t a, int64_t c, int64_t m, int64_t M)
{   
    double rndnum = 0;
    //zwar r_0 und damit ein Glied mehr, aber nur N Glieder ausgeben
    for(int step = 0; step < M; step++) 
    {
        if(rnew == 0)
        {
                rold = r0;
        }
        else 
        {
                rold = rnew;
        }
        rnew = (a*rold+c)%m;  

        // für Floating-Point-Number [0,1[
        rndnum = ( (double) rnew/( (double) m )) ;

        //wird hoffentlich nie true 
        //debug: wird's auch bei den Parametern nicht, aber falls doch: lieber rausnehmen
        if(rndnum == 1)
        {
            cout << "Zufallszahl darf nicht 1 sein!" << endl;
            continue;
        }
        saveVec(step) = rndnum;
    }
}

//-------------------------------------------------------------------------------
//                             CLASS LinConGen
//                       PUBLIC IMPLEMENTATION SECTION
//-------------------------------------------------------------------------------

//Init
LinConGen::LinConGen(int savedvalues, int bins, double lengthbins)
{
    //hier stehen dann evtl noch ein paar Restnullen
    //durch den Cut bei den if-Abfragen, aber das kostet nicht so viel Aufwand
    saveVec = V::Zero(savedvalues);
    //hist = V::Zero(bins);
    rold = 0;
    rnew = 0;
    r0 = 0;
    a = 0;
    c = 0;
    m = 0;
    M = 0;
    N = savedvalues;
    totalbins = bins;
    lengthpbin = lengthbins;
}

/*void LinConGen::sethistogram(int totalbins, double lengthpbin) 
{
    VectorXd vec(1,5,8);
    for (int i = 0; i < totalbins; i++)
    {
        hist(vec(i))++;
    }
}*/

void LinConGen::calc(int64_t r0,int64_t a, int64_t c, int64_t m)
{
    // Max bei 10^4 oder N<m-1
    if(N > (m-1) )
    {
        M = (m-1) ;
    }
    else
    {
        M = N;
    }
    calccongruentgen(r0,a,c,m,M);
    //sethistogram(totalbins, lengthpbin);
}

void LinConGen::save(string name)
{
    ofstream fout;
    fout.open(name+".dat");
    for ( int step = 0; step < M; step++ )
    {
        fout << saveVec(step) << endl;
    }
    fout.close();
}

void LinConGen::reset()
{
    saveVec.setZero();
    //hist.setZero();
    rnew = 0;
    rold = 0;
    a = 0;
    c = 0;
    m = 0;
    M = 0;
}
