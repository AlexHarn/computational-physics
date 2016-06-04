#include <eigen3/Eigen/Core>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
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
LinConGen::LinConGen(int savedvalues, int bins)
{
    //hier stehen dann evtl noch ein paar Restnullen
    //durch den Cut bei den if-Abfragen, aber das kostet nicht so viel Aufwand
    saveVec = V::Zero(savedvalues);
    hist = V::Zero(bins);
    x = V::Zero(2);
    y = V::Zero(2);
    rold = 0;
    rnew = 0;
    r0 = 0;
    a = 0;
    c = 0;
    m = 0;
    M = 0;
    N = savedvalues;
}

void LinConGen::sethistogram() 
{
    //Eventuell geht das effektiver. Läuft aber!
    for (int i = 0; i < M; i++)
    {
        if(saveVec(i) < 0.1 && saveVec(i) >= 0)
            hist(0)++;
        else if(saveVec(i) < 0.2 && saveVec(i) >= 0.1)
            hist(1)++;
        else if(saveVec(i) < 0.3 && saveVec(i) >= 0.2)
            hist(2)++;
        else if(saveVec(i) < 0.4 && saveVec(i) >= 0.3)
            hist(3)++;
        else if(saveVec(i) < 0.5 && saveVec(i) >= 0.4)
            hist(4)++;
        else if(saveVec(i) < 0.6 && saveVec(i) >= 0.5)
            hist(5)++;
        else if(saveVec(i) < 0.7 && saveVec(i) >= 0.6)
            hist(6)++;
        else if(saveVec(i) < 0.8 && saveVec(i) >= 0.7)
            hist(7)++;
        else if(saveVec(i) < 0.9 && saveVec(i) >= 0.8)
            hist(8)++;
        else if(saveVec(i) < 1.0 && saveVec(i) >= 0.9)
            hist(9)++;
        else
            cout << "Hier ist was schief gelaufen!" << endl;
    }
}

void LinConGen::boxmulleralg()
{   
    /* Zufälliger Seed: */
    srand (time(NULL));

    // Nehme aus gleichverteilten Zahlen aus Generator irgendeine Zahl im savevec
    x(0) = saveVec( rand() % M );
    x(1) = saveVec( rand() % M );

    // Nach 10.11 , S. 176
    y(0) = sqrt(-2*log(x(0)))*cos(2*M_PI*x(1));
    y(1) = sqrt(-2*log(x(0)))*sin(2*M_PI*x(1));
}

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
    sethistogram();
}

void LinConGen::save(string name)
{
    ofstream fout;
    fout.open(name+".dat");
    fout << "#RandomNumber" << endl;
    for ( int step = 0; step < M; step++ )
    {
        fout << saveVec(step) << endl;
    }
    fout.close();
    fout.open(name+"_histogram.dat");
    fout << "#Histogrammshäufigkeiten" << endl;
    for ( int step = 0; step < hist.size(); step++ )
    {
        fout << hist(step) << endl;
    }
    fout.close();
    fout.open(name+"_pair.dat");
    fout << "#x" << "\t" << "y" << endl;
    for(int i = 0; i < M; i+=2)
    {
        fout << saveVec(i) << "\t" << saveVec(i+1) << endl;              
    }
}

void LinConGen::reset()
{
    saveVec.setZero();
    hist.setZero();
    x.setZero();
    y.setZero();
    rnew = 0;
    rold = 0;
    a = 0;
    c = 0;
    m = 0;
    M = 0;
}
