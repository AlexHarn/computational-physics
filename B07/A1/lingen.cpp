#include <eigen3/Eigen/Core>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "lingen.h"

using namespace std;
void LinConGen::calcCongruent(int64_t r0, int64_t a, int64_t c, int64_t m, int64_t N)
{
    if ( N > m )
        throw invalid_argument("N > m!");
    V64 r(N);
    this->r = V::Zero(N);
    r(0) = r0;

    for ( int n = 0; n < N; n++ )
        r(n+1) = ( a*r(n) + c )%m;

    this->r = r.cast<double>()/m;
}

//-------------------------------------------------------------------------------
//                             CLASS LinConGen
//                       PUBLIC IMPLEMENTATION SECTION
//-------------------------------------------------------------------------------

LinConGen::LinConGen(int bins)
{
    //hier stehen dann evtl noch ein paar Restnullen
    //durch den Cut bei den if-Abfragen, aber das kostet nicht so viel Aufwand
    x = V::Zero(2);
    y = V::Zero(2);
    this->bins = bins;
}

//void LinConGen::boxmulleralg()
//{
    //[> Zufälliger Seed: <]
    //srand (time(NULL));

    //// Nehme aus gleichverteilten Zahlen aus Generator irgendeine Zahl im savevec
    //x(0) = saveVec( rand() % M );
    //x(1) = saveVec( rand() % M );

    //// Nach 10.11 , S. 176
    //y(0) = sqrt(-2*log(x(0)))*cos(2*M_PI*x(1));
    //y(1) = sqrt(-2*log(x(0)))*sin(2*M_PI*x(1));
//}

void LinConGen::saveHist(string name)
{
    // Histogramm erstellen
    Eigen::VectorXi hist(bins);
    for ( int n = 0; n < r.size(); n++ )
        hist((int) r(n)*bins)++;

    ofstream fout;
    fout.open(name+".dat");
    fout << "#bin\t count" << endl;

    for ( int i = 0; i < bins; i++ )
        fout << i/bins << "\t" << hist(i) << endl;

    fout.close();
}

