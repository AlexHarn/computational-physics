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

LinConGen::LinConGen(int bins)//, bool a2)
{
    //hier stehen dann evtl noch ein paar Restnullen
    //durch den Cut bei den if-Abfragen, aber das kostet nicht so viel Aufwand
    x = V::Zero(2);
    y = V::Zero(2);
    this->bins = bins;

    mean = 0;
    variance = 0;
    Zn = 0;
    Sn = 0;
    /*if(a2)
    {
        calc(1234,16807,0,2147483647);
    }*/
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

// geht effizienter
//Gibt Wahrscheinlichkeit aus, dass Wert in den jew Blöcken des Histogramms liegt
//nicht wirklich die Wahrscheinlichkeit, bei der ein Element aus dem gesamten Zahlenraum gezogen wird
//(ideale) Gleichverteilung--> Wahrscheinlichkeiten für alle Elemente gleich groß; 
double LinConGen::prob(double numbers)
{   
        /*double sum = hist.sum();
        if(i < 0.1 && i >= 0)
            return hist(0)/sum;
        else if(i < 0.2 && i >= 0.1)
            return hist(1)/sum;
        else if(i < 0.3 && i >= 0.2)
            return hist(2)/sum;
        else if(i < 0.4 && i >= 0.3)
            return hist(3)/sum;
        else if(i < 0.5 && i >= 0.4)
            return hist(4)/sum;
        else if(i < 0.6 && i >= 0.5)
            return hist(5)/sum;
        else if(i < 0.7 && i >= 0.6)
            return hist(6)/sum;
        else if(i < 0.8 && i >= 0.7)
            return hist(7)/sum;
        else if(i < 0.9 && i >= 0.8)
            return hist(8)/sum;
        else if(i < 1.0 && i >= 0.9)
            return hist(9)/sum;
        else
        { 
            cout << "Hier ist was schief gelaufen!" << endl;
            return -1;
        }*/
        return 1.0/numbers;
}

/*void LinConGen::clt()
{ 
    // Zufälliger Seed:
    srand (time(NULL));

    //n-te Teilsumme (zufällig) 1000 bis N (100000)
    int numbers = rand() % 10000 + 1000;

    // N sollte groß sein, da N -> infty gegen die Standardnormalvert konvergiert
    for(int i = 0; i<numbers; i++)
    {
        Sn += saveVec(i);
        mean += saveVec(i) * prob( numbers ) ; //ehemals prob(saveVec(i))
    }
    for(int j = 0; j<numbers;j++)
    {
        variance += pow( ( saveVec(j)-mean ) ,2) * prob( numbers ); //ehemals s oben
    }

    // zentraler Grenzwertsatz nach Wikipedia
    Zn = ( Sn - mean )/( sqrt(variance) );

    cout << "Anzahl an Werten " << numbers << endl;
    cout << "Varianz " <<variance << endl;
    cout << "Mittelwert " <<mean << endl;
    cout << "Zufallszahl " << Zn << endl;
}*/
/*
    Mittelwert 0 bei alternierenden Messwerten mit jeweils gleicher Wahrscheinlichkeit
    Standardabweichung = 1, also Varianz =1, (Erwartungswert sei Mittelwert), also Var(X) = sum E(X^2)*P(X=x),,,,
    Nachteile der Methode- Generator a nicht vollständig gleichverteilt 
    (dadurch keine vollständige Gültigkeit der Wahrscheinlichkeit), außerdem kommulieren Fehler des Generators
    zB endliche Länge der Periode, Abwesenheit von Korrelationen zwischen Zahlen,,,
*/

/*void LinConGen::fullreset()
{
    saveVec.setZero();
    hist.setZero();
    x.setZero();
    y.setZero();
    rnew = 0;
    rold = 0;
    r0 = 0;
    a = 0;
    c = 0;
    m = 0;
    M = 0;
    mean = 0;
    variance = 0;
    Sn = 0;
    Zn = 0;
}

void LinConGen::partreset()
{
    x.setZero();
    y.setZero();
    mean = 0;
    variance = 0;
    Sn = 0;
    Zn = 0;
}*/

