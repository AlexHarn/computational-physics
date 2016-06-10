#include <eigen3/Eigen/Core>
#include <cmath>
#include <fstream>
#include <iostream>
#include "lingen.h"
#ifndef M_PI
#define M_PI 3.14159265359
#endif

using namespace std;

void LinConGen::congruent(int64_t r0, unsigned int a, unsigned int c, unsigned int m, unsigned int N=0)
{
    if ( N == 0 )
        N = m - 1;

    V64 r(N+1);
    this->r.resize(N+1);
    r(0) = r0;

    if ( N < m )
    {
        for ( unsigned int n = 0; n < N; n++ )
            r(n+1) = ( a*r(n) + c )%m;
    }
    else // Zahlen wiederholen sich für n >= m, n = m (also r(0)) wird ausgelassen, der Seed gehört nicht zu den Zufallszahlen
    {
        for ( unsigned int n = 0; n < m - 2; n++ )
            r(n+1) = ( a*r(n) + c )%m;
        for ( unsigned int n = m - 1; n < N+1; n++ )
            r(n) = r(( (n+1) % m ) + 1);
    }

    //this->r = ( r/(double)m ).cast<double>();
    this->r = r.cast<double>()/m;
}

void LinConGen::boxMuller()
{
    assert(r.size() > 2);
    dist.resize(r.size()-1);
    double tmp;
    for ( unsigned int n = 0; n < dist.size(); n+=2 )
    {
        tmp = sqrt(-2*log(r(n)));
        dist(n) = tmp*cos(2*M_PI*r(n+1));
        dist(n+1) = tmp*sin(2*M_PI*r(n+1));
    }
}

void LinConGen::centralLimit(int N)
{
    assert(r.size() > N);
    dist.resize((int) r.size()/N);
    for ( unsigned int n = 0; n < dist.size(); n++ )
        dist(n) = r.segment(n*N, N).sum() - N/2;
}

void LinConGen::neumann(unsigned int N)
{
    unsigned int n = 0;
    unsigned int zaehler = 0;
    double x;
    double y;
    tempdist.resize(N);
    while ( zaehler < N )       // dist soll am Ende N Einträge haben
    {
        x = dist(n)+M_PI/2;     // x = Zufallszahl aus verschobener Gaußverteilung (g(x-M_PI/2))
        y = r(n+1)*1.5*exp(-pow((x-M_PI/2), 2)/2)/( sqrt(2*M_PI) );     // y = gleichverteilte Zahl mal 1.5 mal g(x-M_PI/2)
        if ( y < sin(x)/2 )     // test ob y < p(x)
        {
            if ( x >= 0 && x <= M_PI )     // test ob x aus [0, M_PI]
            {
                tempdist(zaehler) = x;
                zaehler++;
            }
        }
        n++;
        assert(n < dist.size());
    }
    dist.resize(N);
    dist = tempdist;
}

void LinConGen::transform(unsigned int N)
{
    dist.resize(N);
    double temp;
    unsigned int n = 0;
    unsigned int zaehler = 0;
    while ( zaehler < N )
    {
        temp = pow(r(n+1), 1./3);
        if ( temp < 1 )
        {
            dist(zaehler) = temp;
            zaehler++;
        }
        n++;
        assert(n < r.size()-1);
    }
}

void LinConGen::save(string name)
{
    ofstream fout;
    fout.open(name+".dat");
    fout << "#r(n), n = Zeilennummer" << endl;
    for ( unsigned int n = 1; n < r.size(); n++ )
        fout << r(n) << endl;
    fout.close();
}

void LinConGen::saveDist(string name)
{
    ofstream fout;
    fout.open(name+".dat");
    fout << "#dist(n), n = Zeilennummer" << endl;
    for ( unsigned int n = 0; n < dist.size(); n++ )
        fout << dist(n) << endl;
    fout.close();
}

void LinConGen::saveHist(string name, int bins)
{
    // Histogramm erstellen
    V64 hist(bins);
    for ( int n = 1; n < r.size(); n++ )
        hist((int) r(n)*bins)++;

    ofstream fout;
    fout.open(name+".dat");
    fout << "#bin\t count" << endl;

    for ( int i = 0; i < bins; i++ )
        fout << i/(float) bins << "\t" << hist(i) << endl;

    fout.close();
}
