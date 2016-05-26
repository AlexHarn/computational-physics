#include <math.h>
#include <eigen3/Eigen/Core>
#include <iostream>
#include "initializeSquareGrid.h"

using namespace Eigen;
using namespace std;

MatrixXd erstellequadrgitter(int AnzTeilchen, double L)
{
    //int AnzTeilchen Anzahl an Teilchen für das Gitter
    //double dichte Dichte der Partikel
    MatrixXd particleinfo(2, AnzTeilchen); //2xAnzTeilchen Array mit einem Ortsvektor pro Spalte

    double N = sqrt(AnzTeilchen);

    // Kord setzen:
    Vector2d tempTeil(0, 0);
    // wie gewünscht nach 1a, zaehler ist Spaltenzähler
    int zaehler = 0;
    for ( int n = 0; n<N; n++ )
    {
        for ( int m=0; m<N ;m++ )
        {
            tempTeil << (1+2*n)*L/8.0, (1+2*m)*L/8.0;
            for ( int komp=0; komp<2; komp++ )
            {
                particleinfo(komp, zaehler) = tempTeil(komp);
            }
            zaehler++;
        }
    }
    return particleinfo;
}
