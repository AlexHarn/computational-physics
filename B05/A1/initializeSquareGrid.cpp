#include "../../Eigen/Eigen/Dense" //Eigen liegt oberhalb der B**-Ordner
#include <math.h>
#include <iostream>
#include "initializeSquareGrid.h"

using namespace Eigen;
using namespace std;

/**
  * Erzeugung eines quadratischen Gitters und glm. Anordnung der Teilchen auf dem Gitter,
  * um kleine Abstände zu vermeiden.
  **/

MatrixXd Quadratgitter::erstellequadrgitter(int AnzTeilchen) {
    //int AnzTeilchen Anzahl an Teilchen für das Gitter
    //double dichte Dichte der Partikel
    MatrixXd particleinfo(2,AnzTeilchen); //2xAnzTeilchen Array mit einem Ortsvektor pro Spalte

    double L = 8;

    // Finde das kleinste perfekte Quadrat größer oder gleich der Anz an Partikeln
    // Bietet Gitterplätze für die Teile
    double AnzQuadrat = sqrt(AnzTeilchen);

    // Kord setzen:
    Vector2d tempTeil(0,0);
    // wie gewünscht nach 1a, zaehler ist Spaltenzähler
    int zaehler = 0;
    for(int n = 0; n<=3; n++)
    {
        for(int m=0; m<=3;m++)
        {
            tempTeil << 1+2*n*1/8.0 * L , 1+2*m*1/8.0*L;
            for(int komp=0;komp<2;komp++)
            {
                particleinfo(komp,zaehler) = tempTeil(komp);
            }
            zaehler++;
        }
    }
    return particleinfo;
}
