#include "../../Eigen/Eigen/Dense" //Eigen liegt oberhalb der B**-Ordner
#include <iostream>
#include "ljforces.h"
#include "periodRB.h"
#include <math.h>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
  * Berechnet die Kräfte auf ein Ensemble von LJ-Teilchen in 2D
  * Gibt die Kräfte auf alle Teilchen in einer bestimmten Ordnung zurück
  * $\sigma$ (Längenskala) = $\epsilon$ (Energieskala) = 1
  **/

MatrixXd LJForces::kraft(MatrixXd particleinfo, double L)
{
    int AnzTeilchen = particleinfo.cols();
    int dim = particleinfo.rows();
    MatrixXd forces(dim, AnzTeilchen); 
    //selbe Dim wie particleinfo, rows: 2


    //Matrixeinträge auf 0 setzen, kann man sicher effizienter schreiben
    for(int teile = 0; teile< AnzTeilchen; teile++) {
        for(int komp = 0; komp<dim; komp++) {
            forces(komp, teile) = 0;
        }
    }
    // Abstandsvektor
    Vector2d dr(0,0);
    for(int TeilA = 0; TeilA < (AnzTeilchen-1); TeilA++) //über alle TeilchenPAARE (<=> Grenzen: Ränder) gehen
    {
        for(int TeilB = (TeilA+2); TeilB < AnzTeilchen; TeilB++) //s. Komm erste Schleife
        {
            for(int komp=0; komp<dim; komp++) //für x und y (geht sicher schöner, wenn man schon mit Vec rechnet)
            {
                // Distanz zwischen zwei Teilchen    
                dr(komp) = particleinfo(komp,TeilA) - particleinfo(komp,TeilB);
                //cout << dr << endl;

                // Period RB beachten
                PeriodRB periodRB;
                dr = periodRB.kurzerWeg(dr, L);
                // dr2 = dr^2
                double dr2 = dr.dot(dr);

                /* 
                 * LJ-Potential:
                 * V(r) = 4*epsilon* ((sigma/r)^12 - (sigma/r)^6), sigma = epsilon = 1 (vgl Aufgabe)
                 * F = -gradV
                 * F(r)_x = 4 * (x/r) * (12*(1/r)^13 - 6*(1/r)^7)
                 *        = x *48 * (1/r)^8 * ((1/r)^6 - 0.5) 
                 * für y analog
                 */

                // Zuerst 1/r^2
                double r2 = 1.0/dr2;
                // Klammerterm nach "48"
                double force = pow(r2, 4) * (pow(r2,3) - 0.5);
                
                // Newtons Dritte (Actio = Reactio) für die beiden Teilchen, also WW
                forces(komp, TeilA) = forces(komp, TeilA) + dr(komp)*force; //das aufaddieren (Reactio)
                forces(komp, TeilB) = forces(komp, TeilB) - dr(komp)*force; //oder abziehen (Actio)
            }
        }
    }
    //Faktor von der Kraft
    forces = forces*48;
    return forces;
}
