#include <eigen3/Eigen/Dense>
#include <iostream>
#include <math.h>
#include <iostream>
#include "ljforces.h"
#include "periodRB.h"

using namespace Eigen;
using namespace std;

void kraft(MatrixXd &forces, MatrixXd particleinfo, double L, bool active, VectorXd &bins, double &V)
{
    int AnzTeilchen = particleinfo.cols();

    // Abstandsvektor
    Vector2d dr(0,0);
    double absr = 0;

    for ( int TeilA = 0; TeilA < (AnzTeilchen-1); TeilA++ ) //über alle TeilchenPAARE (<=> Grenzen: Ränder) gehen
    {
        Vector2d tempA(0, 0);
        tempA = particleinfo.block(0, TeilA, 2, 1); //Gibt die TeilA-te Teilchen-kords aus
        for ( int TeilB = ( TeilA + 1 ); TeilB < AnzTeilchen; TeilB++ ) //s. Komm erste Schleife
        {
            Vector2d tempB(0, 0);
            tempB = particleinfo.block(0, TeilB, 2, 1);

            // Distanz zwischen zwei Teilchen
            dr = tempA - tempB;
            // Period RB beachten, also kürzesten Weg zum nächsten Teil (auch in Nachbarbox)
            dr = kurzerWeg(dr, L);
            // Nimmt den 0-Vektor raus
            if ( dr(0) == 0 && dr(1) == 0 )
            {
                continue;
            }
            /* Potential(Betrag von dr), über alle TeilchenPaare summiert
             * ist Epot, V_LJ(r) = 4*(pow(pow(r,-1),12)-pow(pow(r,-1),6));
             **/

            absr = dr.norm();

            /* ** Histogramm ** */
            if ( active )
            {
                bins((int) ( absr*2*bins.size()/L ))++;
            }

            V += 4*( pow(absr, -12) - pow(absr, -6) );

            /*
             * LJ-Potential:
             * V(r) = 4*epsilon* ((sigma/r)^12 - (sigma/r)^6), sigma = epsilon = 1 (vgl Aufgabe)
             * F = -gradV
             * F(r)_x = 4 * (x/r) * (12*(1/r)^13 - 6*(1/r)^7)
             *        = x *48 * (1/r)^8 * ((1/r)^6 - 0.5)
             * für y analog
             */

            // Zuerst 1/r^2
            double r2 = 1.0/dr.squaredNorm();
            // Klammerterm nach "48"
            double force = pow(r2, 4) * ( pow(r2, 3) - 0.5 );

            //temp 2dim Vektoren für Matrixspalten (im Prinzip Kraftkords der Teile)
            Vector2d tempf1(0, 0);
            tempf1 = forces.block(0, TeilA, 2, 1);
            Vector2d tempf2(0, 0);
            tempf2 = forces.block(0, TeilB, 2, 1);

            // Newtons Dritte (Actio = Reactio) für die beiden Teilchen, also WW
            tempf1 += dr*force; //das aufaddieren (Reactio)
            tempf2 -= dr*force; //oder abziehen (Actio)

            //Zusammenpacken. Da wir vorher mit Vektoren rechnen, muss jede Komp vorher durch sein
            for ( int komp=0; komp<2; komp++ )
            {
                forces(komp, TeilA) = tempf1(komp)*dr(komp);
                forces(komp, TeilB) = tempf2(komp)*dr(komp);
            }
        }
    }
    //Faktor von der Kraft
    forces = forces*48;
}
