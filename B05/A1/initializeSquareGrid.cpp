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

MatrixXd Quadratgitter::erstellequadrgitter(int AnzTeilchen, double dichte) {
    //int AnzTeilchen Anzahl an Teilchen für das Gitter
    //double dichte Dichte der Partikel
    MatrixXd particleinfo(2,AnzTeilchen); //2xAnzTeilchen Array mit einem Ortsvektor pro Spalte

    // Größe der Box (RB: Volumen, bzw. hier Fläche, da 2D); Teilchendichte umgeformt, auf Dim der Länge
    double L = pow((AnzTeilchen/dichte),1/3.0);

    // Finde das kleinste perfekte Quadrat größer oder gleich der Anz an Partikeln
    // Bietet Gitterplätze für die Teile
    double AnzQuadrat = 2;

    while (pow(AnzQuadrat,2) < AnzTeilchen)
    {
       AnzQuadrat = AnzQuadrat + 1;
    }

    // 2D Index zum Zählen der Plätze, Startposition untere linke Kante
    Vector2d aktPos;
    aktPos << 0, 0;
    // Einkommentieren, wenn wie unten gewünscht
    /*Vector2d Schrittweite;
    Schrittweite << 0.5, 0.5;*/

    // Kord setzen:
    Vector2d tempTeil(0,0);
    // wie gewünscht nach 1a, zaehler ist Spaltenzähler
    int zaehler = 0;
    for(int n = 0; n<=3; n++) {
        for(int m=0; m<=3;m++) {
            //45/46 geht leider nicht in einer Zeile mit ().. Mag Eigen nicht..
            tempTeil << 1+2*n , 1+2*m;
            tempTeil = tempTeil *1/8.0 * L;
            for(int komp=0;komp<2;komp++)
            {
                particleinfo(komp,zaehler) = tempTeil(komp);
            }
            zaehler++;
        }
    }
    for(int Teilchen = 0; Teilchen<AnzTeilchen ;Teilchen++)
    {
        //Abwandlung der Aufgabe:
        // Kord setzen: Teilchenpos in der Matrix setzen auf jew Kante
        /*Vector2d tempTeil(0,0);
        tempTeil = particleinfo.block(0,Teilchen,2,1);
        tempTeil = (aktPos + Schrittweite)*(L/AnzQuadrat);
        for(int komp=0;komp<2;komp++)
        {
            particleinfo(komp,Teilchen) = tempTeil(komp);
        }*/


        // Rastern des Pointers
        aktPos(0) = aktPos(0) + 1;
        if ( aktPos(0) == AnzQuadrat )
        {
            aktPos(0) = 0;
            aktPos(1) = aktPos(1) + 1;
        }    
    }
    return particleinfo;
}
