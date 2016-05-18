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
    for(int j = 0; j < AnzTeilchen; j++) { //Anz-Teil-Spalten in der Matrix
        for(int i = 0; i < 2; i++) { //durch Beginn bei 0 passt die Grenze wieder: 2dim Vektor
            particleinfo(i,j) = 0;        //füllt mit Nullen
        }
    }
    // Größe der Box (RB: Volumen, bzw. hier Fläche, da 2D); Teilchendichte umgeformt, auf Dim der Länge
    double L = pow((AnzTeilchen/dichte),0.5);

    // Finde das kleinste perfekte Quadrat größer oder gleich der Anz an Partikeln
    // Bietet Gitterplätze für die Teile
    double AnzQuadrat = 2; //Seite des perfekten Quadrates

    while (pow(AnzQuadrat,2) < AnzTeilchen)
    {
       AnzQuadrat = AnzQuadrat + 1;
    }

    // 2D Index zum Zählen der Plätze, Startposition untere linke Kante
    Vector2d aktPos;
    aktPos << 0, 0;
    Vector2d Schrittweite;
    Schrittweite << 0.5, 0.5;

    for(int Teilchen = 1; Teilchen<AnzTeilchen ;Teilchen++)
    {
        // Kord setzen
        //Addition ist in Eigen etwas schräg. aktPos + (0.5,0.5) fkt nicht
        for(int komp = 0; komp < 2; komp++)
        {
            particleinfo(komp,Teilchen) = (aktPos(komp) + Schrittweite(komp))*(L/AnzQuadrat); //Teilchenpos in der Matrix setzen auf jew Kante
        //Doppelpunkt gibt komplette Zeile,bzw hier Spalte aus (also (x,y) des Teilchens).
        }

        // Rastern des Pointers
        aktPos(0) = aktPos(0) + 1;
        if ( aktPos(0) == AnzQuadrat )
        {
            aktPos(0) = 0;
            aktPos(1) = aktPos(1) + 1;
        }
    }
    //cout << "Matrix mit allen Vektoren: " << particleinfo << endl;
    //cout << "Pointer: " << aktPos << endl;
    //cout << "perfektes Quadrat: " << AnzQuadrat << endl;
    //cout << "Länge: "<< L << endl;
    return particleinfo;
}
