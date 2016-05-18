#include "../../Eigen/Eigen/Dense" //Eigen liegt oberhalb der B**-Ordner
#include <math.h>
#include <iostream>
#include "mdsimulation.h"
#include "initializeSquareGrid.h"
#include <chrono>

using namespace Eigen;

void simulate() {
    // *******************
    // *    Init         *
    // *******************

    // Konfig-Einstellungen
    int AnzTeilchen = 1000;     // Anzahl an Teilchen
    double Dichte = 0.32;       // Dichte der Teilchen
    int masse = 1;              // Masse der Teilchen (1 lt Aufgabe)

    // Simulationsparameter
    double dt = 0.0001;         // Integrationszeit
    double dt2 = dt*dt;         // (Integrationszeit)^2
    double Temperatur = 2.0;    // Temperatur für die Simulation

    int N = 10000;          // Anzahl an Integrationsschritten
    int TastFreq = 10;    // TastFrequenz
    int zaehler = 0;        // TastZähler
    int DruckFreq = 100;    // Druckfrequenz

    // Erstelle Gitter
    MatrixXd = parInfuLaenge(2, AnzTeilchen);
    parInfuLaenge = this->ErstelleQuadrGitter(AnzTeilchen,Dichte);

    //aktueller Unix-Timestamp für rand-Nummern
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    // Setze Anfangsgeschwindigkeit mit zufälligen Nummern
    Vector2d v(0,0);
    MatrixXd vmatrix; //enthält alle Geschwindigkeitsvektoren
    for (int komp = 1; komp < AnzTeilchen, komp++)
    {
        for(int i = 1; i<2; i++)
        {
          //Setze Geschwindigkeit auf rand Zahl
          v(i,komp) = distribution(generator);
        }
        vmatrix << v(komp);
    }


    Vector2d vges; //Schwerpunktsgeschwindigkeit
    //v_i aufsummieren
    for (int index = 1; index < AnzTeilchen; index++) {
      vges += v(index)/AnzTeilchen;
    }

    for (int komp = 1; komp < AnzTeilchen, komp++)
    {
        for(int i = 1; i<2; i++)
        {
          v(i,komp) = v(i,komp) - vges; //Schwerpunktsgeschwindigkeit abziehen, sonst Drift
        }
    }

    //Skalierungsfaktor bestimmen
}
