#include "../../Eigen/Eigen/Dense" //Eigen liegt oberhalb der B**-Ordner
#include <math.h>
#include <iostream>
#include "mdsimulation.h"
#include "initializeSquareGrid.h"
#include "ljforces.h"
#include <chrono>

using namespace Eigen;
using namespace std;

/**
  * simulate(..) macht folgendes:
  * 1. Gitter erstellen (s. initializeSquareGrid.cpp)
  * 2. zufällige Geschwindigkeiten generieren
  * 3. Geschwindigkeiten über feste Temperatur neu skalieren
  * 4.
  **/

void simulate(int AnzTeilchen, double Dichte, double Temperatur) {
    // Erstelle Gitter
    MatrixXd particleInfo(2, AnzTeilchen);
    Quadratgitter quadratgitter;
    particleInfo = quadratgitter.erstellequadrgitter(AnzTeilchen,Dichte);
    double L = pow((AnzTeilchen/Dichte),0.5);

    //aktueller Unix-Timestamp für rand-Nummern
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    // Setze Anfangsgeschwindigkeit mit zufälligen Nummern
    MatrixXd vmatrix(2, AnzTeilchen); //enthält alle Geschwindigkeitsvektoren
    Vector2d vges(0,0); //Schwerpunktsgeschwindigkeit, init mit 0
    double vges2; //Mittlere Geschwindigkeitsquadrate, aufsummiert.
    for (int zeile = 0; zeile < AnzTeilchen; zeile++)
    {
        Vector2d v(0,0);
        for(int i = 0; i<2; i++) 
        /*
        * x und y, hier hätte man sicher auch direkt die Vektoren als ganze addieren können, also Zeile 33 aus der for-schleife, dann "<<" nutzen. Allerdings hat das mit der Dimension 2xAnzTeil-Matrix vs. 2-dim Vektor nicht so einfach geklappt. Das umgeht das Problem.
        */
        {
            //Setze Geschwindigkeit auf rand Zahl (pro Komponente)
            v(i) = distribution(generator);
            vges(i) = (vges(i) + v(i))/AnzTeilchen;
            /*v_i aufsummieren, vges (r.S.) ist jeweils der alte volle v-Vektor, geteilt durch AnzTeilchen wg. Schwerpunktsgeschwindigkeit
            */
            vmatrix(i,zeile) = (v(i)-vges(i)); //Schwerpunktsgeschwindigkeit abziehen, sonst Drift
        }
        //Init: Ekin = 1/2mv^2 = Dim kbT/2 (ÄquipartTheorem) , Dim = 2, kb=1, m=1

        //Skalarprodukt von v => v^2, vges2= v^2/N für Skalierungsfaktor "vskalierungfaktor"
        vges2 += v.dot(v)/AnzTeilchen; 
    }
    //MatrixXd fullInfo(4, AnzTeilchen); //Fügt mir beide Teilmatrizen mit Ort und Geschwindigkeit zusammen, noch unused (evtl später nützlich)
    /*fullInfo << particleInfo,
                vmatrix;
    */

    // Skalierungsfaktor für die Geschwindigkeiten
    double vskalierungfaktor = sqrt(2*Temperatur/vges2);
    for (int zeile = 0; zeile < AnzTeilchen; zeile++)
    {
        for(int i=0; i<2; i++) //für 2 dim
        {
            //Geschwindigkeiten neu skalieren, um bestimmte E bei bestimmter T zu erzeugen.
            vmatrix(i,zeile) = vskalierungfaktor*vmatrix(i,zeile);
        }
    }

    //Jetzt Kraftberechnung (s. ljforces.cpp), ToDo
    LJForces ljforces;
    ljforces.kraft(particleInfo, L);
}

int main() {
    // *******************
    // *    Init         *
    // *******************

    // Konfig-Einstellungen
    int AnzTeilchen = 1000;     // Anzahl an Teilchen
    double Dichte = 0.32;       // Dichte der Teilchen
    //int masse = 1;              // Masse der Teilchen (1 lt Aufgabe)

    // Simulationsparameter
    double dt = 0.0001;         // Integrationszeit
    double dt2 = dt*dt;         // (Integrationszeit)^2
    double Temperatur = 2.0;    // Temperatur für die Simulation

    int N = 10000;          // Anzahl an Integrationsschritten
    int TastFreq = 10;    // TastFrequenz
    int zaehler = 0;        // TastZähler
    int DruckFreq = 100;    // Druckfrequenz

    simulate(AnzTeilchen, Dichte, Temperatur);
    return 0;
}
