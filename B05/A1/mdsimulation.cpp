#include "../../Eigen/Eigen/Dense" //Eigen liegt oberhalb der B**-Ordner
#include <math.h>
#include <iostream>
#include "mdsimulation.h"
#include "initializeSquareGrid.h"
#include "ljforces.h"
#include "periodRB.h"
#include <chrono>

using namespace Eigen;
using namespace std;

/**
  * simulate(..) macht folgendes:
  * 1. Gitter erstellen (s. initializeSquareGrid.cpp)
  * 2. zufällige Geschwindigkeiten generieren
  * 3. Geschwindigkeiten über feste Temperatur neu skalieren
  * 4. Kraftberechnung
  * 5. Integrationsschritt (+Zähler, aktuell auskommentiert)
  **/

void simulate(int AnzTeilchen, double Dichte, double Temperatur, double dt, double dt2, int N, int DruckFreq) {
    // Erstelle Gitter
    MatrixXd particleInfo(2, AnzTeilchen);
    Quadratgitter quadratgitter;
    // Hier Fehler, weil Positionen von Nachbarteilchen gleich (darf nicht sein)!
    particleInfo = quadratgitter.erstellequadrgitter(AnzTeilchen,Dichte);
    double L = pow((AnzTeilchen/Dichte),1.0/3);

    //aktueller Unix-Timestamp für rand-Nummern
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    // Setze Anfangsgeschwindigkeit mit zufälligen Nummern
    MatrixXd vmatrix(2, AnzTeilchen); //enthält alle Geschwindigkeitsvektoren
    Vector2d vges(0,0); //Schwerpunktsgeschwindigkeit, init mit 0, lt. Aufgabe a
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
            v(i) = distribution(generator); // wie in Übung besprochen
            vges(i) = (vges(i) + v(i))/AnzTeilchen; //Berechnung der Schwerpunktsgeschwindigkeit, Aufg. b
            /** TODO Das sollte vermutlich in die .dat-Datei, damit die ÜL zufrieden sind */
            /*v_i aufsummieren, vges (r.S.) ist jeweils der alte volle v-Vektor, geteilt durch AnzTeilchen wg. Schwerpunktsgeschwindigkeit
            */
            //Schwerpunktsgeschwindigkeit abziehen, sonst Drift, in neue Matrix, um alle Vektoren mizunehmen
            vmatrix(i,zeile) = (v(i)-vges(i)); 
        }
        //Init: Ekin = 1/2mv^2 = Dim kbT/2 (ÄquipartTheorem) , Dim = 2, kb=1, m=1

        //Skalarprodukt von v => v^2, vges2= v^2/N für Skalierungsfaktor "vskalierungfaktor"
        vges2 += v.dot(v)/AnzTeilchen; 
    }

    // Skalierungsfaktor für die Geschwindigkeiten (lt. Aufgabe a)
    double vskalierungfaktor = sqrt(2*Temperatur/vges2);
    for (int zeile = 0; zeile < AnzTeilchen; zeile++)
    {
        for(int i=0; i<2; i++) //für 2 dim
        {
            //Geschwindigkeiten neu skalieren, um bestimmte E bei bestimmter T zu erzeugen.
            vmatrix(i,zeile) = vskalierungfaktor*vmatrix(i,zeile);
        }
    }

    //Jetzt Kraftberechnung (s. ljforces.cpp)
    LJForces ljforces;
    MatrixXd forces(2, AnzTeilchen);
    forces = ljforces.kraft(particleInfo, L);

    // Histogramm (fehlt)
    // TODO: Kin, pot Energie berechnen, Temperatur berechnen, bzw. in .dat-File eintragen <-> Plots
    // Nach wie vielen Zeitschritten äquilib. ihr System? 
    // T, g(r) messen für verschiedene T0 (N=16, L = 8 sigma). Phasen erkennen. -< Plots

    /*****************
     * MD-Simulation *
     *****************/

    double zeitzaehler = 0;
    
    for(int step=0; step < N; step++)
    {
        /* ** Erster Integrationsschritt **  */
        
        // Aktualisiere Positionen
        //Ich liebe Eigen :-), das aktualisiert hier alle Kords auf einmal :-D
        particleInfo = particleInfo + dt*vmatrix + 0.5*dt2*forces;    
        
        // RB beachten
        for(int Teil = 0; Teil < AnzTeilchen; Teil++)
        {
                //Init eines Temp-Vektors mit Teilkords
                Vector2d tempRB(0,0);
                tempRB = particleInfo.block(0,Teil,2,1);

                //RB beachten
                PeriodRB periodRB;
                tempRB = periodRB.kurzerWeg(tempRB, L);
        }
        
        // Geschwindigkeiten aktualisieren, alle Kords auf einmal *.*
        vmatrix = vmatrix + 0.5*dt*forces;
    
        /** Neue Kräfte berechnen **/        
        forces = ljforces.kraft(particleInfo, L);
        
        /* ** Zweiter Integrationsschritt ** */

        // Geschwindigkeiten aktualisieren
        vmatrix = vmatrix + 0.5*dt*forces;

        // evtl später Thermostat

        /* ** Zeitschritt ** */
        zeitzaehler = zeitzaehler + dt;

        /* ** Schrittzähler, einkommentieren, wenn gewünscht ** */
        /*if(step % DruckFreq == 0)
        {
            cout << step << endl;
        }*/
    }
}

int main() {
    // *******************
    // *    Init         *
    // *******************

    // Konfig-Einstellungen
    int AnzTeilchen = 16;     // Anzahl an Teilchen
    double Dichte = 0.32;       // Dichte der Teilchen
    //int masse = 1;              // Masse der Teilchen (1 lt Aufgabe)

    // Simulationsparameter
    double dt = 0.0001;         // Integrationszeit
    double dt2 = dt*dt;         // (Integrationszeit)^2
    double Temperatur = 1.0;    // Temperatur für die Simulation, T(0) = 1 lt. Aufgabe

    int N = 10000;          // Anzahl an Integrationsschritten
    int TastFreq = 10;    // TastFrequenz
    int zaehler = 0;        // TastZähler
    int DruckFreq = 100;    // Druckfrequenz

    simulate(AnzTeilchen, Dichte, Temperatur, dt, dt2, N, DruckFreq);
    return 0;
}
