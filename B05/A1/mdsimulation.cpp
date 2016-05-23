#include "../../Eigen/Eigen/Dense" //Eigen liegt oberhalb der B**-Ordner
#include <math.h>
#include <iostream>
#include "mdsimulation.h"
#include "initializeSquareGrid.h"
#include "ljforces.h"
#include "periodRB.h"
#include <chrono>
#include <fstream>
#include <string>

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

void simulate(int AnzTeilchen, double Dichte, double Temperatur, double dt, double dt2, int N, int DruckFreq, string fname) {
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

    for (int zeile = 0; zeile < AnzTeilchen; zeile++)
    {
        Vector2d v(0,0);
        MatrixXd vschwmatrix(2, AnzTeilchen);
        for(int i = 0; i<2; i++) 
        /*
        * x und y, hier hätte man sicher auch direkt die Vektoren als ganze addieren können, also Zeile 33 aus der for-schleife, dann "<<" nutzen. Allerdings hat das mit der Dimension 2xAnzTeil-Matrix vs. 2-dim Vektor nicht so einfach geklappt. Das umgeht das Problem.
        */
        {
            //Setze Geschwindigkeit auf rand Zahl (pro Komponente)
            v(i) = distribution(generator); // wie in Übung besprochen
            vges(i) = (vges(i) + v(i))/AnzTeilchen; //Berechnung der Schwerpunktsgeschwindigkeit, Aufg. b
            /*v_i aufsummieren, vges (r.S.) ist jeweils der alte volle v-Vektor, geteilt durch AnzTeilchen wg. Schwerpunktsgeschwindigkeit
            */
            //Schwerpunktsgeschwindigkeit abziehen, sonst Drift, in neue Matrix, um alle Vektoren mizunehmen
            vmatrix(i,zeile) = (v(i)-vges(i));
            // Für die .dat-Datei Schwerpunktsgeschwindigkeiten speichern
            vschwmatrix(i, zeile) = vges(i);
        }
        //Init: Ekin = 1/2mv^2 = Nf kbT/2 (ÄquipartTheorem) , Nf = 2N-2, kb=1, m=1 (2 Raumdimensionen, 2 FG pro Schwerpunktsimpuls)

        double Ekin = Ekin + 0.5*v.dot(v);
        double T = 2*Ekin/(2*AnzTeilchen-2);
        // Epot wird in ljforces bestimmt
        /* ** Isokinetischer Thermostat ** */

        // Skalierungsfaktor für die Geschwindigkeiten (lt. Aufgabe a) nach 5.40
        double vskalierungfaktor = sqrt(Temperatur/T);

        for(int i=0; i<2; i++) //für 2 dim
        {
            //Geschwindigkeiten neu skalieren, um bestimmte E bei bestimmter T zu erzeugen (s. Kap 5.5.1).
            vmatrix(i,zeile) = vskalierungfaktor*vmatrix(i,zeile);
        }
    }

    //Jetzt Kraftberechnung (s. ljforces.cpp)
    LJForces ljforces;
    MatrixXd forces(2, AnzTeilchen);
    forces = ljforces.kraft(particleInfo, L);

    // Histogramm (fehlt), Epot spinnt noch rum, V wird in ljforces.cpp einfach riesig..
    //Äquiliebr.zeit bestimmen: Dann, wenn Ekin ~= Epot
//erste Idee: umgeht double-Vergleich, indem beide doubles auf int gecastet werden. Dann sind mir Nachkommastellen egal
//geht das schöner?
    // if ((int)Ekin == (int)Epot) 
        //cout << "Äquizeit" << zeitzaehler << endl;
    // TODO: Ekin, Epot, Temperatur in .dat-File eintragen <-> Plots
    // Nach wie vielen Zeitschritten äquilib. ihr System? 
    // TODO g(r) implementieren
    // T, g(r) messen für verschiedene T0 (N=16, L = 8 sigma). Phasen erkennen. -< Plots
    // Aufg e: Was passiert nun bei T = 0.01? Schauen Sie sich
    // die potentielle, kinetische und Gesamtenergie während der Äquilibrierung an und messen
    // Sie g(r) an dem äquilibrierten System. Es kann auch interessant sein sich Snapshots, also
    //Bilder der Konfiguration des Systems, anzuschauen. -> Plots der Energien, Paarkorrfkt.

    /****************************
     * kanonische MD-Simulation *
     ****************************/

    double zeitzaehler = 0;
    
    for(int step=0; step < N; step++)
    {
        /* ** Erster Integrationsschritt **  */
        
        // Aktualisiere Positionen
        //Ich liebe Eigen :-), das aktualisiert hier alle Kords auf einmal :-D
        //Geschwindigkeits-Verlet-Algorithmus nach 4.25, h = dt
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
        // Rest vom v-Verlet
        vmatrix = vmatrix + 0.5*dt*forces;
    
        /** Neue Kräfte berechnen **/        
        forces = ljforces.kraft(particleInfo, L);
        
        /* ** Zweiter Integrationsschritt ** */

        // Geschwindigkeiten aktualisieren, v-Verlet..
        vmatrix = vmatrix + 0.5*dt*forces;

        /* ** Zeitschritt ** */
        zeitzaehler = zeitzaehler + dt;

        /* ** Schrittzähler, einkommentieren, wenn gewünscht ** */
        /*if(step % DruckFreq == 0)
        {
            cout << step << endl;
        }*/
    }
}

void save(double AnzTeilchen, double Dichte, double Temperatur, double dt, double dt2, int N, int DruckFreq, string fname)
{
    ofstream fout;
    fout.open(fname);  
    fout << "Schwerpunktsgeschwindigkeit v_schwpkt(t) / (m/s)" << "\t" << "Temperatur T / K" << "\t" << "E_pot(t)" << "\t" << "E_kin(t)" << "\t" << "Korrelationspaarfkt g(r)";

    simulate(AnzTeilchen, Dichte, Temperatur, dt, dt2, N, DruckFreq, fname);
    //fout << vschwmatrix;
    fout.close();
}


int main() {
    // *******************
    // *    Init         *
    // *******************

    // Konfig-Einstellungen
    int AnzTeilchen = 16;       // Anzahl an Teilchen
    double Dichte = 0.32;       // Dichte der Teilchen
    //int masse = 1;              // Masse der Teilchen (1 lt Aufgabe), in den jew. Berechnung 1 gesetzt
    string fname = "A1.dat";

    // Simulationsparameter
    double dt = 0.01;                       // Integrationszeit, entspricht h im Verlet dt = 0.01
    double dt2 = dt*dt;                     // (Integrationszeit)^2
    double Temperatur = 1.0;                // Temperatur für die Simulation, T(0) = 1 lt. Aufgabe c
    //double Temperatur2 = 0.01;              // Temperatur für die Simulation, T(0) = 0.01 lt. Aufgabe c
    //double Temperatur3 = 100.0;             // Temperatur für die Simulation, T(0) = 100 lt. Aufgabe c

    int N = 10;          // gesamte Simulationszeit (in Anzahl an Integrationsschritten), hier später 10^4 Schritte (0.01->1000000)
    int DruckFreq = 100;    // Druckfrequenz

    save(AnzTeilchen, Dichte, Temperatur, dt, dt2, N, DruckFreq, fname);
    //save(AnzTeilchen, Dichte, Temperatur2, dt, dt2, N, DruckFreq, fname);
    //save(AnzTeilchen, Dichte, Temperatur3, dt, dt2, N, DruckFreq, fname);
    return 0;
}
