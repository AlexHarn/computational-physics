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

void simulate(int AnzTeilchen, double  dichte, double Temperatur, double dt, int N, int L, int DruckFreq, string fname, int bins, double Dr) {

    /****************************
     *      Initialisierung     *
     ****************************/

    // (Integrationszeit)^2
    double dt2 = dt*dt;        
    // Erstelle Gitter
    MatrixXd particleInfo(2, AnzTeilchen);
    Quadratgitter quadratgitter;
    particleInfo = quadratgitter.erstellequadrgitter(AnzTeilchen);

    //aktueller Unix-Timestamp für rand-Nummern
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> distributionR(0.0, 1.0);
    std::uniform_real_distribution<double> distributionPhi(0.0, 2*M_PI);

    // Setze Anfangsgeschwindigkeit mit zufälligen Nummern
    MatrixXd vmatrix(2, AnzTeilchen); //enthält alle Geschwindigkeitsvektoren
    Vector2d vges(0,0); //Schwerpunktsgeschwindigkeit, init mit 0
    double r = 0;
    double phi = 0;
    Vector2d v(0,0);
    double Ekin = 0;

    for (int zeile = 0; zeile < AnzTeilchen; zeile++)
    {
        r = distributionR(generator);
        phi = distributionPhi(generator);
        v << r*cos(phi), r*sin(phi);
        vges += v;
        for (int i = 0; i<2; i++)
        {
            vmatrix(i, zeile) = v(i);
        }
    }
    vges = vges/AnzTeilchen;
    for (int zeile = 0; zeile < AnzTeilchen; zeile++)
    {
        for(int i = 0; i<2; i++)
        {
        /*
        * x und y, hier hätte man sicher auch direkt die Vektoren als ganze addieren können, also Zeile 33 aus der for-schleife, dann "<<" nutzen. Allerdings hat das mit der Dimension 2xAnzTeil-Matrix vs. 2-dim Vektor nicht so einfach geklappt. Das umgeht das Problem.
        */
            //Schwerpunktsgeschwindigkeit abziehen, sonst Drift, in neue Matrix, um alle Vektoren mizunehmen
            vmatrix(i,zeile) = (vmatrix(i,zeile)-vges(i));
        }

        //Init: Ekin = 1/2mv^2 = Nf kbT/2 (ÄquipartTheorem) , Nf = 2N-2, kb=1, m=1 (2 Raumdimensionen, 2 FG pro Schwerpunktsimpuls)
        Ekin = Ekin + 0.5*(vmatrix(0,zeile)*vmatrix(0,zeile)+vmatrix(1,zeile)*vmatrix(1,zeile));
        // Epot wird in ljforces bestimmt
    }
    
    double T = 2*Ekin/(2*AnzTeilchen-2);    

    /* ** Isokinetischer Thermostat ** */
    // Skalierungsfaktor für die Geschwindigkeiten (lt. Aufgabe a) nach 5.40
    double vskalierungfaktor = sqrt(Temperatur/T);
    for (int zeile = 0; zeile < AnzTeilchen; zeile++)
    {
        for(int i=0; i<2; i++) //für 2 dim
        {
            //Geschwindigkeiten neu skalieren, um bestimmte E bei bestimmter T zu erzeugen (s. Kap 5.5.1).
            vmatrix(i,zeile) = vskalierungfaktor*vmatrix(i,zeile);
        }
    }


    /****************************
     *     Kraftberechnung      *
     ****************************/


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
    // T, g(r) messen für verschiedene T0. Phasen erkennen. -< Plots
    // Aufg e: Was passiert nun bei T = 0.01? Schauen Sie sich
    // die potentielle, kinetische und Gesamtenergie während der Äquilibrierung an und messen
    // Sie g(r) an dem äquilibrierten System. Es kann auch interessant sein sich Snapshots, also
    //Bilder der Konfiguration des Systems, anzuschauen. -> Plots der Energien, Paarkorrfkt.

    /****************************
     * kanonische MD-Simulation *
     ****************************/

    double zeitzaehler = 0;
    bool active = false;
    bool mittelung = false;
    double sumanzahl = 0;
    
    for(int step=0; step < N; step++)
    {
        /* ** Erster Integrationsschritt **  */
        
        // Aktualisiere Positionen
        //r(t)..
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
            tempRB = periodRB.umklapp(tempRB, L);

            for(int zeile = 0; zeile < 2; zeile++)
            {
                particleInfo(zeile, Teil) = tempRB(zeile);
            }
        }
        
        // Geschwindigkeiten aktualisieren, alle Kords auf einmal
        // Rest vom v-Verlet
        vmatrix = vmatrix + 0.5*dt*forces;
    
        /** Neue Kräfte berechnen **/        
        forces = ljforces.kraft(particleInfo, L);
        
        /* ** Zweiter Integrationsschritt ** */

        // Geschwindigkeiten aktualisieren, v-Verlet.. TODO Darf das forces sein?
        vmatrix = vmatrix + 0.5*dt*forces;

        /* ** Zeitschritt ** */
        zeitzaehler = zeitzaehler + dt;

        /* ** Schrittzähler, praktisch für Zeitmittelung, einkommentieren, wenn gewünscht ** */
        if(step % DruckFreq == 0)
        {
            cout << step << endl;
            mittelung = true;
        }
        else
        {
            mittelung = false;
        }

        /* ** Histogramm ** */
        VectorXd anzahl(bins);
        for(int part = 0; part < AnzTeilchen; part++)
        {
            for(int veckomp = 0; veckomp < anzahl.size(); veckomp++)
            {   
                for(int komp; komp < 2; komp++)
                {
                    if(forces(komp,part) < (veckomp+1)*L/(2*anzahl.size()))  
                    {   
                        anzahl(veckomp) += 1;
                        continue;
                    }
                }
            }
        }
        // Anzahl der Paare pro Zeitschritt addieren
        sumanzahl += anzahl.sum();

        /* ** Warte Äquilibrierungsphase ab ** */
        if(step = 10)//10000 nach c
        {
            active = true;
        }
        /* Differenzfläche zwischen Kreis mit r = r_i und Kreis mit r = r_a,
         * Delta A = PI*((l*Delta r)^2-(l-1*Delta r)^2) :: S.85 im Skript */
        double DA = M_PI*(pow((bins*Dr),2)-pow(((bins-1)*Dr),2));;
        double g = 0;

        if(active && mittelung)
        {
            //Messe Temperatur und g(r)
            for (int Teil = 0; Teil < AnzTeilchen; Teil++)
            {
                Ekin = Ekin + 0.5*(vmatrix(0,Teil)*vmatrix(0,Teil)+vmatrix(1,Teil)*vmatrix(1,Teil)); 
                // Epot wird in ljforces bestimmt
            }
            T = 2*Ekin/(2*AnzTeilchen-2);  

            //Korrfkt g(r)
            g = sumanzahl/(dichte*AnzTeilchen*DA);
        }

        if(step = N)
        {
            active = false;
        }
    }
}

void save(double AnzTeilchen, double dichte, double Temperatur, double dt, int N, int L, int DruckFreq, string fname, int bins, double Dr)
{
    ofstream fout;
    fout.open(fname);  
    fout << "Schwerpunktsgeschwindigkeit v_schwpkt(t) / (m/s)" << "\t" << "Temperatur T / K" << "\t" << "E_pot(t)" << "\t" << "E_kin(t)" << "\t" << "Korrelationspaarfkt g(r)" << endl;

    simulate(AnzTeilchen, dichte, Temperatur, dt, N, L, DruckFreq, fname, bins, Dr);
    VectorXd anzahl(bins);
    for(int veckomp = 0; veckomp < anzahl.size(); veckomp++)
    {
        fout << (veckomp+0.5)*L/(2*anzahl.size()) << "\t" << anzahl(veckomp) << endl;
    }
    //fout << vschwmatrix;
    fout.close();
}


int main() {
    // *******************
    // *    Init         *
    // *******************

    // Konfig-Einstellungen
    int AnzTeilchen = 16;       // Anzahl an Teilchen
    //int masse = 1;              // Masse der Teilchen (1 lt Aufgabe), in den jew. Berechnung 1 gesetzt
    double dichte = 0.32;
    string fname = "A1.dat";
    double L = 8.0; 
    int bins = 15;
    double Dr = L / 2;

    // Simulationsparameter
    double dt = 0.01;                       // Integrationszeit, entspricht h im Verlet dt = 0.01
    double Temperatur = 1.0;                // Temperatur für die Simulation, T(0) = 1 lt. Aufgabe c
    //double Temperatur2 = 0.01;              // Temperatur für die Simulation, T(0) = 0.01 lt. Aufgabe c
    //double Temperatur3 = 100.0;             // Temperatur für die Simulation, T(0) = 100 lt. Aufgabe c

    int N = 100;          // gesamte Simulationszeit (in Anzahl an Integrationsschritten), hier später 10^4 Schritte (0.01->1000000)
    int DruckFreq = 100;    // Druckfrequenz

    save(AnzTeilchen, dichte, Temperatur, dt, N, L, DruckFreq, fname, bins, Dr);
    //save(AnzTeilchen, dichte, Temperatur2, dt, N, L, DruckFreq, fname, bins, Dr);
    //save(AnzTeilchen, dichte, Temperatur3, dt, N, L, DruckFreq, fname, bins, Dr);
    return 0;
}
