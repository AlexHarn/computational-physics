#include <eigen3/Eigen/Core>
#include <math.h>
#include <iostream>
#include <chrono>
#include <fstream>
#include <string>
#include <random>
#include "initializeSquareGrid.h"
#include "ljforces.h"
#include "periodRB.h"

#define CALI 10

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

void simulate(int AnzTeilchen, double Temperatur, double dt, int N, string fname, int L, int bins)
{
    MatrixXd savedata(4, N); // Schwerpunktsgeschwindigket, Temperatur, potentielle Energie, kinetische Energie
    savedata(1, 0) = Temperatur;

    /****************************
     *      Initialisierung     *
     ****************************/

    // (Integrationszeit)^2
    double dt2 = dt*dt;
    // Erstelle Gitter
    MatrixXd particleInfo(2, AnzTeilchen);
    particleInfo = erstellequadrgitter(AnzTeilchen);

    //aktueller Unix-Timestamp für rand-Nummern
    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 generator(seed);
    uniform_real_distribution<double> distributionR(0.0, 1.0);
    uniform_real_distribution<double> distributionPhi(0.0, 2*M_PI);

    // Setze Anfangsgeschwindigkeit mit zufälligen Nummern
    MatrixXd vmatrix(2, AnzTeilchen); //enthält alle Geschwindigkeitsvektoren
    Vector2d vges(0,0); //Schwerpunktsgeschwindigkeit, init mit 0
    double r = 0;
    double phi = 0;
    Vector2d v(0, 0);

    for ( int zeile = 0; zeile < AnzTeilchen; zeile++ )
    {
        r = distributionR(generator);
        phi = distributionPhi(generator);
        v << r*cos(phi), r*sin(phi);
        vges += v;
        for ( int i = 0; i<2; i++ )
        {
            vmatrix(i, zeile) = v(i);
        }
    }

    vges = vges/AnzTeilchen;
    savedata(0, 0) = 0;
        
    for ( int zeile = 0; zeile < AnzTeilchen; zeile++ )
    {
        for ( int i = 0; i<2; i++ )
        {
            //Schwerpunktsgeschwindigkeit abziehen, sonst Drift, in neue Matrix, um alle Vektoren mizunehmen
            vmatrix(i, zeile) = vmatrix(i, zeile) - vges(i);
        }

        //Init: Ekin = 1/2mv^2 = Nf kbT/2 (ÄquipartTheorem) , Nf = 2N-2, kb=1, m=1 (2 Raumdimensionen, 2 FG pro Schwerpunktsimpuls)
        savedata(3, 0) += 0.5*( pow(vmatrix(0, zeile), 2) + pow(vmatrix(1,zeile), 2) );
    }

    double T = 2*savedata(3, 0)/( 2*AnzTeilchen - 2 );


    // Skalierungsfaktor für die Geschwindigkeiten (lt. Aufgabe a) nach 5.40
    double vskalierungfaktor = sqrt(Temperatur/T);
    for ( int zeile = 0; zeile < AnzTeilchen; zeile++ )
    {
        for ( int i=0; i<2; i++ )
        {
            //Geschwindigkeiten neu skalieren, um bestimmte E bei bestimmter T zu erzeugen (s. Kap 5.5.1).
            vmatrix(i, zeile) = vskalierungfaktor*vmatrix(i, zeile);
        }
    }

    for ( int zeile = 0; zeile < AnzTeilchen; zeile++ )
    {
        //Init: Ekin = 1/2mv^2 = Nf kbT/2 (ÄquipartTheorem) , Nf = 2N-2, kb=1, m=1 (2 Raumdimensionen, 2 FG pro Schwerpunktsimpuls)
        savedata(3, 0) += 0.5*( pow(vmatrix(0, zeile), 2) + pow(vmatrix(1,zeile), 2) );
    }

    /*
     * Histogramm (fehlt), Epot spinnt noch rum, V wird in ljforces.cpp einfach riesig..
     * Äquiliebr.zeit bestimmen: Dann, wenn Ekin ~= const
     * erste Idee: umgeht double-Vergleich, indem beide doubles auf int gecastet werden. Dann sind mir Nachkommastellen egal
     * geht das schöner?
     * if ((int)Ekin == (int)Epot)
     *    cout << "Äquizeit" << zeitzaehler << endl;
     * TODO: Ekin, Epot, Temperatur in .dat-File eintragen <-> Plots
     * Nach wie vielen Zeitschritten äquilib. ihr System?
     * TODO g(r) implementieren
     * T, g(r) messen für verschiedene T0. Phasen erkennen. -< Plots
     * Aufg e: Was passiert nun bei T = 0.01? Schauen Sie sich
     * die potentielle, kinetische und Gesamtenergie während der Äquilibrierung an und messen
     * Sie g(r) an dem äquilibrierten System. Es kann auch interessant sein sich Snapshots, also
     * Bilder der Konfiguration des Systems, anzuschauen. -> Plots der Energien, Paarkorrfkt.
     */

    /****************************
     * kanonische MD-Simulation *
     ****************************/

    double zeitzaehler = 0;
    bool active = false;
    VectorXd hist(bins);

    //Jetzt Kraftberechnung (s. ljforces.cpp)
    MatrixXd forces(2, AnzTeilchen);
    kraft(forces, particleInfo, L, false, hist, savedata(2, 0));

    for ( int step=1; step < N; step++ )
    {
        /* ** Erster Integrationsschritt **  */

        // Aktualisiere Positionen
        //r(t)..
        //Geschwindigkeits-Verlet-Algorithmus nach (4.25), h = dt
        particleInfo = particleInfo + dt*vmatrix + 0.5*dt2*forces;

        // RB beachten
        for (int Teil = 0; Teil < AnzTeilchen; Teil++)
        {
            //Init eines Temp-Vektors mit Teilkords
            Vector2d tempRB(0, 0);
            tempRB = particleInfo.block(0,Teil,2,1);

            //RB beachten
            umklapp(tempRB, L);

            for (int zeile = 0; zeile < 2; zeile++)
            {
                particleInfo(zeile, Teil) = tempRB(zeile);
            }
        }

        // Geschwindigkeiten aktualisieren, alle Kords auf einmal
        // Rest vom v-Verlet
        vmatrix += 0.5*dt*forces;

        /** Neue Kräfte berechnen **/
        kraft(forces, particleInfo, L, active, hist, savedata(2, step));

        vmatrix += 0.5*dt*forces;

        /* ** Zeitschritt ** */
        zeitzaehler += dt;

        /* ** Warte Äquilibrierungsphase ab ** */
        if (step == CALI)
        {
            active = true;
        }

        //Messe Temperatur
        for (int Teil = 0; Teil < AnzTeilchen; Teil++)
        {
            savedata(3, step) += 0.5*(pow(vmatrix(0,Teil), 2)+pow(vmatrix(1,Teil), 2));
        }
        savedata(1, step) = 2*savedata(3, step)/(2*AnzTeilchen-2);

        for ( int zeile = 0; zeile < AnzTeilchen; zeile++ )
        {
            for ( int i = 0; i<2; i++ )
            {
                vges(i) += vmatrix(i, zeile);
            }
        }
        vges = vges/AnzTeilchen;
        savedata(0, step) = vges.norm();
    }
    double Tm = 0;
    for ( int i = CALI; i<N; i++ )
    {
       Tm += savedata(1, i); 
    }
    Tm /= N - CALI;
    /*
     *double DA = 0;
     *for(int i = 0; i< bins;i++)
     *{
     *    DA = DA + M_PI*(pow((i*Dr),2)-pow(((i-1)*Dr),2));
     *}
     *double g = 0;
     */
    ////Korrfkt g(r)
    //g = sumanzahl/(dichte*AnzTeilchen*DA);
    ofstream fout;
    fout.open(fname+"_savedata.dat");
    fout << "#t\tSchwerpunktsgeschwindigkeit\tTemperatur\tE_pot\tE_kin" << endl;
    for ( int i = 0; i<N; i++ )
    {
        fout << i*dt << "\t";
        for ( int j = 0; j<4; j++ )
            fout << savedata(j, i) << "\t";
        fout << endl;
    }
    fout.close();
    
    fout.open(fname+"_hist.dat");
    fout << "#Abstand\tAnzahl" << endl;
    for ( int i = 0; i<bins; i++ )
    {
        fout << (i+0.5)*L/( 2*bins ) << "\t" << hist(i) << endl;
    }
    fout.close();
}

int main()
{
    int AnzTeilchen = 16;         // Anzahl an Teilchen
    double L = 8.0;
    int bins = 15;

    // Simulationsparameter
    double dt = 0.01;                       // Integrationszeit, entspricht h im Verlet dt = 0.01
    //double Temperatur = 0.01;              // Temperatur für die Simulation, T(0) = 0.01 lt. Aufgabe c
    //double Temperatur = 100.0;             // Temperatur für die Simulation, T(0) = 100 lt. Aufgabe c

    int N = 100;            // gesamte Simulationszeit (in Anzahl an Integrationsschritten), hier später 10^4 Schritte (0.01->1000000)

    simulate(AnzTeilchen, 1, dt, N, "test", L, bins);
    return 0;
}
