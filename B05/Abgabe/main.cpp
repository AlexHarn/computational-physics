#include <eigen3/Eigen/Core>
#include <cmath>
#include <iostream>
#include <chrono>
#include <fstream>
#include <string>
#include <random>
#include <math.h>
#include <eigen3/Eigen/Core>
#include <iostream>

using namespace Eigen;
using namespace std;

MatrixXd erstellequadrgitter(int AnzTeilchen, double L)
{
    //int AnzTeilchen Anzahl an Teilchen für das Gitter
    //double dichte Dichte der Partikel
    MatrixXd particleinfo(2, AnzTeilchen); //2xAnzTeilchen Array mit einem Ortsvektor pro Spalte

    double N = sqrt(AnzTeilchen);

    // Kord setzen:
    Vector2d tempTeil(0, 0);
    // wie gewünscht nach 1a, zaehler ist Spaltenzähler
    int zaehler = 0;
    for ( int n = 0; n<N; n++ )
    {
        for ( int m=0; m<N ;m++ )
        {
            tempTeil << (1+2*n)*L/8.0, (1+2*m)*L/8.0;
            for ( int komp=0; komp<2; komp++ )
            {
                particleinfo(komp, zaehler) = tempTeil(komp);
            }
            zaehler++;
        }
    }
    return particleinfo;
}

void umklapp(Vector2d &tempRB, double L)
{
    for ( int i = 0; i<2; i++ )
    {
        if ( tempRB(i) < 0 )
        {
            tempRB(i) += L;
        }
        if ( tempRB(i) > L )
        {
            tempRB(i) -= L;
        }
    }
}

void kurzerWeg(Vector2d &dr, double L)
{
    /*  Nicht hübsch aber effektiv. Frage in einem Kreis alle 9 in Frage kommenden Positionen ab.
        Reihenfolge (gegen den UZS):
        9 8 7
        2 1 6
        3 4 5
    */
    const double rc = L/2.0;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(0) = dr(0) - L;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(1) = dr(1)-L;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(0) = dr(0)+L;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(0) = dr(0)+L;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(1) = dr(1)+L;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(1) = dr(1)+L;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(0) = dr(0)-L;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(0) = dr(0)-L;
    if (dr.norm() <= rc)
    {
        return;
    }
    // Wenn dr außerhalb liegt, gebe 0 0 zurück.
    dr << 0, 0;
    return;
}

void kraft(MatrixXd &forces, MatrixXd &particleinfo, double L, bool active, VectorXd &bins, double &V)
{
	forces.setZero();
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
            kurzerWeg(dr, L);
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
                forces(komp, TeilA) = tempf1(komp);
                forces(komp, TeilB) = tempf2(komp);
            }
        }
    }
    //Faktor von der Kraft
    forces = forces*48;
}
#define CALI 100

/**
  * simulate(..) macht folgendes:
  * 1. Gitter erstellen (s. initializeSquareGrid.cpp)
  * 2. zufällige Geschwindigkeiten generieren
  * 3. Geschwindigkeiten über feste Temperatur neu skalieren
  * 4. Kraftberechnung
  * 5. Integrationsschritt (+Zähler, aktuell auskommentiert)
  **/

void simulate(int AnzTeilchen, double Temperatur, double dt, int N, string fname, int L, int bins, bool thermo)
{
    MatrixXd savedata(4, N); // Schwerpunktsgeschwindigket, Temperatur, potentielle Energie, kinetische Energie
	savedata.setZero();

    /****************************
     *      Initialisierung     *
     ****************************/

    // (Integrationszeit)^2
    double dt2 = dt*dt;

    // Erstelle Gitter
    MatrixXd particleInfo = erstellequadrgitter(AnzTeilchen, L);

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

    for ( int zeile = 0; zeile < AnzTeilchen; zeile++ )
    {
        for ( int i = 0; i<2; i++ )
        {
            //Schwerpunktsgeschwindigkeit abziehen, sonst Drift, in neue Matrix, um alle Vektoren mizunehmen
            vmatrix(i, zeile) -= vges(i);
        }

        //Init: Ekin = 1/2mv^2 = Nf kbT/2 (ÄquipartTheorem) , Nf = 2N-2, kb=1, m=1 (2 Raumdimensionen, 2 FG pro Schwerpunktsimpuls)
        savedata(3, 0) += 0.5*( pow(vmatrix(0, zeile), 2) + pow(vmatrix(1,zeile), 2) );
    }

    vges << 0,0;
	for ( int zeile = 0; zeile < AnzTeilchen; zeile++ )
    {
        for ( int i = 0; i<2; i++ )
        {
            vges(i) += vmatrix(i, zeile);
        }
    }
    vges = vges/AnzTeilchen;
    savedata(0, 0) = vges.norm();

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

	savedata(3, 0) = 0;
    for ( int zeile = 0; zeile < AnzTeilchen; zeile++ )
    {
        //Init: Ekin = 1/2mv^2 = Nf kbT/2 (ÄquipartTheorem) , Nf = 2N-2, kb=1, m=1 (2 Raumdimensionen, 2 FG pro Schwerpunktsimpuls)
        savedata(3, 0) += 0.5*( pow(vmatrix(0, zeile), 2) + pow(vmatrix(1,zeile), 2) );
    }
    savedata(1, 0) = 2*savedata(3, 0)/( 2*AnzTeilchen - 2 );

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
        particleInfo += dt*vmatrix + 0.5*dt2*forces;


        // Geschwindigkeiten aktualisieren, alle Kords auf einmal
        // Rest vom v-Verlet
        vmatrix += 0.5*dt*forces;

        /** Neue Kräfte berechnen **/
        kraft(forces, particleInfo, L, active, hist, savedata(2, step));

        vmatrix += 0.5*dt*forces;

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
		vges << 0, 0;
        for ( int zeile = 0; zeile < AnzTeilchen; zeile++ )
        {
            for ( int i = 0; i<2; i++ )
            {
                vges(i) += vmatrix(i, zeile);
            }
        }
        vges = vges/AnzTeilchen;
        savedata(0, step) = vges.norm();

        //if ( thermo )
        //{
            //// Skalierungsfaktor für die Geschwindigkeiten (lt. Aufgabe a) nach 5.40
            //double vskalierungfaktor = sqrt(Temperatur/T);
            //for ( int zeile = 0; zeile < AnzTeilchen; zeile++ )
            //{
                //for ( int i=0; i<2; i++ )
                //{
                    ////Geschwindigkeiten neu skalieren, um bestimmte E bei bestimmter T zu erzeugen (s. Kap 5.5.1).
                    //vmatrix(i, zeile) = vskalierungfaktor*vmatrix(i, zeile);
                //}
            //}
        //}
    }

    double Tm = 0;
    for ( int i = CALI; i<N; i++ )
    {
       Tm += savedata(1, i);
    }
    Tm /= N - CALI;

	double DA = 0;
    for (int i = 0; i < bins; i++)
    {
        DA = M_PI * ( pow((i+1.0)*(L/(2.0*bins)), 2)-pow((i)*L/(2.0*bins), 2) );
        hist(i) /= 0.5 * DA * ( N - CALI ) * pow(AnzTeilchen, 2) / pow(L, 2);
    }

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

    fout.open(fname+"_paar.dat");
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
    int bins = 100;

    // Simulationsparameter
    double dt = 0.01;                       // Integrationszeit, entspricht h im Verlet dt = 0.01

    int N = 10000;            // gesamte Simulationszeit (in Anzahl an Integrationsschritten), hier später 10^4 Schritte (0.01->1000000)

    #pragma omp parallel
    simulate(AnzTeilchen, 1, dt, N, "a1", L, bins, false);
    simulate(AnzTeilchen, 0.01, dt, N, "a2", L, bins, false);
    simulate(AnzTeilchen, 100, 0.0001, N*100, "a3", L, bins, false);
    //simulate(AnzTeilchen, 1, dt, N, "a1thermo", L, bins, true);
    //simulate(AnzTeilchen, 0.01, dt, N, "a2thermo", L, bins, true);
    //simulate(AnzTeilchen, 100, 0.0001, N, "a3thermo", L, bins, true);

    return 0;
}
