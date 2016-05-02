//-------------------------------------------------------------------------------
//                          Computational Physics 2016
//                               Blatt 3 Aufgabe 4
//-------------------------------------------------------------------------------
//               Anwendung des Runge Kutta Verfahrens 4. Ordnung
//               auf das Keppler Problem.
//-------------------------------------------------------------------------------
#include <sstream>
#include <cmath>
#include "rungkutt.h"
#include "kepler.h"

using namespace std;

int main()
{
    #pragma omp parallel // wirkt Wunder!
    // ACHTUNG: Je nach Hardware muss man hier ordentlich was anpassen
    // Die Datei für h = 1e-6 ist ca 1,5 GB groß

    // Verschiedene Schrittweiten an Ellipsen testen
    #pragma omp for
    for ( int i = 1;  i < 6; i++ )
    {
        double h = pow(10, -i);
        stringstream ss;
        ss << h;
        doWork(1, h, 10, {1, 0, 0}, {0, 0.5, 0}, "h="+ss.str()+".dat");
    }
    // Damit sich die Ellipse exakt schließt braucht man wirklich seeehr kleine
    // Schrittweiten aber 1e-6 ist schon sehr gut und 1e-5 ganz ok

    ////  Kreis
    doWork(1, 1e-5, 10, {1, 0, 0}, {0, 1, 0}, "kreis.dat");

    // Kleine Anfangsgeschwindigkeit: Problem ist die Divergenz im Ursprung
    doWork(1, 1e-5, 5, {1, 0, 0}, {0, 0, 0}, "kleines_v.dat");


    // Forwärts/Rückwärts tests durchführen
    forwNback(1, 1e-5, 10e4, {1, 0, 0}, {0, 0.5, 0}, "forwback.dat1");

    // alpha = 1.1 und 0.9 ausprobieren
    doWork(1.1, 1e-5, 10, {1, 0, 0}, {0, 0.5, 0}, "alpha=1.1.dat");
    doWork(1.1, 1e-5, 10, {1, 0, 0}, {0, 0.5, 0}, "alpha=0.9.dat");

    return 0;
}

