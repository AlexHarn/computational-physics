//-------------------------------------------------------------------------------
//                      Computational Physics 2016 Blatt 4
//                           Aufgabe 1: Doppelpendel
//-------------------------------------------------------------------------------
//              Anwendung des Runge-Kutta-Verfahren 4. Ordnung auf das
//              nicht lineare, chaotische Problem des Doppelpendels.
//-------------------------------------------------------------------------------
#include <cmath>
#include "Dpendulum.h"

int main()
{
    Dpendulum pendulum;
    pendulum.doEverything(0.1, sqrt(2)*0.1, 0, 0, 1e-5, 3, "pos_sqrt2.dat");
    pendulum.doEverything(0.1, -sqrt(2)*0.1, 0, 0, 1e-5, 3, "neg_sqrt2.dat");
}
