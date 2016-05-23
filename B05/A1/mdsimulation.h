#include <string>

#ifndef MDSimulation_H
#define MDSimulation_H
    //Methoden
class MDSimulation
{
    public:
        // Speichert alle Observablen
        void save(double AnzTeilchen, double Dichte, double Temperatur, double dt, double dt2, int N, int DruckFreq, std::string fname);

        /** MD-Simulation
          * 1. Gitter erstellen (s. initializeSquareGrid.cpp)
          * 2. zufällige Geschwindigkeiten generieren
          * 3. Geschwindigkeiten über feste Temperatur neu skalieren
          * 4. Kraftberechnung
          * 5. Integrationsschritt (+Zähler, aktuell auskommentiert)
          **/
        void simulate(int AnzTeilchen, double Dichte, double Temperatur, double dt, double dt2, int N, int DruckFreq, std::string fname);
};
#endif
