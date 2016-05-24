#include <string>

#ifndef MDSimulation_H
#define MDSimulation_H
    //Methoden
class MDSimulation
{
    public:
        // Speichert alle Observablen
        void save(double AnzTeilchen, double dichte, double Temperatur, double dt, int N, int L, int DruckFreq, std::string fname, int bins, double Dr);

        /** MD-Simulation
          * 1. Gitter erstellen (s. initializeSquareGrid.cpp)
          * 2. zufällige Geschwindigkeiten generieren
          * 3. Geschwindigkeiten über feste Temperatur neu skalieren
          * 4. Kraftberechnung
          * 5. Integrationsschritt (+Zähler, aktuell auskommentiert)
          **/
        void simulate(int AnzTeilchen, double dichte, double Temperatur, double dt, int N, int L, int DruckFreq, std::string fname, int bins, double Dr);
};
#endif
