#include <string>

#ifndef DPENDULUM_H
#define DPENDULUM_H

class Dpendulum
{
        void static f(double /*t*/, double* y, double* dy, double g);
        /* f für Runge-Kutta mit y' = f(y)
         *      y: y vektor (input)
         *      dy: y' vektor (output)
         *      g: Gravitationskonstante (wird später gebunden)
         */
        void calcE();
        /* Berechnet U und T aus y
         */
        void reset();
        /* Setzt das Pendel zurück (löscht die aktuellen Daten)
         */
        const double g = 9.81;
        int N = 0;
        double tN = 0, lh = 0;
        double* y0;
        double** y, ** T, ** U;
        bool ready = false, energy = false, active = false;

    public:
        ~Dpendulum();
        /* Destruktor
         */
        void setInitial(double theta1, double theta2, double ddtTheta1, double ddtTheta2);
        /* Setzt die AB
         */
        void swing(double h, double t);
        /* Berechnet Schwingung mit Schrittweite h bis Zeitpunkt t
         */
        void swing(double h, int n);
        /* Berechnet Schwingung mit Schrittweite h und Schrittzahl n
         */
        void save(std::string fname);
        /* Speichert die Daten in Datei fname
         */
        void doEverything(double theta1, double theta2, double ddtTheta1, double ddtTheta2, double h, double t, std::string fname);
        /* Macht alles in einem Aufruf
         */
};

#endif
