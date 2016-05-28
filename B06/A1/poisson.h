#include <eigen3/Eigen/Core>

#ifndef POISSION_H
#define POISSION_H

typedef Eigen::MatrixXd M;

class PoissionSquare
/* Löst die Poissiongleichung auf einem Quadrat mit gegebenen Maßen,
 * Randbedingungen und gegebener Ladungsverteilung.
 */
{
    private:
        M rho, phi, ex, ey;
        double delta, eps;

        void calcP();
        /* Berechnet das Potential mittels Gauß-Seidel-Iteration
         * und Genauigkeit eps
         */
        void calcE();
        /* Berechnet das E-Feld aus dem Potential
         */

    public:
        PoissionSquare(double lx, double ly, double delta, double eps);
        /* Initiallisiert ein Quadrat [0, lx] x [0, ly] mit Diskretisierung delta
         * und Genauigkeit eps
         */
        void calc();
        /* Wrapper Methode für alle Berechnungen
         */
        void setConstBC(double right, double top, double left, double bottom);
        /* Setzt die RB an den 4 Rändern auf die jeweils gegebene Konstante.
         * Zu Beachten: Die Reihenfolge ist eine Prioritätenliste, die Ecken werden
         * entsprechend überschrieben!
         */
        void addQ(double x, double y, double Q);
        /* Fügt eine Ladung am Ort (x, y) mit Betrag Q
         */
        void save(std::string name);
        /* Speichert die aktuellen Ergebnisse in Dateien mit Namen "[name]_[postfix].dat"
         */
        void reset();
        /* Setzt das Quadrat zurück
         */
};

#endif
