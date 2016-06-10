#include <eigen3/Eigen/Core>

#ifndef LINGEN_H
#define LINGEN_H

typedef Eigen::VectorXd V;
typedef Eigen::Matrix<int64_t, Eigen::Dynamic, 1> V64;

class LinConGen
/* linear kongruenter Generator r_(n+1) = (a*r_n+c)%m
 */
{
    public:
        V r, dist, tempdist;
        /* r: Zuletzt mit boxMuller() generierte Zufallszahlen
         * dist: Zuletzt generierte Verteilung
         */

        void congruent(int64_t r0, unsigned int a, unsigned int c, unsigned int m, unsigned int N);
        /* Generiert N linear kongurente Zufallszahlen
         */
        void boxMuller();
        /* Berechnet eine Gaußverteilung mittels Box-Muller-Algorithmus
         * Erwartungswert 0, Standardabweichung 1
         */
        void centralLimit(int N=12);
        /* Berechnet eine Gaußverteilung mittels zentralem Grenzwertsatz
         * Erwartungswert 0
         * Mit N = 12: Standardabweichung 1
         */
        void neumann(unsigned int N);
        /* Berechnet eine sin(x)/2 - Verteilung nach dem von
         * von Neumannschen Rückweisungsverfahren
         */
        void transform(unsigned int N);
        /* Berechnet eine 3x^2 - Verteilung aus gleichverteilten
         * Zufallszahlen
         */
        void save(std::string name);
        /* Speichert die zuletzt generierten Zufallszahlen
         */
        void saveDist(std::string name);
        /* Speichert die zuletzt generierte Verteilung
         */
        void saveHist(std::string name, int bins);
        /* Speichert das Histogramm zu den zuletzt generierten Zufallszahlen
         * Nur noch enthalten, falls explizit gefordert. Die Histogramme werden aber mit matplotlib aus
         * den Zufallszahlen direkt erstellt.
         */
};

#endif
