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
        V r, dist;

        void congruent(int64_t r0, unsigned int a, unsigned int c, unsigned int m, unsigned int N);
        /* Generiert N linear kongurente Zufallszahlen
         */
        void boxMuller();
        /* Berechnet eine Gaußverteilung mittels Box-Muller-Algorithmus
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
