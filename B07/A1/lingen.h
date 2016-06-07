#include <eigen3/Eigen/Core>

#ifndef LINGEN_H
#define LINGEN_H

typedef Eigen::VectorXd V;
typedef Eigen::Matrix<int64_t, Eigen::Dynamic, 1> V64;

class LinConGen
/* linear kongruenter Generator r_(n+1) = (a*r_n+c)%m
 */
{
    private:
        int bins;
    public:
        V r, x, y;
        // M ist eine Hilfsvariable, die sich ausgehend von N dznamisch auf maximal 10^4 oder kleiner setzt.

        LinConGen(int bins);
        /* Konstruktor
         */
        void calcCongruent(int64_t r0, int64_t a, int64_t c, int64_t m, int64_t N);
        /* Generiert N linear kongurente Zufallszahlen
         */
        void boxmulleralg();

        void saveHist(std::string name);
        /* Speichert das Histogramm zum aktuellen r
         */
};

#endif
