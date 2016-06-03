#include <eigen3/Eigen/Core>

#ifndef LINGEN_H
#define LINGEN_H

typedef Eigen::VectorXd V;

class LinConGen
/* 
 */
{
    private:
        int64_t rold,rnew;

        void calccongruentgen(int64_t r0, int64_t a, int64_t c, int64_t m, int64_t M);
        /* linear kongruenter Generator r_(n+1) = (a*r_n+c)%m
         */

        //void sethistogram(int totalbins, double lengthpbin);

    public:
        V saveVec, hist;
        int64_t a,c,r0,m, N, M;
        int totalbins, lengthpbin;
        // M ist eine Hilfsvariable, die sich ausgehend von N dznamisch auf maximal 10^4 oder kleiner setzt.

        LinConGen(int savedvalues, int bins, double lengthbins);
        /* 
         */
        void calc(int64_t r0, int64_t a, int64_t c, int64_t m);
        /* Wrapper Methode für alle Berechnungen
         */
        void save(std::string name);
        /* Speichert die aktuellen Ergebnisse in Dateien mit Namen "[name]_[postfix].dat"
         */
        void reset();
        /* Setzt zurück
         */
};

#endif
