//-------------------------------------------------------------------------------
//                         Computational Physics 2016
//                               Übungsblatt 7
//-------------------------------------------------------------------------------
//            Implementierung eines linearen kongruenten Generators
//            und einiger Verteilunge.
//-------------------------------------------------------------------------------
#include "lingen.h"

int main()
{
    LinConGen gen;

    //// Aufgabe 1
    //  b)
    //  (i)
    gen.congruent(1234, 20, 120, 6075, 1e5);
    gen.save("1b1");

    //  (ii)
    gen.congruent(1234, 137, 187, 256, 1e5);
    gen.save("1b2");

    //  (iii)
    gen.congruent(123456789, 65539, 0, 2147483648, 1e5);
    gen.save("1b3");

    ////  (iv)
    gen.congruent(1234, 16807, 0, 2147483647, 1e5);
    gen.save("1b4");


    // Aufgabe 2
    //  a)
    gen.boxMuller();
    gen.saveDist("2a");

    // b)
    gen.congruent(1234, 16807, 0, 2147483647, 12e5);
    gen.centralLimit();
    gen.saveDist("2b");

    return 0;
}
