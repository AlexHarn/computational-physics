//-------------------------------------------------------------------------------
//                          Computational Physics 2016
//                                Übungsblatt 6
//-------------------------------------------------------------------------------
//       Lösung der Poission Gleichung auf einem Quadrat mit konstanten
//       Randbedingungen und verschiedenen Ladungsverteilungen.
//-------------------------------------------------------------------------------
#include "poisson.h"

int main()
{
    PoissionRect square(1, 1, 0.05, 1e-5);

    // a)
    square.calc();
    square.save("a");

    // b)
    square.setConstBC(1, 0, 0, 0);
    square.calc();
    square.save("b");

    // c)
    square.reset();
    square.addQ(0.5, 0.5, 1);
    square.calc();
    square.save("c");

    // e)
    //  Dipol
    square.reset();
    square.addQ(0.25, 0.25, 1);
    square.addQ(0.75, 0.75, -1);
    square.calc();
    square.save("e-dipol");
    //  Quadrupol
    square.reset();
    square.addQ(0.25, 0.25, 1);
    square.addQ(0.25, 0.75, -1);
    square.addQ(0.75, 0.75, 1);
    square.addQ(0.75, 0.25, -1);
    square.calc();
    square.save("e-quadrupol");

    return 0;
}
