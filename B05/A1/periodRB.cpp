#include <eigen3/Eigen/Dense>
#include "periodRB.h"

using namespace Eigen;
void umklapp(Vector2d &tempRB, double L)
{
    for ( int i = 0; i<2; i++ )
    {
        if ( tempRB(i) < -L )
        {
            tempRB(i) += L;
        }
        if ( tempRB(i) > L )
        {
            tempRB(i) -= L;
        }
    }
}


void kurzerWeg(Vector2d &dr, double L)
{
    /*  Nicht hübsch aber effektiv. Frage in einem Kreis alle 9 in Frage kommenden Positionen ab.
        Reihenfolge (gegen den UZS):
        9 8 7
        2 1 6
        3 4 5
    */
    const double rc = L/2.0;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(0) = dr(0) - L;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(1) = dr(1)-L;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(0) = dr(0)+L;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(0) = dr(0)+L;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(1) = dr(1)+L;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(1) = dr(1)+L;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(0) = dr(0)-L;
    if (dr.norm() <= rc)
    {
        return;
    }

    dr(0) = dr(0)-L;
    if (dr.norm() <= rc)
    {
        return;
    }
    // Wenn dr außerhalb liegt, gebe 0 0 zurück.
    dr << 0, 0;
    return;
}
