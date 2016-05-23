#include "../../Eigen/Eigen/Dense" //Eigen liegt oberhalb der B**-Ordner
#include "periodRB.h"

using namespace Eigen;
/**
  * Klasse untersucht, ob der kürzeste Weg zur nächsten Box oder in eigener Box liegt. Wird in ljforces.cpp und
  * mdsimulation.cpp aufgerufen.
  */
Vector2d PeriodRB::kurzerWeg(Vector2d dr, double L)
{

               // 1. Halbe Boxlänge in jede Richtung
                double hL = L/2.0;    

                /* 2. Abstandsvec sollte in [-hL_x, hL_x] & [-hL_y, hL_y] liegen
                 * (Weg zur Nachbarbox kürzer als in eigener Box)
                 * Sonst kürzeren Weg nehmen
                 */
                for(int dim = 0; dim < 2; dim++) {
                    if(dr(dim) > hL)
                    {
                        dr(dim) = dr(dim)-L;
                    }
                    else if (dr(dim) < -hL)
                    {
                        dr(dim) = dr(dim)+L;
                    }
                }   

                return dr;
}
