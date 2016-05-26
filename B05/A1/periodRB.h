#include <eigen3/Eigen/Core>
#ifndef PERIODRB_H
#define PERIODRB_H
void umklapp(Eigen::Vector2d &tempRB, double L);
/*
 * Klappt den Vektor tempRB um.
 */
void kurzerWeg(Eigen::Vector2d &dr, double L);
/**
  * untersucht, ob der kürzeste Weg zur nächsten Box oder in eigener Box liegt. Wird in ljforces.cpp und
  * mdsimulation.cpp aufgerufen.
  */
#endif
