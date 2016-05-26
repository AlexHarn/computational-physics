#include <eigen3/Eigen/Core>
#ifndef LJFORCES_H
#define LJFORCES_H
void kraft(Eigen::MatrixXd &forces, Eigen::MatrixXd &particleinfo, double L, bool active, Eigen::VectorXd &bins, double &V);
/*
  * Berechnet die Kräfte auf ein Ensemble von LJ-Teilchen in 2D
  * Gibt die Kräfte auf alle Teilchen in einer bestimmten Ordnung zurück
  * sigma (Längenskala) = epsilon (Energieskala) = 1
  * Berechnet Dir auch Epot in jedem Zeitschritt
  * Ach ja und das Histogramm liefern wir auch gleich mti!
  */
#endif
