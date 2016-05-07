#include <functional>

#ifndef RUNGKUTT_H
#define RUNGKUTT_H

void rungkutt(std::function<void(double t, double* y, double* out)> f, const double N, const double h, const double t_0, double** y, const int d);
/* Führt das Runge-Kutta-Verfahren 4. Ordnung durch. Neue Implementierung ohne std:vector und in
 * allgemeiner Formulierung wie in der Vorlesung.
 *      DGL: y'(t) = f(t, y(t))
 *      N: Schrittzahl
 *      h: Schrittweite
 *      t_0: Startzeit
 *      y: y Vektoren zu allen Zeiten. Bei Aufruf müssen Starbedingungen gesetzt sein
 *      d: Länge der y Vektoren (!!! Nicht Anzahl der Raumdimensionen !!!)
 */
#endif
