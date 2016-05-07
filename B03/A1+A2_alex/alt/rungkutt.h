#include <functional>
#include <vector>

using namespace std;

#ifndef RUNGKUTT_H
#define RUNGKUTT_H

void rungkutt(function<void(vector<double> &r, vector<double> &force)> F, double t_0, double t_N, double h, vector<vector<double>> &r, vector<vector<double>> &v);
/* Führt Runge-Kutta 4. Ordnung durch.
 *      F: Kraftfeld.
 *          r: Ortsvektor (input)
 *          force: Kraft am Ort r (output)
 *      t_0: Startzeit
 *      t_N: Endzeit
 *      h: Schrittweite
 *      r: Orte zu allen Zeiten
 *      Muss nicht mit richtiger Größe übergeben werden!
 *          r(0): Startwert
 *      v: Geschwindigkeiten zu allen Zeiten
 *      Muss nicht mit richtiger Größe übergeben werden!
 *          v(0): Startwert
 */

#endif
