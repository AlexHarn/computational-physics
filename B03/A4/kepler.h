#include <vector>
#include <functional>

using namespace std; // Unsauber aber hier egal

#ifndef KEPLER_H
#define KEPLER_H
//-------------------------------------------------------------------------------
//                              Globale Parameter
//-------------------------------------------------------------------------------

const double G = 1; // Gravitationskonstante
const double m = 1; // Teilchenmasse

//-------------------------------------------------------------------------------

void energy(function<double(vector<double> &r)> V, vector<vector<double>> &r, vector<vector<double>> &v, vector<double> &U, vector<double> &T);
/* Berechnet kinetische Energie T und potentielle energie U
 *      F: Kraftfeld
 *      V: Potential
 *      r: Ortsvektoren
 *      v: Geschwindigkeitsvektoren
 *      U: Potentielle Energie an Orten r mit v (output)
 *      T: Kinetische Energie an Orten r mit v (output)
 */

void savedat(string dest, function<double(vector<double> &r)> V,  double h, vector<vector<double>> &r, vector<vector<double>> &v);
/* Speichermethode
 */

void doWork(double alpha, double h, int t, vector<double> r_0, vector<double> v_0, string fname);
/* Führt einen Durchgang mit gegebenen Parametern durch und speichert das
 * Ergebnis.
 *      alpha: alpha im Potential/Kraftfeld
 *      h: Schrittweite
 *      t: Zu simmulierende Gesamtzeit
 *      r_0: Startwert für den Ort
 *      v_0: Startwert für die Geschwindigkeit
 *      fnmae: Name der Zieldatei
 */

void forwNback(double alpha, double h, int N, vector<double> r_0, vector<double> v_0, string fname);
/* Führt für n = 1e3, 2e3, ..., N Vorwärts-Rückwärts-Tests mit gegebenen Parametern durch und speichert das
 * Ergebnis.
 *      alpha: alpha im Potential/Kraftfeld
 *      h: Schrittweite
 *      N: maximal Durchzuführende Schritte
 *      r_0: Startwert für den Ort
 *      v_0: Startwert für die Geschwindigkeit
 *      fnmae: Name der Zieldatei
 */

void F(vector<double> &r, vector<double> &force, double alpha);
/* Kraftfeld des Keppler Problems
 *      r: Ortsvektor (input)
 *      force: Kraft am Ort r (output)
 */

double vabs(vector<double> &v);
/* Berechnet |v|
 */

double V(vector<double> &r, double alpha);
/* Potential am Ort r
 * Ja man könnte das auch einfach differenzieren und sich F sparen, wenn man
 * beides aber analytisch kennt sollte man das auch lieber nutzen.
 */

void cross(vector<double> a, vector<double> b, vector<double> &result);
/* Berechnet Kreuzprodukt result = a x b
 * NUR FÜR D=3 !!!
 */

void angm(vector<vector<double>> &r, vector<vector<double>> &v, vector<double> &L);
/* Berechnet den Betrag des Drehimpulses L = | r x p |
 */

void area(vector<vector<double>> &r, vector<double> &A);
/* Berechnet für jeden Zeitpunkt die INSGESAMT überstrichene Fläche. Indem die
 * Teilflächen durch Dreiecke genähert werden, dessen Flächen nach Satz des
 * Herons berechnet werden.
 *      r: r[t]
 *      A: A[t] (Insgesamt)
 */

void lenzrunge(vector<vector<double>> &r, vector<vector<double>> &v, vector<vector<double>> &LR);
/* Berechnet den Lenz-Runge Vektor für alle Zeiten.
 *      r: r[t]
 *      v: v[t]
 *      LR: Lenz-Runge Vektoren (output)
 */
#endif
