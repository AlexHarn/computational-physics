//-------------------------------------------------------------------------------
//                          Computational Physics 2016
//                           Blatt 3 Aufgaben 1 und 2
//-------------------------------------------------------------------------------
//        Implementierung des Runge-Kutta Verfahrens 4. Ordnung und Test
//        am Beispiel des harmonischen Oszillators.
//-------------------------------------------------------------------------------
#include <fstream>
#include <sstream>
#include <iostream>
#include <functional>
#include <vector>
#include <cmath>
#include "rungkutt.h"

using namespace std;

void F(vector<double> &r, vector<double> &force);
/* Kraftfeld des harmonischen Oszillators
 *      r: Ortsvektor (input)
 *      force: Kraft am Ort r (output)
 */

void energy(function<void(vector<double> &r, vector<double> &force)> F, vector<vector<double>> &r, vector<vector<double>> &v, vector<double> &U, vector<double> &T);
/* Berechnet kinetische Energie T und potentielle energie U *      F: Kraftfeld
 *      r: Ortsvektoren
 *      v: Geschwindigkeitsvektoren
 *      U: Potentielle Energie an Orten r mit v (output)
 *      T: Kinetische Energie an Orten r mit v (output)
 */

void savedat(string dest, double h, vector<vector<double>> &r, vector<vector<double>> &v, vector<double> &U, vector<double> &T);
/* Speichermethode
 */

void doWork(double h, int Ts, vector<double> r_0, vector<double> v_0, string fname);
/* Führt einen Durchgang mit gegebenen Parametern durch und speichert das
 * Ergebnis.
 *      h: Schrittweite
 *      Ts: Anzahl an zu berechnenden Perioden
 *      r_0: Startwert für den Ort
 *      v_0: Startwert für die Geschwindigkeit
 *      fnmae: Name der Zieldatei
 */

int main()
{
    double h = 1e-4;
    #pragma omp parallel

    // 3D v_0 = 0
    doWork(h, 2, {1, 1, 1}, {0, 0, 0}, "v_0=0.dat");

    // v nicht parallel zu r
    doWork(h, 2, {1, 0} , {0, 1}, "v_senkrecht_r.dat");

    // Verschiedene Schrittweiten testen
    //h = 1;
    //while ( h > 1e-5 )
    //{
        //stringstream ss;
        //ss << h;
        //cout << ss.str() << " ..." << endl;
        //doWork(h, 5, {1}, {0}, "h="+ss.str()+".dat");
        //h *= 0.1;
    //}
    return 0;
}

void F(vector<double>& r, vector<double>& force)
{
    //force.resize(r.size());
    for ( unsigned i = 0; i < r.size(); i++ )
    {
        force[i] = - r[i];
    }
}

void F1(double *r,  double *force, double d)
{
    for ( unsigned i = 0; i < d; i++ )
    {
        force[i] = - r[i];
    }
}

void energy(vector<vector<double>> &r, vector<vector<double>> &v, vector<double> &U, vector<double> &T)
{
    double s = r.size();
    double d = r[0].size();
    U.resize(s, 0);
    T.resize(s, 0);
    vector<double> force(d, 0);

    //#pragma omp parallel for
    for ( int t = 0; t < s; t++ )
    {
       for ( int i = 0; i < d; i++  )
       {
           U[t] += 0.5*r[t][i]*r[t][i];
           T[t] += 0.5*v[t][i]*v[t][i];
       }
    }
}

void savedat(string dest, double h, vector<vector<double>> &r, vector<vector<double>> &v, vector<double> &U, vector<double> &T)
{
    ofstream out;
    out.open(dest);
    out.precision(10); // gespeicherte Nachkommastellen ehöhen!!! (Zu viel ist aber auch nicht gut sonst braucht python ewig)

    for ( unsigned t = 0; t < r.size(); t++ )
    {
        out << t*h << "\t";
        for ( unsigned i = 0; i <  r[0].size(); i++ )
        {
            out << r[t][i] << "\t";
        }
        for ( unsigned i = 0; i <  r[0].size(); i++ )
        {
            out << v[t][i];
                out << "\t";
        }
        out << U[t] << "\t" << T[t] << endl;
    }
    out.close();
}

void doWork(double h, int Ts, vector<double> r_0, vector<double> v_0, string fname)
{
    vector<double> U, T;
    vector<vector<double>> r, v;

    int N = ceil( Ts*2*3.1415926/h );
    int d = r_0.size();
    double** y = new double*[N];
    for ( int i = 0; i<N; i++ )
        y[i] = new double[d*2];

    // Anfangsbedingungen setzen
    for ( int i = 0; i<d; i++ )
    {
        y[0][i] = r_0[i];
        y[0][i+d] = v_0[i];
    }

    auto f = [d](double t, double* y, double* out)
    {
        for ( int i = 0; i<d; i++ )
        {
            out[i] = y[i+d];
        }
        F1(y, out+d, d);
    };
    // Runge Kutter 4. Ordnung durchführen
    rungkutt(f, N, h, 0, y, d*2);

    r.resize(N, vector<double>(d, 0));
    v.resize(N, vector<double>(d, 0));
    for ( int i = 0; i<N; i++ )
    {
        for ( int j = 0; j<d; j++ )
        {
            r[i][j] = y[i][j];
            v[i][j] = y[i][j+d];
        }
    }

    // Energie berechnen
    energy(r, v, U, T);

    // Ergebnis speichern
    savedat(fname, h, r, v, U, T);
}
