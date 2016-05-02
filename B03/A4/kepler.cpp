#include <vector>
#include <functional>
#include <fstream>
#include <sstream>
#include <cmath>
#include "kepler.h"
#include "rungkutt.h"

using namespace std;

void energy(function<double(vector<double> &r)> V, vector<vector<double>> &r, vector<vector<double>> &v, vector<double> &U, vector<double> &T)
{
    double s = r.size();
    double d = r[0].size();
    U.resize(s, 0);
    T.resize(s, 0);

    //#pragma omp parallel for
    for ( int t = 0; t < s; t++ )
    {
       U[t] = V(r[t]);
       for ( int i = 0; i < d; i++  )
       {
           T[t] += 0.5*v[t][i]*v[t][i];
       }
    }
}

void savedat(string dest, function<double(vector<double> &r)> V,  double h, vector<vector<double>> &r, vector<vector<double>> &v)
{
    vector<double> U, T, L, A;
    vector<vector<double>> LR;
    ofstream out;
    out.open(dest);
    out.precision(10); // gespeicherte Nachkommastellen ehöhen!!! (Zu viel ist aber auch nicht gut sonst braucht python ewig)

    // Energie berechnen
    energy(V, r, v, U, T);

    // Drehimpuls berechnen
    angm(r, v, L);

    // A(t) berechnen
    area(r, A);

    // Lenz-Runge Vektor berechnen
    lenzrunge(r, v, LR);

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
        for ( unsigned i = 0; i <  r[0].size(); i++ )
        {
            out << LR[t][i];
                out << "\t";
        }
        out << L[t] << "\t" << A[t] << "\t" << U[t] << "\t" << T[t] << endl;
    }
    out.close();
}

void doWork(double alpha, double h, int t, vector<double> r_0, vector<double> v_0, string fname)
{
    vector<vector<double>> r, v;
    // Anfangsbedingungen setzen
    r.push_back(r_0);
    v.push_back(v_0);

    // alpha binden
    auto Fa = bind(F, placeholders::_1, placeholders::_2, alpha);
    auto Va = bind(V, placeholders::_1, alpha);

    // Runge Kutter 4. Ordnung durchführen
    rungkutt(Fa, 0, t, h, r, v);

    // Ergebnis weiter bearbeiten und speichern
    savedat(fname, Va,  h, r, v);
}

vector<double> inverse(vector<double> v)
{
    vector<double> r;
    for ( auto vi : v )
        r.push_back(-vi);
    return r;
}

void forwNback(double alpha, double h, int N, vector<double> r_0, vector<double> v_0, string fname)
{
    vector<vector<double>> r, dr, v;

    // alpha binden
    auto Fa = bind(F, placeholders::_1, placeholders::_2, alpha);

    //#pragma omp for
    for ( int n = 1000; n<=N; n+=1000 )
    {
        // Anfangsbedingungen setzen
        r.push_back(r_0);
        v.push_back(v_0);

        // Runge Kutter 4. Ordnung durchführen, umdrehen und nochmal
        rungkutt(Fa, 0, n*h, h, r, v);
        r.erase(r.begin()+1, r.end());
        v.erase(v.begin()+1, v.end());
        v.back() = inverse(v.back());
        rungkutt(Fa, 0, n*h, h, r, v);
        dr.push_back(r.back());

        // Fehler berechnen
        for ( unsigned i = 0; i<r_0.size(); i++)
        {
            dr.back()[i] = r_0[i] - dr.back()[i];
        }
    }

    // Ergebnis speichern
    ofstream out;
    out.precision(10);
    out.open(fname);
    double s;
    for ( unsigned i = 0; i<dr.size(); i++ )
    {
        s = 0;
        for ( unsigned j = 0; j<dr[0].size(); j++ )
        {
            s += dr[i][j];
        }
        out << (i+1)*1000 << "\t" << abs(s) << endl;
    }
}

double vabs(vector<double> &v)
{
    double R = 0;
    for ( auto vi : v )
    {
        R += vi*vi;
    }
    return sqrt(R);
}

double V(vector<double> &r, double alpha)
{
    return -m*G/pow(vabs(r), alpha);
}

void F(vector<double> &r, vector<double> &force, double alpha)
{
    //force.resize(r.size()); Darum muss sich rungkutt ab jetzt kümmern, könnte sonst
    //Probleme mit den ks geben

    double R = vabs(r);
    for ( unsigned i = 0; i<r.size(); i++ )
    {
        force[i] = -G*r[i]/pow(R, alpha+2);
    }
}

void cross(vector<double> a, vector<double> b, vector<double> &result)
{
    // Hier ist auch wichtig, dass a und b NICHT per reference übergeben
    // werden! sonst gibt bspw. cross(a, b, b) Probleme!
    // Weiß nicht ob der compiler gut genug is um das nicht zig mal auszuführen
    // in den schleifen
    //if ( a.size() != 3 )
        //throw invalid_argument("Nur für d=3 implementiert!");
    //if ( a.size() != b.size() )
        //throw invalid_argument("Dimensionen stimmen nicht überein!");

    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}

void angm(vector<vector<double>> &r, vector<vector<double>> &v, vector<double> &L)
{
    L.resize(r.size());
    vector<double> l(r[0].size(), 0);
    for ( unsigned t = 0; t<r.size(); t++ )
    {
        cross(r[t], v[t], l);
        L[t] = m*vabs(l);
    }
}

void area(vector<vector<double>> &r, vector<double> &A)
{
    vector<double> av(r[0].size(), 0);
    double s, a, b, c;
    c = vabs(r[0]);

    A.resize(r.size());
    A[0] = 0;
   for ( unsigned t = 1; t<r.size(); t++ )
   {
       for ( unsigned i = 0; i<r[0].size(); i++)
           av[i] = r[t][i] - r[t-1][i];
       a = vabs(av);
       b = vabs(r[t]);
       s = ( a + b + c )/2;
       A[t] = A[t-1] + sqrt(( s - a )*( s - b )*( s - c ));
       c = b;
   }
}

void lenzrunge(vector<vector<double>> &r, vector<vector<double>> &v, vector<vector<double>> &LR)
{
    vector<double> vrv(r[0].size(), 0);
    double R, pre;
    pre = m/G;
    LR.resize(r.size());
    for ( unsigned t = 0; t<r.size(); t++ )
    {
        cross(r[t], v[t], vrv);
        cross(v[t], vrv, vrv);
        R = vabs(r[t]);
        LR[t].resize(r[0].size());
        for ( unsigned i = 0; i<r[0].size(); i++ )
            LR[t][i] = pre*vrv[i] - r[t][i]/R;
    }
}
