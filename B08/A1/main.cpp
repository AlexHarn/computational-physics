#include <fstream>
#include <cmath>
#include <iostream>
#include <functional>

using namespace std;

// 1a) x_n+1 = r*x_n*( 1 - x_n )
// 1b) x_n+1 = r*x_n - x_n^3

void iterate(function<double (double, double)> f, double x0, double rmax, double rstep, int icali, int branches,  string fname)
{
    ofstream fout;
    fout.open(fname);
    fout << "#r, x" << endl;
    double x;
    for ( double r = 0; r <= rmax; r+=rstep )
    {
        x = x0;
        for ( int i = 0; i < icali; i++ )
        {
            x = f(r, x);
        }
        for ( int i = 0; i < pow(2, branches); i++ ) // maximal branches Verzweigungen werden geplottet
        {
            fout << r << "\t" << x << endl;
            x = f(r, x);
        }
    }
    fout.close();
}

int main()
{
    // Aufgabe 1a
    iterate([](double r, double x){ return r*x*( 1 - x ); }, 0.5, 4.5, 0.001, 1000, 5, "1a.dat");

    // Aufgabe 1b
    iterate([](double r, double x){ return r*x - pow(x, 3); }, 0.5, 3.5, 0.001, 1000, 5, "1b-pos.dat");
    iterate([](double r, double x){ return r*x - pow(x, 3); }, -0.5, 3.5, 0.001, 1000, 5, "1b-neg.dat");

    return 0;
}
