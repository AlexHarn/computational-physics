#include <fstream>
#include <cmath>
#include <functional>
#include <vector>
#include <iostream>
#include <string>
#include <assert.h>

using namespace std;

double f2n(double r, double x, unsigned int n)
{
    for ( int i = 0; i < pow(2, n); i++ )
        x = r*x*( 1 - x );
    return x;
}

double df2n(double r, double x, unsigned int n)
{
    double y = 0;
    for ( int i = 0; i < pow(2, n); i++ )
    {
        y = x*( 1 - x ) + r*( 1 - 2*x )*y;
        x = r*x*( 1 - x );
    }
    return -y;
}

double regulaFalsi(function<double (double)> f, double x, double y, double eps)
{
    double z, fx, fy, fz;
    do
    {
        fx = f(x);
        fy = f(y);
        z = ( x*fy - y*fx )/( fy - fx );
        fz = f(z);

        assert(fx*fy<0);
        if ( fx*fz > 0 )
            x = z;
        else
            y = z;
    } while ( f(z) > eps );
    return z;
}

double newton(function<double (double)> f, function<double (double)> df, double x, double eps)
{
    double fx;
    do
    {
        fx = f(x);
        x = x - fx/df(x);
    } while ( fx > eps );
    return x;
}

int main()
{
    auto g = [](double r, unsigned n) { return 0.5 - f2n(r, 0.5, n); };
    ofstream fout;
    fout.precision(7);
    cout.precision(7);

    // a
    for ( int n = 0; n < 4; n++ )
    {
        fout.open("2a_n="+to_string(n)+".dat");
        fout << "#r g(r)" << endl;
        for ( double r = 0; r <= 3.5699; r+= 1e-4 )
            fout << r << "\t" << g(r, n) << endl;
        fout.close();
    }

    // Schranken
    vector<double> x = { 3., 3.4, 3.554 };
    vector<double> y = { 3.5, 3.6, 3.556 };
    
    // a nur in den Schranken plotten
    for ( int n = 1; n < 4; n++ )
    {
        fout.open("2a_n="+to_string(n)+"_schranken.dat");
        fout << "#r g(r)" << endl;
        for ( double r = x[n-1]; r <= y[n-1]; r+= 1e-5 )
            fout << r << "\t" << g(r, n) << endl;
        fout.close();
    }

    // b
    vector<double> r(3);
    for ( unsigned i = 0; i < 3; i++ )
    {
        r[i] = regulaFalsi(bind(g, placeholders::_1, i+1), x[i], y[i], 1e-6);
        cout << "[Regula Falsi] R_" << i << " = " << r[i] << endl;
    }

    // c
    cout << "[Regula Falsi] genäherte Feigenbaumkonstante: " << ( r[1] - r[0] )/( r[2] - r[1] ) << endl;

    // d
    r.resize(15);
    r[0] = 2;
    r[1] = 1 + sqrt(5);
    fout.open("d.dat1");
    fout << "#R_i für i>2, delta_i" << endl;
    double delta = 5;
    for ( unsigned i = 2; i < r.size(); i++ )
    {
        r[i] = newton(bind(g, placeholders::_1, i+1), bind(df2n, placeholders::_1, 0.5, i+1), r[i-1] + ( r[i-1] - r[i-2] )/delta, 1e-6);
        delta = ( r[i-1] - r[i-2] )/( r[i] - r[i-1] );
        fout << r[i] << "\t" << delta << endl;
    }
    fout.close();
    cout << "[Newton] genäherte Feigenbaumkonstante: " << delta << endl;
    return 0;
}
