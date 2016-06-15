#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;

// 1a) x_n+1 = r*x_n*( 1 - x_n )
// 1b) x_n+1 = r*x_n - x_n^3

int main()
{
    // Aufgabe 1a
    double x0 = 0.5; // Startwert
    int iMax = 1000; // Iterationsschritte bevor geplottet wird
    double r0 = 2.4;
    double rMax = 3.65;
    double x = x0;
    ofstream fout;
    fout.open("1a.dat");
    fout << "#r, x" << endl;
    for ( double r = r0; r <= rMax; r+=0.001 )
    {
        x = x0;
        for ( int i = 0; i < iMax; i++ )
        {
            x = r*x*( 1 - x );
        }
        for ( int i = 0; i < pow(2, 5); i++ ) // maximal 5 Verzweigungen werden geplottet
        {
            fout << r << "\t" << x << endl;
            x = r*x*( 1 - x );
        }
    }
    fout.close();

    // Aufgabe 1b
    x0 = 1; // Startwert
    iMax = 1000; // Iterationsschritte bevor geplottet wird
    r0 = 1.9;
    rMax = 2.5;
    fout.open("1b.dat");
    fout << "#r, x" << endl;
    for ( double r = r0; r <= rMax; r+=0.001 )
    {
        x = x0;
        for ( int i = 0; i < iMax; i++ )
        {
            x = r*x - pow(x, 3);
        }
        for ( int i = 0; i < pow(2, 5); i++ ) // maximal 5 Verzweigungen werden geplottet
        {
            fout << r << "\t" << x << endl;
            x = r*x - pow(x, 3);
        }
    }
    fout.close();

    return 0;
}