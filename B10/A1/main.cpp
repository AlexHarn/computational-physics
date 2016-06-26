#include <random>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

vector<double> ransin(int N)
{
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0, M_PI);

    int count = 0;
    vector<double> r(N);
    while ( count < N )
    {
        r[count] = dist(mt);
        if ( dist(mt)/( 2*M_PI ) < sin(r[count])/2 )
            count++;
    }
    return r;
}

int main()
{
    // a)
    vector<double> theta = ransin(1e6);
    ofstream fout("a.dat");
    for ( auto &t : theta )
        fout << t << endl;
    fout.close();

    // b)
    vector<double> omega(theta.size());
    double s = 0;
    fout.open("b.dat");
    for ( unsigned i = 0; i < omega.size(); i++ )
    {
        omega[i] = 0.5*( 3 * pow(cos(theta[i]), 2) - 1 );
        s += omega[i];
        fout << omega[i] << endl;
    }
    cout << "Omega gemittelt: "<< s/omega.size() << endl;
    fout.close();
    return 0;
}
