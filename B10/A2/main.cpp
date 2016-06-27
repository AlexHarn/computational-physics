#include <random>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

double montePi(int N)
{
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0, 1);
    double s = 0, x, y;
    #pragma omp parallel for reduction(+:s)
    for ( int i = 0; i < N; i++ )
    {
        x = dist(mt);
        y = dist(mt);
        if ( x*x + y*y < 1 )
            s++;
    }
    return 4*s/N;
}

int main()
{
    // a)
    //cout << "Pi = " << montePi(1e8) << endl;

    // b)
    vector<vector<double>> peace(100, vector<double>(2, 0));
    #pragma omp parallel for
    for ( unsigned i = 0; i < peace.size(); i++ )
    {
        peace[i][0] = pow(10, 6.0*i/(peace.size() - 1) + 1);
        peace[i][1] = montePi(peace[i][0]);
    }
    ofstream fout("b_errors.dat");
    fout << "#N\tError" << endl;
    for ( auto &p : peace )
        fout << p[0] << "\t" << abs(M_PI - p[1]) << endl;
    fout.close();

    vector<double> pis(1000);
    #pragma omp parallel for
    for ( unsigned i = 0; i < pis.size(); i++ )
        pis[i] = montePi(1000);
    fout.open("b_hist.dat");
    for ( auto &pi : pis )
        fout << pi << endl;
    fout.close();

    return 0;
}
