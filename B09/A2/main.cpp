#include <Eigen/Eigenvalues>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <iostream>

using namespace std;
using namespace Eigen;

MatrixXd Hb(const double L, const double delta, const double lam)
{
    const int N = 2*L/delta + 1;
    const double d2 = delta*delta;
    assert(round(L/delta) == L/delta);

    MatrixXd H(N, N);
    H(0, 0) = 2/d2 + d2*pow(L/delta, 2) + lam*d2*d2*pow(L/delta, 4);
    H(0, 1) = -1/d2;
    for ( int i = 1; i < N-1; i++ )
    {
        double n = i - L/delta;
        H(i, i-1) = -1/d2;
        H(i, i+1) = H(i, i-1);
        H(i, i) = 2/d2 + d2*n*n + lam*d2*d2*pow(n, 4);
    }
    H(N-1, N-2) = -1/d2;
    H(N-1, N-1) = H(0, 0); //2/d2 + d2*pow(L/delta, 2) + lam*d2*d2*pow(L/delta, 4);
    return H;
}

MatrixXd Hc(const double lam, const int N=50)
{
    MatrixXd H(N, N);
    for ( int i = 4; i < N; i++ )
    {
        H(i, i-4) = lam*sqrt(i*( i - 1 )*( i - 2 )*( i - 3 ));
        H(i-4, i) = H(i, i-4);
    }
    for ( int i = 2; i < N; i++ )
    {
        H(i, i-2) = lam*sqrt(i*( i - 1 ))*( 4*i - 2 );
        H(i-2, i) = H(i, i-2);
    }
    for ( int i = 0; i < N; i++ )
    {
        H(i, i) = 8*i + 4 + lam*( 6*i*i + 6*i + 3 );
    }
    return H/4;
}

void work(MatrixXd H, string save, int saveNum=10)
{
    VectorXd evs = H.eigenvalues().real();
    sort(evs.data(), evs.data()+evs.size());

    ofstream fout(save);
    fout.precision(20);
    fout << "# Eigenvalues" << endl;
    for ( int i = 0; i < saveNum; i++ )
        fout << evs(i) << endl;
    fout.close();
}

int main()
{
    // 2b)
    work(Hb(10, 0.1, 0), "2b0.dat");
    work(Hb(10, 0.1, 0.2), "2b0.2.dat");

    // 2c)
    work(Hc(0), "2c0.dat");
    work(Hc(0.2), "2c0.2.dat");

    return 0;
}
