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

void work(const Ref<const MatrixXd> H, string save)
{
    VectorXd evs = H.eigenvalues().real();
    sort(evs.data(), evs.data()+evs.size());
    evs = evs.head(10);

    ofstream fout(save);
    fout.precision(20);
    fout << "# Eigenvalues" << endl;
    for ( int i = 0; i < evs.size(); i++ )
        fout << evs(i) << endl;
    fout.close();
}

int main()
{
    // 2b)
    //  lam = 0
    work(Hb(10, 0.1, 0), "2b0.dat");

    cout << "superstrange" << endl;

    //  lam = 0.2
    work(Hb(10, 0.1, 0.2), "2b0.2.dat");

    return 0;
}
