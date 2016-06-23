#include <Eigen/Eigenvalues>
#include <fstream>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;

MatrixXd xH(const double L, const double delta, const double lam)
{
    //assert(round(2*L/delta) == 2*L/delta);
    const int N = round(2*L/delta + 1);
    const double d2 = delta*delta;

    MatrixXd H(N, N);
    H.setZero();
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
    H(N-1, N-1) = H(0, 0);
    return H;
}

MatrixXd eH(const double lam, const int N=50)
{
    MatrixXd H(N, N);
    H.setZero();
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

void work(Ref<const MatrixXd> H, string save, int saveNum=10)
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
    #pragma omp parallel
    {
        // 2b)
        work(xH(10, 0.1, 0), "2b_l=0.dat");
        work(xH(10, 0.1, 0.2), "2b_l=0.2.dat");

        // 2c)
        work(eH(0), "2c_l=0.dat");
        work(eH(0.2), "2c_l=0.2.dat");
    }

    // 2d)
    #define N 41
    vector<VectorXd> evs(N);
    #pragma omp parallel for
    for ( int i = 0; i<N; i++ )
    {
        evs[i] = xH(10, 20./( i + 9 ), 0.2).eigenvalues().real();
        sort(evs[i].data(), evs[i].data()+evs[i].size());
    }
    ofstream fout("2d_x.dat");
    fout.precision(20);
    fout << "#N\tEigenwerte" << endl;
    for ( unsigned i = 0; i<evs.size(); i++ )
    {
        fout << i + 10 << "\t";
        for ( int j = 0; j < 10; j++ )
            fout << evs[i](j) << "\t";
        fout << endl;
    }
    fout.close();

    evs.resize(N);
    #pragma omp parallel for
    for ( int n = 10; n < 10+N; n++ )
    {
        evs[n-10] = eH(0.2, n).eigenvalues().real();
        sort(evs[n-10].data(), evs[n-10].data()+evs[n-10].size());
    }
    fout.open("2d_e.dat");
    fout << "#N\tEigenwerte" << endl;
    for ( int n = 10; n < 10+N; n++ )
    {
        fout << n << "\t";
        for ( int i = 0; i < 10; i++ )
            fout << evs[n-10](i) << "\t";
        fout << endl;
    }
    return 0;
}
