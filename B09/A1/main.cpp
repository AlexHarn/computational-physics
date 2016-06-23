#include <Eigen/Eigenvalues>
#include <assert.h>
#include <algorithm>
#include <fstream>

using namespace std;
using namespace Eigen;

VectorXd powerIteration(Ref<MatrixXd> A, const Ref<const VectorXd> v0, const double eps=1e-6)
{
    assert(A.rows() == A.cols() && A.rows() == v0.size());
    VectorXd r = VectorXd(A.rows());
    for ( int i = 0; i < r.size(); i++ )
    {
        VectorXd w(v0);
        double ev;
        do
        {
            ev = r(i);
            w = A*w;
            w /= w.norm();
            r(i) = w.transpose()*A*w;
        } while ( abs(r(i) - ev) > eps );
        A -= r(i)*w*w.transpose();
    }
    return r;
}

int main()
{
    Matrix4d A;
    A << 1, -2, 2, 4,
         -2, 3, -1, 0,
         2, -1, 6, 3,
         4, 0, 3, 5;

    // a)
    Vector4d a = A.eigenvalues().real();
    sort(a.data(), a.data()+a.size(), greater<double>());

    // b)
    Vector4d b = powerIteration(A, Vector4d({1, 2, 3, 4}));

    // save
    ofstream fout("A1.dat");
    fout.precision(20);
    fout << "# Eigen, powerIteration" << endl;
    for ( int i = 0; i < 4; i++ )
        fout << a(i) << "\t" << b(i) << endl;
    fout.close();
    return 0;
}
