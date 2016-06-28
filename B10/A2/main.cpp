#include <random>
#include <functional>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

double montePi(const int N)
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

double monteEllipse(const double a, const double b, const int N, function<double (const double x, const double y)> f = [](double , double ) { return 1; })
{
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(-1, 1);
    double s = 0, a2 = a*a, b2 = b*b, x, y;
    #pragma omp parallel for reduction(+:s)
    for ( int i = 0; i < N; i++ )
    {
        x = a*dist(mt);
        y = b*dist(mt);
        if ( x*x/a2 + y*y/b2 < 1 )
            s += f(x, y);
    }
    return s/N*4*a*b;
}

int main()
{
    ofstream fout;
    cout.precision(6);
    fout.precision(6);
    // a)
    cout << "a) Pi = " << montePi(1e8) << endl;

    // b)
    vector<vector<double>> peace(100, vector<double>(2, 0));
    #pragma omp parallel for
    for ( unsigned i = 0; i < peace.size(); i++ )
    {
        peace[i][0] = pow(10, 6.0*i/(peace.size() - 1) + 1);
        peace[i][1] = montePi(peace[i][0]);
    }
    fout.open("b_errors.dat");
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

    // c)
    vector<vector<double>> area(10, vector<double>(10, 0));
    #pragma omp parallel for
    for ( int a = 1; a <= 10; a++ )
        for ( int b = 1; b <= 10; b++ )
            area[a-1][b-1] = monteEllipse(a, b, 1e6);
    fout.open("c.dat");
    for ( int a = 1; a <= 10; a++ )
        for ( int b = 1; b <= 10; b++ )
            fout << a << "\t" << b << "\t" << area[a-1][b-1] << endl;
    fout.close();

    // d)
    cout << "d) I = " << monteEllipse(sqrt(2), 1, 1e8, [](double x, double ) { return exp(-x*x); }) << endl;

    return 0;
}
