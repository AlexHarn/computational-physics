#include <random>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

double simulate(unsigned N, double H)
{
    double p = exp(-2*abs(H)), s, m=0;
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0, 1);
    if ( dist(mt) < 0.5)
        s = -1;
    else
        s = 1;

    for ( int i = 0; i < N; i++)
    {
        if ( s*H > 0 ) // --> Delta E > 0
        {
            if ( dist(mt) < p )
            {
                s *= -1;
            }
        }
        else
        {
            s *= -1;
        }
        m += s;
    }

    return m/N;
}

int main()
{
    ofstream fout("b.dat");
    fout << "#H\tm" << endl;
    int N = 1e4;
    double HMIN = -5, HMAX = 5;
    vector<double> savedata(N);
    #pragma omp parallel for
    for ( int i = 0; i < N; i++)
        savedata[i] = simulate(1e5, (1.0*i)/N*(HMAX-HMIN)+HMIN);
    for ( int i = 0; i < N; i++)
        fout << (1.0*i)/N*(HMAX-HMIN)+HMIN << "\t" << savedata[i] << endl;

    return 0;
}