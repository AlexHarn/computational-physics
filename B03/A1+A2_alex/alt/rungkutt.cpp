#include <iostream>
#include <functional>
#include <vector>
#include <math.h> // für std::ceil

using namespace std;

void rungkutt(function<void(vector<double> &r, vector<double> &force)> F, double t_0, double t_N, double h, vector<vector<double>> &r, vector<vector<double>> &v)
{
    // Größen der Ausgabevektoren vorher festlegen, sollte dringend gemacht werden, da .push_back
    // sehr sehr viel langsamer ist
    int N = ceil( ( t_N - t_0 )/h + 1); // Zaunpfahlproblem

    int d = r[0].size(); // Anzahl der Dimensionen
    r.resize(N, vector<double>(d, 0));
    v.resize(N, vector<double>(d, 0));

    vector<double> k1(d*2, 0), k2(d*2, 0), k3(d*2, 0), k4(d*2, 0), temp(d*2, 0); // ks für Geschwindigkeits- und Ortsbestimmung (darum d*2)

    for ( double t = 0; t < N - 1; t++ )
    {
        F(r[t], k1);
        for ( int i = 0; i < d; i++ )
        {
            k1[i+d] = v[t][i];
        }

        for ( int i = 0; i < d; i++ )
        {
            k2[i+d] = v[t][i]*k1[i]*h*0.5;
            temp[i] = r[t][i] + 0.5*k1[i]*h;
        }
        F(temp, k2);

        for ( int i = 0; i < d; i++ )
        {
            k3[i+d] = v[t][i]*k2[i]*h*0.5;
            temp[i] = r[t][i] + 0.5*k2[i]*h;
        }
        F(temp, k3);

        for ( int i = 0; i < d; i++ )
        {
            k4[i+d] = v[t][i]*k3[i]*h;
            temp[i] = r[t][i] + k3[i]*h;
        }
        F(temp, k4);

        for ( int i = 0; i < d; i++ )
        {
            v[t+1][i] = v[t][i] + h/6*( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] );
            r[t+1][i] = r[t][i] + h*( k1[i+d] + 2*k2[i+d] + 2*k3[i+d] + k4[i+d]);
        }
    }
}
