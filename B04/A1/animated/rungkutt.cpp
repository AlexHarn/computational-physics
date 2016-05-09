#include <functional>
#include <cmath>

using namespace std;

void rungkutt(function<void(double t, double* y, double* out)> f, const double N, const double h, const double t_0, double** y, const int d)
{
	double* k1 = new double[d];
	double* k2 = new double[d];
	double* k3 = new double[d];
	double* k4 = new double[d];
	double* temp = new double[d];

	for ( int t = 0; t < N-1; t++ )
	{
		f(t_0 + t*h, y[t], k1);

		for ( int i = 0; i < d; i++ )
		{
			temp[i] = y[t][i] + 0.5*h*k1[i];
		}
		f(t_0 + h/2 + t*h, temp, k2);

		for ( int i = 0; i < d; i++ )
		{
			temp[i] = y[t][i] + 0.5*h*k2[i];
		}
		f(t_0 + h/2 + t*h, temp, k3);

		for ( int i = 0; i < d; i++ )
		{
			temp[i] = y[t][i] + h*k3[i];
		}
		f(t_0 + h + t*h, temp, k4);

		for ( int i = 0; i < d; i++ )
		{
			y[t+1][i] = y[t][i] + h/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
		}
	}

	delete []k1;
	delete []k2;
	delete []k3;
	delete []k4;
	delete []temp;
}
