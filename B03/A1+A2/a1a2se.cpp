#include <fstream>
#include <iostream>
#include <vector>
using namespace std;


double** rk4(int N, double x0, double v0, double (*F)(double r), double t0, double tmax){
	double** y = 0;
	y = new double*[2];
	for(int i = 0; i < 2; i++){
		y[i] = new double[N];
	}
	double k1[2];
	double k2[2];
	double k3[2];
	double k4[2];
	double h = (tmax - t0)/N;
	y[0][0] = x0;
	y[1][0] = v0;
	
	for(int i = 0; i < N-1; i++){
		k1[0] = h*y[1][i];
		k1[1] = h*F(y[0][i]);
		
		k2[0] = h*(y[1][i]+0.5*k1[1]);
		k2[1] = h*F(y[0][i]+0.5*k1[0]);
		
		k3[0] = h*(y[1][i]+0.5*k2[1]);
		k3[1] = h*F(y[0][i]+0.5*k2[0]);
		
		k4[0] = h*(y[1][i]+k3[1]);
		k4[1] = h*F(y[0][i]+k3[0]);
		
		y[0][i+1] = y[0][i] + 1.0/6 * (k1[0]+2*k2[0]+2*k3[0]+k4[0]);
		y[1][i+1] = y[1][i] + 1.0/6 * (k1[1]+2*k2[1]+2*k3[1]+k4[1]);
	}
	
	return y;
}

double F(double r){
	return -r;
}


int main(){
	
	double** y  = rk4(1000, 1, 0, &F, 0, 10);
	
    ofstream data;
    data.open("data2.dat", ios::out);
    data << "#Zeit\tOrt\tGeschwindigkeit" << endl;
    for (int i = 0; i < 1000; i++) {
		data << i << "\t" << y[0][i] << "\t" << y[1][i] << endl;
	}
    data.close();

	return 0;
}