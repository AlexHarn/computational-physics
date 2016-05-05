#include <fstream>
#include <iostream>
#include <vector>
using namespace std;


double*** rk4(int N, double x0[], double v0[], double (*F)(double r), double t0, double tmax, int d){ //N Stützstellen, x0,v0,F(x), t0 tmax intervall
	double*** y = 0;
    y = new double**[2];
    for (int i = 0; i < 4; ++i) { //x+v+Energie in Epot und Ek
        y[i] = new double*[d];
        for (int j = 0; j < d; ++j) {
            y[i][j] = new double[N];
         }
    }

    double k1[2][d];
    double k2[2][d];
    double k3[2][d];
    double k4[2][d];    

	double h = (tmax - t0)/N;

    for(int i = 0; i < d; i++){
        y[0][i][0] = x0[i];
    	y[1][i][0] = v0[i];
    }
	
	for(int i = 0; i < N-1; i++){  // alle yN durchlaufen
        for(int j = 0; j < d; j++){  // alle Dimensionen durchlaufen
            k1[0][j] = h*y[1][j][i]; // s. f(y) Blatt 3. Absatz
		    k1[1][j] = h*F(y[0][j][i]);
        }

        for(int j = 0; j < d; j++){
    		k2[0][j] = h*(y[1][j][i]+0.5*k1[1][j]);
    		k2[1][j] = h*F(y[0][j][i]+0.5*k1[0][j]);
        }

        for(int j = 0; j < d; j++){
		    k3[0][j] = h*(y[1][j][i]+0.5*k2[1][j]);
		    k3[1][j] = h*F(y[0][j][i]+0.5*k2[0][j]);
        }

        for(int j = 0; j < d; j++){		
		    k4[0][j] = h*(y[1][j][i]+k3[1][j]);
		    k4[1][j] = h*F(y[0][j][i]+k3[0][j]);		
        }

        for(int j = 0; j < d; j++){
		    y[0][j][i+1] = y[0][j][i] + 1.0/6 * (k1[0][j]+2*k2[0][j]+2*k3[0][j]+k4[0][j]);
		    y[1][j][i+1] = y[1][j][i] + 1.0/6 * (k1[1][j]+2*k2[1][j]+2*k3[1][j]+k4[1][j]);
        }

        for(int j = 0; j < d; j++){		
		    y[2][j][i+1] = 0.5*y[1][j][i]*y[1][j][i];  //Ekin = 0.5*m*v**2, m=1
		    y[3][j][i+1] = 0.5*y[0][j][i]*y[0][j][i];		//Epot = 0.5*k*x**2, k=1
        }
	}
	
	return y;
}

double F(double r){
	return -r;
}


int main(){
    int d = 3;
	double x0[d] = {1, 1, 1};
	double v0[d] = {0, 0, 0};
	double*** y  = rk4(2000, x0, v0, &F, 0, 30, d);
	
    ofstream data;
    data.open("data2a1.dat", ios::out);
    data << "#Zeit\tx\ty\tz\tvx\tvy\tvz\tekin\tepot" << endl;
    for (int i = 0; i < 2000; i++) {
    	data << i << "\t" << y[0][0][i] << "\t" << y[0][1][i] << "\t" << y[0][2][i] << "\t" << y[1][0][i] << "\t" << y[1][1][i] << "\t" << y[1][2][i] << "\t" << y[2][0][i] << "\t" << y[3][0][i] << "\t" << endl;
	}
    data.close();

	double x1[3] = {1, 0};
	double v1[3] = {0, 1};

	double*** y1  = rk4(2000, x0, v0, &F, 0, 30, d);
	
    data.open("data2a2.dat", ios::out);
    data << "#Zeit\tx\ty\tz\tvx\tvy\tvz\tekin\tepot" << endl;
    for (int i = 0; i < 2000; i++) {
    	data << i << "\t" << y1[0][0][i] << "\t" << y1[0][1][i] << "\t" << y1[0][2][i] << "\t" << y1[1][0][i] << "\t" << y1[1][1][i] << "\t" << y1[1][2][i] << "\t" << y1[2][0][i] << "\t" << y1[3][0][i] << endl;
	}
    data.close();

	return 0;
}
