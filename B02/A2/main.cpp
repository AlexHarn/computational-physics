#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

double a = 1;

double mittelWuerfel(double (*f)(double mp, double x, double y, double z), double mp, double a, double b, int N) {
    double sum = 0;
    double h = (b-a) / N;
    for (int i = 0; i < N; i++) { //Unterteilen in kleine Würfel
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                if ((mp == a+h/2+i*h) && (a+h/2+j*h == 0) && (a+h/2+k*h == 0)) {                // Division durch 0, eine Stützstelle wird rausgenommen
                    cout << "Singularitaet! Messpunkt " << mp << " ist ungenau!" << endl;
                } else {
                    sum += f(mp, a+h/2+i*h, a+h/2+j*h, a+h/2+k*h);      // erst z-, dann y-, dann x-Integration
                }
            }
        }
    }
    return pow(h, 3) * sum;
}


double f(double mp, double x, double y, double z) {             // mp = Messpunkt, probl wenn mp=x und y=z=0
    return pow(pow(mp - x, 2) + pow(y, 2) + pow(z, 2), -0.5);
}

double f2(double mp, double x, double y, double z) {             // mp=x/a
    return x * pow(pow(mp - x, 2) + pow(y, 2) + pow(z, 2), -0.5);
}

int main() {
    ofstream data;
    data.open("data2a+b.dat", ios::out);
    data << "x-Werte\tFunktionswerte" << endl;
    cout << "Aufgabenteil a)" << endl;    
    for (int n = 11; n <= 80; n++) {
        data << 0.1*n << "\t" << mittelWuerfel(&f, 0.1*n, -a, a, 100) << endl; // phi ist dimlos über *4pieps_0/rho_0, x_0=x/a, int von -a,a
        cout << (n-10.0)/70.0*100 << "\%" << endl; //Erwartung: 1/x-Abfall (-Grad phi = E), passt, da Werte gg 0 laufen. 
    }
    cout << "Aufgabenteil b)" << endl;
    for (int n = 0; n <= 10; n++) {
        data << 0.1*n << "\t" << mittelWuerfel(&f, 0.1*n, -a, a, 100) << endl; // phi ist dimlos über *4pieps_0/rho_0, x_0=x/a, int von -a,a
        cout << n/10.0*100 << "\%" << endl; //Erwartung: 1/x-Abfall (-Grad phi = E), passt, da Werte gg 0 laufen. Glaubt mir 
    }
    data.close();
    cout << "Aufgabenteil c)" << endl;
    data.open("data2c.dat", ios::out);
    for (int n = 0; n <= 80; n++) {
        data << 0.1*n << "\t" << mittelWuerfel(&f2, 0.1*n, -a, a, 100) << endl; // x/a fkt , grauer Bereich ist Würfel
        cout << n/80.0*100 << "\%" << endl;
    }
    data.close();
    return 0;
}
