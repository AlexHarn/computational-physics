/* ----------------------------------------------------------------------------- 
 *        Computational Physics 2016 Blatt 01 Aufgabe 1: Drehmomente
 *                            !!! C++11 Standard !!!
 * -----------------------------------------------------------------------------          
 *        3 verschiedene Integrationsroutinen. Die Verfahren sind nach der
 *        Trapezregel, Mittelpunktsregel und Simpsonregel implementiert.
 *        Dabei sind zwei Integrale I1 und I2 um die Implementierungen
 *        zu testen.
 * ---------------------------------------------------------------------------*/
#include <iostream>
#include <cmath>        
using namespace std;

// Implementierung der Trapezregel
double trapezregel(double (*f)(double x), double a, double b, int N) {
    double sum = 0;     // Summenwert sum und
    double h;           // Abstand h initialisieren
    h = (b-a) / N;   // h berechnen
    for (int i = 1; i < N; i++)     // Mittelteil
        sum += f(a+i*h);
    sum += 0.5 * (f(a)+f(b));       // Anfang und Ende
    return h * sum;
}

double mittelpunktsregel(double (*f)(double x), double a, double b, int N) {
    double sum = 0;
    double h;
    h = (b-a) / N;
    sum += f(a+h/2);    // Anfang
    for (int i = 1; i < N-1; i++)   // Mittelteil
        sum += f(a+h/2+i*h);
    sum += f(b-h/2);    // Ende
    return h * sum;
}

double simpsonregel(double (*f)(double x), double a, double b, int N) {
    if(N % 2){      // Die Schrittzahl muss gerade sein.
        cout << "N ist nicht gerade! N wird jetzt um 1 erhoeht." << endl;
        N+=1;
    }
    double sum = 0;
    double h;
    h = (b-a) / N;
    sum += 1.0/3.0 * f(a);        // Anfang
    for (int i = 1; i < N; i+=2)      // Mittelteil
        sum += 4.0/3.0 * f(a+i*h);      // dabei werden die Elemente abwechselnd mit 4/3 und 2/3 gewichtet
    for (int i = 2; i < N; i+=2)
        sum += 2.0/3.0 * f(a+i*h);
    sum += 1.0/3.0 * f(b);      // Ende
    return h * sum;
}

double I1(double x) {           // Erste Testfunktion (siehe Aufgabe 3)
    return exp(-x)/x;
}

double I2(double x) {           // Zweite Testfunktion (siehe Aufgabe 3)
    return x * sin(1/x);
}

int main() {
    // Test der Implementierungen
    cout << "Erstes Integral:" << endl;
    cout << "Trapezregel: I_1 = " << trapezregel(&I1, 1, 100, 10000) << endl;
    cout << "Mittelpunktsregel: I_1 = " << mittelpunktsregel(&I1, 1, 100, 10000) << endl;
    cout << "Simpsonregel: I_1 = " << simpsonregel(&I1, 1, 100, 10000) << endl << endl;

    cout << "Zweites Integral:" << endl;
    cout << "Trapezregel: I_1 = " << trapezregel(&I2, 0, 1, 10000) << endl;
    cout << "Mittelpunktsregel: I_1 = " << mittelpunktsregel(&I2, 0, 1, 10000) << endl;
    cout << "Simpsonregel: I_1 = " << simpsonregel(&I2, 0, 1, 10001) << endl << endl;
    
    return 0;
}
