#include <iostream>
#include <cmath>        
using namespace std;

double delta = 0.0001;

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

double iterativSimpson(double (*f)(double x), double eps, double a, double b, int N) {
    double simpsonNeu = 0; // Eine Art iterative Simpson-Variable, welche (vgl do-while-schleife unten) die Simpson-Regel für kleinere kritische Integralgrenzen (untere Grenze) erneut aufruft. Dabei wird immer die                    
                           // Mitte des Intervalls als neue Grenze gesetzt. simpsonNeu ist dann der neu berechnete Integralwert
    double simpsonAlt = 0; //alter Integralwert
    double c = (a+b)/2; // Setze untere Grenze auf Mitte des Intervalls
    do {
        simpsonAlt = simpsonNeu;
        simpsonNeu = simpsonAlt + simpsonregel(f, c, b, N); // Neuer Integralwert wird mit alten Integralwert * Schrittweite summiert.
        b = c; // Setze untere Grenze auf Mitte des Intervalls (gleiches Spiel wie oben für neuen Wert)
        c = (a+b)/2;
    } while(abs((simpsonNeu-simpsonAlt)/simpsonNeu) > eps); //Mache solange, bis Fehler bei 10^(-5) (bzw 10^(-6), vgl Kommentar unten). Dabei nähert dieser Wert sich nur dem rel. Fehler von 10^(-5), da der nächste Schritt über 
                                                            //10^(-5) hinaus ginge.
    return simpsonNeu;
}


double I1(double x) {           // Erste Teilaufgabe (außerhalb Singularität)
    return exp(x)/x;
}

double I11(double x) {           // Erste Teilaufgabe (innerhalb  Singularität)
    return delta + 0.5*pow(delta, 2)*x+1.0/6*pow(delta, 3)*pow(x, 2); // Taylor
}

double I2(double x) {           // umgestellt siehe Zusammenfassung.pdf
    return exp(-x)/sqrt(x) + pow(x, -1.5) * exp(-1/x);
}

int main() {
    cout.precision(10);
    // Test der Implementierungen
    cout << "Erstes Integral:" << endl;
    cout << "Simpsonregel: I_1 = " << simpsonregel(&I1, -1, -delta, 10000) + simpsonregel(&I11, -1, 1, 10000) + simpsonregel(&I1, delta, 1, 10000) << endl;
    cout << "Zweites Integral:" << endl;
    cout << "Iterativ: I_2 = " << iterativSimpson(&I2, pow(10, -6), 0, 1, 100000) << endl; // Programm hört bei 2*10^(-5) auf, nächster Schritt wäre <10^(-5), daher nach unten etwas großzügiger sein
    return 0;
}
