#include <iostream>
#include <cmath>
#include <functional>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
using namespace std;
// rungkutt.h
class rungkutt {
rungkutt(std::function<void(double t, double* y, double* out)> f, const double N, const double h, const double t_0, double** y, const int d);
    /* Führt das Runge-Kutta-Verfahren 4. Ordnung durch. Neue Implementierung ohne std:vector und in
     * allgemeiner Formulierung wie in der Vorlesung.
     *      DGL: y'(t) = f(t, y(t))
     *      N: Schrittzahl
     *      h: Schrittweite
     *      t_0: Startzeit
     *      y: y Vektoren zu allen Zeiten. Bei Aufruf müssen Starbedingungen gesetzt sein
     *      d: Länge der y Vektoren (!!! Nicht Anzahl der Raumdimensionen !!!)
     */
};
// Header Dpendulum.h
class Dpendulum
{
        void static f(double /*t*/, double* y, double* dy, double g);
        /* f für Runge-Kutta mit y' = f(y)
         *      y: y vektor (input)
         *      dy: y' vektor (output)
         *      g: Gravitationskonstante (wird später gebunden)
         */
        //void calcE();
        /* Berechnet U und T aus y
         */
        void reset();
        /* Setzt das Pendel zurück (löscht die aktuellen Daten)
         */
        const double g = 9.81;
        int N = 0;
        double tN = 0, lh = 0;
        double* y0;
        double** y, ** T, ** U;
        bool ready = false, active = false; //energy = false

    public:
        ~Dpendulum();
        /* Destruktor
         */
        void setInitial(double theta1, double theta2, double ddtTheta1, double ddtTheta2);
        /* Setzt die AB
         */
        void swing(double h, double t);
        /* Berechnet Schwingung mit Schrittweite h bis Zeitpunkt t
         */
        void swing(double h, int n);
        /* Berechnet Schwingung mit Schrittweite h und Schrittzahl n
         */
        void save(std::string fname);
        /* Speichert die Daten in Datei fname
         */
        void doEverything(double theta1, double theta2, double ddtTheta1, double ddtTheta2, double h, double t, std::string fname);
        /* Macht alles in einem Aufruf
         */
		void teilC(double Energie, double h, double t, std::string fname);
		/* Geändertes save, mehrere Aufrufe mit anderen Startwerten aber gleicher Energie
		*/
		void saveC(std::string fname);
		/* siehe oben
		*/
};

//rungkutt.cpp
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

// dependulum.cpp
void Dpendulum::f(double /*t*/, double* y, double* dy, double g)
{
    // einige mehrfach benötigten Terme vorher berechnen um Rechenzeit zu sparen
    const double cost2t1 = cos(y[1] - y[0]);
    const double sint2t1 = sin(y[1] - y[0]);
    const double sint1 = sin(y[0]);
    const double sint2 = sin(y[1]);
    const double pre = 1/( 1 - 0.5*cost2t1*cost2t1 );

    // Alles einsetzen
    dy[0] = y[2]; // d/dt theta_1 = theta_3
    dy[1] = y[3]; // d/dt theta_2 = theta_4
    dy[2] = pre*( 0.5*g*sint2*cost2t1 + 0.5*y[2]*y[2]*sint2t1*cost2t1 - g*sint1 + 0.5*y[3]*y[3]*sint2t1 ); // vergleiche Zettel
    dy[3] = pre*(g*sint1*cost2t1 - 0.5*y[3]*y[3]*sint2t1*cost2t1 - g*sint2 - 0.5*y[2]*y[2]*sint2t1 ); // vergleiche Zettel

    // In Kleinwinkelnäherung
    //dy[2] = g*(y[1] - 2*y[0]);
    //dy[3] = 2*g*(y[0] - y[1]);
}

/*void Dpendulum::calcE()			// Energieberechnung bei Aufgabe 2 nicht nötig
{
    if ( energy )
        throw "Energie bereits berechnet!";

    const double m = 1;

    T = new double*[this->N];
    U = new double*[this->N];

    #pragma omp parallel for
    for ( int i = 0; i<N; i++ )
    {
        T[i] = new double[2];
        U[i] = new double[2];

        T[i][0] = 0.5*m*y[i][2]*y[i][2];
        U[i][0] = 2*m*g*(1-cos(y[i][0]));

        T[i][1] = 0.5*m*( y[i][2]*y[i][2] + y[i][3]*y[i][3] + 2*y[i][2]*y[i][3]*cos(y[i][0] - y[i][1]) );
        U[i][1] = m*g*(1-cos(y[i][1]));
    }
    energy = true;
}*/

void Dpendulum::reset()
{
    if ( !active )
    {
        if ( ready )
        {
            delete[] y0;
            ready = false;
        }
        return;
    }
    for ( int i = 0; i<N; i++ )
    {
        delete[] y[i];
        /*if ( energy )
        {
            delete[] T[i];
            delete[] U[i];
        }*/
    }
    delete[] y;
    /*if ( energy )
    {
        delete[] T;
        delete[] U;
    }*/
    N = 0;
    tN = 0;
    active = false;
    ready = false;
    //energy = false;
}

// PUBLIC:
Dpendulum::~Dpendulum()
{
    this->reset();
}

void Dpendulum::setInitial(double theta1, double theta2, double ddtTheta1, double ddtTheta2)
{
    if ( ready )
        throw "Bereits readyiallisiert";
    y0 = new double[4] { theta1, theta2, ddtTheta1, ddtTheta2 };
    ready = true;
}

void Dpendulum::swing(double h, double t)
{
    this->swing(h, (int) ceil( t/h )+1);
}

void Dpendulum::swing(double h, int n)
{
    using namespace std::placeholders;

    if ( !ready )
        throw "Noch nicht readyiallisiert!";

    if ( tN > 0 )
    {
        double** ny = new double*[N+n];
        std::copy(y, y + N, ny);
        delete[] y;
        y = ny;
        for ( int i = N; i<N+n; i++ )
            y[i] = new double[4];
        rungkutt(std::bind(Dpendulum::f, _1, _2, _3, g), n, h, tN, y+N-1, 4);
    }
    else
    {
        y = new double*[n];
        for ( int i = 1; i<n; i++ )
            y[i] = new double[4];
        y[0] = y0;
        rungkutt(std::bind(Dpendulum::f, _1, _2, _3, g), n, h, 0, y, 4);
    }
    N += n;
    tN += h*n;
    lh = h;
    active = true;
}

void Dpendulum::save(string fname)
{
    if ( !active )
        throw "Noch nichts berechnet!";

    //this->calcE();
    ofstream fout;
    fout.open(fname);
    for ( int i = 0; i<N; i++ )
    {
        fout << i*lh << "\t";
        for ( int j = 0; j<4; j++ ) {
            fout << y[i][j] << "\t";
		}
        fout << endl;
        //fout << T[i][0] << "\t" << T[i][1] << "\t";	//	Energieberechnung bei Aufgabe 2 nicht nötig
        //fout << U[i][0] << "\t" << U[i][1] << "\t";
        //fout << T[i][0]+T[i][1]+U[i][0]+U[i][1] << endl;
	}
    fout.close();
    return;
}


void Dpendulum::saveC(string fname)
{
    if ( !active )
        throw "Noch nichts berechnet!";

    ofstream fout;
    fout.open(fname/*, ios_base::app*/);
    for ( int i = 1; i<N; i++ )
    {
		if ( (y[i][1] * y[i-1][1] < 0) && (y[i][2] * cos(y[i][0]) + y[i][3] * cos(y[i][1]) > 0) ) {			// Aufgabe 2c
			fout << (y[i][2] + y[i-1][2])/2 << "\t" << (y[i][0] + y[i-1][0])/2 << endl;
		}
	}
    fout.close();
    return;
}

void Dpendulum::doEverything(double theta1, double theta2, double ddtTheta1, double ddtTheta2, double h, double t, std::string fname)
{
    this->reset();
    this->setInitial(theta1, theta2, ddtTheta1, ddtTheta2);
    this->swing(h, t);
    this->save(fname);
}

void Dpendulum::teilC(double E, double h, double t, std::string fname)
{
/*	this->reset();											// Das war die usprüngliche Implementierung mit nur 4 Startwerten
    this->setInitial(acos((2*g-E)/(2*g)), 0, 0, 0);
    this->swing(h, t);
    this->saveC(fname+"_1.dat");

	this->reset();
    this->setInitial(0, acos((g-E)/g), 0, 0);
    this->swing(h, t);
    this->saveC(fname+"_2.dat");

	this->reset();
    this->setInitial(0, 0, sqrt(E), 0);
    this->swing(h, t);
    this->saveC(fname+"_3.dat");

	this->reset();
    this->setInitial(0, 0, 0, sqrt(2*E));
    this->swing(h, t);
    this->saveC(fname+"_4.dat");*/

    ofstream fout;
    // fout.open(fname);
    srand(time(NULL));
    for(int i = 0; i < 80; i++){
        this->reset();
        double Etest = 0;
        double start[4] = {0, 0, 0, 0};
        do {
            start[0] = (rand() % 10000)/10000.0;
            start[1] = (rand() % 10000)/10000.0;
            start[2] = (rand() % 10000)/10000.0;
            start[3] = 0;
            Etest = 0.5*start[2]*start[2] + 2*g*(1-cos(start[0])) + 0.5*start[2]*start[2] + g*(1-cos(start[1]));
        } while(Etest > E);
        start[3] = sqrt(start[2]*start[2]*cos(start[0]-start[1])*cos(start[0]-start[1]) + 2*(E-Etest) ) - start[2] * cos(start[0]-start[1]);
        this->setInitial(start[0], start[1], start[2], start[3]);
        this->swing(h, t);
        this->saveC(fname+"_"+to_string(i)+".dat");
		cout << (i+1)/80.0 * 100 << "%" << endl;
    }
}
//a2.cpp mit main
int main()
{
    Dpendulum pendulum;
	cout << "Erstelle a1.dat" << endl;
    pendulum.doEverything(0.1, sqrt(2)*0.1, 0, 0, 1e-5, 3, "a1.dat");
	cout << "Erstelle a2.dat" << endl;
    pendulum.doEverything(0, 0, 0, 4.472, 1e-5, 3, "a2.dat");
	cout << "Erstelle a3.dat" << endl;
    pendulum.doEverything(0, 0, 0, 11.832, 1e-5, 3, "a3.dat");

    double eps = pow(10,-2);

	cout << "Erstelle b2.dat" << endl;
    pendulum.doEverything(eps, 0, 0, 4.472, 1e-5, 3, "b2.dat");
	cout << "Erstelle b3.dat" << endl;
    pendulum.doEverything(eps, 0, 0, 11.832, 1e-5, 3, "b3.dat");

	cout << "Erstelle c1.dat" << endl;
	pendulum.teilC(3, 1e-5, 50, "c1");
	cout << "Erstelle c2.dat" << endl;
	pendulum.teilC(10, 1e-5, 50, "c2");
	cout << "Erstelle c3.dat" << endl;
	pendulum.teilC(20, 1e-5, 50, "c3");
    return 0;
}
