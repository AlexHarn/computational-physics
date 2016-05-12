//-------------------------------------------------------------------------------
//                      Computational Physics 2016 Blatt 4
//                           Aufgabe 1: Doppelpendel
//-------------------------------------------------------------------------------
//              Anwendung des Runge-Kutta-Verfahren 4. Ordnung auf das
//              nicht lineare, chaotische Problem des Doppelpendels.
//-------------------------------------------------------------------------------



//-------------------------------------------------------------------------------
//                               FILE: rungkutt.h
//-------------------------------------------------------------------------------
#include <functional>

#ifndef RUNGKUTT_H
#define RUNGKUTT_H

void rungkutt(std::function<void(double t, double* y, double* out)> f, const double N, const double h, const double t_0, double** y, const int d);
/* Führt das Runge-Kutta-Verfahren 4. Ordnung durch. Neue Implementierung ohne std:vector und in
 * allgemeiner Formulierung wie in der Vorlesung.
 *      DGL: y'(t) = f(t, y(t))
 *      N: Schrittzahl
 *      h: Schrittweite
 *      t_0: Startzeit
 *      y: y Vektoren zu allen Zeiten. Bei Aufruf müssen Starbedingungen gesetzt sein
 *      d: Länge der y Vektoren (!!! Nicht Anzahl der Raumdimensionen !!!)
 */
#endif
//-------------------------------------------------------------------------------
//                              FILE: Dpendulum.h
//-------------------------------------------------------------------------------
#include <string>

#ifndef DPENDULUM_H
#define DPENDULUM_H

class Dpendulum
{
        void static f(double /*t*/, double* y, double* dy, double g);
        /* f für Runge-Kutta mit y' = f(y)
         *      y: y vektor (input)
         *      dy: y' vektor (output)
         *      g: Gravitationskonstante (wird später gebunden)
         */
        void calcE();
        /* Berechnet U und T aus y
         */
        void calcCartesian();
        /* Berechnet die kartesischen Koordinaten für alle y
         */
        void cartesian(int n, double* xy);
        /* Berechnet die kartesischen Koordinaten für y_n
         */
        void reset();
        /* Setzt das Pendel zurück (löscht die aktuellen Daten)
         */
        const double g = 9.81;
        int N = 0;
        double tN = 0, lh = 0;
        double* y0;
        double** y, ** xy, ** T, ** U;
        bool ready = false, energy = false, cart = false, active = false;

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
        void getXY(double *xy);
        /* Gibt die aktuellen XY-Koordinaten aus
         */
        void save(std::string fname);
        /* Speichert die Daten in Datei fname
         */
        void doEverything(double theta1, double theta2, double ddtTheta1, double ddtTheta2, double h, double t, std::string fname);
        /* Macht alles in einem Aufruf
         */
};

#endif
//-------------------------------------------------------------------------------
//                              FILE: rungkutt.cpp
//-------------------------------------------------------------------------------
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
//-------------------------------------------------------------------------------
//                             FILE: Dpendulum.cpp
//-------------------------------------------------------------------------------
#include <functional>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

// PRIVATE:
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

void Dpendulum::calcE()
{
    if ( energy )
        throw "Energie bereits berechnet!";

    const double m = 1;

    T = new double*[N];
    U = new double*[N];

    #pragma omp parallel for
    for ( int i = 0; i<N; i++ )
    {
        T[i] = new double[2];
        U[i] = new double[2];

        T[i][0] = 0.5*m*y[i][2]*y[i][2];
        U[i][0] = 2*m*g*(1-cos(y[i][0])); // Ursprung so verschoben, dass U = 0 für theta = 0

        T[i][1] = 0.5*m*( y[i][2]*y[i][2] + y[i][3]*y[i][3] + 2*y[i][2]*y[i][3]*cos(y[i][0] - y[i][1]) );
        U[i][1] = m*g*(1-cos(y[i][1]));
    }
    energy = true;
}

void Dpendulum::calcCartesian()
{
    if ( cart )
        throw "Kartesische Koordinaten bereits berechnet!";

    xy = new double*[N];

    #pragma omp parallel for
    for ( int i = 0; i<N; i++ )
    {
        xy[i] = new double[4];
        cartesian(i, xy[i]);
    }
    cart = true;
}

void Dpendulum::cartesian(int n, double* xy)
{
    xy[0] = sin(y[n][0]);
    xy[1] = cos(y[n][0]);
    xy[2] = xy[0] + sin(y[n][1]);
    xy[3] = xy[1] + cos(y[n][1]);
}

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
        if ( energy )
        {
            delete[] T[i];
            delete[] U[i];
        }
        if ( cart )
            delete[] xy[i];
    }
    delete[] y;
    if ( energy )
    {
        delete[] T;
        delete[] U;
    }
    if ( cart )
       delete[] xy;
    N = 0;
    tN = 0;
    active = false;
    ready = false;
    energy = false;
    cart = false;
}

// PUBLIC:
Dpendulum::~Dpendulum()
{
    reset();
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
    swing(h, (int) ceil( t/h )+1);
}

void Dpendulum::swing(double h, int n)
{
    using namespace std::placeholders;

    if ( !ready )
        throw "Noch nicht readyiallisiert!";

    if ( tN > 0 )
    {
        double** ny = new double*[N+n];
        copy(y, y + N, ny);
        delete[] y;
        y = ny;
        for ( int i = N; i<N+n; i++ )
            y[i] = new double[4];
        rungkutt(bind(Dpendulum::f, _1, _2, _3, g), n, h, tN, y+N-1, 4);
    }
    else
    {
        y = new double*[n];
        for ( int i = 1; i<n; i++ )
            y[i] = new double[4];
        y[0] = y0;
        rungkutt(bind(Dpendulum::f, _1, _2, _3, g), n, h, 0, y, 4);
    }
    N += n;
    tN += h*n;
    lh = h;
    active = true;
}

void Dpendulum::getXY(double *xy)
{
    cartesian(N-1, xy);
}

void Dpendulum::save(string fname)
{
    if ( !active )
        throw "Noch nichts berechnet!";
    calcE();
    calcCartesian();
    ofstream fout;
    fout.open(fname);
    for ( int i = 0; i<N; i++ )
    {
        fout << i*lh << "\t";
        for ( int j = 0; j<4; j++ )
            fout << y[i][j] << "\t";
        for ( int j = 0; j<4; j++ )
            fout << xy[i][j] << "\t";
        fout << T[i][0] << "\t" << T[i][1] << "\t";
        fout << U[i][0] << "\t" << U[i][1] << "\t";
        fout << T[i][0]+T[i][1]+U[i][0]+U[i][1] << endl;
    }
    fout.close();
    return;
}

void Dpendulum::doEverything(double theta1, double theta2, double ddtTheta1, double ddtTheta2, double h, double t, std::string fname)
{
    reset();
    setInitial(theta1, theta2, ddtTheta1, ddtTheta2);
    swing(h, t);
    save(fname);
}
//-------------------------------------------------------------------------------
//                                FILE: a1.cpp
//-------------------------------------------------------------------------------
#include <cmath>

int main()
{
    Dpendulum pendulum;
    pendulum.doEverything(0.1, sqrt(2)*0.1, 0, 0, 1e-5, 3, "pos_sqrt2.dat");
    pendulum.doEverything(0.1, -sqrt(2)*0.1, 0, 0, 1e-5, 3, "neg_sqrt2.dat");
}
