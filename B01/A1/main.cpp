/* -----------------------------------------------------------------------------
 *        Computational Physics 2016 Blatt 01 Aufgabe 1: Drehmomente
 *                            !!! C++11 Standard !!!
 * -----------------------------------------------------------------------------
 *        2N+1 x 2N+1 fixierte Magnetische Momente. Berechnet wird das
 *        Drehmoment auf das Moment in der Mitte in Abhängigkeit zu dessen
 *        Winkel mit der y-Achse im ferro- (i) und antiferromagnetischen (ii)
 *        Fall.
 * ---------------------------------------------------------------------------*/
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES

using namespace std;

// globale Parameter
const double a = 0.01; // Gitterkonstante
const double M = 1; // Betrag der Momente

double energy(double theta, int k, int l, bool anti = false)
{
    /*
     * Berechnet die Energie von n_00 im Feld von n_kl
     * theta: Winkel zwischen m=n_00 und y-achse
     * k, l: Position von m_kl
     * anti: false für ferro- true für antiferromagnetischen Fall
     */
    // Energie im eigenen Feld ist 0
    if ( k==l && k==0 )
        return 0;
    // ferro/antiferro wert für n_kl
    int af = anti ? -1 : 1;
    // R Betrag
    double R = a*sqrt(pow(k, 2) + pow(l, 2));
    // m in kartesischen Koordinaten
    double m_x = M*sin(theta);
    double m_y = M*cos(theta);

    // Skalarprodukt R_norm*m berechnen
    double Rm = a*( m_x*k + m_y*l )/R;
    // Skalarprodukt R_norm*n berechnen
    double Rn = M*af*a*l/R;
    // Skalarprodukt n*m
    double mn = af*M*m_y;

    return 1e-7*pow(R, -3)*( -3*Rm*Rn + mn );
}

double sumE(double theta, int N, bool anti)
{
    /*
     * Berechnet die Gesamtenergie von m_00
     */
    double E = 0;
    if ( !anti )
    {
        for ( int k = -N; k <= N; k++ )
            for ( int l = -N; l <= N; l++ )
                E += energy(theta, k, l);
    } else
    {
        for ( int k = -N; k <= N; k++ )
            for ( int l = -N; l <= N; l++ )
                E += energy(theta, k, l, (bool)( ( l + k ) % 2 ));
    }
    return E;
}

double dE_dTheta(double h, double (*E)(double, int, bool), double theta,  int N, bool anti)
{
   /*
    * Berechnet die numerische Ableitung von E nach theta
    */
   return ( E(theta + h, N, anti) - E(theta - h, N, anti) )/( 2*h );
}

double T(double theta, int N, bool anti)
{
    /*
     * Berechnet T = m x B(0)
     */
    // m_00 in kartesischen Koordinaten
    double m_x = M*sin(theta);
    double m_y = M*cos(theta);
    // Variablen initiallisieren
    double B_x, B_y, factor;
    B_x = B_y = 0;
    int af = 1;

    // alle durchgehen
    for ( int k = -N; k<=N; k++ )
        for ( int l = -N; l<=N; l++ )
        {
            // Fallunterscheidung
            if (anti && ( ( l + k ) % 2 )==0)
                af = -1;
            else
                af = 1;
            // es wirkt kein Drehmoment durch das eigene Feld
            if ( k==0 && k==l )
                continue;
            // Abstandsbetrag
            double R = a*sqrt(pow(k, 2) + pow(l, 2));
            // Skalarprodukt R*n berechnen
            double Rn = a*l*af*M;
            /* Drehmomente aufsummieren, dabei ist
             * n = (0, af*M)^T
             * R = a*(k, l)^T
             */
            factor = 1e-7*pow(R, -5);
            B_x += factor*( 3*a*k*Rn /* - 0 */ );
            B_y += factor*( 3*a*l*Rn - af*M*pow(R, 2) );
        }
    /* T zeig in z-Richtung. Man könnte die Berechnung von T auch
     * direkt in die Schleife einarbeiten und B_x, B_y nicht getrennt speichern
     * aber ich finde es so übersichtlicher hier spielt Effizienz jetzt keine
     * große Rolle.
     */
    return m_x*B_y-m_y*B_x;
}

int main()
{
    vector<int> Ns; // Verschiedene Werte für N
    Ns.insert(Ns.end(), { 2, 5, 10 });
    double h = M_PI/180; // Schrittweite für theta
    for ( auto &N : Ns )
    {
        // output file vorbereiten
        ofstream datafile;
        datafile.open("N="+to_string(N)+".dat");
        datafile << "Theta\tE_f\tT_f_abl\tT_f_kreuz\tE_af\tT_af_abl\tT_af_kreuz" << endl;
        // gesamten Raumwinkel für theta mit Schrittweite h durchgehen
        for ( double theta = 0; theta<2*M_PI; theta += h )
        {
            // ferromagnetischer Fall
            datafile << 180/M_PI*theta << "\t"
                     << sumE(theta, N, false) << "\t"
                     << abs(dE_dTheta(1e-7, &sumE, theta, N, false)) << "\t"
                     << abs(T(theta, N, false)) << "\t";
            // antiferromagnetischer Fall
            datafile << sumE(theta, N, true) << "\t"
                     << abs(dE_dTheta(1e-7, &sumE, theta, N, true)) << "\t"
                     << abs(T(theta, N, true))
                     << endl;
        }
        datafile.close();
    }
    return 0;
}
