/* ----------------------------------------------------------------------------- 
 *        Computational Physics 2016 Blatt 01 Aufgabe 1: Drehmomente
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
     * theta: Winkel zwischen m_00 und y-achse
     * k, l: Position von m_kl
     * anti: false für ferro- true für antiferromagnetischen Fall
     */
    int af = anti ? -1 : 1;
    if ( k==l && k==0 )
        return 0;
    // Vorfaktor
    double factor = pow(10,-7)/( pow(pow(k, 2) + pow(l, 2),3/2)*pow(a, 3) );
    // Winkel zwischen y-Achse und R
    double phi=0;
    if ( k>0 )
        phi = M_PI/2 - atan(l/k);
    else if ( k<0 && l>0 )
        phi = 1.5*M_PI - atan(l/k);
    else if ( k<0 && l<=0 )
        phi = - M_PI/2 - atan(l/k);
    else if ( k==0 && l>0 )
        phi = 0;
    else 
        phi = M_PI;
    // Skalarprodukt R_norm*m berechnen
    double Rm = M*cos(theta - phi);
    // Skalarprodukt R_norm*n berechnen
    double Rn = af*M*cos(phi);
    // Skalarprodukt n*m
    double mn = af*pow(M, 2)*cos(theta);
    return factor*( -3*Rm*Rn + mn );
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
        datafile.open("a1_N="+to_string(N)+".dat");
        datafile << "Theta\tE_f\tT_f_abl\tT_f_dreh\tE_af\tT_af_abl\tT_af_dreh" << endl;
        // Gesamtenergie initiallisieren
        double E = 0; 
        // gesamten Raumwinkel für theta mit Schrittweite h durchgehen 
        for ( double theta = 0; theta<2*M_PI; theta += h )
        {
            E = 0;
            /*
             * ferromagnetischer Fall 
             * Ja man könnte beides in eine Doppelschleife packen,
             * so war es aber einfacher zu entwickeln und ist überischtlicher
             * wenn auch etwas ineffizienter.
             */
            for ( int k = -N; k <= N; k++ ) 
                for ( int l = -N; l <= N; l++ ) 
                    E += energy(theta, k, l);
            datafile << 180/M_PI*theta << "\t" << E << "\t\t\t";  
            /*
             * antiferromagnetischer Fall
             */
            E = 0;
            for ( int k = -N; k <= N; k++ ) 
                for ( int l = -N; l <= N; l++ ) 
                        E += energy(theta, k, l, (bool)( l+k % 2 ));
            datafile  << E;     

            datafile << endl;
        }
        datafile.close();
    }
    return 0;
}
