#include <random>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

#ifndef M_PI
#define M_PI 3.1415
#endif

using namespace std;

vector<double> ransin(int N)
{
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0, M_PI);

    int count = 0;
    vector<double> r(N);
    while ( count < N )
    {
        r[count] = dist(mt);
        if ( dist(mt)/( 2*M_PI ) < sin(r[count])/2 )
            count++;
    }
    return r;
}

double tLife(double tau)
{
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0, 1);
    
    double temp;
    do
    {
        temp = dist(mt);
    } while(temp == 0);
    return -tau*log(temp);
}

void simulate(vector<double> tp, double tau, double alpha, int N)
{
    vector<double> theta = ransin(N), magnetization(tp.size(), 0.0), phi1(tp.size(), 0.0), phi2(tp.size(), 0.0);
    // theta: startwerte, magnetization: naja.. formel 6, phi1: phi(0, tp) für alle tp's speichern, phi2: siehe phi1
    double t, intSize = tLife(tau), omega;
    // t: Zeit, intSize: Größe des Intervalls bzw Lebensdauer, omega: Formel 1
    int count, count2;
    // count: Zählt für phi1 in welchem Interval man sich befindet, count2: Zählt für phi2
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0, 2*M_PI);
    // gleichverteilte Zufallszahl für die Berechnung von dem neuen Theta
    
    // #pragma omp parallel for  <-- Das macht noch Probleme, nochmal angucken!
    for ( int i = 0; i < N; i++) // alle Teilchen durchgehen
    {
        count = 0; // alles auf Null setzen
        count2 = 0;
        fill(phi1.begin(), phi1.end(), 0);
        fill(phi2.begin(), phi2.end(), 0);
        t = 0;
        
        do
        {
            omega = 0.5*( 3 * pow(cos(theta[i]), 2) - 1 ); // brauch man noch häufiger
            
            // phi(0, tp) berechnen [phi1 genannt]
            if (t < tp.back()-intSize)
            {
                if ( t > tp[count] ) // True --> ins nächste Interval gehen
                {
                    phi1[count] += omega*(intSize - ( t - tp[count] )); // Restteil berechnen
                    while ( (t > tp[count]) && (count < tp.size()-1) ) // solange count vergrößern, bis der nächste Messpunkt größer als das aktuelle t ist
                    {
                        phi1[count+1] = phi1[count]; // Intervall (0, tp[count+1]) hat das Intervall (0, tp[count]) als Teilintervall
                        if ( t > tp[count+1] ) // Restteil berechnen (siehe oben)
                        {
                            phi1[count+1] += omega*(tp[count+1] - tp[count] );
                        }
                        else
                        {
                            phi1[count+1] += omega*(t - tp[count+1] );
                        }
                        count++;
                    }
                }
                else // wenn t < tp[count] ist --> phi einfach integrieren
                {
                    phi1[count] += omega*intSize;
                }
            }
            
            // phi(tp, 2tp) berechnen [phi2 genannt]
            if ( (t > tp[count2]) && (t < 2*tp[count2]) )
            {
                if ( t - intSize > tp[count2] ) // testen ob ich das erste mal größer als tp bin
                {
                    phi2[count2] += omega*intSize;
                }
                else // wenn ja, dann hier rein --> nur Teilintervall integrieren
                {
                    phi2[count2] = omega*( t - tp[count2] );
                }
            }
            else if ( t > tp[count2] ) // True, wenn t ebenfalls größer als 2tp ist --> das phi2 ist fertig, nächstes phi2 berechnen
            {
                phi2[count2] += omega*(intSize - ( t - 2*tp[count2] )); // Restteil, alles analog wie bei phi1
                while ( t > tp[count2] && (count2 < tp.size()-1) )
                {
                    count2++;
                    phi2[count2] = phi2[count2-1]+phi1[count2-1]-phi1[count2];
                    if ( t > tp[count2] )
                    {
                        phi2[count2] += omega*( 2*tp[count2] - 2*tp[count2-1] );
                    }
                    else
                    {
                        phi2[count2] += omega*( t - 2*tp[count2] );
                    }
                }
            }

            // theta springt
            theta[i] = acos( cos(alpha)*cos(theta[i]) - sin(alpha)*sin(theta[i])*cos(dist(mt)) );
            t += intSize; // Zeit erhöhen
            intSize = tLife(tau); // nächstes intSize (Intervall Größe) berechnen
        } while ( t <= 2*tp.back() ); // <= 2tp.back(), da für phi2 auch bis 2tp integriert werden muss
        for ( int j = 0; j < tp.size(); j++)
        {
            magnetization[j] += cos(phi1[j] - phi2[j]); // Mittelung der Magnetisierung über alle Teilchen (1/N kommt bei der Ausgabe)
        }
    }
    ofstream fout("c.dat");
    fout << "#tp\tM" << endl;
    for ( int i = 0; i < tp.size(); i++)
        fout << tp[i] << "\t" << magnetization[i]/N << endl;
    fout.close();
}

int main()
{
    // a)
    vector<double> theta = ransin(1e6);
    ofstream fout("a.dat");
    for ( auto &t : theta )
        fout << t << endl;
    fout.close();

    // b)
    vector<double> omega(theta.size());
    double s = 0;
    fout.open("b.dat");
    for ( unsigned i = 0; i < omega.size(); i++ )
    {
        omega[i] = 0.5*( 3 * pow(cos(theta[i]), 2) - 1 );
        s += omega[i];
        fout << omega[i] << endl;
    }
    cout << "Omega gemittelt: "<< s/omega.size() << endl;
    fout.close();

    // c)
    vector<double> tp(30);
    for ( int i = 0; i < tp.size(); i++)
        tp[i] = pow(10, 3.0*i/(tp.size() - 1) - 1); // logarithmische Skalierung der Messpunkte von 10^-1 bis 10^2
    simulate(tp, 0.75, 10*M_PI/180, 1e4);

    return 0;
}
