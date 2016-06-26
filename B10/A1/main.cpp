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
    double t, intSize, omega;
    // t: Zeit, intSize: Größe des Intervalls bzw Lebensdauer, omega: Formel 1
    int count, count2;
    // count: Zählt für phi1 in welchem Interval man sich befindet, count2: Zählt für phi2
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0, 2*M_PI);
    // gleichverteilte Zufallszahl für die Berechnung von dem neuen Theta
    
    for ( int i = 0; i < N; i++) // alle Teilchen durchgehen
    {
        count = 0; // alles auf Null setzen
        count2 = 0;
        fill(phi1.begin(), phi1.end(), 0);
        fill(phi2.begin(), phi2.end(), 0);
        t = 0;
        
        do
        {
            intSize = tLife(tau);
            omega = 0.5*( 3 * pow(cos(theta[i]), 2) - 1 );
            if ( t > tp[count] ) // True --> ins nächste Interval gehen
            {
                phi1[count] += omega*(intSize - ( t - tp[count] )); // Restteil berechnen
                phi2[count] += omega*( t - tp[count] ); // phi2 bekommt die andere Hälfte vom neuen Interval
                if ( count < tp.size()-1 )
                {
                    phi1[count+1] = phi1[count]; // Interval [0, tp_n+1] ist Teilinterval von [0, tp]
                    phi1[count+1] += omega*( t - tp[count] );
                    count++;
                    while ( (t > tp[count]) && (count < tp.size()-1) )
                    {
                        phi1[count+1] = phi1[count];
                        count++;
                    }
                }
            }
            else
            {
                phi1[count] += omega*intSize;
            }
            
            if ( (t > tp[count2]) && (t < 2*tp[count2]) )
            {
                phi2[count2] += omega*intSize;
            }
            else if ( t > tp[count2] )
            {
                phi2[count2] += omega*(intSize - ( t - tp[count2] ));
                if ( count2 < tp.size()-1 )
                {
                    count2++;
                    if ( t > tp[count2] )
                    {
                        phi2[count2] = phi2[count2-1]+phi1[count2-1]-phi1[count2];
                        phi2[count2] += omega*( t - tp[count2-1] );
                    }
                    while ( t > tp[count2] && (count2 < tp.size()-1) )
                    {
                        phi2[count2+1] = phi2[count2];
                        count2++;
                    }
                }
            }

            
            t += intSize;
            theta[i] = acos( cos(alpha)*cos(theta[i]) - sin(alpha)*sin(theta[i])*cos(dist(mt)) );
            
            
        } while ( t <= 2*tp.back() );
        for ( int j = 0; j < tp.size(); j++)
        {
            magnetization[j] += cos(phi1[j] - phi2[j])/N;
        }
    }
    
    ofstream fout("c.dat");
    fout << "#tp\tM" << endl;
    for ( int i = 0; i < tp.size(); i++)
        fout << tp[i] << "\t" << magnetization[i] << endl;
    fout.close();
}

int main()
{/*
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
   */ 
    // c)
    vector<double> tp(30);
    for ( int i = 0; i < tp.size(); i++)
        tp[i] = pow(10, 3.0*i/(tp.size() - 1) - 1);
    simulate(tp, 0.75, 10*M_PI/180, 1e3);
    // simulate(tp, 0.05, 10*M_PI/180, 1e3); // das sieht super aus, also scheint der Fehler in einer der nervigen verschachtelten if Abfragen zu liegen.. irgendwo überlappen sich noch die Intervalle..
    return 0;
}
