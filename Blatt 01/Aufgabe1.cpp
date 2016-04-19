#include<iostream>
#define _USE_MATH_DEFINES
#include <cmath> 
#include <fstream>
using namespace std;

const double a = 0.01; 
const double M = 1;

double energy(double theta, int k, int l) {
  // Vorfaktor
  double factor = pow(10,-7)/(pow(pow(k, 2)+pow(l, 2),3/2)*pow(a, 3));
  // Winkel zwischen y-Achse und R
  double phi=0;
  if( k>0 )
     phi = M_PI/2 - atan(l/k);
  else if( k<0 && l>0 )
      phi = 1.5*M_PI - atan(l/k);
  else if( k<0 && l<=0 )
      phi = - M_PI/2 - atan(l/k);
  else if (k == 0 && l > 0)
      phi = 0;
  else 
      phi = M_PI;
 //  cout << phi/M_PI * 180 << endl;
  // Skalarprodukt R_norm*m berechnen
  double Rm = M*cos(theta - phi);
  // Skalarprodukt R_norm*n berechnen
  double Rn = M*cos(phi);
  // Skalarprodukt n*m
  double mn = pow(M, 2)*cos(theta);
  return factor*(-3*Rm*Rn+mn);
}

int main() { 
  int N = 2;
  double h = M_PI/8; // Schrittweite für theta
  // open filestream
  ofstream datafile;
  datafile.open("dat.txt");
  double E = 0;
  for(double theta = 0; theta < 2*M_PI; theta+=h) {
      E = 0;
      for(int k = -N; k <= N; k++) 
         for(int l = -N; l <= N; l++) 
            if(!(k == 0 && l == 0))
                E += energy(theta, k, l);
      datafile << theta << "\t" << E  << endl;     
  }
  datafile.close();
  return 0;
}
