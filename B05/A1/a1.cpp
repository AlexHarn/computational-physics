#include <math.h>
#include <iostream>
#include <string> 
#include "../../Eigen/Eigen/Dense" //Eigen liegt oberhalb der B**-Ordner
#include "initializeSquareGrid.h" 

using namespace Eigen;
using namespace std;

#ifndef MDSimulation_H
#define MDSimulation_H
class MDSimulation {    
    //const double m = 1; //Masse
    //const double N = 6000; //Teilanzahl
    public:
        //~MDSimulation();
        /*void MachAlles(Vector2d r, Vector2d v);
        double ljpotential(double r);*/
};
#endif

/*double MDSimulation::ljpotential(double r) {
    return 4*(pow(pow(r,-1),12)-pow(pow(r,-1),6));
}

void MDSimulation::MachAlles(Vector2d r, Vector2d v, Vector2d F) {
    cout << to_string(this->ljpotential(5)) << endl;
    cout << "Done." << endl;   
}*/

int main() {
    int AnzTeilchen = 1000;
    double dichte = 40;
    QuadrGitter quadrGitter;
    quadrGitter.ErstelleQuadrGitter(AnzTeilchen, dichte);
    return 0; //nicht zwingend notwendig, aber was soll's
}   
