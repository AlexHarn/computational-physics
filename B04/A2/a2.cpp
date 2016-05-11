#include <iostream>
#include <cmath>
#include "Dpendulum.h"

using namespace std;

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
