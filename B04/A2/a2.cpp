#include <cmath>
#include "Dpendulum.h"

int main()
{
    Dpendulum pendulum;
    pendulum.doEverything(0.1, sqrt(2)*0.1, 0, 0, 1e-5, 3, "a1.dat");
    pendulum.doEverything(0, 0, 0, 4.472, 1e-5, 3, "a2.dat");
    pendulum.doEverything(0, 0, 0, 11.832, 1e-5, 3, "a3.dat");

    double eps = pow(10,-2);

    pendulum.doEverything(eps, 0, 0, 4.472, 1e-5, 3, "b2.dat");
    pendulum.doEverything(eps, 0, 0, 11.832, 1e-5, 3, "b3.dat");
    return 0;
}
