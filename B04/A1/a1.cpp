#include <cmath>
#include "Dpendulum.h"

int main()
{
    Dpendulum pendulum;
    pendulum.doEverything(0.1, sqrt(2)*0.1, 0, 0, 1e-5, 3, "test.dat");
}
