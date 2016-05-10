#include <functional>
#include <cmath>
#include "../A1/Dpendulum.h"


int main(){
    Dpendulum pendulum;
    pendulum.doEverything(0.1, sqrt(2)*0.1, 0, 0, 1e-5, 3, "test2ap1.dat");
    pendulum.doEverything(0.1, -sqrt(2)*0.1, 0, 0, 1e-5, 3, "test2ap2.dat");
    pendulum.doEverything(0.1, sqrt(2)*0.1, 0, 0, 1e-3, 3, "test2bp1.dat");//störung aus 2b für beide AB
    pendulum.doEverything(0.1, -sqrt(2)*0.1, 0, 0, 1e-3, 3, "test2bp2.dat");
    return 0;
}
