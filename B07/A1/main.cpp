//-------------------------------------------------------------------------------
//                          Computational Physics 2016
//                                Übungsblatt 7
//-------------------------------------------------------------------------------
//            Implementierung eines linearen kongruenten Generators
//-------------------------------------------------------------------------------
#include "lingen.h"

int main()
{
    LinConGen gen(10);

    // b)
    //(i)
    gen.calcCongruent(1234, 20, 120, 6075, 6074);
    gen.saveHist("bi_hist");

    ////(ii)
    //gen.calc(1234,137,187,256);
    //gen.save("bii");
    //gen.reset();

    ////(iii)
    //gen.calc(123456789,65539,0,2147483648);
    //gen.save("biii");
    //gen.reset();

    ////(iv)
    ////7^5 = 16807, 2^(31) -1 = 2147483647
    //gen.calc(1234,16807,0,2147483647);
    //gen.save("biv");
    //gen.reset();

    ////A2 - ist Generator b(iv)
    //gen.calc(1234,16807,0,2147483647);
    //gen.boxmulleralg();
    //gen.save("2a");
    //gen.reset();

    return 0;
}
