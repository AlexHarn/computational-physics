#include <eigen3/Eigen/Core>
#include <cmath>
#include <fstream>
#include <iostream>
#include "poisson.h"

using namespace std;
//-------------------------------------------------------------------------------
//                             CLASS PoissionSquare
//                        PRIVATE IMPLEMENTATION SECTION
//-------------------------------------------------------------------------------
void PoissionSquare::calcP()
{
    double max;
    double temp = 0;
    do
    {
        max = 0;
        for ( int i = 1; i< phi.rows() - 1; i++ )
        {
            for ( int j = 1; j < phi.cols() - 1;  j++ )
            {
                temp = 0.25*( phi(i + 1, j) + phi(i - 1, j) + phi(i, j + 1) + phi(i, j - 1)) - 0.25*delta*delta*rho(i, j);
                if ( fabs(temp - phi(i, j)) > max )
                    max = fabs(temp);
                phi(i, j) = temp;
            }
        }
    } while ( max >= eps );
}
void PoissionSquare::calcE()
{
	for ( int x = 1; x < phi.rows() - 1; x++ )
	{
		for (int y = 1; y < phi.cols() - 1; y++)
		{
			ex(x, y) = ( phi(x + 1, y) - phi(x - 1, y) )/( 2*delta );
			ey(x, y) = ( phi(x, y + 1) - phi(x, y - 1) )/( 2*delta );
		}
	}
}
//-------------------------------------------------------------------------------
//                             CLASS PoissionSquare
//                       PUBLIC IMPLEMENTATION SECTION
//-------------------------------------------------------------------------------
PoissionSquare::PoissionSquare(double lx, double ly, double delta, double eps)
{
    double J = lx/delta, L = ly/delta;
    if ( floor(J) != J || floor(L) != L )
        throw invalid_argument("Die gegebene Abmessungen lassen sich nicht mit gegebenem delta diskretisieren!");
    this->delta = delta;
    this->eps = eps;
    rho = M::Zero(J, L);
    phi = M::Zero(J, L);
    ex = M::Zero(J, L);
    ey = M::Zero(J, L);
}

void PoissionSquare::calc()
{
    calcP();
    calcE();
}

void PoissionSquare::setConstBC(double top, double bottom, double right, double left)
{
    for ( int i = 0; i < phi.cols(); i++ )
        phi(0, i) = left;
    for ( int i = 0; i < phi.cols(); i++ )
        phi(phi.rows() - 1, i) = right;
    for ( int i = 0; i < phi.rows(); i++ )
        phi(i, 0) = bottom;
    for ( int i = 0; i < phi.rows(); i++ )
        phi(i, phi.cols() - 1) = top;
}

void PoissionSquare::addQ(double x, double y, double Q)
{
    rho((int) ( x/delta ), (int) ( y/delta )) = Q;
}

void PoissionSquare::save(string name)
{
    ofstream fout;
    fout.open(name+"_phi.dat");
    for ( int j = 0; j < phi.cols(); j++ )
    {
        for ( int i = 0; i < phi.rows(); i++ )
            fout << phi(i, j) << "\t";
        fout << endl;
    }
    fout.close();
    fout.open(name+"_ex.dat");
    for ( int j = 0; j < phi.cols(); j++ )
    {
        for ( int i = 0; i < phi.rows(); i++ )
            fout << ex(i, j) << "\t";
        fout << endl;
    }
    fout.close();
    fout.open(name+"_ey.dat");
    for ( int j = 0; j < phi.cols(); j++ )
    {
        for ( int i = 0; i < phi.rows(); i++ )
            fout << ey(i, j) << "\t";
        fout << endl;
    }
    fout.close();
}

void PoissionSquare::reset()
{
    rho.setZero();
    phi.setZero();
    ex.setZero();
    ey.setZero();
}
