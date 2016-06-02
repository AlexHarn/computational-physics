//-------------------------------------------------------------------------------
//                          Computational Physics 2016
//                                Übungsblatt 6
//-------------------------------------------------------------------------------
//       Lösung der Poission Gleichung auf einem Quadrat mit konstanten
//       Randbedingungen und verschiedenen Ladungsverteilungen.
//-------------------------------------------------------------------------------

#include <eigen3/Eigen/Core>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

typedef Eigen::MatrixXd M;

class PoissionRect
/* Löst die Poissiongleichung auf einem Quadrat mit gegebenen Maßen,
 * Randbedingungen und gegebener Ladungsverteilung.
 */
{
    private:
        M rho, phi, ex, ey;
        double delta, eps, influ;

        void calcP();
        /* Berechnet das Potential mittels Gauß-Seidel-Iteration
         * und Genauigkeit eps
         */
        void calcE();
        /* Berechnet das E-Feld aus dem Potential
         */
        void calcInflu();
        /* Berechnet die influenzierte Ladung auf dem Rand
         */

    public:
        PoissionRect(double lx, double ly, double delta, double eps);
        /* Initiallisiert ein Quadrat [0, lx] x [0, ly] mit Diskretisierung delta
         * und Genauigkeit eps
         */
        void calc();
        /* Wrapper Methode für alle Berechnungen
         */
        void setConstBC(double right, double top, double left, double bottom);
        /* Setzt die RB an den 4 Rändern auf die jeweils gegebene Konstante.
         * Zu Beachten: Die Reihenfolge ist eine Prioritätenliste, die Ecken werden
         * entsprechend überschrieben!
         */
        void addQ(double x, double y, double Q);
        /* Fügt eine Ladung am Ort (x, y) mit Betrag Q
         */
        void save(std::string name);
        /* Speichert die aktuellen Ergebnisse in Dateien mit Namen "[name]_[postfix].dat"
         */
        void reset();
        /* Setzt das Quadrat zurück
         */
};

using namespace std;
//-------------------------------------------------------------------------------
//                             CLASS PoissionRect
//                        PRIVATE IMPLEMENTATION SECTION
//-------------------------------------------------------------------------------
void PoissionRect::calcP()
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

void PoissionRect::calcE()
{
    #pragma omp parallel

    // Feld innerhalb des Gebietes mit zweiseitigem Differenzenqoutienten berechnen
    #pragma omp for schedule(static, 1) collapse(2)
	for ( int x = 1; x < phi.rows() - 1; x++ )
	{
		for ( int y = 1; y < phi.cols() - 1; y++ )
		{
			ex(x, y) = ( phi(x + 1, y) - phi(x - 1, y) )/( 2*delta );
			ey(x, y) = ( phi(x, y + 1) - phi(x, y - 1) )/( 2*delta );
		}
	}

    /*
     * Auf dem Rand gehts nur mit dem einseitigen Differenzenqoutienten
     * auf den Seiten ist nur die Normalkomponente von 0 verschieden
     */
    // Oben und Unten
    #pragma omp for
    for ( int x = 0; x < phi.rows(); x++ )
    {
        ey(x, phi.cols() - 1) = ( phi(x, phi.cols() - 1) - phi(x, phi.cols() - 2) )/delta;
        ey(x, 0) = ( phi(x, 1) - phi(x, 0) )/delta;
    }
    // Links und rechts
    #pragma omp for
    for ( int y = 0; y < phi.cols(); y++ )
    {
        ex(0, y) = ( phi(1, y) - phi(0, y) )/delta;
        ex(phi.cols() - 1, y) = ( phi(phi.cols() - 1,y) - phi(phi.cols() - 2, y) )/delta;
    }
}

void PoissionRect::calcInflu()
// Die Methode kann man mit Eigen bestimmt irgendwie eleganter implementieren
{
    influ = 0;
    /*
     * Da die Klasse nur rechteckige Gebiete behandelt kann einfach über alle 4 Seiten summiert werden
     * Vorgehensweise für sigma:
     * Alle Seiten im mathematisch positivem sinn durchnummerieren, oben start. -E_n ist dann jeweils:
     * Oben: -ey, Unten: +ey, Links: ex, Rechts: -ex
     * Vorgehensweise für die Influenzladung:
     * Integral wird zur Summe über sigma * delta, dabei vorzeichen für mathematisch positive
     * Richtung beachten.
     */

    // Oben und unten
    #pragma parallel for reduction(+:influ)
    for ( unsigned x = 0; x<ex.rows(); x++ )
    {
        influ -= ey(x, ey.cols() - 1);
        influ += ey(x, 0);
    }
    // Links und rechts
    #pragma parallel for reduction(+:influ)
    for ( unsigned y = 0; y<ex.cols(); y++ )
    {
        influ += ex(0, y);
        influ -= ex(ex.rows() - 1, y);
    }   

    //Linienintegral
    influ *= delta;
}

//-------------------------------------------------------------------------------
//                             CLASS PoissionRect
//                       PUBLIC IMPLEMENTATION SECTION
//-------------------------------------------------------------------------------
PoissionRect::PoissionRect(double lx, double ly, double delta, double eps)
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

void PoissionRect::calc()
{
    calcP();
    calcE();
    calcInflu();
}

void PoissionRect::setConstBC(double top, double bottom, double right, double left)
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

void PoissionRect::addQ(double x, double y, double Q)
{
    rho((int) ( x/delta ), (int) ( y/delta )) = Q;
}

void PoissionRect::save(string name)
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
    cout << "Influenzierte Ladung für "+name+": " << influ << endl;
}

void PoissionRect::reset()
{
    rho.setZero();
    phi.setZero();
    ex.setZero();
    ey.setZero();
    influ = 0;
}

int main()
{
    PoissionRect square(1, 1, 0.05, 1e-5);

    // a)
    square.calc();
    square.save("a");

    //// b)
    square.setConstBC(1, 0, 0, 0);
    square.calc();
    square.save("b");

    // c)
    square.reset();
    square.addQ(0.5, 0.5, 1);
    square.calc();
    square.save("c");

    // e)
    //  Dipol
    square.reset();
    square.addQ(0.25, 0.25, 1);
    square.addQ(0.75, 0.75, -1);
    square.calc();
    square.save("e-dipol");
    //  Quadrupol
    square.reset();
    square.addQ(0.25, 0.25, 1);
    square.addQ(0.25, 0.75, -1);
    square.addQ(0.75, 0.75, 1);
    square.addQ(0.75, 0.25, -1);
    square.calc();
    square.save("e-quadrupol");

    return 0;
}
