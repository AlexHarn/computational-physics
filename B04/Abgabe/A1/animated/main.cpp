#include <cmath>
#include <thread>
#include <iostream>
#include <SFML/Graphics.hpp>
#include "Dpendulum.h"

using namespace sf;

Dpendulum pendulum;
bool active = true, drawing = false;
double** cords = NULL;
double* offset = new double[2] {250, 250};
double scale = 100;

//   Hinweis: Weil ich keine Zeit hatte vernünftige Zeitdarstellung zu implementieren, muss je nach
//   CPU die Schrittweite so eingestellt werden, dass die Animation in einer vernünftigen
//   Geschwindigkeit läuft. Das ist natürlich unnötig genau für
//   schnellere CPUs...

void work()
{
    // Zwei pointer, die immer umgebogen werden, um zu vermeiden, dass der mainthread "unfertige"
    // Werte liest (weiß nicht ob das wirklich nötig ist)
    double* temp;
    double* xy = new double[4];
    double* xy1 = new double[4];
    while ( active )
    {
        pendulum.swing(1e-6, 1000);
        temp = xy;
        pendulum.getXY(temp);
        //std::cout << temp[3] << std::endl;

        // Werte an Fenster anpassen
        for ( int i = 0; i<4; i++ )
            temp[i] *= scale;
        temp[0] += offset[0];
        temp[1] += offset[1];
        temp[2] += offset[0];
        temp[3] += offset[1];

        cords = &temp;
        xy = xy1;
        xy1 = temp;
    }
}

int main()
{
    // Pendel initiallisieren
    pendulum.setInitial(3.141, 3.141, 0, 0); // Chaotisch
    //pendulum.setInitial(0.1, sqrt(2)*0.1, 0, 0); // Paralele Mode
    //pendulum.setInitial(0.1, -sqrt(2)*0.1, 0, 0); // Antiparalele Mode

    // Fenster vorbereiten
    RenderWindow window(VideoMode(500, 500), "Doppelpendel");
    CircleShape m1(10), m2(10);
    m1.setFillColor(Color::Green);
    m2.setFillColor(Color::Red);
    m1.setOrigin(10, 10);
    m2.setOrigin(10, 10);
    sf::VertexArray lines(sf::LinesStrip, 3);
    lines[0].position = sf::Vector2f(offset[0], offset[1]);
    window.draw(lines);

    // Workthread starten
    std::thread worker(work);
    // Auf Workthread warten
    while ( cords==NULL )
        continue;

    while ( window.isOpen() )
    {
        Event event;
        while ( window.pollEvent(event) )
            if ( event.type == Event::Closed )
                window.close();

        window.clear();

        lines[1].position = sf::Vector2f((*cords)[0], (*cords)[1]);
        lines[2].position = sf::Vector2f((*cords)[2], (*cords)[3]);
        m1.setPosition((*cords)[0], (*cords)[1]);
        m2.setPosition((*cords)[2], (*cords)[3]);

        window.draw(lines);
        window.draw(m1);
        window.draw(m2);
        window.display();
    }
    active = false;
    worker.join();

    return 0;
}
