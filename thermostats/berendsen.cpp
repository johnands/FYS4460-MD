#include "berendsen.h"
#include "atom.h"
#include "../system.h"
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

Berendsen::Berendsen(System &system, double temperatureHeatBath, double tau) :
        Thermostat(system, temperatureHeatBath, tau) {
}


void Berendsen::applyThermostat(double temperature) {
    // rescale velocities of all atoms to keep a constant temperature

    // velocity scale factor
    double gamma = sqrt(1 + m_dtTau*(m_temperatureHeatBath/temperature - 1));
    //cout << m_dtTau << endl;

    // rescale all velocitites
    for (Atom *atom : m_system.atoms()) {
        atom->velocity *= gamma;
    }

}
