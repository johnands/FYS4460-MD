#include "berendsen.h"
#include "atom.h"
#include "../system.h"
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

Berendsen::Berendsen(System &system, double temperature, double dtTau) :
        Thermostat(system) {
    m_temperatureHeatBath = temperature;
    m_dtTau = dtTau;
}


void Berendsen::applyThermostat(double temperature) {
    // rescale velocities of all atoms to keep a constant temperature

    // velocity scale factor
    double gamma = sqrt(1 + m_dtTau*(m_temperatureHeatBath/temperature - 1));

    cout << "gamma = " << gamma << endl;
    cout << "m_dtTAu = " << m_dtTau << endl;
    cout << "heat bath = " << m_temperatureHeatBath << endl;
    cout << "temp = " << temperature << endl;

    // rescale all velocitites
    for (Atom *atom : m_system.atoms()) {
        atom->velocity *= gamma;
    }

}
