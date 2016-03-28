#include "andersen.h"
#include "atom.h"
#include "../system.h"
#include <cmath>
#include <iostream>
#include "../math/random.h"

using std::cout;
using std::endl;

Andersen::Andersen(System &system, double temperatureHeatBath, double tau) :
        Thermostat(system, temperatureHeatBath, tau) {
}


void Andersen::applyThermostat(double temperature) {
    // simulate collisions between atoms in system and atoms from
    // heat bath

    for (int i=0; i < m_system.atoms().size(); i++) {
        Atom *atom1 = m_system.atoms()[i];

        for (int j=i+1; j < m_system.atoms().size(); j++) {
            Atom *atom2 = m_system.atoms()[j];

            vec3 dr = atom1->position - atom2->position;
            //cout << "radius" << m_radiusOfAtom << endl;

            if (dr.lengthSquared() < 1.42) {

                if (Random::nextDouble() < m_dtTau) {
                    atom1->resetVelocityMaxwellian(m_temperatureHeatBath);
                }

                if (Random::nextDouble() < m_dtTau) {
                    atom2->resetVelocityMaxwellian(m_temperatureHeatBath);
                }
            }
        }
    }
}
