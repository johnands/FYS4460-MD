#include "spheres.h"
#include "../atom.h"
#include "../math/random.h"
#include <iostream>

Spheres::Spheres(System &system) :
    Porosities(system) {
    makePores();
}

void Spheres::makePores() {

    for (int i=0; i < m_numberOfPores; i++) {

        // center of random sphere
        vec3 sphereCenter;
        sphereCenter[0] = Random::nextDouble()*m_system.systemSize().x();
        sphereCenter[1] = Random::nextDouble()*m_system.systemSize().y();
        sphereCenter[2] = Random::nextDouble()*m_system.systemSize().z();

        // random radius 2nm-3nm
        double sphereRadius = m_radius0 + Random::nextDouble()*(m_radius1 - m_radius0);

        // freeze all atoms inside sphere
        for (Atom *atom : m_system.atoms()) {
            vec3 dr = atom->position - sphereCenter;
            double distanceFromCenter = dr.length();

            if (distanceFromCenter < sphereRadius) {
                atom->setMovingAtom(false);
                atom->setName("NM");
            }
        }
    }
    computePorosity();
    halfDensity();
}


