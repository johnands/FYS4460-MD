#include "porosities.h"
#include "../atom.h"
#include <iostream>

Porosities::Porosities(System *system) {
    m_system = system;
}


void Porosities::computePorosity() {

    int movingAtoms = 0;
    for (Atom *atom : m_system->atoms()) {
        if (atom->movingAtom()) {
            movingAtoms++;
        }
    }

    m_numberOfMovingAtoms = movingAtoms;
    m_porosity = (double) movingAtoms / m_system->atoms().size();
    m_system->setPorosity(m_porosity);
    std::cout << "Porosity: " << m_porosity << std::endl;
}


void Porosities::halfDensity() {

    // remove half of the fluid atoms
    bool everySecondAtom = true;
    for (int i=0; i < m_system->atoms().size(); i++) {
        if (m_system->atoms()[i]->movingAtom()) {
            if (everySecondAtom) {
                m_system->removeAtom(i);
                i--;
                everySecondAtom = false;
            }
            else {
                everySecondAtom = true;
            }
        }
    }
    std::cout << "N with half-density " << m_system->atoms().size() << std::endl;
}
