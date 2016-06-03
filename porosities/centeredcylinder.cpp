#include "centeredcylinder.h"
#include "../atom.h"
#include "../system.h"
#include <iostream>

CenteredCylinder::CenteredCylinder(System *system) :
    Porosities(system) {
    makePores();
    std::cout << "poreradius: " << m_poreRadius << std::endl;
}

void CenteredCylinder::makePores() {

    for (Atom *atom : m_system->atoms()) {
        double yFromCenter = atom->position.y() - m_system->systemSizeHalf().y();
        double zFromCenter = atom->position.z() - m_system->systemSizeHalf().z();
        double distanceFromCenter = yFromCenter*yFromCenter + zFromCenter*zFromCenter;

        if (distanceFromCenter > m_poreRadius) {
            atom->setMovingAtom(false);
            atom->setName("NM");
        }
    }
    computePorosity();
    halfDensity();
}

