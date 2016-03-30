#include "centeredcylinder.h"
#include "../atom.h"
#include "../system.h"

CenteredCylinder::CenteredCylinder(System *system) :
    Porosities(system) {
    makePores();
}

void CenteredCylinder::makePores() {

    for (Atom *atom : m_system->atoms()) {
        double xFromCenter = atom->position.x() - m_system->systemSizeHalf().x();
        double yFromCenter = atom->position.y() - m_system->systemSizeHalf().y();
        double distanceFromCenter = xFromCenter*xFromCenter + yFromCenter*yFromCenter;

        if (distanceFromCenter > m_poreRadius) {
            atom->setMovingAtom(false);
            atom->setName("NM");
        }
    }
    computePorosity();
    halfDensity();
}

