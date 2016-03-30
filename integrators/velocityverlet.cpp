#include "velocityverlet.h"
#include "../system.h"
#include <iostream>

using std::cout;
using std::endl;


VelocityVerlet::VelocityVerlet(System *system) :
    Integrator(system) {

}

void VelocityVerlet::integrate(double dt)
{
    double dtHalf = dt/2.0;

    // calculate forces on initial state
    if (m_firstStep) {
        m_system->calculateForces();
        m_firstStep = false;
    }

    for (Atom *atom : m_system->atoms()) {
        if (atom->movingAtom()) {
            atom->velocity += atom->force*dtHalf / atom->mass();
            atom->position += atom->velocity*dt;
        }
    }

    m_system->applyPeriodicBoundaryConditions();
    m_system->calculateForces();

    for (Atom *atom : m_system->atoms()) {
        if (atom->movingAtom()) {
            atom->velocity += atom->force*dtHalf / atom->mass();
        }
    }


}
