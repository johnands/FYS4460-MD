#include "velocityverlet.h"
#include "../system.h"
#include <iostream>

using std::cout;
using std::endl;

void VelocityVerlet::integrate(System *system, double dt)
{
    double dtHalf = dt/2.0;

    // calculate forces on initial state
    if (m_firstStep) {
        system->calculateForces();
        m_firstStep = false;
    }

    for (Atom *atom : system->atoms()) {
        atom->velocity += atom->force*dtHalf / atom->mass();
        atom->position += atom->velocity*dt;
    }


    system->applyPeriodicBoundaryConditions();
    system->calculateForces();

    for (Atom *atom : system->atoms()) {
        atom->velocity += atom->force*dtHalf / atom->mass();
    }


}
