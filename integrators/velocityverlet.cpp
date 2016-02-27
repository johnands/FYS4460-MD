#include "velocityverlet.h"
#include "../system.h"
#include "../atom.h"
#include <iostream>

using std::cout;
using std::endl;

void VelocityVerlet::integrate(System *system, double dt)
{
    double dtHalf = dt/2.0;

    // calculate forces on initial state
    if (system->steps() == 0) {
        system->calculateForces();
    }

    vec3 v_half;
    for(Atom *atom : system->atoms()) {
        v_half = atom->velocity + atom->force*dtHalf / atom->mass();
        atom->position += v_half*dt;
    }

    system->applyPeriodicBoundaryConditions();
    system->calculateForces();

    for(Atom *atom : system->atoms()) {
        atom->velocity = v_half + atom->force*dt / atom->mass();
    }



}
