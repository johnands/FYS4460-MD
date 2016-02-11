#include "velocityverlet.h"
#include "../system.h"
#include "../atom.h"

void VelocityVerlet::integrate(System *system, double dt)
{
    // calculate forces on initial state
    if (system->m_steps == 1) {system->calculateForces();}

    vec3 v_half;
    for(Atom *atom : system->atoms()) {
        v_half = atom->velocity + atom->force*dt / (2*atom->mass());
        atom->position += v_half*dt;
        system->calculateForces();
        atom->velocity = v_half + atom->force*dt / atom->mass();
    }

    system->applyPeriodicBoundaryConditions();

}
