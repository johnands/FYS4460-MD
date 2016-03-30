#include "eulercromer.h"
#include "../system.h"

EulerCromer::EulerCromer(System *system) :
    Integrator(system) {

}


void EulerCromer::integrate(double dt)
{
    m_system->calculateForces();
    for (Atom *atom : m_system->atoms()) {
        atom->velocity += atom->force*dt / atom->mass();
        atom->position += atom->velocity*dt;
    }

    m_system->applyPeriodicBoundaryConditions();
}
