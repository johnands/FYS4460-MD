#include "system.h"
#include <iostream>
#include <fstream>
#include "integrators/integrator.h"
#include "potentials/potential.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

using std::cout;
using std::endl;

System::System()
{

}

/*System::~System()
{
    delete m_potential;
    delete m_integrator;
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}*/

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention

    // check if atoms is outside of box [0, b*Nx] x [0, b*Ny] x [0, b*Nz]
    for (Atom *atom : m_atoms) {

        // is position less than zero?
        if (atom->position.x() < 0) { atom->position[0] += m_systemSize.x(); }
        if (atom->position.y() < 0) { atom->position[1] += m_systemSize.y(); }
        if (atom->position.z() < 0) { atom->position[2] += m_systemSize.z(); }

        // is position more than b*N?
        if (atom->position.x() > m_systemSize.x()) { atom->position[0] -= m_systemSize.x(); }
        if (atom->position.y() > m_systemSize.y()) { atom->position[1] -= m_systemSize.y(); }
        if (atom->position.z() > m_systemSize.z()) { atom->position[2] -= m_systemSize.z(); }
    }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero

    vec3 totalMomentum;

    // find total momentum
    for (Atom *atom : m_atoms) {
        totalMomentum += atom->mass() * atom->velocity;
    }

    // remove momentum equally on each atom
    int N = m_atoms.size();
    for (Atom *atom : m_atoms) {
        atom->velocity -= totalMomentum / (atom->mass() * (double) N);
    }
}

void System::resetForcesOnAllAtoms() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature, double mass,
                              bool BoltzmannDist, double maxMinVelocity) {

    double lc = latticeConstant;
    double size = numberOfUnitCellsEachDimension*lc;
    setSystemSize(vec3(size, size, size));
    setSystemSizeHalf(vec3(size/2.0, size/2.0, size/2.0));

    // create FCC lattice
    int index = 0;
    for (int i=0; i < numberOfUnitCellsEachDimension; i++) {
        for (int j=0; j < numberOfUnitCellsEachDimension; j++) {
            for (int k=0; k < numberOfUnitCellsEachDimension; k++) {
                // find position of all four atoms in each cell

                Atom *atom1 = new Atom(mass);
                atom1->position.set(lc*i, lc*j, lc*k);
                atom1->storeInitialPosition();
                if (BoltzmannDist) { atom1->resetVelocityMaxwellian(temperature); }
                else               { atom1->resetVelocityUniform(maxMinVelocity); }
                m_atoms.push_back(atom1);
                atom1->setIndex(index);
                index++;

                Atom *atom2 = new Atom(mass);
                atom2->position.set(lc*(i+0.5), lc*(j+0.5), lc*k);
                atom2->storeInitialPosition();
                if (BoltzmannDist) { atom2->resetVelocityMaxwellian(temperature); }
                else               { atom2->resetVelocityUniform(maxMinVelocity); }
                m_atoms.push_back(atom2);
                atom2->setIndex(index);
                index++;

                Atom *atom3 = new Atom(mass);
                atom3->position.set(lc*i, lc*(j+0.5), lc*(k+0.5));
                atom3->storeInitialPosition();
                if (BoltzmannDist) { atom3->resetVelocityMaxwellian(temperature); }
                else               { atom3->resetVelocityUniform(maxMinVelocity); }
                m_atoms.push_back(atom3);
                atom3->setIndex(index);
                index++;

                Atom *atom4 = new Atom(mass);
                atom4->position.set(lc*(i+0.5), lc*j, lc*(k+0.5));
                atom4->storeInitialPosition();
                if (BoltzmannDist) { atom4->resetVelocityMaxwellian(temperature); }
                else               { atom4->resetVelocityUniform(maxMinVelocity); }
                m_atoms.push_back(atom4);
                atom4->setIndex(index);
                index++;
            }
        }
    }
}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    m_potential->setPotentialEnergy(0.0);
    m_potential->setPressure(0.0);
    m_potential->calculateForces();
}

void System::step(double dt) {
    m_integrator->integrate(this, dt);
    m_steps++;
    m_time += dt;
}
