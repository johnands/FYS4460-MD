#include "system.h"
#include <iostream>
#include <fstream>
#include "integrators/integrator.h"
#include "potentials/potential.h"
#include "thermostats/thermostat.h"
#include "porosities/porosities.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "io.h"
#include "math/random.h"
#include <string>

using std::cout;
using std::endl;

System::System() {

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

    // check if atoms is outside of box [0, b*Nx] x [0, b*Ny] x [0, b*Nz]
    for (Atom *atom : m_atoms) {
        if (atom->movingAtom()) {

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

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature,
                              double mass,
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

void System::runSimulation() {

    m_statisticsSampler =  new StatisticsSampler(this);
    IO movie; // To write the state to file

    if (getMakeXYZ()) {
        movie.open(getXYZName());
        // sample initial state
        movie.saveState(this);
    }
    double dt = getTimeStep();

    // sample initial state
    m_statisticsSampler->sample(0);

    cout << "Timestep Time Temperature Pressure Density KineticEnergy PotentialEnergy TotalEnergy" << endl;

    clock_t start, finish;
    start = clock();
    for (int timeStep=0; timeStep < getNumberOfTimeSteps(); timeStep++) {
        step(dt);
        m_statisticsSampler->sample(timeStep);
        if (m_radialDistribution && m_steps > m_thermalization) {
            m_statisticsSampler->sampleRadialDistribution(50);
        }
        if ( !(timeStep % 100) ) {
            // print sample every 100 timesteps
            cout << steps() << "      " << time() << "    "
                 << m_statisticsSampler->temperature() << "    "
                 << m_statisticsSampler->pressure() << "  "  << m_statisticsSampler->density() << "    "
                 << m_statisticsSampler->kineticEnergy() << "     "
                 << m_statisticsSampler->potentialEnergy() << "      "
                 << m_statisticsSampler->totalEnergy() << endl;
        }
        /*if ( m_steps == 1002) {
            { movie.saveState(this); }
        }*/
        /*if ( !(timeStep % 10) || timeStep == 1 ) {
            if (getMakeXYZ()) { movie.saveState(this); }
        }*/
        if ( m_steps > m_thermalization && !(timeStep % 10) ) {
            if (m_makeXYZ) { movie.saveState(this); }
        }
    }
    finish = clock();
    cout << "Time elapsed: " << ((finish-start)/CLOCKS_PER_SEC) << endl;

    if (getMakeXYZ()) { movie.close(); }
}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    m_potential->setPotentialEnergy(0.0);
    m_potential->setPressure(0.0);
    m_potential->calculateForces();
}

void System::step(double dt) {
    m_integrator->integrate(dt);
    if (getUseThermostat() && m_steps < m_thermalization) {
        m_thermostat->applyThermostat(m_statisticsSampler->temperature());
    }
    /*else if (m_steps == 2*m_thermalization) {
        m_pores->makePores();
    }
    else if (m_steps > 2*m_thermalization && m_steps < (2*m_thermalization + m_thermalization)) {
        m_thermostat->applyThermostat(m_statisticsSampler->temperature());
    }*/
    m_steps++;
    m_time += dt;
}

void System::readFromStateFile(const char *filename, double mass, double latticeConstant,
                               int numberOfUnitCellsEachDimension) {

    double lc = latticeConstant;
    double size = numberOfUnitCellsEachDimension*lc;
    setSystemSize(vec3(size, size, size));
    setSystemSizeHalf(vec3(size/2.0, size/2.0, size/2.0));

    // open file
    std::ifstream input;
    input.open(filename, std::ios::in);

    if (input.is_open()) { cout << "yes" << endl; }

    // skip first two lines
    std::string dummyLine;
    std::getline(input, dummyLine);

    // process file
    std::string name;
    double x, y, z;
    double vx, vy, vz;
    int index = 0;
    for ( std::string line; std::getline(input, line); ) {
            input >> name >> x >> y >> z >> vx >> vy >> vz;

            // convert positions to MD units
            x = UnitConverter::lengthFromAngstroms(x);
            y = UnitConverter::lengthFromAngstroms(y);
            z = UnitConverter::lengthFromAngstroms(z);

            // set atom in grid
            Atom *atom = new Atom(mass);
            atom->position.set(x, y, z);
            atom->velocity.set(vx, vy, vz);
            atom->storeInitialPosition();
            m_atoms.push_back(atom);
            atom->setIndex(index);
            index++;
    }
}

void System::removeAtom(int index) {
    m_atoms.erase(m_atoms.begin() + index);
    for (int i=index; i < m_atoms.size(); i++) {
        m_atoms[i]->adjustIndex(-1);
    }
}
