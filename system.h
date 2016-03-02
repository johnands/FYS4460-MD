#ifndef SYSTEM_H
#define SYSTEM_H
#include "atom.h"
#include "math/vec3.h"
#include <vector>

class Potential; class Integrator;
using std::vector;

class System
{
private:
    vec3 m_systemSize;
    vec3 m_systemSizeHalf;
    bool m_periodicBoundaries;
    vector<Atom*> m_atoms;
    Potential* m_potential = nullptr;
    Integrator* m_integrator = nullptr;
    double m_time = 0;
    int m_steps = 0;


public:
    System();
    //~System();
    void resetForcesOnAllAtoms();
    void createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature, double mass,
                          bool BoltzmannDist, double maxMinVelocity);
    void applyPeriodicBoundaryConditions();
    void removeTotalMomentum();
    void calculateForces();
    void step(double dt);


    // Setters and getters
    vector<Atom *>& atoms() { return m_atoms; } // Returns a reference to the std::vector of atom pointers
    vec3 systemSize() { return m_systemSize; }
    void setSystemSize(vec3 systemSize) { m_systemSize = systemSize; }
    vec3 systemSizeHalf() { return m_systemSizeHalf; }
    void setSystemSizeHalf(vec3 systemSizeHalf) { m_systemSizeHalf = systemSizeHalf; }
    void setPeriodicBoundaries(bool periodicBoundaries) { m_periodicBoundaries = periodicBoundaries; }
    bool getPeriodicBoundaries() { return m_periodicBoundaries; }

    Potential *potential() { return m_potential; }
    void setPotential(Potential *potential) { m_potential = potential; }
    double time() { return m_time; }
    void setTime(double time) { m_time = time; }
    Integrator *integrator() { return m_integrator; }
    void setIntegrator(Integrator *integrator) { m_integrator = integrator; }
    int steps() { return m_steps; }
    void setSteps(int steps) { m_steps = steps; }
    double volume() { return m_systemSize[0]*m_systemSize[1]*m_systemSize[2]; }
};
#endif
