#ifndef SYSTEM_H
#define SYSTEM_H
#include "atom.h"
#include "math/vec3.h"
#include <vector>

class Potential; class Integrator;
class Thermostat; class StatisticsSampler;
class Porosities;

using std::vector;

class System
{
private:
    vec3 m_systemSize;
    vec3 m_systemSizeHalf;
    bool m_periodicBoundaries = true;
    bool m_makeXYZ = false;
    bool m_radialDistribution = false;
    bool m_useThermostat = false;
    bool m_usePores = false;
    bool m_writeSample = false;
    bool m_useExternalForce = false;
    int m_thermalization;
    const char *m_XYZName;
    double m_porosity;
    double m_temperature;
    vector<Atom*> m_atoms;
    Potential* m_potential = nullptr;
    Integrator* m_integrator = nullptr;
    Thermostat* m_thermostat = nullptr;
    Porosities* m_pores = nullptr;
    StatisticsSampler* m_statisticsSampler = nullptr;
    double m_time = 0;
    int m_steps = 0;
    int m_numberOfTimeSteps;
    double m_timeStep;


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
    void runSimulation();
    void readFromStateFile(const char *filename, double mass, double latticeConstant,
                           int numberOfUnitCellsEachDimension);
    void removeAtom(int index);


    // Setters and getters
    vector<Atom *>& atoms() { return m_atoms; } // Returns a reference to the std::vector of atom pointers
    void setAtoms(vector<Atom *>& atoms) { m_atoms = atoms; }
    vec3 systemSize() { return m_systemSize; }
    void setSystemSize(vec3 systemSize) { m_systemSize = systemSize; }
    vec3 systemSizeHalf() { return m_systemSizeHalf; }
    void setSystemSizeHalf(vec3 systemSizeHalf) { m_systemSizeHalf = systemSizeHalf; }
    void setPeriodicBoundaries(bool periodicBoundaries) { m_periodicBoundaries = periodicBoundaries; }
    bool getPeriodicBoundaries() { return m_periodicBoundaries; }
    void setNumberOfTimeSteps(int numberOfTimeSteps) {m_numberOfTimeSteps = numberOfTimeSteps; }
    int getNumberOfTimeSteps() { return m_numberOfTimeSteps; }
    void setMakeXYZ(bool makeXYZ) {m_makeXYZ = makeXYZ; }
    bool getMakeXYZ() { return m_makeXYZ; }
    void setRadialDistribution(bool radialDistribution) {m_radialDistribution = radialDistribution; }
    bool getRadialDistribution() { return m_radialDistribution; }
    void setUseThermoStat(bool useThermostat) { m_useThermostat = useThermostat; }
    bool getUseThermostat() { return m_useThermostat; }
    void setWriteSample(bool writeSample) { m_writeSample = writeSample; }
    bool getWriteSample() { return m_writeSample; }
    void setXYZName(const char *filename) { m_XYZName = filename; }
    const char *getXYZName() { return m_XYZName; }
    void setUseExternalForce(bool useExternalForce) { m_useExternalForce = useExternalForce; }
    bool getUseExternalForce() { return m_useExternalForce; }
    void setTemperature(double temperature) { m_temperature = temperature; }
    double getTemperature() { return m_temperature; }
    void setThermalization(int thermalization) { m_thermalization = thermalization; }
    int getThermalization() { return m_thermalization; }

    Potential *potential() { return m_potential; }
    void setPotential(Potential *potential) { m_potential = potential; }
    double time() { return m_time; }
    void setTime(double time) { m_time = time; }
    void setTimeStep(double timeStep) { m_timeStep = timeStep; }
    double getTimeStep() { return m_timeStep; }
    Integrator *integrator() { return m_integrator; }
    void setIntegrator(Integrator *integrator) { m_integrator = integrator; }
    StatisticsSampler *statisticsSampler() { return m_statisticsSampler; }
    void setStatisticsSampler(StatisticsSampler *statisticsSampler) { m_statisticsSampler = statisticsSampler; }
    Thermostat *getThermostat() { return m_thermostat; }
    void setThermostat(Thermostat *thermostat) { m_thermostat = thermostat; }
    Porosities *getPores() { return m_pores; }
    void setPores(Porosities *porosities) { m_pores = porosities; }
    double getPorosity() { return m_porosity; }
    void setPorosity(double porosity) { m_porosity = porosity; }
    int steps() { return m_steps; }
    void setSteps(int steps) { m_steps = steps; }
    double volume() { return m_systemSize[0]*m_systemSize[1]*m_systemSize[2]; }
};
#endif
