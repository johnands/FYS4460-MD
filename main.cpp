#include "math/random.h"
#include "potentials/lennardjones.h"
#include "potentials/lennardjonescelllist.h"
#include "integrators/eulercromer.h"
#include "integrators/velocityverlet.h"
#include "thermostats/berendsen.h"
#include "thermostats/andersen.h"
#include "porosities/centeredcylinder.h"
#include "porosities/spheres.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include "celllist.h"
#include <iostream>
#include <time.h>

using namespace std;

int main(int numberOfArguments, char **argumentList)
{
    int numberOfUnitCells = 20;
    //double initialTemperature = UnitConverter::temperatureFromSI(300);  // measured in Kelvin
    double initialTemperature = 1.05;
    double latticeConstant    = UnitConverter::lengthFromAngstroms(5.72); // measured in angstroms

    double mass = UnitConverter::massFromSI(6.63352088e-26); // mass of Argon atom

    double dt = 0.01; //UnitConverter::timeFromSI(5e-14); // Measured in seconds

    cout << "One unit of length is "      << UnitConverter::lengthToSI(1.0)      << " meters"        << endl;
    cout << "One unit of velocity is "    << UnitConverter::velocityToSI(1.0)    << " meters/second" << endl;
    cout << "One unit of time is "        << UnitConverter::timeToSI(1.0)        << " seconds"       << endl;
    cout << "One unit of mass is "        << UnitConverter::massToSI(1.0)        << " kg"            << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K"             << endl;
    cout << "One unit of energy is "      << UnitConverter::energyToSI(1.0)      << " J"             << endl;

    bool BoltzmannDist = true;             // initial velocities given by Boltzmann distribution
    double maxMinVelocity = 0.5;           // uniformly distributed velocities [-v, v]
    double tau = 10*dt;

    System *system = new System();
    system->createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature,
                            mass, BoltzmannDist, maxMinVelocity);
    system->setPores(new Spheres(system));

    //system->setPotential(new LennardJones(system, 1.0, 1.0));
    system->setPotential(new LennardJonesCellList(system, 1.0, 1.0, 2.5));
    //system->setIntegrator(new EulerCromer(system));
    system->setTimeStep(dt);
    system->setIntegrator(new VelocityVerlet(system));
    system->setThermostat(new Berendsen(system, initialTemperature, tau));
    system->setUseThermoStat(false);

    system->setNumberOfTimeSteps(1001);
    system->setTemperature(initialTemperature);
    system->removeTotalMomentum();

    //system->setPeriodicBoundaries(true);
    system->setWriteSample(true);
    system->setRadialDistribution(false);
    system->setMakeXYZ(false);
    system->setXYZName("test.xyz");

    system->runSimulation();

    return 0;
}
