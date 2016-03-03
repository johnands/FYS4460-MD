#include "math/random.h"
#include "potentials/lennardjones.h"
#include "potentials/lennardjonescelllist.h"
#include "integrators/eulercromer.h"
#include "integrators/velocityverlet.h"
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
    int numberOfUnitCells = 8;
    double initialTemperature = UnitConverter::temperatureFromSI(100);  // measured in Kelvin
    double latticeConstant    = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms
    cout << latticeConstant << endl;

    // If a first argument is provided, it is the number of unit cells
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);

    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));

    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3) latticeConstant= UnitConverter::lengthFromAngstroms(atof(argumentList[3]));

    double mass = UnitConverter::massFromSI(6.63352088e-26); // mass of Argon atom

    double dt = 0.01; //UnitConverter::timeFromSI(5e-14); // Measured in seconds

    cout << "One unit of length is "      << UnitConverter::lengthToSI(1.0)      << " meters"        << endl;
    cout << "One unit of velocity is "    << UnitConverter::velocityToSI(1.0)    << " meters/second" << endl;
    cout << "One unit of time is "        << UnitConverter::timeToSI(1.0)        << " seconds"       << endl;
    cout << "One unit of mass is "        << UnitConverter::massToSI(1.0)        << " kg"            << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K"             << endl;
    cout << "One unit of energy is "      << UnitConverter::energyToSI(1.0)      << " J"             << endl;

    bool BoltzmannDist = false;             // initial velocities given by Boltzmann distribution
    double maxMinVelocity = 0.5;           // uniformly distributed velocities [-v, v]
    bool periodicBoundaries = true;
    bool writeSampleToFile = false;

    System system;
    system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature, mass, BoltzmannDist, maxMinVelocity);
    //system.setPotential(new LennardJones(system, 1.0, 1.0));
    system.setPotential(new LennardJonesCellList(system, 1.0, 1.0, 2.5)); // You must insert correct parameters here
    //system.setIntegrator(new EulerCromer());
    system.setIntegrator(new VelocityVerlet());
    system.removeTotalMomentum();
    system.setPeriodicBoundaries(periodicBoundaries);

    StatisticsSampler statisticsSampler(system, writeSampleToFile);
    IO movie; // To write the state to file
    movie.open("animateVelocityDist.xyz");

    cout << "Timestep Time Temperature Pressure Density KineticEnergy PotentialEnergy TotalEnergy" << endl;

    // sample initial state
    movie.saveState(&system);

    clock_t start, finish;
    start = clock();
    for (int timestep=0; timestep < 501; timestep++) {
        system.step(dt);
        statisticsSampler.sample(timestep);
        //statisticsSampler.sampleRadialDistribution(50);
        if( !(timestep % 100) ) {
            // print the timestep every 100 timesteps
            cout << system.steps() << "      " << system.time() << "    " << statisticsSampler.temperature() << "    "
                 << statisticsSampler.pressure() << "  "  << statisticsSampler.density() << "    "
                 << statisticsSampler.kineticEnergy() << "     " << statisticsSampler.potentialEnergy() << "      "
                 << statisticsSampler.totalEnergy() << endl;
        }
        movie.saveState(&system);
    }
    finish = clock();
    cout << "Time elapsed: " << ((finish-start)/CLOCKS_PER_SEC) << endl;

    movie.close();

    return 0;
}
