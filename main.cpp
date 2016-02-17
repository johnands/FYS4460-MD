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

using namespace std;

int main(int numberOfArguments, char **argumentList)
{

    /*double a = 7.9;
    int b = 12;
    int c = 24;
    int d = a/b*c;
    cout << d << endl;*/


    int numberOfUnitCells = 2;
    double initialTemperature = UnitConverter::temperatureFromSI(100.0);  // measured in Kelvin
    double latticeConstant    = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms

    // If a first argument is provided, it is the number of unit cells
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);

    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));

    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3) latticeConstant= UnitConverter::lengthFromAngstroms(atof(argumentList[3]));

    double mass = UnitConverter::massFromSI(6.63352088e-26); // mass of Argon atom

    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds

    cout << "One unit of length is "      << UnitConverter::lengthToSI(1.0)      << " meters"        << endl;
    cout << "One unit of velocity is "    << UnitConverter::velocityToSI(1.0)    << " meters/second" << endl;
    cout << "One unit of time is "        << UnitConverter::timeToSI(1.0)        << " seconds"       << endl;
    cout << "One unit of mass is "        << UnitConverter::massToSI(1.0)        << " kg"            << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K"             << endl;
    cout << "One unit of energy is "      << UnitConverter::energyToSI(1.0)      << " J"             << endl;

    System system;
    system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature, mass);
    system.setPotential(new LennardJonesCellList(1.0, 1.0, 3.0)); // You must insert correct parameters here
    system.setIntegrator(new EulerCromer());
    system.removeTotalMomentum();

    StatisticsSampler statisticsSampler;
    IO movie; // To write the state to file
    movie.open("movie.xyz");

    cout << "Timestep Time Temperature KineticEnergy PotentialEnergy TotalEnergy" << endl;
    for (int timestep=0; timestep<1000; timestep++) {
        system.step(dt);
        statisticsSampler.sample(system);
        if( !(timestep % 100) ) {
            // Print the timestep every 100 timesteps
            cout << system.steps() << "      " << system.time() << "      " << statisticsSampler.temperature() << "      " << statisticsSampler.kineticEnergy() << "      " << statisticsSampler.potentialEnergy() << "      " << statisticsSampler.totalEnergy() << endl;
        }
        movie.saveState(&system);
    }

    movie.close();

    return 0;
}
