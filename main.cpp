#include "math/random.h"
#include "potentials/lennardjones.h"
#include "potentials/lennardjonescelllist.h"
#include "potentials/neuralnetwork.h"
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
#include <armadillo>
#include <fstream>

using namespace std;

int main(int numberOfArguments, char **argumentList)
{
    int numberOfUnitCells = 10;
    //double initialTemperature = UnitConverter::temperatureFromSI();  // measured in Kelvin
    double initialTemperature = 1.5;
    //double latticeConstant    = UnitConverter::lengthFromAngstroms(5.72); // measured in angstroms
    double latticeConstant    = UnitConverter::lengthFromAngstroms(5.26);
    cout << "lattice: " << latticeConstant << endl;

    double mass = UnitConverter::massFromSI(6.63352088e-26); // mass of Argon atom

    double dt = 0.01; //UnitConverter::timeFromSI(5e-14); // Measured in seconds

    cout << "One unit of length is "      << UnitConverter::lengthToSI(1.0)      << " meters"        << endl;
    cout << "One unit of velocity is "    << UnitConverter::velocityToSI(1.0)    << " meters/second" << endl;
    cout << "One unit of time is "        << UnitConverter::timeToSI(1.0)        << " seconds"       << endl;
    cout << "One unit of mass is "        << UnitConverter::massToSI(1.0)        << " kg"            << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K"             << endl;
    cout << "One unit of energy is "      << UnitConverter::energyToSI(1.0)      << " J"             << endl;

    bool BoltzmannDist = true;             // initial velocities given by Boltzmann distribution
    double maxMinVelocity = 1.5;           // uniformly distributed velocities [-v, v]
    double tau = 10*dt;

    bool readFromFile = false;
    bool usePores = false;

    System *system = new System();
    if (readFromFile) {
        system->readFromStateFile("thermalizedFluidNc20T15Nt1001.xyz", mass, latticeConstant,
                                      numberOfUnitCells);
    }
    else {
        system->createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature,
                                mass, BoltzmannDist, maxMinVelocity);
    }

    system->setPores(new CenteredCylinder(system, usePores));

    //system->setPotential(new LennardJones(system, 1.0, 1.0));
    system->setPotential(new LennardJonesCellList(system, 1.0, 1.0, 2.5, 3.0));
    //system->setPotential(new NeuralNetwork(system, "../TensorFlow/TrainingData/04.11-14.58.07/graph.dat", 2.5, 3.0));
    //system->setIntegrator(new EulerCromer(system));
    system->setTimeStep(dt);
    system->setIntegrator(new VelocityVerlet(system));
    system->setThermostat(new Berendsen(system, initialTemperature, tau));
    system->setUseThermoStat(false);
    system->setThermalization(0);

    system->setNumberOfTimeSteps(501);
    system->setTemperature(initialTemperature);
    system->removeTotalMomentum();

    //system->setPeriodicBoundaries(true);
    system->setUseExternalForce(false);
    system->setWriteSample(false);
    system->setRadialDistribution(false);
    system->setMakeXYZ(false);
    system->setXYZName("testCell.xyz");

    system->runSimulation();

    // test if network represent LJ well
    /*NeuralNetwork *networkz = new NeuralNetwork(system, "../TensorFlow/TrainingData/04.11-14.58.07/graph.dat", 2.5, 3.0);
    ofstream outFile;
    outFile.open("../TensorFlow/networkTest.dat", ios::out);
    arma::vec vector = arma::linspace<arma::vec>(0.8, 2.5, 1000);
    for (int i=0; i < arma::size(vector)(0); i++) {
        double energy = networkz->network(vector(i));
        double derivative = networkz->backPropagation();
        outFile << energy << " " << derivative << endl;
    }*/

    return 0;
}
