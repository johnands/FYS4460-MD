#include "examples.h"
#include "math/random.h"
#include "potentials/lennardjones.h"
#include "potentials/lennardjonescelllist.h"
#include "potentials/neuralnetwork.h"
#include "potentials/tensorflownetwork.h"
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
#include <string>

using std::cout;
using std::endl;


int Examples::lennardJonesFCC() {

    int numberOfUnitCells = 12;
    double initialTemperature = 1.0;  // measured in Kelvin
    double latticeConstant    = UnitConverter::lengthFromAngstroms(5.26);
    double mass = UnitConverter::massFromSI(6.63352088e-26); // mass of Argon atom
    double dt = 0.01;  // Measured in seconds

    bool BoltzmannDist = true;             // initial velocities given by Boltzmann distribution
    double maxMinVelocity = 1.5;           // uniformly distributed velocities [-v, v]
    double tau = 10*dt;                    // relaxation time for Berendsen and Andersen

    bool usePores = false;

    System *system = new System();
    system->createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature,
                             mass, BoltzmannDist, maxMinVelocity);
    system->setPores(new CenteredCylinder(system, usePores));

    system->setPotential(new LennardJones(system, 1.0, 1.0));
    system->setTimeStep(dt);
    system->setIntegrator(new VelocityVerlet(system));
    system->setThermostat(new Berendsen(system, initialTemperature, tau));
    system->setUseThermoStat(false);
    system->setThermalization(0);

    system->setNumberOfTimeSteps(501);
    system->setTemperature(initialTemperature);
    system->removeTotalMomentum();

    system->setUseExternalForce(false);
    system->setWriteSample(false);
    system->setRadialDistribution(false);
    system->setMakeXYZ(false);
    system->setXYZName("LJFCC.xyz");

    system->runSimulation();
}


int Examples::lennardJonesFCCCellList() {

    int numberOfUnitCells = 12;
    double initialTemperature = 1.0;  // measured in Kelvin
    double latticeConstant    = UnitConverter::lengthFromAngstroms(5.26);
    double mass = UnitConverter::massFromSI(6.63352088e-26); // mass of Argon atom
    double dt = 0.01;  // Measured in seconds

    bool BoltzmannDist = true;             // initial velocities given by Boltzmann distribution
    double maxMinVelocity = 1.5;           // uniformly distributed velocities [-v, v]
    double tau = 10*dt;                    // relaxation time for Berendsen and Andersen

    bool usePores = false;

    System *system = new System();
    system->createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature,
                             mass, BoltzmannDist, maxMinVelocity);
    system->setPores(new CenteredCylinder(system, usePores));

    system->setPotential(new LennardJonesCellList(system, 1.0, 1.0, 2.5, 3.0));
    system->setTimeStep(dt);
    system->setIntegrator(new VelocityVerlet(system));
    system->setThermostat(new Berendsen(system, initialTemperature, tau));
    system->setUseThermoStat(false);
    system->setThermalization(0);

    system->setNumberOfTimeSteps(501);
    system->setTemperature(initialTemperature);
    system->removeTotalMomentum();

    system->setUseExternalForce(false);
    system->setWriteSample(false);
    system->setRadialDistribution(false);
    system->setMakeXYZ(false);
    system->setXYZName("LJNeighbourFCC.xyz");

    system->runSimulation();
}


int Examples::lennardJonesFCCNeuralNetwork() {

    int numberOfUnitCells = 12;
    double initialTemperature = 1.0;  // measured in Kelvin
    double latticeConstant    = UnitConverter::lengthFromAngstroms(5.26);
    double mass = UnitConverter::massFromSI(6.63352088e-26); // mass of Argon atom
    double dt = 0.01;  // Measured in seconds

    bool BoltzmannDist = true;             // initial velocities given by Boltzmann distribution
    double maxMinVelocity = 1.5;           // uniformly distributed velocities [-v, v]
    double tau = 10*dt;                    // relaxation time for Berendsen and Andersen

    bool usePores = false;

    System *system = new System();
    system->createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature,
                             mass, BoltzmannDist, maxMinVelocity);
    system->setPores(new CenteredCylinder(system, usePores));

    system->setPotential(new NeuralNetwork(system, "../TensorFlow/TrainingData/28.11-17.40.35/graph.dat", 2.5, 3.0));
    system->setTimeStep(dt);
    system->setIntegrator(new VelocityVerlet(system));
    system->setThermostat(new Berendsen(system, initialTemperature, tau));
    system->setUseThermoStat(false);
    system->setThermalization(0);

    system->setNumberOfTimeSteps(501);
    system->setTemperature(initialTemperature);
    system->removeTotalMomentum();

    system->setUseExternalForce(false);
    system->setWriteSample(false);
    system->setRadialDistribution(false);
    system->setMakeXYZ(false);
    system->setXYZName("LJNeighbourFCC.xyz");

    system->runSimulation();

}


int Examples::lennardJonesFCCTensorFlow() {


    int numberOfUnitCells = 12;
    double initialTemperature = 1.0;  // measured in Kelvin
    double latticeConstant    = UnitConverter::lengthFromAngstroms(5.26);
    double mass = UnitConverter::massFromSI(6.63352088e-26); // mass of Argon atom
    double dt = 0.01;  // Measured in seconds

    bool BoltzmannDist = true;             // initial velocities given by Boltzmann distribution
    double maxMinVelocity = 1.5;           // uniformly distributed velocities [-v, v]
    double tau = 10*dt;                    // relaxation time for Berendsen and Andersen
    string graphFileName ("../TensorFlow/TrainingData/18.11-16.12.57/frozen_graph.pb");

    bool usePores = false;

    System *system = new System();
    system->createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature,
                             mass, BoltzmannDist, maxMinVelocity);
    system->setPores(new CenteredCylinder(system, usePores));

    system->setPotential(new TensorFlowNetwork(system, graphFileName, 2.5, 3.0));
    system->setTimeStep(dt);
    system->setIntegrator(new VelocityVerlet(system));
    system->setThermostat(new Berendsen(system, initialTemperature, tau));
    system->setUseThermoStat(false);
    system->setThermalization(0);

    system->setNumberOfTimeSteps(501);
    system->setTemperature(initialTemperature);
    system->removeTotalMomentum();

    system->setUseExternalForce(false);
    system->setWriteSample(false);
    system->setRadialDistribution(false);
    system->setMakeXYZ(false);
    system->setXYZName("LJNeighbourFCC.xyz");

    return system->runSimulation();
}


int Examples::lennardJonesLiquid() {

    int numberOfUnitCells = 12;
    double initialTemperature = 0.851;  // measured in Kelvin
    double latticeConstant    = UnitConverter::lengthFromAngstroms(5.72);
    double mass = UnitConverter::massFromSI(6.63352088e-26); // mass of Argon atom
    double dt = 0.01;  // Measured in seconds

    bool BoltzmannDist = true;             // initial velocities given by Boltzmann distribution
    double maxMinVelocity = 1.5;           // uniformly distributed velocities [-v, v]
    double tau = 10*dt;                    // relaxation time for Berendsen and Andersen

    bool usePores = false;

    System *system = new System();
    system->createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature,
                             mass, BoltzmannDist, maxMinVelocity);
    system->setPores(new CenteredCylinder(system, usePores));

    system->setPotential(new LennardJonesCellList(system, 1.0, 1.0, 2.5, 3.0));
    system->setTimeStep(dt);
    system->setIntegrator(new VelocityVerlet(system));
    system->setThermostat(new Berendsen(system, initialTemperature, tau));
    system->setUseThermoStat(false);
    system->setThermalization(0);

    system->setNumberOfTimeSteps(501);
    system->setTemperature(initialTemperature);
    system->removeTotalMomentum();

    system->setUseExternalForce(false);
    system->setWriteSample(false);
    system->setRadialDistribution(false);
    system->setMakeXYZ(false);
    system->setXYZName("LJNeighbourFCC.xyz");

    return system->runSimulation();
}

int Examples::loadFromFile() {

    int numberOfUnitCells = 12;
    double initialTemperature = 0.851;  // measured in Kelvin
    double latticeConstant    = UnitConverter::lengthFromAngstroms(5.72);
    double mass = UnitConverter::massFromSI(6.63352088e-26); // mass of Argon atom
    double dt = 0.01;  // Measured in seconds

    bool BoltzmannDist = true;             // initial velocities given by Boltzmann distribution
    double maxMinVelocity = 1.5;           // uniformly distributed velocities [-v, v]
    double tau = 10*dt;                    // relaxation time for Berendsen and Andersen

    bool usePores = false;

    System *system = new System();
    system->readFromStateFile("thermalizedFluidNc20T15Nt1001.xyz", mass, latticeConstant,
                              numberOfUnitCells);
    system->setPores(new CenteredCylinder(system, usePores));

    system->setPotential(new LennardJonesCellList(system, 1.0, 1.0, 2.5, 3.0));
    system->setTimeStep(dt);
    system->setIntegrator(new VelocityVerlet(system));
    system->setThermostat(new Berendsen(system, initialTemperature, tau));
    system->setUseThermoStat(false);
    system->setThermalization(0);

    system->setNumberOfTimeSteps(501);
    system->setTemperature(initialTemperature);
    system->removeTotalMomentum();

    system->setUseExternalForce(false);
    system->setWriteSample(false);
    system->setRadialDistribution(false);
    system->setMakeXYZ(false);
    system->setXYZName("LJNeighbourFCC.xyz");

    return system->runSimulation();
}


int Examples::lennardJonesFCCNanoCylinder() {

    int numberOfUnitCells = 12;
    double initialTemperature = 0.851;  // measured in Kelvin
    double latticeConstant    = UnitConverter::lengthFromAngstroms(5.72);
    double mass = UnitConverter::massFromSI(6.63352088e-26); // mass of Argon atom
    double dt = 0.01;  // Measured in seconds

    bool BoltzmannDist = true;             // initial velocities given by Boltzmann distribution
    double maxMinVelocity = 1.5;           // uniformly distributed velocities [-v, v]
    double tau = 10*dt;                    // relaxation time for Berendsen and Andersen

    bool usePores = true;

    System *system = new System();
    system->createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature,
                             mass, BoltzmannDist, maxMinVelocity);
    system->setPores(new CenteredCylinder(system, usePores));

    system->setPotential(new LennardJonesCellList(system, 1.0, 1.0, 2.5, 3.0));
    system->setTimeStep(dt);
    system->setIntegrator(new VelocityVerlet(system));
    system->setThermostat(new Berendsen(system, initialTemperature, tau));
    system->setUseThermoStat(false);
    system->setThermalization(0);

    system->setNumberOfTimeSteps(501);
    system->setTemperature(initialTemperature);
    system->removeTotalMomentum();

    system->setUseExternalForce(false);
    system->setWriteSample(false);
    system->setRadialDistribution(false);
    system->setMakeXYZ(false);
    system->setXYZName("LJNeighbourFCC.xyz");

    return system->runSimulation();
}


int Examples::lennardJonesFCCNanoSpheres() {

    int numberOfUnitCells = 12;
    double initialTemperature = 0.851;  // measured in Kelvin
    double latticeConstant    = UnitConverter::lengthFromAngstroms(5.72);
    double mass = UnitConverter::massFromSI(6.63352088e-26); // mass of Argon atom
    double dt = 0.01;  // Measured in seconds

    bool BoltzmannDist = true;             // initial velocities given by Boltzmann distribution
    double maxMinVelocity = 1.5;           // uniformly distributed velocities [-v, v]
    double tau = 10*dt;                    // relaxation time for Berendsen and Andersen

    bool usePores = true;

    System *system = new System();
    system->createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature,
                             mass, BoltzmannDist, maxMinVelocity);
    system->setPores(new Spheres(system, usePores));

    system->setPotential(new LennardJonesCellList(system, 1.0, 1.0, 2.5, 3.0));
    system->setTimeStep(dt);
    system->setIntegrator(new VelocityVerlet(system));
    system->setThermostat(new Berendsen(system, initialTemperature, tau));
    system->setUseThermoStat(false);
    system->setThermalization(0);

    system->setNumberOfTimeSteps(501);
    system->setTemperature(initialTemperature);
    system->removeTotalMomentum();

    system->setUseExternalForce(false);
    system->setWriteSample(false);
    system->setRadialDistribution(false);
    system->setMakeXYZ(false);
    system->setXYZName("LJNeighbourFCC.xyz");

    return system->runSimulation();
}


int Examples::computeTemperatureFluctuations() {

}


int Examples::computeRadialDistributionFunction() {

}


int Examples::compareNeuralNetworkError() {

    // check if error of NN has same shape as in python
    int numberOfPoints = 500;
    arma::vec distances = arma::linspace<arma::vec>(0.8, 2.5, numberOfPoints);
    NeuralNetwork *networkPotential = new NeuralNetwork(system, "../TensorFlow/TrainingData/28.11-17.40.35/graph.dat", 2.5, 3.0);
    for (int i=0; i < numberOfPoints; i++) {
        double energy = networkPotential->network(distances(i));

    }
}
