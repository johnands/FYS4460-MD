#include "math/random.h"
#include "examples.h"
#include "potentials/lennardjones.h"
#include "potentials/lennardjonescelllist.h"
#include "potentials/neuralnetwork.h"
#include "potentials/manyneighbournn.h"
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
#include <armadillo>
#include <vector>
#include <iomanip>
#include <math.h>
#include <utility>
#include <tuple>
#include <map>

using std::cout;
using std::endl;

int main() {

    cout << "One unit of length is "      << UnitConverter::lengthToSI(1.0)      << " meters"        << endl;
    cout << "One unit of velocity is "    << UnitConverter::velocityToSI(1.0)    << " meters/second" << endl;
    cout << "One unit of time is "        << UnitConverter::timeToSI(1.0)        << " seconds"       << endl;
    cout << "One unit of mass is "        << UnitConverter::massToSI(1.0)        << " kg"            << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K"             << endl;
    cout << "One unit of energy is "      << UnitConverter::energyToSI(1.0)      << " J"             << endl;

    return Examples::lennardJonesFCC();
    //return Examples::lennardJonesFCCCellList();
    return Examples::lennardJonesFCCNeuralNetwork();
    //return Examples::lennardJonesFCCManyNeighbourNeuralNetwork();
    //return Examples::lennardJonesFCCTensorFlow();
    //return Examples::lennardJonesLiquid();
    //return Examples::loadFromFile();
    //return Examples::lennardJonesFCCNanoCylinder();
    //return Examples::lennardJonesFCCNanoSpheres();
    //return Examples::computeTemperatureFluctuations();
    //return Examples::computeRadialDistributionFunction();
    //return Examples::compareNeuralNetworkError();
    //return Examples::compareManyNeighbourNeuralNetworkError();
    //return Examples::testBackpropagation();
}
