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
    //return Examples::lennardJonesFCCNeuralNetwork();
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

    /*std::vector<std::vector<double>> yes(5);
    yes[0].push_back(0);
    cout << yes[0][0] << endl;
    exit(1);*/


    /*arma::mat inputVector(2,3);
    inputVector(0,0) = 2.0;
    inputVector(0,1) = 0.5;
    inputVector(0,2) = 1;
    inputVector.row(1) = inputVector.row(0);
    cout << inputVector << endl;
    double sumYes = 5 + arma::accu(inputVector);
    cout << sumYes << endl;*/

    /*arma::mat no(1,0);
    cout << no*4 << endl;
    cout << arma::accu(4 + no) << endl;
    exit(1);

    arma::mat vector3(3,3, arma::fill::zeros);
    vector3(0,0) = 2;
    vector3(0,1) = 3;
    vector3(0,2) = 2;
    arma::mat inputVectorz = inputVector.head_cols(2);
    cout << "vector3" << arma::size(vector3.row(1)) << endl;

    for (int i=0; i < 5; i++) {
        inputVector += inputVector;
        cout << inputVector << endl;
    }
    cout << inputVectorz << endl;

    arma::rowvec vector2(3);
    vector2(0) = 0;
    vector2(1) = 0.5;
    vector2(2) = 1;
    cout << arma::accu(vector2) << endl;

    arma::rowvec vector1(3);
    vector1(0) = 2;
    vector1(1) = 3;
    vector1(2) = 2;
    cout << vector2 % vector3 << endl;

    arma::arma_version ver;
    std::cout << "ARMA version: "<< ver.as_string() << std::endl;


    arma::mat A(1, 500, arma::fill::randu);
    arma::mat B(1, 1, arma::fill::randu);

    clock_t start, finish;
    start = clock();
    for (int i=0; i < int(1e6); i++) {
        arma::mat product = A % B;
    }
    finish = clock();
    cout << "Time elapsed: " << ((finish-start)/CLOCKS_PER_SEC) << endl;


    arma::mat C(1, 500, arma::fill::randu);
    double D = Random::nextDouble();

    start = clock();
    for (int i=0; i < int(1e6); i++) {
        arma::mat product = C*D;
    }
    finish = clock();
    cout << "Time elapsed: " << ((finish-start)/CLOCKS_PER_SEC) << endl;*/


    /*start = clock();
    for (int i=0; i < int(1e7); i++) {
        arma::mat product(1,500);
        for (int j=0; j < 500; j++) {
            product(0,j) = A(0,j) * B(0,0);
        }
    }
    finish = clock();
    cout << "Time elapsed: " << ((finish-start)/CLOCKS_PER_SEC) << endl;*/

    /*start = clock();
    for (int i=0; i < int(1e8); i++) {
        arma::mat product(1,500);
        for (int j=0; j < 500; j++) {
            product(0,j) = A(0,j) * B(0,j);
        }
    }
    finish = clock();
    cout << "Time elapsed: " << ((finish-start)/CLOCKS_PER_SEC) << endl;*/

    /*arma::rowvec C(500, arma::fill::randu);
    arma::rowvec D(500, arma::fill::randu);

    start = clock();
    for (int i=0; i < 100000000; i++) {
        arma::rowvec product = C % D;
    }
    finish = clock();
    cout << "Time elapsed: " << ((finish-start)/CLOCKS_PER_SEC) << endl;

    arma::mat A(1, 50000, arma::fill::randu);
    arma::mat B(1, 500, arma::fill::randu);

    clock_t start, finish;
    start = clock();
    double sum1;
    for (int i=0; i < 100000000; i++) {
        sum1 = arma::accu(A);
    }
    finish = clock();
    cout << "Time elapsed: " << ((finish-start)/CLOCKS_PER_SEC) << endl;

    arma::rowvec C(50000, arma::fill::randu);
    arma::rowvec D(500, arma::fill::randu);

    start = clock();
    double sum2;
    for (int i=0; i < 100000000; i++) {
        sum2 = arma::accu(C);
    }
    finish = clock();
    cout << "Time elapsed: " << ((finish-start)/CLOCKS_PER_SEC) << endl;*/

}
