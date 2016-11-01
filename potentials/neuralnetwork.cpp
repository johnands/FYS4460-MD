#include "neuralnetwork.h"
#include "celllist.h"
#include "../math/activationfunctions.h"

NeuralNetwork::NeuralNetwork(System *system, const char *filename, double cutOffDistance) : Potential(system) {

    m_filename = filename;
    readFromFile();

    m_cellList = new CellList(system, cutOffDistance);
    m_cellList->setupCells();
    m_cellList->updateCells();
    m_updateLists = 0;
}

void NeuralNetwork::readFromFile() {

    std::ifstream input;
    input.open(m_filename, std::ios::in);

    // check if file successfully opened
    if ( !input.is_open() ) std::cout << "File is not opened" << std::endl;

    // process first line
    std::string activation;
    input >> m_nLayers >> m_nNodes >> activation;
    std::cout << m_nLayers << " " << m_nNodes << " " << activation << std::endl;

    // skip a blank line
    std::string dummyLine;
    std::getline(input, dummyLine);

    // process file
    // store all weights in a temporary vector
    // that will be reshaped later
    std::vector<arma::mat> weightsTemp;
    for ( std::string line; std::getline(input, line); ) {
        //std::cout << line << std::endl;

        if ( line.empty() )
        { std::cout << "yes" << std::endl; break;}


        // store all weights in a vector
        double buffer;                  // have a buffer string
        std::stringstream ss(line);     // insert the string into a stream

        // while there are new weights on current line, add them to vector
        arma::mat matrix(1,m_nNodes);
        int i = 0;
        while ( ss >> buffer ) {
            matrix(0,i) = buffer;
            i++;
        }
        weightsTemp.push_back(matrix);
    }

    // can put all biases in vector directly
    // no need for temporary vector
    for ( std::string line; std::getline(input, line); ) {

        // store all weights in vector
        double buffer;                  // have a buffer string
        std::stringstream ss(line);     // insert the string into a stream

        // while there are new weights on current line, add them to vector
        arma::mat matrix(1,m_nNodes);
        int i = 0;
        while ( ss >> buffer ) {
            matrix(0,i) = buffer;
            i++;
        }
        m_biases.push_back(matrix);
    }

    // write out all weights and biases
    for (const auto i : weightsTemp)
        std::cout << i << std::endl;
    std::cout << std::endl;
    for (const auto i : m_biases)
        std::cout << i << std::endl;

    // resize weights and biases matrices to correct shapes
    m_weights.resize(m_nLayers+1);

    // first hidden layer
    m_weights[0]  = weightsTemp[0];

    // following hidden layers
    for (int i=0; i < m_nLayers-1; i++) {
        m_weights[i+1] = weightsTemp[i*m_nNodes+1];
        for (int j=1; j < m_nNodes; j++) {
            m_weights[i+1] = arma::join_cols(m_weights[i+1], weightsTemp[i+1+j]);
        }
    }

    // output layer
    arma::mat outputLayer = weightsTemp.back();
    m_weights[m_nLayers] = arma::reshape(outputLayer, m_nNodes, 1);

    // reshape bias of output node
    m_biases[m_nLayers].shed_cols(1,m_nNodes-1);

    // write out entire system for comparison
    for (const auto i : m_weights)
        std::cout << i << std::endl;

    for (const auto i : m_biases)
        std::cout << i << std::endl;
}

double NeuralNetwork::network(double dataPoint) {
    // the data needs to be a 1x1 armadillo matrix
    // maybe more than one data point can be processed simultaneously?

    // convert data point to 1x1 matrix
    arma::mat data(1,1); data(0,0) = dataPoint;

    // send data through network
    // use relu as activation except for output layer
    m_preActivations.resize(m_nLayers+2);

    m_activations.resize(m_nLayers+1);
    m_preActivations[0] = data*m_weights[0] + m_biases[0];
    m_activations[0] = ActivationFunctions::sigmoid(m_preActivations[0]);
    for (int i=1; i < m_nLayers; i++) {
        // relu for all hidden layers except last one
        if (i < m_nLayers-1) {
            m_preActivations[i] = m_activations[i-1]*m_weights[i] + m_biases[i];
            m_activations[i] = ActivationFunctions::sigmoid(m_preActivations[i]);
        }

        // sigmoid on last hidden layer
        else {
            m_preActivations[i] = m_activations[i-1]*m_weights[i] + m_biases[i];
            m_activations[i] = ActivationFunctions::sigmoid(m_preActivations[i]);
        }
    }
    // no activation function for output layer (i.e. linear or identity activation function)
    m_preActivations[m_nLayers] = m_activations[m_nLayers-1]*m_weights[m_nLayers] + m_biases[m_nLayers];
    m_activations[m_nLayers] = m_preActivations[m_nLayers];

    //std::cout << activations[m_nLayers](0,0) << std::endl;
    return m_activations[m_nLayers](0,0);
}


double NeuralNetwork::backPropagation() {
    // find derivate of output w.r.t. intput, i.e. dE/dr_ij
    // need to find the "error" terms for all the nodes in all the layers

    // the derivative of the output neuron's activation function w.r.t.
    // its input is propagated backwards.
    // the output activation function is f(x) = x, so this is 1
    m_derivatives.resize(m_nLayers+1);
    arma::mat output(1,1); output.fill(1);
    m_derivatives[m_nLayers] = output;

    // we can thus compute the error vectors for the other layers
    for (int i=m_nLayers-1; i >= 0; i--) {
        m_derivatives[i] = ( m_derivatives[i+1]*m_weights[i+1].t() ) %
                           ActivationFunctions::sigmoidDerivative(m_preActivations[i]);
        std::cout << arma::size(m_derivatives[i]) << std::endl;
    }

    std::cout << arma::size(m_derivatives[0]) << std::endl;
    return m_derivatives[0](0,0);
}


void NeuralNetwork::calculateForces() {

    double potentialEnergy = 0;
    double pressure = 0;

    // update lists every 5th time step
    if (m_updateLists==5) {
        m_cellList->clearCells();
        m_cellList->updateCells();
        m_updateLists = 0;
    }

    m_updateLists++;

    for (int i=0; i < m_system->atoms().size(); i++) {
        Atom *atom1 = m_system->atoms()[i];
        vec3 dr, forceOnAtom;

        // loop over all atoms in atom1's neighbour list
        for (int j=0; j < m_cellList->getNeighbours()[i].size(); j++) {

            Atom *atom2 = m_cellList->getNeighbours()[i][j];

            // calculate distance vector
            dr = atom1->position - atom2->position;

            // make sure we're using shortest distance component-wise (periodic boundary conditions)
            for (int dim=0; dim < 3; dim++) {
                if      (dr[dim] > m_system->systemSizeHalf()[dim])  { dr[dim] -= m_system->systemSize()[dim]; }
                else if (dr[dim] < -m_system->systemSizeHalf()[dim]) { dr[dim] += m_system->systemSize()[dim]; }
            }

            // calculate potential energy with ANN using distance as input
            double distance = dr.length();
            double energy = network(distance);
            potentialEnergy += energy;

            // calculate force: backpropagation...
            double dEdr = backPropagation();
            double drInverse = 1.0/distance;
            forceOnAtom = -dEdr*drInverse*dr;
            //std::cout << forceOnAtom << std::endl;

            // add contribution to force on atom i and j
            atom1->force += forceOnAtom;
            atom2->force -= forceOnAtom;   // Newton's third law

            if (m_system->getUseExternalForce()) {
                atom1->force[0] += 1.0;
            }

            // dot product of Fij and dr
            pressure += forceOnAtom.dot(dr);
        }
    }

    m_potentialEnergy = potentialEnergy;
    m_pressure = m_inverseVolume*pressure;
}









