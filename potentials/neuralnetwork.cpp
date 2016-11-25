#include "neuralnetwork.h"
#include "celllist.h"
#include "../math/activationfunctions.h"

NeuralNetwork::NeuralNetwork(System *system, const char *filename,
                             double rCut, double neighbourCut) :
               Potential(system) {

    m_filename = filename;
    readFromFile();

    m_cellList = new CellList(system, rCut, neighbourCut);
    m_cellList->setupCells();
    m_cellList->updateCells();

    m_rCut = rCut;
    m_rCutSquared = rCut*rCut;
    m_neighbourCut = neighbourCut;

    double r2 = 1.0 / m_rCutSquared;
    m_potentialCut = 4*r2*r2*r2*(r2*r2*r2 - 1);

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

    // set sizes
    m_preActivations.resize(m_nLayers+2);
    m_activations.resize(m_nLayers+2);
    m_derivatives.resize(m_nLayers+2);

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
    /*for (const auto i : weightsTemp)
        std::cout << i << std::endl;
    std::cout << std::endl;
    for (const auto i : m_biases)
        std::cout << i << std::endl;*/

    // resize weights and biases matrices to correct shapes
    m_weights.resize(m_nLayers+1);

    // first hidden layer
    m_weights[0]  = weightsTemp[0];

    // following hidden layers
    for (int i=0; i < m_nLayers-1; i++) {
        m_weights[i+1] = weightsTemp[i*m_nNodes+1];
        for (int j=1; j < m_nNodes; j++) {
            m_weights[i+1] = arma::join_cols(m_weights[i+1], weightsTemp[i*nodes+1+j]);
        }
    }

    // output layer
    arma::mat outputLayer = weightsTemp.back();
    m_weights[m_nLayers] = arma::reshape(outputLayer, m_nNodes, 1);

    // reshape bias of output node
    m_biases[m_nLayers].shed_cols(1,m_nNodes-1);

    m_weightsTransposed.resize(m_nLayers+1);
    // obtained transposed matrices
    for (int i=0; i < m_weights.size(); i++)
        m_weightsTransposed[i] = m_weights[i].t();

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

    // linear activation for input layer
    m_preActivations[0] = data;
    m_activations[0] = m_preActivations[0];

    // hidden layers
    for (int i=0; i < m_nLayers; i++) {
        // weights and biases starts at first hidden layer:
        // weights[0] are the weights connecting input layer to first hidden layer
        m_preActivations[i+1] = m_activations[i]*m_weights[i] + m_biases[i];
        m_activations[i+1] = ActivationFunctions::sigmoid(m_preActivations[i+1]);
    }

    // linear activation for output layer
    m_preActivations[m_nLayers+1] = m_activations[m_nLayers]*m_weights[m_nLayers] + m_biases[m_nLayers];
    m_activations[m_nLayers+1] = m_preActivations[m_nLayers+1];

    // return activation of output neuron, which is a 1x1-matrix
    return m_activations[m_nLayers+1](0,0);
}


double NeuralNetwork::backPropagation() {
    // find derivate of output w.r.t. intput, i.e. dE/dr_ij
    // need to find the "error" terms for all the nodes in all the layers

    // the derivative of the output neuron's activation function w.r.t.
    // its input is propagated backwards.
    // the output activation function is f(x) = x, so this is 1
    arma::mat output(1,1); output.fill(1);
    m_derivatives[m_nLayers+1] = output;

    // we can thus compute the error vectors for the other layers
    for (int i=m_nLayers; i > 0; i--) {
        m_derivatives[i] = ( m_derivatives[i+1]*m_weightsTransposed[i] ) %
                           ActivationFunctions::sigmoidDerivative(m_preActivations[i]);
    }

    // linear activation function for input neuron
    m_derivatives[0] = m_derivatives[1]*m_weightsTransposed[0];

    return m_derivatives[0](0,0);
}


void NeuralNetwork::calculateForces() {

    double potentialEnergy = 0;
    double pressure = 0;

    // update lists every 5th time step
    if (m_updateLists >= 20) {
        m_cellList->clearCells();
        m_cellList->updateCells();
        m_updateLists = 0;
    }

    m_updateLists++;

    for (int i=0; i < m_system->atoms().size(); i++) {
        Atom *atom1 = m_system->atoms()[i];

        // loop over all atoms in atom1's neighbour list
        for (int j=0; j < m_cellList->getNeighbours()[i].size(); j++) {

            Atom *atom2 = m_cellList->getNeighbours()[i][j];

            // calculate distance vector
            vec3 dr = atom1->position - atom2->position;

            // make sure we're using shortest distance component-wise (periodic boundary conditions)
            for (int dim=0; dim < 3; dim++) {
                if      (dr[dim] > m_system->systemSizeHalf()[dim])  { dr[dim] -= m_system->systemSize()[dim]; }
                else if (dr[dim] < -m_system->systemSizeHalf()[dim]) { dr[dim] += m_system->systemSize()[dim]; }
            }

            double dr2 = dr.lengthSquared();

            // calculate potential energy with ANN using distance as input if
            // distance is shorter than cutoff
            vec3 forceOnAtom;
            if (dr2 <= m_rCutSquared) {
                double distance = dr.length();
                double energy = network(distance);
                potentialEnergy += energy;

                // calculate force: backpropagation...
                double dEdr = backPropagation();
                double drInverse = 1.0 / distance;
                forceOnAtom = -dEdr*drInverse*dr;
            }
            else {
                // if not, zero force and energy
                forceOnAtom = dr*0;
            }

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


void NeuralNetwork::clearVectors() {

    for (int i=0; i < m_derivatives.size(); i++) {
        m_derivatives[i].zeros();
        m_preActivations[i].zeros();
        m_activations[i].zeros();
    }
}







