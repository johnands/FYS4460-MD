#include "manyneighbournn.h"
#include "celllist.h"
#include "../math/activationfunctions.h"

using std::cout;
using std::endl;

ManyNeighbourNN::ManyNeighbourNN(System *system, const char *filename,
                                 double rCut, double neighbourCut,
                                 int numberOfNeighbours) :
                 NeuralNetwork(system, filename, rCut, neighbourCut, numberOfNeighbours) {

}


double ManyNeighbourNN::network(arma::mat inputVector) {
    // input vector is a 1xinputs vector

    cout << "yes" << endl;
    // linear activation for input layer
    m_preActivations[0] = inputVector;
    m_activations[0] = m_preActivations[0];
    cout << "shape: " << arma::size(m_activations[0]) << endl;

    // hidden layers
    for (int i=0; i < m_nLayers; i++) {
        // weights and biases starts at first hidden layer:
        // weights[0] are the weights connecting input layer to first hidden layer
        m_preActivations[i+1] = m_activations[i]*m_weights[i] + m_biases[i];
        m_activations[i+1] = ActivationFunctions::sigmoid(m_preActivations[i+1]);
        cout << "shape: " << arma::size(m_activations[i+1]) << endl;
    }

    // linear activation for output layer
    m_preActivations[m_nLayers+1] = m_activations[m_nLayers]*m_weights[m_nLayers] + m_biases[m_nLayers];
    m_activations[m_nLayers+1] = m_preActivations[m_nLayers+1];
    cout << "shape: " << arma::size(m_activations[m_nLayers+1]) << endl;

    // return activation of output neuron, which is a 1x1-matrix
    return m_activations[m_nLayers+1](0,0);
}


arma::mat ManyNeighbourNN::backPropagation() {
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

    return m_derivatives[0];
}


void ManyNeighbourNN::calculateForces() {

    double potentialEnergy = 0;
    double pressure = 0;

    // update lists every 20th time step
    if (m_updateLists >= 20) {
        m_cellList->clearCells();
        m_cellList->updateCells();
        m_updateLists = 0;
    }
    m_updateLists++;

    arma::mat distanceNeighbours(1, m_numberOfNeighbours);
    std::vector<std::vector<double> > drNeighbours(m_numberOfNeighbours, std::vector<double>(3));
    std::vector<Atom *> atomNeighbours(m_numberOfNeighbours);
    for (int i=0; i < m_system->atoms().size(); i++) {
        Atom *atom1 = m_system->atoms()[i];

        // loop over all atoms in atom1's neighbour list
        int neighbourCount = 0;
        for (int j=0; j < m_cellList->getNeighbours()[i].size(); j++) {
            //std::cout << m_cellList->getNeighbours()[i].size() << std::endl;

            Atom *atom2 = m_cellList->getNeighbours()[i][j];

            // calculate distance vector
            double dr[] = {atom1->position[0] - atom2->position[0],
                           atom1->position[1] - atom2->position[1],
                           atom1->position[2] - atom2->position[2]};

            // make sure we're using shortest distance component-wise (periodic boundary conditions)
            for (int dim=0; dim < 3; dim++) {
                if      (dr[dim] > m_system->systemSizeHalf()[dim])  { dr[dim] -= m_system->systemSize()[dim]; }
                else if (dr[dim] < -m_system->systemSizeHalf()[dim]) { dr[dim] += m_system->systemSize()[dim]; }
            }

            double dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

            // add the 20 first neighbours to input vector
            // regardless of cut etc.
            if (neighbourCount < m_numberOfNeighbours) {
                double distance = sqrt(dr2);
                distanceNeighbours(1, neighbourCount) = distance;
                drNeighbours[neighbourCount][0] = dr[0];
                drNeighbours[neighbourCount][1] = dr[1];
                drNeighbours[neighbourCount][2] = dr[2];
                atomNeighbours[neighbourCount] = atom2;
                neighbourCount++;
            }
            else {
                // if not, zero force and energy
                vec3 forceOnAtom;
                atom2->force += forceOnAtom;
            }
        }

        // after making input vector, compute energy contribution from all neighbours
        double energy = network(distanceNeighbours);
        potentialEnergy += energy;

        // find dEdr for all neighbouring atoms
        arma::mat dEdr = backPropagation();

        // loop through the atoms that has been sent through the network
        // and calculate forces
        for (int k=0; k < m_numberOfNeighbours; k++) {

            vec3 forceOnAtom;
            double drInverse = 1.0 / distanceNeighbours(1,k);
            forceOnAtom[0] = -dEdr(1,k)*drInverse*drNeighbours[k][0];
            forceOnAtom[1] = -dEdr(1,k)*drInverse*drNeighbours[k][1];
            forceOnAtom[2] = -dEdr(1,k)*drInverse*drNeighbours[k][2];

            atom1->force += forceOnAtom;
            atomNeighbours[k]->force -= forceOnAtom;   // Newton's third law

            // dot product of Fij and dr
            pressure += forceOnAtom[0]*drNeighbours[k][0] +
                        forceOnAtom[1]*drNeighbours[k][1] +
                        forceOnAtom[2]*drNeighbours[k][2];
        }
    }

    m_potentialEnergy = potentialEnergy;
    m_pressure = m_inverseVolume*pressure;
}
