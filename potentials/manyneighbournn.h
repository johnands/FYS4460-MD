#ifndef MANYNEIGHBOURNN_H
#define MANYNEIGHBOURNN_H
#include "neuralnetwork.h"

class ManyNeighbourNN : public NeuralNetwork {

public:
    ManyNeighbourNN(System *system, const char *filename,
                    double rCut, double neighbourCut);
    arma::mat network(arma::mat inputVector);
    arma::mat backPropagation();
    void calculateForces();
};

#endif // MANYNEIGHBOURNN_H
