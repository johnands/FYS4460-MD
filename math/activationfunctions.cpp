#include "activationfunctions.h"
#include <cmath>

namespace ActivationFunctions {

    arma::mat relu(arma::mat matrix) {
        // loop through vector of Wx + b and apply the relu function
        // i.e. replacing all negative elements with zeros
 /*       for (auto i : matrix[1]) {
            if (i < 0)
                i = 0;
        }
        return matrix;*/
    }


    arma::mat sigmoid(arma::mat matrix) {

        return 1./(1 + arma::exp(-matrix));
    }

    arma::mat reluDerivative(arma::mat matrix) {
        // the derivative of relu is 1 for positive numbers
        // and zero for negative numbers

        for (int i=0; i < arma::size(matrix)[1]; i++) {
            if (matrix(0,i) > 0)
                matrix(0,i) = 1;
            else
                matrix(0,i) = 0;
        }
    }

    arma::mat sigmoidDerivative(arma::mat matrix) {

        return sigmoid(matrix) % (1 - sigmoid(matrix));
        //return arma::exp(-matrix)/((1 + arma::exp(-matrix))*(1 + arma::exp(-matrix)));

    }
}

