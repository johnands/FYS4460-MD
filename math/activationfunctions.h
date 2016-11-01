#include <armadillo>

namespace ActivationFunctions {
    arma::mat relu(arma::mat matrix);
    arma::mat sigmoid(arma::mat matrix);
    arma::mat reluDerivative(arma::mat matrix);
    arma::mat sigmoidDerivative(arma::mat matrix);
}
