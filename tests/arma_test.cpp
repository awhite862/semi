#include <armadillo>

/** \test Validate that arma addition compliles and works correctly. **/
int run_addition_test() {
    arma::mat A(2,2), B(2,2), C(2,2);
    A.fill(1.0); B.fill(1.0); C.fill(2.0);
    arma::mat test = A + B;
    if ((C - test).max() > 1e-12) return 1;
    return 0;
}

/** \test Validate that arma multiplication compliles and works correctly. **/
int run_multiplication_test() {
    arma::mat A(2,2), I(2,2);
    A.fill(2.0); I.eye();
    arma::mat test = I * A;
    if ((test - A).max() > 1e-12) return 1;
    return 0;
}

/** \test Validate that arma eigensolver compiles and works correctly/ **/
int run_eigensolver_test() {
    arma::mat A(2,2); 
    A.fill(0.0); A(0,0) = 1.0; A(1,1) = 2.0;
    arma::vec E;
    arma::mat C;
    arma::eig_sym(E, C, A);
    if (std::abs(E(0) - 1.0) > 1e-12) return 1;
    if (std::abs(E(1) - 2.0) > 1e-12) return 1;
    return 0;
}

int main() {
    return 

    run_addition_test() |
    run_multiplication_test() |
    run_eigensolver_test() |

    0;
}
