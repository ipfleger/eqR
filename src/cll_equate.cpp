#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>
#include <vector>

// [[Rcpp::plugins(cpp11)]]

//' Evaluate the Exponential of a Polynomial (Vectorized C++ Version)
 //'
 //' @description
 //' This is the C++ equivalent of the `exp_polynomial` R function. It calculates
 //' the value of exp(P(x)) where P(x) is a polynomial. It is optimized to
 //' evaluate the polynomial for a vector of x values.
 //'
 //' @param x_vec A numeric vector of values at which to evaluate the function.
 //' @param params A numeric vector of polynomial coefficients (e.g., c0, c1, c2, ...).
 //' @return A numeric vector containing the results.
 //' @export
 // [[Rcpp::export]]
 Rcpp::NumericVector exp_polynomial_cpp(Rcpp::NumericVector x_vec, Rcpp::NumericVector params) {
    int n = x_vec.size();
    int p = params.size();
    Rcpp::NumericVector result(n);

    for (int i = 0; i < n; ++i) {
       double x = x_vec[i];
       double sum_poly = 0.0;
       for (int j = 0; j < p; ++j) {
          sum_poly += params[j] * std::pow(x, j);
       }
       result[i] = std::exp(sum_poly);
    }
    return result;
 }

//' Bivariate Exponential Polynomial (Vectorized C++ Version)
 //'
 //' @description
 //' Calculates the value of the bivariate log-linear density for a grid of X and Y scores.
 //'
 //' @param x_vec A numeric vector for the first variable.
 //' @param y_vec A numeric vector for the second variable.
 //' @param bivar_params A numeric vector of the bivariate polynomial coefficients.
 //' @param x_powers The powers of x for each coefficient.
 //' @param y_powers The powers of y for each coefficient.
 //' @return A numeric vector containing the evaluated density values.
 //' @export
 // [[Rcpp::export]]
 Rcpp::NumericVector biv_exp_polynomial_vec_cpp(Rcpp::NumericVector x_vec, Rcpp::NumericVector y_vec, Rcpp::NumericVector bivar_params, Rcpp::IntegerVector x_powers, Rcpp::IntegerVector y_powers) {
    int n = x_vec.size();
    if (y_vec.size() != n) {
       Rcpp::stop("x_vec and y_vec must have the same length.");
    }
    int num_params = bivar_params.size();
    Rcpp::NumericVector result(n);

    for (int i = 0; i < n; ++i) {
       double x = x_vec[i];
       double y = y_vec[i];
       double sum_poly = 0.0;
       for (int j = 0; j < num_params; ++j) {
          sum_poly += bivar_params[j] * std::pow(x, x_powers[j]) * std::pow(y, y_powers[j]);
       }
       result[i] = std::exp(sum_poly);
    }
    return result;
 }


//' Integrand for Marginal CDF Calculation (C++ Version)
 //'
 //' @description
 //' This function prepares the values for numerical integration to find the marginal
 //' density.
 //'
 //' @param x A numeric vector of x scores.
 //' @param y_nodes A numeric vector of quadrature nodes for the y-dimension.
 //' @param bivar_params A numeric vector of bivariate polynomial coefficients.
 //' @param x_powers Integer vector of x powers for the polynomial.
 //' @param y_powers Integer vector of y powers for the polynomial.
 //' @return A matrix where each column corresponds to an x value and each row
 //'         to a y_node, containing the bivariate density values.
 //' @export
 // [[Rcpp::export]]
 Rcpp::NumericMatrix integrand_fully_vectorized_cpp(Rcpp::NumericVector x, Rcpp::NumericVector y_nodes, Rcpp::NumericVector bivar_params, Rcpp::IntegerVector x_powers, Rcpp::IntegerVector y_powers) {
    int nx = x.size();
    int ny = y_nodes.size();
    Rcpp::NumericMatrix result(ny, nx);

    for (int i = 0; i < nx; ++i) {
       Rcpp::NumericVector x_col(ny, x[i]);
       Rcpp::NumericVector density_col = biv_exp_polynomial_vec_cpp(x_col, y_nodes, bivar_params, x_powers, y_powers);
       for (int j = 0; j < ny; ++j) {
          result(j, i) = density_col[j];
       }
    }
    return result;
 }
