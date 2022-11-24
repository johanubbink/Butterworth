#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>

#include "butterworth.h"

namespace filter
{
Butterworth::Butterworth(double wc, double dt, int n)
{
  std::cout << "Creating a filter with: \n";
  std::cout << "    wc: " << wc << "\n";
  std::cout << "    dt: " << dt << "\n";
  std::cout << "    n : " << n << "\n";

  // calculate the a coefficients
  Eigen::VectorXd a_coeff = Butterworth::calculate_a_coeff(wc, n);
  std::cout << "a_coeff : \n" << a_coeff << "\n";

  // calculate the b coefficients
  // scaled for the cutoff frequency
  Eigen::VectorXd b_coeff(1);
  b_coeff << std::pow(wc,n);
  std::cout << "b_coeff : \n" << b_coeff << "\n";
}

/**
 * @brief Calculates coefficients for the butterworth tf
 * 
 * based on https://en.wikipedia.org/wiki/Butterworth_filter
 * recursive formulation.
 * 
 * @param wc 
 * @param n 
 * @return Eigen::VectorXd 
 */
Eigen::VectorXd Butterworth::calculate_a_coeff(double wc, int n)
{
  // vector to store the coefficients
  Eigen::VectorXd a_coeff(n+1);

  // first coefficient is always 1
  a_coeff[0] = 1;

  double gamma =  (M_PI)/(2*n);

  // calculate the coefficients
  for (int i = 0; i < n; i++)
  {
    std::cout << i << "\n";
    a_coeff[i+1] = (std::cos(i*gamma)*a_coeff[i])/(std::sin((i+1)*gamma));
  }

  // change normalized coefficients for cutoff frequency
  for (int i = 0; i < n+1; i++)
    a_coeff[i] = a_coeff[i]*std::pow(wc,i);

  return a_coeff;
}

}  // namespace filter