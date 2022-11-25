#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>

#include "butterworth.h"

namespace filter
{
/**
 * @brief Construct a new Butterworth:: Butterworth object
 *
 * @param wc Cutoff frequency of filter [rad/s]
 * @param dt Sampling time of filter
 * @param n  Order of butterworth filter
 * @param size Number of input channels
 */
Butterworth::Butterworth(double wc, double dt, int n, int size)
{
  std::cout << "Creating a filter with: \n";
  std::cout << "    wc: " << wc << "\n";
  std::cout << "    dt: " << dt << "\n";
  std::cout << "    n : " << n << "\n";

  // calculate the a coefficients
  Eigen::VectorXd a_coeff = Butterworth::calculate_a_coeff(wc, n);

  // calculate the b coefficients
  // scaled for the cutoff frequency
  double b_coeff = std::pow(wc, n);

  // convert transfer function to continuous state space
  ContinuousSS cont_ss = Butterworth::tf2ss(a_coeff, b_coeff);

  // std::cout << cont_ss.Ac << "\n";
  discrete_sys = Butterworth::continuous2discrete(cont_ss, dt);
  // std::cout << discrete_sys.Ad << "\n";
  // std::cout << discrete_sys.Bd << "\n";
  // initialize the states
  state_x = Eigen::VectorXd::Zero(n, size);

  // std::cout << state_x << '\n';
}

/**
 * @brief Calculates coefficients for the butterworth tf
 *
 * based on https://en.wikipedia.org/wiki/Butterworth_filter
 * recursive formulation.
 *
 * @param wc Cutoff frequency of filter [rad/s]
 * @param n Order of butterworth filter
 * @return Eigen::VectorXd Coefficients of transfer function
 */
Eigen::VectorXd Butterworth::calculate_a_coeff(double wc, int n)
{
  // vector to store the coefficients
  Eigen::VectorXd a_coeff(n + 1);

  // first coefficient is always 1
  a_coeff[0] = 1;

  double gamma = (M_PI) / (2 * n);

  // calculate the coefficients
  for (int i = 0; i < n; i++)
  {
    a_coeff[i + 1] = (std::cos(i * gamma) * a_coeff[i]) / (std::sin((i + 1) * gamma));
  }

  // change normalized coefficients for cutoff frequency
  for (int i = 0; i < n + 1; i++)
    a_coeff[i] = a_coeff[i] * std::pow(wc, i);

  return a_coeff;
}

/**
 * @brief Convert transfer function to continuous state space
 *
 * @param a_coeff Coefficients of the denominator
 * @param b_coeff Coefficients of the numerator
 */
ContinuousSS Butterworth::tf2ss(Eigen::VectorXd a_coeff, double b_coeff)
{
  std::cout << "Converting transfer function to ss \n";

  // Determine size of the system
  int n = a_coeff.size() - 1;

  // Create the A matrix (control form)
  Eigen::MatrixXd Ac = Eigen::MatrixXd::Zero(n, n);
  // set bottom row to coefficients
  Ac.block(n - 1, 0, 1, n) = -a_coeff.segment(1, n).reverse().transpose();
  // set upper block to identity
  Ac.block(0, 1, n - 1, n - 1) = Eigen::MatrixXd::Identity(n - 1, n - 1);

  // Create the B vector
  Eigen::VectorXd Bc = Eigen::VectorXd::Zero(n, 1);
  Bc[n - 1] = b_coeff;

  // Create the C vector
  Eigen::RowVectorXd C = Eigen::RowVectorXd::Zero(1, n);
  C[0] = 1;

  // Store within a structure
  ContinuousSS cont_system;
  cont_system.Ac = Ac;
  cont_system.Bc = Bc;
  cont_system.Cc = C;

  return cont_system;
}

/**
 * @brief Converts continuous system into a discrete system
 *
 * @param cont_sys A continuous state space model
 * @param dt Sampling time of discrete system
 * @return DiscreteSS An equivelent discrete system
 */
DiscreteSS Butterworth::continuous2discrete(ContinuousSS cont_sys, double dt)
{
  std::cout << "Convert continuous system to a discrete system \n";

  // Create a discrete system
  DiscreteSS disc_sys;

  int n = cont_sys.Bc.size();
  Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(n, n);

  // Determine discrete a matrix
  // use (approximate) matrix exponential
  disc_sys.Ad = eye + cont_sys.Ac * dt + cont_sys.Ac * cont_sys.Ac * dt * dt / 2 +
                cont_sys.Ac * cont_sys.Ac * cont_sys.Ac * dt * dt * dt / (2 * 3);

  // determine discrete B
  disc_sys.Bd = cont_sys.Ac.colPivHouseholderQr().solve((disc_sys.Ad - eye) * cont_sys.Bc);

  // C matrix stays the same
  disc_sys.Cd = cont_sys.Cc;

  return disc_sys;
}

/**
 * @brief Apply input to filter and receive output
 *
 * @param u Input applied
 * @return std::vector<double> Ouput of filter
 */
std::vector<double> Butterworth::step(std::vector<double> u)
{
  // convert input to eigen vector
  Eigen::VectorXd u_vec = Eigen::VectorXd::Map(u.data(), u.size());

  // apply the input
  state_x = discrete_sys.Ad * state_x + discrete_sys.Bd * u_vec;
  Eigen::RowVectorXd y_out = discrete_sys.Cd * state_x;

  // convert from eigen back to vector
  std::vector<double> y (y_out.data(), y_out.size() + y_out.data());
 
  return y;
}
}  // namespace filter