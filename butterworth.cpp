/**
 * @file butterworth.cpp
 * @author Johan Ubbink (johan.ubbink@kuleuven.be)
 * @brief Basic implementation of a butterworth filter
 * @version 0.1
 * @date 2022-11-25
 *
 *
 */

#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>

#include "butterworth.h"

namespace butter
{
  Butterworth::Butterworth(double wc, double dt, int n, int size)
  {
    // calculate the a coefficients
    Eigen::VectorXd a_coeff = Butterworth::calculate_a_coeff(wc, n);

    // calculate the b coefficients
    // scaled for the cutoff frequency
    double b_coeff = std::pow(wc, n);

    // convert transfer function to continuous state space
    ContinuousSS cont_ss = Butterworth::tf2ss(a_coeff, b_coeff);

    // convert to a discrete system
    discrete_sys = Butterworth::continuous2discrete(cont_ss, dt);

    // initialize the states
    state_x = Eigen::MatrixXd::Zero(n, size);
  }

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

  ContinuousSS Butterworth::tf2ss(Eigen::VectorXd a_coeff, double b_coeff)
  {
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

  DiscreteSS Butterworth::continuous2discrete(ContinuousSS cont_sys, double dt)
  {
    // Create a discrete system
    DiscreteSS disc_sys;

    int n = cont_sys.Bc.size();
    Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(n, n);

    // Determine discrete a matrix
    disc_sys.Ad = Butterworth::expm(cont_sys.Ac * dt);
    // determine discrete B
    disc_sys.Bd = cont_sys.Ac.colPivHouseholderQr().solve((disc_sys.Ad - eye) * cont_sys.Bc);

    // C matrix stays the same
    disc_sys.Cd = cont_sys.Cc;

    return disc_sys;
  }

  std::vector<double> Butterworth::step(const std::vector<double> &u)
  {
    // convert input to eigen vector
    Eigen::RowVectorXd u_vec = Eigen::RowVectorXd::Map(u.data(), u.size());

    // apply the input
    state_x = discrete_sys.Ad * state_x + discrete_sys.Bd * u_vec;
    Eigen::RowVectorXd y_out = discrete_sys.Cd * state_x;

    // convert from eigen back to vector
    std::vector<double> y(y_out.data(), y_out.size() + y_out.data());

    return y;
  }

  Eigen::MatrixXd Butterworth::expm(Eigen::MatrixXd A)
  {
    int n = A.cols();

    Eigen::MatrixXd expA_k = Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd accumulator = Eigen::MatrixXd::Identity(n, n);
    for (int i = 1; i < 10; i++)
    {
      accumulator = accumulator * (A / i);
      expA_k = expA_k + accumulator;
    }

    return expA_k;
  }
} // namespace filter