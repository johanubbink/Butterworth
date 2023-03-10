#ifndef BUTTERWORTH_H
#define BUTTERWORTH_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <vector>

namespace filter
{
  /**
   * @brief Struct to represent continuous state space model.
   *
   */
  struct ContinuousSS
  {
    Eigen::MatrixXd Ac;
    Eigen::VectorXd Bc;
    Eigen::RowVectorXd Cc;
  };

  struct DiscreteSS
  {
    Eigen::MatrixXd Ad;
    Eigen::VectorXd Bd;
    Eigen::RowVectorXd Cd;
  };

  class Butterworth
  {
  public:
    /**
     * @brief Construct a new Butterworth:: Butterworth object
     *
     * @param wc Cutoff frequency of filter [rad/s]
     * @param dt Sampling time of filter
     * @param n  Order of butterworth filter
     * @param size Number of input channels
     */
    Butterworth(double wc, double dt, int n, int size);

    /**
     * @brief Apply input to filter and receive output
     *
     * @param u Input applied
     * @return std::vector<double> Ouput of filter
     */
    std::vector<double> step(const std::vector<double> &u);

  private:
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
    Eigen::VectorXd calculate_a_coeff(double wc, int n);

    /**
     * @brief Convert transfer function to continuous state space
     *
     * @param a_coeff Coefficients of the denominator
     * @param b_coeff Coefficients of the numerator
     */
    ContinuousSS tf2ss(Eigen::VectorXd a_coeff, double b_coeff);

    /**
     * @brief Converts continuous system into a discrete system
     *
     * @param cont_sys A continuous state space model
     * @param dt Sampling time of discrete system
     * @return DiscreteSS An equivelent discrete system
     */
    DiscreteSS continuous2discrete(ContinuousSS cont_sys, double dt);

        DiscreteSS discrete_sys;
    Eigen::MatrixXd state_x;

    /**
     * @brief Calculates (approximate) matrix exponential
     * https://en.wikipedia.org/wiki/Matrix_exponential
     * Based on the power series
     * @param A Input matrix
     * @return Eigen::MatrixXd Exponential of the matrix
     */
    Eigen::MatrixXd expm(Eigen::MatrixXd A);
  };

} // namespace filter

#endif