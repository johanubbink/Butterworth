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
  Butterworth(double wc, double dt, int n, int size);
  std::vector<double> step(std::vector<double> u);

private:
  Eigen::VectorXd calculate_a_coeff(double wc, int n);
  ContinuousSS tf2ss(Eigen::VectorXd a_coeff, double b_coeff);
  DiscreteSS continuous2discrete(ContinuousSS cont_sys, double dt);
  DiscreteSS discrete_sys;
  Eigen::MatrixXd state_x;
  Eigen::MatrixXd expm(Eigen::MatrixXd A);
};

}  // namespace filter

#endif