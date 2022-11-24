#ifndef BUTTERWORTH_H
#define BUTTERWORTH_H

#include <eigen3/Eigen/Dense>

namespace filter
{
class Butterworth
{
public:
  Butterworth(double wc, double dt, int n);
private:
  Eigen::VectorXd calculate_a_coeff(double wc, int n);
};

}  // namespace filter

#endif