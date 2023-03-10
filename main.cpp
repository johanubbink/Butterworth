#include <iostream>
#include <vector>
#include "butterworth.h"

int main(int, char **)
{
  // create a butterworth filter with cutoff frequency of 50 rad/s
  // and sampling time of 1/200 s
  // und 2 input channels filtered in parallel
  butter::Butterworth filter{50, 1.0 / 200, 4, 2};

  // apply input to filter
  std::vector<double> u{4, 3};

  // receive output
  std::vector<double> y = filter.step(u);

  // print output
  std::cout << y[0] << "\n";
  std::cout << y[1] << "\n";
}
