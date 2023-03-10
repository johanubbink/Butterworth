# Butterworth filter

This library implements a Butterworth filter in c++. The library allows for specifying the cutoff frequency, the sampling time and the order of the filter. The filter can be used for filtering of signals with multiple input channels in parallel.

## Dependencies
- eigen3


## Installation

Make sure you have the dependencies installed.
```bash
sudo apt install libeigen3-dev
```

Install the code with:
```bash
git clone git@github.com:johanubbink/Butterworth.git
cd butterworth
mkdir build
cd build
cmake ..
make
sudo make install
```


## Usage

```c++

#include <iostream>
#include <vector>
#include "butterworth.h"

int main(int, char **)
{
  // create a butterworth filter with 
  //    cutoff frequency of 50 rad/s
  //    sampling time of 1/200 s
  //    2 input channels filtered in parallel
  butter::Butterworth filter {50, 1.0 / 200, 4, 2};

  // apply input to filter
  std::vector<double> u{4, 3};

  // receive output
  std::vector<double> y = filter.step(u);

  // print output
  std::cout << y[0] << "\n";
  std::cout << y[1] << "\n";
}
    
```

## Implementation details

The following steps are used to implement the filter:
-  The coefficients of the transfer function are calculated using the recursive formula [recursive formula](https://en.wikipedia.org/wiki/Butterworth_filter#Normalized_Butterworth_polynomials).
-  The continuous time transfer function is transformed into a continuous time state space model in the control canonical form.
- The continuous time state space model is converted to a discrete time state space model using the matrix exponential.
- This discrete time state space model is used to implement the filter.
