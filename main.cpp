#include <iostream>
#include <vector>
#include "butterworth.h"


int main(int, char**) {

    // create a filter
    filter::Butterworth filter = filter::Butterworth(50, 1.0/200, 3, 1);

    std::vector<double> u {1};

    std::vector<double> y = filter.step(u);

    std::cout << y[0] << "\n";


}
