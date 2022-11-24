#include <iostream>

#include "butterworth.h"


int main(int, char**) {

    // create a filter
    filter::Butterworth filter = filter::Butterworth(50, 1/200, 3, 1);



}
