cmake_minimum_required(VERSION 3.0.0)
project(butterworth-filter VERSION 0.1.0)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include(CTest)
enable_testing()

add_executable(butterworth-filter main.cpp butterworth.cpp butterworth.h)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
