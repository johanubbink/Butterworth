cmake_minimum_required(VERSION 3.0.0)
project(butterworth-filter VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 11)


find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include(CTest)
enable_testing()

include_directories(${EIGEN3_INCLUDE_DIR})

add_library(butterworth butterworth.cpp butterworth.h)
target_link_libraries(butterworth Eigen3::Eigen)

add_executable(butterworth-example main.cpp)
target_link_libraries(butterworth-example butterworth)

# install 
install(TARGETS butterworth DESTINATION lib)
install(FILES butterworth.h DESTINATION include)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)


