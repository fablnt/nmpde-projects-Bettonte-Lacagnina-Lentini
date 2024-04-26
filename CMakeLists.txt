cmake_minimum_required(VERSION 3.12.0)
project(FisherKolmogorov LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(main3D src/main.cpp src/FisherKolmogorov.cpp)
deal_ii_setup_target(main3D)
add_executable(main1D src/main1D.cpp src/FisherKolmogorov1D.cpp)
deal_ii_setup_target(main1D)
