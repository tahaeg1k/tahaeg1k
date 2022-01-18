cmake_minimum_required(VERSION 3.21)
project(base)

set(CMAKE_CXX_STANDARD 14)

include_directories(~/CLionProjects/tahaeg1k/options-pricing/include)

add_library(base OptionsPricingModel.cpp)

add_executable(base OptionsPricingModel.cpp)
