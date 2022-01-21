cmake_minimum_required(VERSION 3.21)

project(pricing DESCRIPTION "CMakeLists for pricing")

set(CMAKE_CXX_STANDARD 14)

include_directories(include)
include_directories(base)
include_directories(black-scholes)
include_directories(binomial)
include_directories(monte-carlo)

add_library(lib
        include/Cdf.h include/OptionsPricingModel.h
        base/OptionsPricingModel.cpp
        black-scholes/BlackScholes.h black-scholes/BlackScholes.cpp
        binomial/Binomial.h binomial/Binomial.cpp
        monte-carlo/MonteCarlo.h monte-carlo/MonteCarlo.cpp)

target_include_directories(lib PUBLIC include)
target_include_directories(lib PUBLIC base)
target_include_directories(lib PUBLIC black-scholes)
target_include_directories(lib PUBLIC binomial)
target_include_directories(lib PUBLIC monte-carlo)

add_executable(pricing_test main_bs.cpp)
target_link_libraries(pricing_test lib)

