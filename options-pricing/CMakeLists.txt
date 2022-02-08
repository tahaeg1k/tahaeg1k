cmake_minimum_required(VERSION 3.21)

project(pricing DESCRIPTION "CMakeLists for pricing")

set(CMAKE_CXX_STANDARD 14)

include_directories(include)
include_directories(base)
include_directories(black-scholes)
include_directories(monte-carlo)

add_library(lib
        include/OptionsPricingModel.h
        include/OptionsPricingModel.cpp
        black-scholes/BlackScholes.h black-scholes/BlackScholes.cpp
        monte-carlo/MonteCarlo.h monte-carlo/MonteCarlo.cpp include/random_singleton.cpp include/random_singleton.h)

target_include_directories(lib PUBLIC include)
target_include_directories(lib PUBLIC base)
target_include_directories(lib PUBLIC black-scholes)
target_include_directories(lib PUBLIC monte-carlo)

add_executable(random_test main_rd.cpp)
add_executable(pricing_test main_bs.cpp)
add_executable(asian_test main_asian.cpp)

target_link_libraries(random_test lib)
target_link_libraries(pricing_test lib)
target_link_libraries(asian_test lib)

