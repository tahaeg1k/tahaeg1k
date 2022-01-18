cmake_minimum_required(VERSION 3.21)
project(blackscholes)

set(CMAKE_CXX_STANDARD 14)

include_directories(~/CLionProjects/tahaeg1k/options-pricing/include)
include_directories(~/CLionProjects/tahaeg1k/options-pricing/base)

#add_subdirectory(/Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/base)

add_executable(blackscholes BlackScholes.h BlackScholes.cpp main.cpp)

target_link_libraries(blackscholes PUBLIC base)