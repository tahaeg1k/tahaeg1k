cmake_minimum_required(VERSION 3.21)
project(blackscholes)

set(CMAKE_CXX_STANDARD 14)

# Add the given directories to the ones the compiler uses to search for include files
include_directories(~/CLionProjects/tahaeg1k/options-pricing/include)
include_directories(~/CLionProjects/tahaeg1k/options-pricing/base)

#add_subdirectory(/Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/base)
#find_library(BASE_LIB base_library)

add_executable(blackscholes BlackScholes.h BlackScholes.cpp main.cpp)
target_link_libraries(blackscholes PUBLIC include_library)

# Target link should come after add_executable :
#target_link_libraries(blackscholes PUBLIC "${BASE_LIB}")
#target_link_libraries(blackscholes PUBLIC base)