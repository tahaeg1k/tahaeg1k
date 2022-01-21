cmake_minimum_required(VERSION 3.21)
project(blackscholes)

set(CMAKE_CXX_STANDARD 14)

# Add the given directories to the ones the compiler uses to search for include files
include_directories(~/CLionProjects/tahaeg1k/options-pricing/include)
include_directories(~/CLionProjects/tahaeg1k/options-pricing/base)

#find_library(BASE_LIB base_library)
link_directories(~/CLionProjects/tahaeg1k/options-pricing/)
find_library (
        INC_LIB
        NAMES include_library # what to look for
        HINTS "~/CLionProjects/tahaeg1k/options-pricing/base/" # where to look
)
message(STATUS "INC_LIB: [${INC_LIB}]")

add_executable(blackscholes BlackScholes.h BlackScholes.cpp main.cpp)



# check if we found the library

if(INC_LIB)
    target_link_libraries(blackscholes PUBLIC ${INC_LIB})
else()
    message(Cannot link black scholes to INC_LIB)
endif()



#target_link_libraries(blackscholes PUBLIC include_library)

# Target link should come after add_executable :
#target_link_libraries(blackscholes PUBLIC "${BASE_LIB}")
#target_link_libraries(blackscholes PUBLIC base)