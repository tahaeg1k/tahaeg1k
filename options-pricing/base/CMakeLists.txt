cmake_minimum_required(VERSION 3.21)
project(base DESCRIPTION "CMakeLists for the base folder")

set(CMAKE_CXX_STANDARD 14)

# Add the given directories to the ones the compiler uses to search for include files
include_directories(~/CLionProjects/tahaeg1k/options-pricing/include)

# This is used to locate ourselves :
message("This is the base source directory : ${base_SOURCE_DIR}")
# Sets a list of necessary headers
set(INCLUDE_HEADER_LIST "~/CLionProjects/tahaeg1k/options-pricing/include/Cdf.h"
        "~/CLionProjects/tahaeg1k/options-pricing/include/OptionsPricingModel.h")

# Library target to be constructed from the source files listed in the command invocation.
# For now we set it to STATIC ( to be reviewed later). Any targets that link to base_library will need
# to have the include directory (role of PUBLIC)

add_library(include_library STATIC OptionsPricingModel.cpp ${INCLUDE_HEADER_LIST})

# We need include directory, and users of our library will need it too
#target_include_directories(include_library PUBLIC ~/CLionProjects/tahaeg1k/options-pricing/include)
message(STATUS "Done include_library ${include_library} here!")
find_library(INC_LIB include_library)
message("found ${INC_LIB}")


add_library(base_library STATIC OptionsPricingModel.cpp)
target_include_directories(base_library PUBLIC ~/CLionProjects/tahaeg1k/options-pricing/include)
message(STATUS "Done base_library ${base_library} here!")

# CMake uses this compile source files (.cpp), the headers are untouched.
add_executable(base OptionsPricingModel.cpp)
#target_link_libraries(base PUBLIC include_library)
target_compile_options(base PUBLIC -Wall -Wextra -Wpedantic)

# IDEs should put the headers in a nice place
source_group(
        TREE "~/CLionProjects/tahaeg1k/options-pricing/include"
        PREFIX "HEADER_FILES"
        FILES ${INCLUDE_HEADER_LIST}
)