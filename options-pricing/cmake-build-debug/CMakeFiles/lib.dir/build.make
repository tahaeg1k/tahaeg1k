# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/lib.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/lib.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/lib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/lib.dir/flags.make

CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.o: CMakeFiles/lib.dir/flags.make
CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.o: ../base/OptionsPricingModel.cpp
CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.o: CMakeFiles/lib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.o -MF CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.o.d -o CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.o -c /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/base/OptionsPricingModel.cpp

CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/base/OptionsPricingModel.cpp > CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.i

CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/base/OptionsPricingModel.cpp -o CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.s

CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.o: CMakeFiles/lib.dir/flags.make
CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.o: ../black-scholes/BlackScholes.cpp
CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.o: CMakeFiles/lib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.o -MF CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.o.d -o CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.o -c /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/black-scholes/BlackScholes.cpp

CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/black-scholes/BlackScholes.cpp > CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.i

CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/black-scholes/BlackScholes.cpp -o CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.s

CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.o: CMakeFiles/lib.dir/flags.make
CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.o: ../financial-products/EuropeanVanilla.cpp
CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.o: CMakeFiles/lib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.o -MF CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.o.d -o CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.o -c /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/financial-products/EuropeanVanilla.cpp

CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/financial-products/EuropeanVanilla.cpp > CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.i

CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/financial-products/EuropeanVanilla.cpp -o CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.s

CMakeFiles/lib.dir/binomial/Binomial.cpp.o: CMakeFiles/lib.dir/flags.make
CMakeFiles/lib.dir/binomial/Binomial.cpp.o: ../binomial/Binomial.cpp
CMakeFiles/lib.dir/binomial/Binomial.cpp.o: CMakeFiles/lib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/lib.dir/binomial/Binomial.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lib.dir/binomial/Binomial.cpp.o -MF CMakeFiles/lib.dir/binomial/Binomial.cpp.o.d -o CMakeFiles/lib.dir/binomial/Binomial.cpp.o -c /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/binomial/Binomial.cpp

CMakeFiles/lib.dir/binomial/Binomial.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lib.dir/binomial/Binomial.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/binomial/Binomial.cpp > CMakeFiles/lib.dir/binomial/Binomial.cpp.i

CMakeFiles/lib.dir/binomial/Binomial.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lib.dir/binomial/Binomial.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/binomial/Binomial.cpp -o CMakeFiles/lib.dir/binomial/Binomial.cpp.s

CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.o: CMakeFiles/lib.dir/flags.make
CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.o: ../monte-carlo/MonteCarlo.cpp
CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.o: CMakeFiles/lib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.o -MF CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.o.d -o CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.o -c /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/monte-carlo/MonteCarlo.cpp

CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/monte-carlo/MonteCarlo.cpp > CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.i

CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/monte-carlo/MonteCarlo.cpp -o CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.s

CMakeFiles/lib.dir/include/random_singleton.cpp.o: CMakeFiles/lib.dir/flags.make
CMakeFiles/lib.dir/include/random_singleton.cpp.o: ../include/random_singleton.cpp
CMakeFiles/lib.dir/include/random_singleton.cpp.o: CMakeFiles/lib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/lib.dir/include/random_singleton.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lib.dir/include/random_singleton.cpp.o -MF CMakeFiles/lib.dir/include/random_singleton.cpp.o.d -o CMakeFiles/lib.dir/include/random_singleton.cpp.o -c /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/include/random_singleton.cpp

CMakeFiles/lib.dir/include/random_singleton.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lib.dir/include/random_singleton.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/include/random_singleton.cpp > CMakeFiles/lib.dir/include/random_singleton.cpp.i

CMakeFiles/lib.dir/include/random_singleton.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lib.dir/include/random_singleton.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/include/random_singleton.cpp -o CMakeFiles/lib.dir/include/random_singleton.cpp.s

# Object files for target lib
lib_OBJECTS = \
"CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.o" \
"CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.o" \
"CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.o" \
"CMakeFiles/lib.dir/binomial/Binomial.cpp.o" \
"CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.o" \
"CMakeFiles/lib.dir/include/random_singleton.cpp.o"

# External object files for target lib
lib_EXTERNAL_OBJECTS =

liblib.a: CMakeFiles/lib.dir/base/OptionsPricingModel.cpp.o
liblib.a: CMakeFiles/lib.dir/black-scholes/BlackScholes.cpp.o
liblib.a: CMakeFiles/lib.dir/financial-products/EuropeanVanilla.cpp.o
liblib.a: CMakeFiles/lib.dir/binomial/Binomial.cpp.o
liblib.a: CMakeFiles/lib.dir/monte-carlo/MonteCarlo.cpp.o
liblib.a: CMakeFiles/lib.dir/include/random_singleton.cpp.o
liblib.a: CMakeFiles/lib.dir/build.make
liblib.a: CMakeFiles/lib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX static library liblib.a"
	$(CMAKE_COMMAND) -P CMakeFiles/lib.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/lib.dir/build: liblib.a
.PHONY : CMakeFiles/lib.dir/build

CMakeFiles/lib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/lib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/lib.dir/clean

CMakeFiles/lib.dir/depend:
	cd /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/cmake-build-debug /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/cmake-build-debug /Users/mohamedtahabennani/CLionProjects/tahaeg1k/options-pricing/cmake-build-debug/CMakeFiles/lib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/lib.dir/depend

