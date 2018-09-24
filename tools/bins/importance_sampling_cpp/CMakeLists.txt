cmake_minimum_required(VERSION 3.1)

# Configure the project
#----------------------
project(importance_sampling)
set(CMAKE_CXX_STANDARD 11)
set(EXECUTABLE_OUTPUT_PATH ../)
set(LIBRARY_OUTPUT_PATH lib/)

# Check required libraries
#----------------------
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)

if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

find_package( Boost 1.47 COMPONENTS random)
if (NOT Boost_FOUND)
    message(STATUS "C++ Boost::random is not installed. Please install it.")
endif ()

# Generate library reconstruction
#--------------------------------
file( GLOB LIB_SOURCES lib/*.cpp)
file( GLOB LIB_HEADERS lib/*.hpp)
add_library( reconstruction ${LIB_SOURCES} ${LIB_HEADERS} )

# Generate file list for main
#-------------------------------
file( GLOB source_files src/* )

# Include libraries
#------------------
include_directories(lib/)
link_directories(lib/)


# Declare executables
#--------------------
add_executable( random_sampling_std src/random_sampling_std.cpp )


# Link libraries
#---------------
target_link_libraries( random_sampling_std reconstruction )
