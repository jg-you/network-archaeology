cmake_minimum_required(VERSION 3.1)

# Configure the project
#----------------------
project(importance_sampling)
set(CMAKE_CXX_STANDARD 11)
set(EXECUTABLE_OUTPUT_PATH ../)
set(LIBRARY_OUTPUT_PATH lib/)

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
