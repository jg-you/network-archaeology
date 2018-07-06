include_directories(${BOOST_INCLUDEDIR})

add_library(growth_models 
            unstructured_models/delayed.cpp
            structured_models/gn.cpp
            structured_models/generalized_gn.cpp
            structured_models/bianconi_barabasi.cpp)
add_executable(growth growth_main.cpp)
target_link_libraries(growth growth_models)
target_link_libraries(growth ${Boost_LIBRARIES})
set_target_properties(growth PROPERTIES RUNTIME_OUTPUT_DIRECTORY ../ )

