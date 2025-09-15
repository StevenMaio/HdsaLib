# Create an interface library for HDSALibROL
add_library(HDSALibROL INTERFACE)

# Specify include directories for the interface library
target_include_directories(HDSALibROL INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/rol>
    ${Trilinos_INCLUDE_DIRS}
)

# Optionally, link against Trilinos libraries if needed
# target_link_libraries(HDSALibROL INTERFACE ${Trilinos_LIBRARIES})