cmake_minimum_required(VERSION 3.6.1)

project(Poisson_Test)

add_executable(Poisson_Test Poisson_Test.cpp)

target_link_libraries(Poisson_Test OedLib)
target_link_libraries(Poisson_Test OedLibImpl)

#add_subdirectory(../../ oed-lib/)  # Include the src directory
add_subdirectory(../impl oed-lib-impl/)  # Include the src directory
