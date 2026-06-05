cmake_minimum_required(VERSION 3.6.1)

project(Poisson_Test)

add_executable(Poisson_OED Poisson_OED.cpp)
add_executable(Poisson_Inversion Poisson_Inversion.cpp)

target_link_libraries(Poisson_OED OedLib)
target_link_libraries(Poisson_OED OedLibImpl)

target_link_libraries(Poisson_Inversion OedLib)
target_link_libraries(Poisson_Inversion OedLibImpl)

#add_subdirectory(../../ oed-lib/)  # Include the src directory
add_subdirectory(../impl oed-lib-impl/)  # Include the src directory
