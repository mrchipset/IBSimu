cmake_minimum_required(VERSION 3.0.0)
project(Vlasov2D VERSION 0.1.0)

find_package(PkgConfig)
pkg_check_modules(IBSimu ibsimu-1.0.6dev)

include_directories(${IBSimu_INCLUDE_DIRS})
LINK_DIRECTORIES(/usr/local/lib)
add_executable(Vlasov2D main.cpp)
target_link_libraries(Vlasov2D ${IBSimu_LIBRARIES})

