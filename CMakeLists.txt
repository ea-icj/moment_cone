# See https://learn.scientific-python.org/development/guides/packaging-compiled/
cmake_minimum_required(VERSION 3.15...3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 17)

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_set_of_vector MODULE src/moment_cone/_set_of_vector.cpp)
install(TARGETS _set_of_vector DESTINATION ${SKBUILD_PROJECT_NAME})
