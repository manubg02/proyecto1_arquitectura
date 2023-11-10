cmake_minimum_required(VERSION 3.26)
project(proyecto1_arquitectura)

set(CMAKE_CXX_STANDARD 17) # Configura el estándar de C++ que desees utilizar, como 17 o 20

add_executable(fluid
        fluid/fluid.cpp
        sim/progargs.hpp
        sim/progargs.cpp
        sim/grid.hpp
        sim/grid.cpp
        sim/block.cpp
        sim/block.hpp)