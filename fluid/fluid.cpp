//
// Created by manu on 15/10/23.
//

#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include "../sim/progargs.hpp"
#include "../sim/grid.hpp"
#include <typeinfo>
int main(int argc, char *argv[]) {

    ProgArgs progArgs(argc, argv);


    // Verificar si hubo errores al analizar los argumentos
    if (progArgs.getErrorCode() != 0) {
        std::cerr << "Error: " << progArgs.getErrorMessage() << std::endl;
        return progArgs.getErrorCode();
    }

    ParticleArray particles;
    guardar_particulas(particles, argv[2]);


    return 0;
}