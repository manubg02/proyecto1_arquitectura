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
    try{
        guardar_particulas(particles, argv[2]);
    }catch (const std::runtime_error &e){
        std::cerr << "Error: " << e.what() << std::endl;
        return -5;
    }

    for (auto i = 0u; i < particles.px.size(); ++i) {
        std::cout << "Particula " << i+1 << " : ("
                  << particles.px[i] << ", "
                  << particles.py[i] << ", "
                  << particles.pz[i] << ", "
                  << particles.hvx[i] << ", "
                  << particles.hvy[i] << ", "
                  << particles.hvz[i] << ", "
                  << particles.vx[i] << ", "
                  << particles.vy[i] << ", "
                  << particles.vz[i] << ")"
                  << " (" << particles.i[i] << ", "
                  << particles.j[i] << ", "
                  << particles.k[i] << ")" <<std::endl;
    }


    return 0;
}