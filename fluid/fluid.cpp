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

int main(int argc, char *argv[]) {

    ProgArgs progArgs(argc, argv);

    // Verificar si hubo errores al analizar los argumentos
    if (progArgs.getErrorCode() != 0) {
        std::cerr << "Error: " << progArgs.getErrorMessage() << std::endl;
        return progArgs.getErrorCode();
    }

    std::ifstream file("../small.fld", std::ios::binary);

    if (file) {
        // Lee los primeros 4 bytes como un entero de 32 bits
        int32_t intValue;
        file.read(reinterpret_cast<char *>(&intValue), sizeof(intValue));

        // Interpreta el entero como un valor de coma flotante de 32 bits
        float ppm = *reinterpret_cast<float *>(&intValue);

        // Imprime el valor en la consola
        std::cout << "Particulas por metro: " << ppm << std::endl;

        //file.seekg(, std::ios::cur);
        file.read(reinterpret_cast<char *>(&intValue), sizeof(intValue));
        int np = *reinterpret_cast<int *>(&intValue);
        std::cout << "Numero de particulas: " << np << std::endl;


        int count = 0; // Inicializa count antes del bucle
        while (count < np) {
            std::cout << "Particula " << count  << ": (";
            for (int i = 0; i < 8; i++) {
                file.read(reinterpret_cast<char *>(&intValue), sizeof(intValue));
                float particula = *reinterpret_cast<float *>(&intValue);
                std::cout << std::fixed << std::setprecision(8) << particula << ", ";


            }
            file.read(reinterpret_cast<char *>(&intValue), sizeof(intValue));
            float lastParticula = *reinterpret_cast<float *>(&intValue);
            std::cout << std::fixed << std::setprecision(8) << lastParticula << ")" << std::endl;
            count += 1;
        }

        return 0;
    }
}