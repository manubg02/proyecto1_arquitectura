//
// Created by manu on 15/10/23.
//

#include <iostream>
#include <string>
#include "../sim/progargs.hpp"
#include "../sim/grid.hpp"

int main(int argc, char *argv[]) {

    ProgArgs progArgs(argc, argv);


    // Verificar si hubo errores al analizar los argumentos
    if (progArgs.getErrorCode() != 0) {
        std::cerr << "Error: " << progArgs.getErrorMessage() << std::endl;
        return progArgs.getErrorCode();
    }

    int timeSteps = progArgs.getTimeSteps();
    std::string inputFile = progArgs.getInputFile();
    std::string outputFile = progArgs.getOutputFile();

    Grid myGrid(timeSteps, inputFile);


    return 0;
}