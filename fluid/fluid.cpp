//
// Created by manu on 15/10/23.
//

#include <iostream>
#include <string>
#include "../sim/progargs.hpp"
#include "../sim/grid.hpp"

// Función principal
int main(int argc, char *argv[]) {
    // Analiza los argumentos de línea de comandos
    ProgArgs progArgs(argc, argv);

    // Verifica si hubo errores al analizar los argumentos
    if (progArgs.getErrorCode() != 0) {
        std::cerr << "Error: " << progArgs.getErrorMessage() << std::endl;
        return progArgs.getErrorCode();
    }

    // Obtiene la cantidad de pasos de tiempo, el archivo de entrada y el archivo de salida desde los argumentos
    int timeSteps = progArgs.getTimeSteps();
    const std::string& inputFile = progArgs.getInputFile();
    const std::string& outputFile = progArgs.getOutputFile();

    // Crea una instancia del objeto Grid con los parámetros proporcionados
    Grid myGrid(timeSteps, inputFile, outputFile);

    return 0;
}