//
// Created by manu on 8/11/23.
//
#include "progargs.hpp"
#include <iostream>
#include <cctype>
#include <fstream>

ProgArgs::ProgArgs(int argc, char *argv[]) {
    parseArguments(argc, argv);
}

void ProgArgs::parseArguments(int argc, char *argv[]) {
    int32_t intValue;
    // Verificar el número correcto de argumentos
    if (argc != 4) {
        errorMessage = "Invalid number of arguments.";
        errorCode = -1;
        return;
    }

    // Obtener el número de pasos de tiempo
    try {
        timeSteps = std::stoi(argv[1]);

        // Verifica si cada carácter de la cadena es un dígito
        for (char c : std::string(argv[1])) {
            if (!std::isdigit(c)) {
                throw std::invalid_argument("");
            }
        }
    } catch (const std::invalid_argument &) {
        errorMessage = "Time steps must be numeric.";
        errorCode = -1;
        return;
    }

    // Verificar si el primer argumento es un número entero positivo
    if (timeSteps <= 0) {
        errorMessage = "Invalid number of time steps.";
        errorCode = -2;
        return;
    }

    // Obtener los nombres de los archivos de entrada y salida
    inputFile = argv[2];
    std::ifstream file(inputFile, std::ios::binary);
    if (!file.is_open()) {
        errorMessage = "Cannot open " + inputFile + " for reading.";
        errorCode = -3;
        return;
    }
    if (file){
        file.seekg(4, std::ios::cur);
        file.read(reinterpret_cast<char *>(&intValue), sizeof(intValue));
        int np = *reinterpret_cast<int *>(&intValue);
        if (np <= 0){
            errorMessage = "Invalid number of particles: 0.";
            errorCode = -5;
            return ;
        }

    }

    outputFile = argv[3];

    errorCode = 0;
}

int ProgArgs::getTimeSteps() const {
    return timeSteps;
}

const std::string &ProgArgs::getInputFile() const {
    return inputFile;
}

const std::string &ProgArgs::getOutputFile() const {
    return outputFile;
}

int ProgArgs::getErrorCode() const {
    return errorCode;
}

const std::string &ProgArgs::getErrorMessage() const {
    return errorMessage;
}