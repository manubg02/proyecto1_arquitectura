//
// Created by manu on 8/11/23.
//
#include "progargs.hpp"
#include <iostream>
#include <cctype>
#include <fstream>

ProgArgs::ProgArgs(int argc, char *argv[]) : timeSteps(0), errorCode(0) {
    // Inicia el proceso de análisis de argumentos al construir una instancia de ProgArgs
    parseArguments(argc, argv);
}

template <typename T>
T leer_binario(std::ifstream& infile) {
    // Lee un valor binario de un archivo y lo devuelve
    T value;
    infile.read(reinterpret_cast<char*>(&value), sizeof(T));
    return value;
}

void ProgArgs::parseArguments(int argc, char *argv[]) {
    // Verifica el número correcto de argumentos
    if (argc != 4) {
        errorMessage = "Invalid number of arguments.";
        errorCode = -1;
        return;
    }

    // Verifica si el primer carácter es un dígito o el signo negativo
    if (!std::isdigit(argv[1][0]) && argv[1][0] != '-') {
        errorMessage = "Time steps must be numeric.";
        errorCode = -1;
        return;
    }

    // Verifica si el resto de la cadena contiene solo dígitos
    for (size_t i = 1; i < std::string(argv[1]).size(); ++i) {
        if (!std::isdigit(argv[1][i])) {
            errorMessage = "Time steps must be numeric.";
            errorCode = -1;
            return;
        }
    }

    // Obtiene el número de pasos de tiempo
    try {
        timeSteps = std::stoi(argv[1]);
    } catch (const std::invalid_argument &) {
        errorMessage = "Time steps must be numeric.";
        errorCode = -1;
        return;
    }

    // Verifica si el número de pasos de tiempo es un entero positivo
    if (timeSteps <= 0) {
        errorMessage = "Invalid number of time steps.";
        errorCode = -2;
        return;
    }

    // Obtiene los nombres de los archivos de entrada y salida
    inputFile = argv[2];
    std::ifstream file(inputFile, std::ios::binary);

    // Verifica si el archivo de entrada se puede abrir para lectura
    if (!file.is_open()) {
        errorMessage = "Cannot open " + inputFile + " for reading.";
        errorCode = -3;
        return;
    }

    // Salta la información del encabezado del archivo de entrada
    file.seekg(4, std::ios::cur);
    int np = static_cast<int>(leer_binario<int>(file));

    // Verifica si el número de partículas es válido
    if (np <= 0) {
        errorMessage = "Invalid number of particles: 0.";
        errorCode = -5;
        return;
    }

    int count = 0;
    float aux;

    // Lee las partículas del archivo para contarlas
    while (file.read(reinterpret_cast<char*>(&aux), sizeof(aux))) {
        file.seekg(-4, std::ios::cur);
        file.seekg(36, std::ios::cur);
        count += 1;
    }

    // Limpia cualquier flag de error en el archivo
    file.clear();

    // Verifica si el número de partículas coincide con la información del encabezado
    if (count != np) {
        errorMessage = "Number of particles mismatch. Header: " + std::to_string(np) + ", Found: " + std::to_string(count);
        errorCode = -5;
        return;
    }

    // Obtiene el nombre del archivo de salida
    outputFile = argv[3];

    // No hay errores en los argumentos
    errorCode = 0;
}

// Obtiene el número de iteraciones
int ProgArgs::getTimeSteps() const {
    return timeSteps;
}

// Obtiene el fichero de entrada
const std::string &ProgArgs::getInputFile() const {
    return inputFile;
}

// Obtiene el fichero de salida
const std::string &ProgArgs::getOutputFile() const {
    return outputFile;
}

// Obtiene el código de error
int ProgArgs::getErrorCode() const {
    return errorCode;
}

// Obtiene el mensaje de error
const std::string &ProgArgs::getErrorMessage() const {
    return errorMessage;
}