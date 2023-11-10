//
// Created by manu on 8/11/23.
//

#ifndef PROYECTO1_ARQUITECTURA_PROGARGS_HPP
#define PROYECTO1_ARQUITECTURA_PROGARGS_HPP

#include <string>

class ProgArgs {
public:
    ProgArgs(int argc, char *argv[]);

    int getTimeSteps() const;
    const std::string &getInputFile() const;
    const std::string &getOutputFile() const;
    int getErrorCode() const;
    const std::string &getErrorMessage() const;

private:
    void parseArguments(int argc, char *argv[]);

    int timeSteps;
    std::string inputFile;
    std::string outputFile;
    int errorCode;
    std::string errorMessage;
};

#endif //PROYECTO1_ARQUITECTURA_PROGARGS_HPP
