#ifndef PROYECTO1_ARQUITECTURA_GRID_HPP
#define PROYECTO1_ARQUITECTURA_GRID_HPP

#include <iostream>
#include <fstream>
#include <vector> // Incluye la librería vector para usar std::vector

// Estructura de las partículas
struct ParticleArray {
    std::vector<float> px, py, pz;
    std::vector<float> hvx, hvy, hvz;
    std::vector<float> vx, vy, vz;
};

// Prototipo de la función para guardar partículas en el grid
void guardar_particulas(ParticleArray& particles, const std::string& filename);

#endif //PROYECTO1_ARQUITECTURA_GRID_HPP