#ifndef PROYECTO1_ARQUITECTURA_GRID_HPP
#define PROYECTO1_ARQUITECTURA_GRID_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <numbers>

// Constantes
const float multiplicador_radio = 1.695;
const float densidad_fluido = 1000;
const float presion_rigidez = 3.0;
const float colision_rigidez = 30000;
const float amortiguamiento = 128.0;
const float viscosidad = 0.4;
const float tamaño_particula = 0.0002;
const float paso_tiempo = 0.001;

// Constantes vectoriales
const std::vector<float> g = {0.0, -9.8, 0.0};
const std::vector<float> bmax = {0.065, 0.1, 0.065};
const std::vector<float> bmin = {-0.065, -0.08, -0.065};

// Estructura de las partículas
struct ParticleArray {
    std::vector<float> px, py, pz;   // Cambiados a double
    std::vector<float> hvx, hvy, hvz; // Cambiados a double
    std::vector<float> vx, vy, vz;   // Cambiados a double
};

// Prototipo de la función para guardar partículas en el grid
void guardar_particulas(ParticleArray& particles, const std::string& filename);

#endif //PROYECTO1_ARQUITECTURA_GRID_HPP