#ifndef PROYECTO1_ARQUITECTURA_GRID_HPP
#define PROYECTO1_ARQUITECTURA_GRID_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <numbers>

// Constantes
const double multiplicador_radio = 1.695;
const double densidad_fluido = 1000;
const double presion_rigidez = 3.0;
const double colision_rigidez = 30000;
const double amortiguamiento = 128.0;
const double viscosidad = 0.4;
const double tamaño_particula = 0.0002;
const double paso_tiempo = 0.001;

// Constantes vectoriales
const std::vector<double> g = {0.0, -9.8, 0.0};
const std::vector<double> bmax = {0.065, 0.1, 0.065};
const std::vector<double> bmin = {-0.065, -0.08, -0.065};

// Estructura de las partículas
struct ParticleArray {
    std::vector<double> px, py, pz;   // Cambiados a double
    std::vector<double> hvx, hvy, hvz; // Cambiados a double
    std::vector<double> vx, vy, vz;   // Cambiados a double
    std::vector<int> i, j, k; //bloque
};

// Prototipo de la función para guardar partículas en el grid
void guardar_particulas(ParticleArray& particles, const std::string& filename);

#endif //PROYECTO1_ARQUITECTURA_GRID_HPP