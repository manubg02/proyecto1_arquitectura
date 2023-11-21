#ifndef PROYECTO1_ARQUITECTURA_GRID_HPP
#define PROYECTO1_ARQUITECTURA_GRID_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <numbers>
#include <string>
#include <typeinfo>
#include <cmath>
#include <sstream>

//Declaracion general para los magic number
#define NUMBER_315 (315)
#define NUMBER_64 (64)
#define NUMBER_9 (9)
#define NUMBER_15 (15)
#define NUMBER_05 (0.5)
#define NUMBER_6 (6)
#define NUMBER_45 (45)

// Constantes
const double multiplicador_radio = 1.695;
const double densidad_fluido = 1000;
const double presion_rigidez = 3;
const double colision_rigidez = 30000;
const double amortiguamiento = 128;
const double viscosidad = 0.4;
const double tamaño_particula = 0.0002;
const double paso_tiempo = 0.001;

// Constantes vectoriales
const std::vector<double> gravedad = {0.0, -9.8, 0.0};
const std::vector<double> bmax = {0.065, 0.1, 0.065};
const std::vector<double> bmin = {-0.065, -0.08, -0.065};

// Estructura de las partículas
struct ParticleArray {
    std::vector<double> px, py, pz;   // Posiciones de cada particula
    std::vector<double> hvx, hvy, hvz; // Gradiente de velocidad de cada particula
    std::vector<double> vx, vy, vz;   // Velocidades de cada particula
    std::vector<double> density; // Densidades de cada particula
    std::vector<double> ax, ay, az; // Aceleraciones de cada particula
    std::vector<int> i, j, k; //bloque  i[0],  i[1], ...
};

struct block{
    std::vector<int> adjacent_blocks; //Se incluye a el mismo
    std::vector<int> index_particle_block;
};


class Grid{
public:
    Grid(int timeSteps, std::string const& file_input, std::string const& file_output);
    void get_parameters(const std::string& file_input);
    void grid_properties();
    void print_grid();
    void meter_particulas(const std::string& file_input, ParticleArray& particles);
    void particle_block(ParticleArray& particles);
    void crear_bloques();
    [[nodiscard]] std::vector<int> calc_adjacent_blocks(int bloque_x, int bloque_y, int bloque_z ) const;
    [[nodiscard]] int operation_block(int operador_x, int operador_y, int operador_z) const;
    void block_type_check(int bloque_x, int bloque_y, int bloque_z, int bloque_id);

    void simulacion(ParticleArray& particles);
    void inicializacion_aceleracion_densidad(ParticleArray& particles);
    void actualizar_ac_den(ParticleArray& particles);
    void incremento_densidad(int i, int j, ParticleArray& particles);
    static double calcular_modulo(int i, int j, ParticleArray& particles);
    [[nodiscard]] double transformacion_densidad(int i, ParticleArray& particles);



    void parameter_values() const;
    void constant_values_calculations();

    void block_particle_assignment(ParticleArray& particles, int index) const;

private:
    int timeSteps;
    double ppm;
    int np;
    double m; //masa de particulas
    double h; //longitud de suavizado
    double nx, ny, nz; //numero de bloques por coordenadas
    double sx, sy, sz; //size del bloque
    std::vector<block> grid_block;
    std::vector<int> nx_equal_0, ny_equal_0, nz_equal_0, nx_less_1, ny_less_1, nz_less_1;
};


// Prototipo de la función para guardar partículas en el grid
void guardar_particulas(ParticleArray& particles, const std::string& filename);

#endif //PROYECTO1_ARQUITECTURA_GRID_HPP