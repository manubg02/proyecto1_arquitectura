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
#include "block.hpp"

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
const double t_particula = 0.0002;
const double paso_tiempo = 0.001;
const double comparar_colision = 1e-10;

// Constantes vectoriales
const std::vector<double> gravedad = {0.0, -9.8, 0.0};
const std::vector<double> bmax = {0.065, 0.1, 0.065};
const std::vector<double> bmin = {-0.065, -0.08, -0.065};

struct calc{
    double h_cuadrado;
    double trans_densidad;
    double h_sexta;
    double calc1;
    double calc2;
    double calc3;
};

// Estructura de las partículas
struct ParticleArray {
    std::vector<double> px, py, pz;   // Posiciones de cada particula
    std::vector<double> hvx, hvy, hvz; // Gradiente de velocidad de cada particula
    std::vector<double> vx, vy, vz;   // Velocidades de cada particula
    std::vector<double> density; // Densidades de cada particula
    std::vector<double> ax, ay, az; // Aceleraciones de cada particula
    std::vector<int> i, j, k; //bloque  i[0],  i[1], ...
};



class Grid{
public:
    Grid(int timeSteps, std::string const& file_input, std::string const& file_output);
    void get_parameters(const std::string& file_input);
    void constantes();
    void grid_properties();
    void print_grid();
    void meter_particulas(const std::string& file_input, ParticleArray& particles);
    void particulas_bloque(ParticleArray& particles);
    void crear_bloques();
    [[nodiscard]] std::vector<int> bloques_adyacentes(int bloque_x, int bloque_y, int bloque_z ) const;
    [[nodiscard]] int calcular_id_bloque(int operador_x, int operador_y, int operador_z) const;
    void bloque_extremo(int bloque_x, int bloque_y, int bloque_z, int bloque_id);

    void simulacion(ParticleArray& particles);
    void reposicionamiento_particulas(ParticleArray& particles);
    void inicializacion_aceleracion_densidad(ParticleArray& particles);
    void actualizar_ac_den(ParticleArray& particles);
    void incremento_densidad(int i, int j, ParticleArray& particles);
    static double calcular_modulo_cuadrado(int i, int j, ParticleArray& particles);
    void transformacion_densidad(int i, ParticleArray& particles);
    void actualizar_aceleracion(int i, int j, ParticleArray& particles);

    void colisiones(ParticleArray& particles);
    void bucle_bloque_x0(const std::vector<int>& block_list, ParticleArray& particles);
    void bucle_bloque_y0(const std::vector<int>& block_list, ParticleArray& particles);
    void bucle_bloque_z0(const std::vector<int>& block_list, ParticleArray& particles);
    void bucle_bloque_xmenos1(const std::vector<int> &block_list, ParticleArray& particles);
    void bucle_bloque_ymenos1(const std::vector<int> &block_list, ParticleArray& particles);
    void bucle_bloque_zmenos1(const std::vector<int> &block_list, ParticleArray& particles);
    static void colisiones_particulas_eje_x(int id, ParticleArray& particles);
    static void colisiones_particulas_eje_y(int id, ParticleArray& particles);
    static void colisiones_particulas_eje_z(int id, ParticleArray& particles);
    static void colisiones_particulas_eje_xmenos1(int id, ParticleArray& particles);
    static void colisiones_particulas_eje_ymenos1(int id, ParticleArray& particles);
    static void colisiones_particulas_eje_zmenos1(int id, ParticleArray& particles);

    void movimiento_particulas(ParticleArray& particles);
    static void act_posicion(int i, ParticleArray& particles);
    static void act_velocidad(int i, ParticleArray& particles);
    static void act_gradiente(int i, ParticleArray& particles);

    void limites_recinto(ParticleArray& particles);
    void limite_bloque_x0(const std::vector<int> &block_list, ParticleArray& particles);
    void limite_bloque_y0(const std::vector<int> &block_list, ParticleArray& particles);
    void limite_bloque_z0(const std::vector<int> &block_list, ParticleArray& particles);
    void limite_bloque_xmenos1(const std::vector<int> &block_list, ParticleArray& particles);
    void limite_bloque_ymenos1(const std::vector<int> &block_list, ParticleArray& particles);
    void limite_bloque_zmenos1(const std::vector<int> &block_list, ParticleArray& particles);
    static void limite_particulas_eje_x(int id, ParticleArray& particles);
    static void limite_particulas_eje_y(int id, ParticleArray& particles);
    static void limite_particulas_eje_z(int id, ParticleArray& particles);
    static void limite_particulas_eje_xmenos1(int id, ParticleArray& particles);
    static void limite_particulas_eje_ymenos1(int id, ParticleArray& particles);
    static void limite_particulas_eje_zmenos1(int id, ParticleArray& particles);

    void escribir_informacion(ParticleArray& particles, std::string const& file_output);


private:
    int timeSteps;
    double ppm;
    int np;
    double m; //masa de particulas
    double h; //longitud de suavizado
    double nx, ny, nz; //numero de bloques por coordenadas
    double sx, sy, sz; //size del bloque
    std::vector<bloque> grid_block;
    std::vector<int> nx_0, ny_0, nz_0, nx_menos1, ny_menos1, nz_menos1;

    calc calculos{};
};


#endif //PROYECTO1_ARQUITECTURA_GRID_HPP