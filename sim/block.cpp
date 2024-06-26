//
// Created by manu on 8/11/23.
//

#include "block.hpp"
#include "grid.hpp"

// Calcula el identificador único de un bloque en función de sus coordenadas.
int Grid::calcular_id_bloque(int operador_x, int operador_y, int operador_z) const {
    return operador_x * static_cast<int>(nz * ny) + operador_y * static_cast<int>(nz) + operador_z;
}

// Crea los bloques del grid y establece sus bloques adyacentes.
void Grid::crear_bloques() {
    int const num_block = static_cast<int>(nx * nz * ny);
    bloques_grid.resize(num_block);

    for (int bloque_x = 0; bloque_x < nx; bloque_x++) {
        for (int bloque_y = 0; bloque_y < ny; bloque_y++) {
            for (int bloque_z = 0; bloque_z < nz; bloque_z++) {
                int const bloque_id = calcular_id_bloque(bloque_x, bloque_y, bloque_z);

                // Establece los bloques adyacentes al bloque actual.
                bloques_grid[bloque_id].adjacent_blocks = bloques_adyacentes(bloque_x, bloque_y, bloque_z);

                // Identifica y almacena bloques en los límites del grid.
                bloque_extremo(bloque_x, bloque_y, bloque_z, bloque_id);
            }
        }
    }
}

// Identifica y almacena bloques en los límites del grid.
void Grid::bloque_extremo(int bloque_x, int bloque_y, int bloque_z, int bloque_id) {
    if (bloque_x == 0) { nx_0.push_back(bloque_id); }
    if (bloque_y == 0) { ny_0.push_back(bloque_id); }
    if (bloque_z == 0) { nz_0.push_back(bloque_id); }

    if (bloque_x == nx - 1) { nx_menos1.push_back(bloque_id); }
    if (bloque_y == ny - 1) { ny_menos1.push_back(bloque_id); }
    if (bloque_z == nz - 1) { nz_menos1.push_back(bloque_id); }
}

// Devuelve una lista de bloques adyacentes a un bloque dado por sus coordenadas.
std::vector<int> Grid::bloques_adyacentes(int bloque_x, int bloque_y, int bloque_z) const {
    std::vector<int> bloques_contiguos;
    bloques_contiguos.reserve(26);

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                //if (dx == 0 && dy == 0 && dz == 0) continue;  // Excluir el bloque actual
                int const contiguo_x = bloque_x + dx;
                int const contiguo_y = bloque_y + dy;
                int const contiguo_z = bloque_z + dz;

                if (contiguo_x >= 0 && contiguo_y >= 0 && contiguo_z >= 0 &&
                    contiguo_x < nx && contiguo_y < ny && contiguo_z < nz) {
                    bloques_contiguos.push_back(calcular_id_bloque(contiguo_x, contiguo_y, contiguo_z));
                }
            }
        }
    }

    return bloques_contiguos;
}

// Asigna cada partícula a su respectivo bloque en el grid.
void Grid::particulas_bloque(ParticleArray& particles) {
    int index = 0;
    while (index < np){
        particles.i[index] = floor((particles.px[index] - bmin[0]) / sx);
        particles.j[index] = floor((particles.py[index] - bmin[1]) / sy);
        particles.k[index] = floor((particles.pz[index] - bmin[1]) / sy);

        // Ajusta las coordenadas de la partícula si se encuentran fuera de los límites del grid.
        if (particles.i[index] > nx - 1) {
            particles.i[index] = static_cast<int>(nx) - 1;
        }

        if (particles.j[index] > ny - 1) {
            particles.j[index] = static_cast<int>(ny) - 1;
        }

        if (particles.k[index] > nz - 1) {
            particles.k[index] = static_cast<int>(nz) - 1;
        }

        if (particles.i[index] < 0) {
            particles.i[index] = 0;
        }
        if (particles.j[index] < 0) {
            particles.j[index] = 0;
        }
        if (particles.k[index] < 0) {
            particles.k[index] = 0;
        }

        // Calcula el identificador del bloque al que pertenece la partícula y lo almacena en el bloque correspondiente.
        int const bloque_particula = calcular_id_bloque(particles.i[index], particles.j[index], particles.k[index]);
        bloques_grid[bloque_particula].index_particle_block.push_back(index);

        index += 1;
    }
}