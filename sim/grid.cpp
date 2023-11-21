//
// Created by manu on 8/11/23.
//

#include "grid.hpp"
#include <stdexcept>
#include <sstream>
#include <typeinfo>



template <typename T>
T read_binary_value(std::ifstream& infile) {
    T value;
    infile.read(reinterpret_cast<char*>(&value), sizeof(T));
    return value;
}



Grid::Grid(int time, std::string const& file_input, std::string const& file_output)
        : timeSteps(time), ppm(0), np(0), m(0), h(0), nx(0), ny(0), nz(0), sx(0), sy(0), sz(0){
    get_parameters(file_input);
    constantes();
    grid_properties();
    print_grid();
    crear_bloques();
    ParticleArray particles;
    meter_particulas(file_input, particles);

    simulacion(particles);




        }


void Grid::meter_particulas(const std::string& filename, ParticleArray& particles) {
    std::ifstream file(filename, std::ios::binary);
    file.seekg(8, std::ios::cur);
    particles.px.resize(np);
    particles.py.resize(np);
    particles.pz.resize(np);
    particles.hvx.resize(np);
    particles.hvy.resize(np);
    particles.hvz.resize(np);
    particles.vx.resize(np);
    particles.vy.resize(np);
    particles.vz.resize(np);
    particles.i.resize(np);
    particles.j.resize(np);
    particles.k.resize(np);

    int count = 0;

    while (count < np) {
        particles.px[count] = static_cast<double>(read_binary_value<float>(file));
        particles.i[count] = (particles.px[count] - bmin[0]) / sx;

        particles.py[count] = static_cast<double>(read_binary_value<float>(file));
        particles.j[count] = (particles.py[count] - bmin[1]) / sy;
        particles.pz[count] = static_cast<double>(read_binary_value<float>(file));
        particles.k[count] = (particles.pz[count] - bmin[2]) / sz;

        particles.hvx[count] = static_cast<double>(read_binary_value<float>(file));
        particles.hvy[count] = static_cast<double>(read_binary_value<float>(file));
        particles.hvz[count] = static_cast<double>(read_binary_value<float>(file));
        particles.vx[count] = static_cast<double>(read_binary_value<float>(file));
        particles.vy[count] = static_cast<double>(read_binary_value<float>(file));
        particles.vz[count] = static_cast<double>(read_binary_value<float>(file));

        count += 1;
    }
    particle_block(particles);
}





void Grid::particle_block(ParticleArray& particles) {
    int index = 0;
    while (index < np){
        particles.i[index] = floor((particles.px[index] - bmin[0]) / sx);
        particles.j[index] = floor((particles.py[index] - bmin[1]) / sy);
        particles.k[index] = floor((particles.pz[index] - bmin[1]) / sy);

        if (particles.i[index] > nx - 1) {
            particles.i[index] = nx - 1;
        }

        if (particles.j[index] > ny - 1) {
            particles.j[index] = ny - 1;
        }

        if (particles.k[index] > nz - 1) {
            particles.k[index] = nz - 1;
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

        int const bloque_particula = operation_block(particles.i[index], particles.j[index], particles.k[index]);
        grid_block[bloque_particula].index_particle_block.push_back(index);

        index += 1;

    }
}

void Grid::get_parameters(const std::string& filename){
    std::ifstream file(filename, std::ios::binary);
    std::ostringstream errorMessage;
    if (file){
        ppm = static_cast<double>(read_binary_value<float>(file));
        np = static_cast<double>(read_binary_value<float>(file));
    }
}

void Grid::constantes() {
    calculos.h_cuadrado = pow(h, 2);
    calculos.h_sexta = pow(h, NUMBER_6);
    calculos.trans_densidad = (NUMBER_315/(NUMBER_64 * M_PI * pow(h,NUMBER_9))) * m;


    double const aux_calc = 1/(M_PI * calculos.h_sexta);
    calculos.calc1 = NUMBER_15 * 3 * m * presion_rigidez * aux_calc * NUMBER_05;
    calculos.calc2 = 2 * densidad_fluido;
    calculos.calc3 = NUMBER_45 * viscosidad * m * aux_calc;
}

void Grid::grid_properties(){
    m = densidad_fluido/std::pow(ppm, 3); // Masa de una particula
    h = multiplicador_radio /ppm; // Longitud de suavizado

    // Numero de bloques de cada dimension
    nx = std::floor((bmax[0] - bmin[0])/h);
    ny = std::floor((bmax[1] - bmin[1])/h);
    nz = std::floor((bmax[2] - bmin[2])/h);

    // Tamaño de los bloques
    sx = (bmax[0] - bmin[0])/nx;
    sy = (bmax[1] - bmin[1])/ny;
    sz = (bmax[2] - bmin[2])/nz;
}

void Grid::print_grid() {
    std::cout << "Numero de particulas: " << np << std::endl;
    std::cout << "Particulas por metro: " << ppm << std::endl;
    std::cout << "Smoothing length: " << h << std::endl;
    std::cout<< "Masa particula: " << m <<std::endl;
    std::cout<< "Tamaño del grid: " << nx << " x " << ny << " x " << nz << std::endl;
    std::cout<< "Numero de bloques: " << nx*ny*nz << std::endl;
    std::cout << "Tamaño de bloque: " << sx << " x " << sy << " x " << sz << std::endl;
}

int Grid::operation_block(int operador_x, int operador_y, int operador_z) const  {
    return operador_x * static_cast<int>(nz*ny) + operador_y * static_cast<int>(nz) + operador_z;
}

void Grid::crear_bloques() {
    int const num_block = static_cast<int>(nx * nz * ny);
    grid_block.resize(num_block);
    for (int bloque_x = 0; bloque_x < nx; bloque_x++) {
        for (int bloque_y = 0; bloque_y < ny; bloque_y++) {
            for (int bloque_z = 0; bloque_z < nz; bloque_z++) {
                int const bloque_id = operation_block(bloque_x, bloque_y, bloque_z);
                //En las contiguos tambien se incluye a si mismo
                grid_block[bloque_id].adjacent_blocks = calc_adjacent_blocks(bloque_x, bloque_y, bloque_z);
                block_type_check(bloque_x, bloque_y, bloque_z, bloque_id);
            }
        }
    }
}

std::vector<int> Grid::calc_adjacent_blocks(int bloque_x, int bloque_y, int bloque_z) const {
    std::vector<int> bloques_contiguos;
    bloques_contiguos.reserve(26);

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                if (dx == 0 && dy == 0 && dz == 0) continue;  // Excluir el bloque actual
                int const contiguo_x = bloque_x + dx;
                int const contiguo_y = bloque_y + dy;
                int const contiguo_z = bloque_z + dz;

                if (contiguo_x >= 0 && contiguo_y >= 0 && contiguo_z >= 0 &&
                    contiguo_x < nx && contiguo_y < ny && contiguo_z < nz) {
                    bloques_contiguos.push_back(operation_block(contiguo_x, contiguo_y, contiguo_z));
                }
            }
        }
    }

    return bloques_contiguos;
}

void Grid::block_type_check(int bloque_x, int bloque_y, int bloque_z, int bloque_id){
    if (bloque_x == 0){nx_0.push_back(bloque_id);}
    if (bloque_y == 0){ny_0.push_back(bloque_id);}
    if (bloque_z == 0){nz_0.push_back(bloque_id);}

    if (bloque_x == nx - 1){nx_menos1.push_back(bloque_id);}
    if (bloque_y == ny- 1){ny_menos1.push_back(bloque_id);}
    if (bloque_z == nz - 1){nz_menos1.push_back(bloque_id);}
}

void Grid::simulacion(ParticleArray& particles){
    //Do while para evitar recalcular la posicion de la particula en el bloque
    //en la primera iteracion
    inicializacion_aceleracion_densidad(particles);
    actualizar_ac_den(particles);
    colisiones(particles);
    // particle_move();
    // iterations_limits_recient();

    for(int time = 1; time < timeSteps; time++){
        inicializacion_aceleracion_densidad(particles);
        // recalculate_particle_position_in_block();
        // act_aceleration_density();
        // colisions();
        // particle_move();
        // iterations_limits_recient();
    }
}

//Transformacion de la densidad y aceleracion
void Grid::inicializacion_aceleracion_densidad(ParticleArray& particles){
    for (int i = 0; i<np; i++){
        particles.ax[i] = gravedad[0];
        particles.ay[i] = gravedad[1];
        particles.az[i] = gravedad[2];
        particles.density[i] = 0;
    }
}

void Grid::actualizar_ac_den(ParticleArray& particles){
    for (int i=0; i<np; i++){
        int const bloque = operation_block(particles.i[i], particles.j[i], particles.k[i]);
        for (int const bloque_adyacente : grid_block[bloque].adjacent_blocks){
            for (int j : grid_block[bloque_adyacente].index_particle_block){
                if (i < j){
                    incremento_densidad(i, j, particles);
                }
            }
        }
        particles.density[i] = transformacion_densidad(i, particles);
    }

    for (int i=0; i<np; i++){
        int const bloque = operation_block(particles.i[i], particles.j[i], particles.k[i]);
        for (int const bloque_adyacente : grid_block[bloque].adjacent_blocks){
            for (int j : grid_block[bloque_adyacente].index_particle_block){
                if (i < j){
                    actualizar_aceleracion(i, j, particles);
                }
            }
        }
    }
}

void Grid::incremento_densidad(int i, int j, ParticleArray& particles){
    double modulo_cuadrado = calcular_modulo(i, j, particles);
    if (modulo_cuadrado < pow(h, 2)){
        double const incremento = pow((pow(h, 2) - modulo_cuadrado), 3);
        particles.density[i] += incremento;
        particles.density[j] += incremento;
    }
}
double Grid::calcular_modulo(int i, int j, ParticleArray& particles){
    //No hacemos la raiz cuadrada porque al hacer el cuadrado del modulo, la raiz se va con el cuadrado
    //es una operacion que nos ahorramos
    double cuadrado_x = pow((particles.px[i] - particles.px[j]), 2);
    double cuadrado_y = pow((particles.py[i] - particles.py[j]), 2);
    double cuadrado_z = pow((particles.pz[i] - particles.pz[j]), 2);

    double const modulo_cuadrado = cuadrado_x + cuadrado_y + cuadrado_z;
    return modulo_cuadrado;
}

double Grid::transformacion_densidad(int i, ParticleArray& particles){
    particles.density[i] += (particles.density[i] + pow(h, 6)) * calculos.trans_densidad;
}

void Grid::actualizar_aceleracion(int i, int j, ParticleArray& particles){
    double modulo_cuadrado = calcular_modulo(i, j, particles);
    if (modulo_cuadrado < pow(h, 2)){
        double const dij = sqrt(std::max(modulo_cuadrado, 1e-12));
        double delta_aij_x;
        double delta_aij_y;
        double delta_aij_z;

        double const parte1 = calculos.calc1 * pow((h-dij),2) * (particles.density[i] + particles.density[j] - calculos.calc2) * (1/dij);
        double auxiliar = 1 / (particles.density[i] * particles.density[j]);

        delta_aij_x = ((particles.px[i]-particles.px[j]) * parte1 + (particles.vx[j]-particles.vx[i]) * calculos.calc3) * auxiliar;
        delta_aij_y = ((particles.py[i]-particles.py[j]) * parte1 + (particles.vy[j]-particles.vy[i]) * calculos.calc3) * auxiliar;
        delta_aij_z = ((particles.pz[i]-particles.pz[j]) * parte1 + (particles.vz[j]-particles.vz[i]) * calculos.calc3) * auxiliar;

        particles.ax[i] += delta_aij_x;
        particles.ay[i] += delta_aij_y;
        particles.az[i] += delta_aij_z;
        particles.ax[j] -= delta_aij_x;
        particles.ay[j] -= delta_aij_y;
        particles.az[j] -= delta_aij_z;
    }
}

// Colisiones
void Grid::colisiones(ParticleArray& particles){
    bucle_aux_block_0(nx_0, 0);
    bucle_aux_block_0(ny_0, 1);
    bucle_aux_block_0(nz_0, 2);

    bucle_aux_block_n_less_1(nx_menos1, 0);
    bucle_aux_block_n_less_1(ny_menos1, 1);
    bucle_aux_block_n_less_1(nz_menos1, 2);

}


