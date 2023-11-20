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
    grid_properties();
    print_grid();
    crear_bloques();
    ParticleArray particles;
    meter_particulas(file_input, particles);

    simulacion();




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
    if (bloque_x == 0){nx_equal_0.push_back(bloque_id);}
    if (bloque_y == 0){ny_equal_0.push_back(bloque_id);}
    if (bloque_z == 0){nz_equal_0.push_back(bloque_id);}

    if (bloque_x == nx - 1){nx_less_1.push_back(bloque_id);}
    if (bloque_y == ny- 1){ny_less_1.push_back(bloque_id);}
    if (bloque_z == nz - 1){nz_less_1.push_back(bloque_id);}
}





void guardar_particulas(ParticleArray& particles, const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    std::ostringstream errorMessage;
    if (file) {

        float ppm;
        file.read(reinterpret_cast<char*>(&ppm), sizeof(ppm));

        int32_t intValue;
        file.read(reinterpret_cast<char*>(&intValue), sizeof(intValue));
        int np = *reinterpret_cast<int*>(&intValue);


        if (np <= 0){
            errorMessage << "Invalid number of particles: "<<np<<".";
            throw std::runtime_error(errorMessage.str());
        }

        std::streampos initialPos = file.tellg();



        int count = 0;
        float aux;
        while (file.read(reinterpret_cast<char*>(&aux), sizeof(aux))) {
            file.seekg(-4, std::ios::cur);
            file.seekg(36, std::ios::cur);

            count += 1;
        }

        file.clear();  // Limpia cualquier flag de error
        file.seekg(initialPos);


        if (count != np) {
            std::ostringstream errorMessage;
            errorMessage << "Number of particles mismatch. Header: " << np << ", Found: " << count;

            throw std::runtime_error(errorMessage.str());
        }





        // Redimensiona los vectores según el número total de partículas


        }

        // Cerrar el archivo después de usarlo
        file.close();



    } else {
        std::cerr << "Error opening the file." << std::endl;
    }
}