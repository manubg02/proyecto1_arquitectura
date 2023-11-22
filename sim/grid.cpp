//
// Created by manu on 8/11/23.
//

#include "grid.hpp"
#include <stdexcept>
#include <sstream>

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
    escribir_informacion(particles, file_output);
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
    particles.density.resize(np);
    particles.ax.resize(np);
    particles.ay.resize(np);
    particles.az.resize(np);

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
        np = static_cast<int>(read_binary_value<int>(file));
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
    movimiento_particulas(particles);
    limites_recinto(particles);

    for(int time = 1; time < timeSteps; time++){
        inicializacion_aceleracion_densidad(particles);
        reposicionamiento_particulas(particles);
        actualizar_ac_den(particles);
        colisiones(particles);
        movimiento_particulas(particles);
        limites_recinto(particles);
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

void Grid::reposicionamiento_particulas(ParticleArray& particles){
    for (block& block: grid_block){
        block.index_particle_block.clear();
    }
    particle_block(particles);
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
        transformacion_densidad(i, particles);
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
    if (modulo_cuadrado < calculos.h_cuadrado){
        double const incremento = pow((calculos.h_cuadrado - modulo_cuadrado), 3);
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

void Grid::transformacion_densidad(int i, ParticleArray& particles){
    particles.density[i] += (particles.density[i] + pow(h, 6)) * calculos.trans_densidad;
}

void Grid::actualizar_aceleracion(int i, int j, ParticleArray& particles){
    double modulo_cuadrado = calcular_modulo(i, j, particles);
    if (modulo_cuadrado < calculos.h_cuadrado){
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
    bucle_bloque_x0(nx_0, particles);
    bucle_bloque_y0(ny_0, particles);
    bucle_bloque_z0(nz_0, particles);

    bucle_bloque_xmenos1(nx_menos1, particles);
    bucle_bloque_ymenos1(ny_menos1, particles);
    bucle_bloque_zmenos1(nz_menos1, particles);

}

void Grid::bucle_bloque_x0(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        block const block = grid_block[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_x(particle_id, particles);
        }
    }
}

void Grid::bucle_bloque_y0(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        block const block = grid_block[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_y(particle_id, particles);
        }
    }
}

void Grid::bucle_bloque_z0(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        block const block = grid_block[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_z(particle_id, particles);
        }
    }
}

void Grid::colisiones_particulas_eje_x(int id, ParticleArray& particles){
    double const limite_eje_x = particles.px[id] + particles.hvx[id] * paso_tiempo;
    double const incremento_x = t_particula - (limite_eje_x - bmin[0]);
    if (incremento_x > comparar_colision){
        double const aux = colision_rigidez * incremento_x - amortiguamiento * particles.vx[id];
        particles.ax[id] += aux;
    }
}

void Grid::colisiones_particulas_eje_y(int id, ParticleArray& particles){
    double const limite_eje_y = particles.py[id] + particles.hvy[id] * paso_tiempo;
    double const incremento_y = t_particula - (limite_eje_y - bmin[1]);
    if (incremento_y > comparar_colision){
        double const aux = colision_rigidez * incremento_y - amortiguamiento * particles.vy[id];
        particles.ay[id] += aux;
    }
}

void Grid::colisiones_particulas_eje_z(int id, ParticleArray& particles){
    double const limite_eje_z = particles.pz[id] + particles.hvz[id] * paso_tiempo;
    double const incremento_z = t_particula - (limite_eje_z - bmin[2]);
    if (incremento_z > comparar_colision){
        double const aux = colision_rigidez * incremento_z - amortiguamiento * particles.vz[id];
        particles.az[id] += aux;
    }
}

void Grid::bucle_bloque_xmenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        block const block = grid_block[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_xmenos1(particle_id, particles);
        }
    }
}
void Grid::bucle_bloque_ymenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        block const block = grid_block[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_ymenos1(particle_id, particles);
        }
    }
}
void Grid::bucle_bloque_zmenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        block const block = grid_block[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_zmenos1(particle_id, particles);
        }
    }
}

void Grid::colisiones_particulas_eje_xmenos1(int id, ParticleArray& particles){
    double const limite_eje_x = particles.px[id] + particles.hvx[id] * paso_tiempo;
    double const incremento_x = t_particula - (bmax[0] - limite_eje_x);
    if (incremento_x > comparar_colision){
        double const aux = colision_rigidez * incremento_x + amortiguamiento * particles.vx[id];
        particles.ax[id] -= aux;
    }
}
void Grid::colisiones_particulas_eje_ymenos1(int id, ParticleArray& particles){
    double const limite_eje_y = particles.py[id] + particles.hvy[id] * paso_tiempo;
    double const incremento_y = t_particula - (bmax[1] - limite_eje_y);
    if (incremento_y > comparar_colision){
        double const aux = colision_rigidez * incremento_y + amortiguamiento * particles.vy[id];
        particles.ay[id] -= aux;
    }
}
void Grid::colisiones_particulas_eje_zmenos1(int id, ParticleArray& particles){
    double const limite_eje_z = particles.pz[id] + particles.hvz[id] * paso_tiempo;
    double const incremento_z = t_particula - (bmax[2] - limite_eje_z);
    if (incremento_z > comparar_colision){
        double const aux = colision_rigidez * incremento_z + amortiguamiento * particles.vz[id];
        particles.az[id] -= aux;
    }
}

//Movimiento
void Grid::movimiento_particulas(ParticleArray& particles){
    for (int i = 0; i<np; i++){
        act_posicion(i, particles);
        act_velocidad(i, particles);
        act_gradiente(i, particles);
    }
}

void Grid::act_posicion(int i, ParticleArray& particles){
    particles.px[i] += particles.hvx[i] * paso_tiempo + particles.ax[i] * pow(paso_tiempo, 2);
    particles.py[i] += particles.hvy[i] * paso_tiempo + particles.ay[i] * pow(paso_tiempo, 2);
    particles.pz[i] += particles.hvz[i] * paso_tiempo + particles.az[i] * pow(paso_tiempo, 2);
}

void Grid::act_velocidad(int i, ParticleArray& particles){
    particles.vx[i] = particles.hvx[i] + (particles.ax[i] * paso_tiempo) * NUMBER_05;
    particles.vy[i] = particles.hvy[i] + (particles.ay[i] * paso_tiempo) * NUMBER_05;
    particles.vz[i] = particles.hvz[i] + (particles.az[i] * paso_tiempo) * NUMBER_05;
}

void Grid::act_gradiente(int i, ParticleArray& particles){
    particles.hvx[i] += particles.ax[i] * paso_tiempo;
    particles.hvy[i] += particles.ay[i] * paso_tiempo;
    particles.hvz[i] += particles.az[i] * paso_tiempo;

}

// Limites Recinto

void Grid::limites_recinto(ParticleArray& particles){
    limite_bloque_x0(nx_0, particles);
    limite_bloque_y0(ny_0, particles);
    limite_bloque_z0(nz_0, particles);

    limite_bloque_xmenos1(nx_menos1, particles);
    limite_bloque_ymenos1(ny_menos1, particles);
    limite_bloque_zmenos1(nz_menos1, particles);
}

void Grid::limite_bloque_x0(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        block const block = grid_block[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_x(particle_id, particles);
        }
    }
}

void Grid::limite_bloque_y0(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        block const block = grid_block[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_y(particle_id, particles);
        }
    }
}

void Grid::limite_bloque_z0(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        block const block = grid_block[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_z(particle_id, particles);
        }
    }
}

void Grid::limite_particulas_eje_x(int id, ParticleArray &particles) {
    double const d_x = particles.px[id] - bmin[0];
    if (d_x < 0){
        particles.px[id] = bmin[0] - d_x;
        particles.vx[id] = -particles.vx[id];
        particles.hvx[id] = -particles.hvx[id];
    }
}
void Grid::limite_particulas_eje_y(int id, ParticleArray &particles) {
    double const d_y = particles.py[id] - bmin[1];
    if (d_y < 0){
        particles.py[id] = bmin[1] - d_y;
        particles.vy[id] = -particles.vy[id];
        particles.hvy[id] = -particles.hvy[id];
    }
}
void Grid::limite_particulas_eje_z(int id, ParticleArray &particles) {
    double const d_z = particles.pz[id] - bmin[2];
    if (d_z < 0){
        particles.pz[id] = bmin[2] - d_z;
        particles.vz[id] = -particles.vz[id];
        particles.hvz[id] = -particles.hvz[id];
    }
}

void Grid::limite_bloque_xmenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        block const block = grid_block[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_xmenos1(particle_id, particles);
        }
    }
}
void Grid::limite_bloque_ymenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        block const block = grid_block[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_ymenos1(particle_id, particles);
        }
    }
}
void Grid::limite_bloque_zmenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        block const block = grid_block[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_zmenos1(particle_id, particles);
        }
    }
}

void Grid::limite_particulas_eje_xmenos1(int id, ParticleArray& particles){
    double const d_x = bmax[0] - particles.px[id];
    if (d_x<0){
        particles.px[id] = bmax[0] + d_x;
        particles.vx[id] = -particles.vx[id];
        particles.hvx[id] = -particles.hvx[id];
    }
}
void Grid::limite_particulas_eje_ymenos1(int id, ParticleArray& particles){
    double const d_y = bmax[1] - particles.py[id];
    if (d_y<0){
        particles.py[id] = bmax[1] + d_y;
        particles.vy[id] = -particles.vy[id];
        particles.hvy[id] = -particles.hvy[id];
    }
}
void Grid::limite_particulas_eje_zmenos1(int id, ParticleArray& particles){
    double const d_z= bmax[2] - particles.pz[id];
    if (d_z<0){
        particles.pz[id] = bmax[2] + d_z;
        particles.vz[id] = -particles.vz[id];
        particles.hvz[id] = -particles.hvz[id];
    }
}

void Grid::escribir_informacion(ParticleArray &particles, const std::string &file_output) {
    std::ofstream archivo(file_output, std::ios::binary);

    if (!archivo.is_open()) {
        std::cerr << "Error: No se pudo abrir el archivo para escritura." << std::endl;
        return;
    }

    float ppm_ = static_cast<float>(ppm);
    archivo.write(reinterpret_cast<const char*>(&ppm_), sizeof(float));
    archivo.write(reinterpret_cast<const char*>(&np), sizeof(int));

    for (int i=0; i<np; i++){
        float px_float = static_cast<float>(particles.px[i]);
        float py_float = static_cast<float>(particles.py[i]);
        float pz_float = static_cast<float>(particles.pz[i]);

        float hvx_float = static_cast<float>(particles.hvx[i]);
        float hvy_float = static_cast<float>(particles.hvy[i]);
        float hvz_float = static_cast<float>(particles.hvz[i]);

        float vx_float = static_cast<float>(particles.vx[i]);
        float vy_float = static_cast<float>(particles.vy[i]);
        float vz_float = static_cast<float>(particles.vz[i]);
        archivo.write(reinterpret_cast<const char*>(&px_float), sizeof(float));
        archivo.write(reinterpret_cast<const char*>(&py_float), sizeof(float));
        archivo.write(reinterpret_cast<const char*>(&pz_float), sizeof(float));
        archivo.write(reinterpret_cast<const char*>(&hvx_float), sizeof(float));
        archivo.write(reinterpret_cast<const char*>(&hvy_float), sizeof(float));
        archivo.write(reinterpret_cast<const char*>(&hvz_float), sizeof(float));
        archivo.write(reinterpret_cast<const char*>(&vx_float), sizeof(float));
        archivo.write(reinterpret_cast<const char*>(&vy_float), sizeof(float));
        archivo.write(reinterpret_cast<const char*>(&vz_float), sizeof(float));
    }

    archivo.close();
}