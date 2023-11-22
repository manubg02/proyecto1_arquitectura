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
    grid_properties();
    constantes();
    print_grid();
    crear_bloques();
    ParticleArray particles;
    meter_particulas(file_input, particles);
    simulacion(particles);
    escribir_informacion(particles, file_output);
        }


void Grid::meter_particulas(const std::string& filename, ParticleArray& particles) {
    std::ifstream file(filename, std::ios::binary);
    if (file){
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
            particles.i[count] = static_cast<int>((particles.px[count] - bmin[0]) / sx);

            particles.py[count] = static_cast<double>(read_binary_value<float>(file));
            particles.j[count] = static_cast<int>((particles.py[count] - bmin[1]) / sy);
            particles.pz[count] = static_cast<double>(read_binary_value<float>(file));
            particles.k[count] = static_cast<int>((particles.pz[count] - bmin[2]) / sz);

            particles.hvx[count] = static_cast<double>(read_binary_value<float>(file));
            particles.hvy[count] = static_cast<double>(read_binary_value<float>(file));
            particles.hvz[count] = static_cast<double>(read_binary_value<float>(file));
            particles.vx[count] = static_cast<double>(read_binary_value<float>(file));
            particles.vy[count] = static_cast<double>(read_binary_value<float>(file));
            particles.vz[count] = static_cast<double>(read_binary_value<float>(file));

            count += 1;
        }
        particulas_bloque(particles);
        file.close();
    }

}

void Grid::get_parameters(const std::string& filename){
    std::ifstream file(filename, std::ios::binary);
    std::ostringstream errorMessage;
    if (file){
        ppm = static_cast<double>(read_binary_value<float>(file));
        np = static_cast<int>(read_binary_value<int>(file));
        file.close();
    }

}

void Grid::constantes() {
    calculos.h_cuadrado = pow(h, 2);
    calculos.h_sexta = pow(h, NUMERO6);
    calculos.trans_densidad = (NUMERO315 / (NUMERO64 * M_PI * pow(h, NUMERO9))) * m;


    double const aux_calc = 1/(M_PI * calculos.h_sexta);
    calculos.calc1 = NUMERO15 * 3 * m * presion_rigidez * aux_calc * MEDIO;
    calculos.calc2 = 2 * densidad_fluido;
    calculos.calc3 = NUMERO45 * viscosidad * m * aux_calc;
}

void Grid::grid_properties(){
    m = densidad_fluido/pow(ppm, 3); // Masa de una particula
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

void Grid::print_grid() const {
    std::cout << "Numero de particulas: " << np << std::endl;
    std::cout << "Particulas por metro: " << ppm << std::endl;
    std::cout << "Smoothing length: " << h << std::endl;
    std::cout<< "Masa particula: " << m <<std::endl;
    std::cout<< "Tamaño del grid: " << nx << " x " << ny << " x " << nz << std::endl;
    std::cout<< "Numero de bloques: " << nx*ny*nz << std::endl;
    std::cout << "Tamaño de bloque: " << sx << " x " << sy << " x " << sz << std::endl;
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
void Grid::inicializacion_aceleracion_densidad(ParticleArray& particles) const{
    for (int i = 0; i<np; i++){
        particles.ax[i] = gravedad[0];
        particles.ay[i] = gravedad[1];
        particles.az[i] = gravedad[2];
        particles.density[i] = 0;
    }
}

void Grid::reposicionamiento_particulas(ParticleArray& particles){
    for (bloque& bloque: bloques_grid){
        bloque.index_particle_block.clear();
    }
    particulas_bloque(particles);
}

void Grid::actualizar_ac_den(ParticleArray& particles){
    for (int i=0; i<np; i++){
        int const bloque = calcular_id_bloque(particles.i[i], particles.j[i], particles.k[i]);
        for (int const bloque_adyacente : bloques_grid[bloque].adjacent_blocks){
            for (int j : bloques_grid[bloque_adyacente].index_particle_block){
                if (i < j){
                    incremento_densidad(i, j, particles);

                }
            }
        }

        transformacion_densidad(i, particles);
    }


    for (int i=0; i<np; i++){
        int const bloque = calcular_id_bloque(particles.i[i], particles.j[i], particles.k[i]);
        for (int const bloque_adyacente : bloques_grid[bloque].adjacent_blocks){
            for (int j : bloques_grid[bloque_adyacente].index_particle_block){
                if (i < j){
                    actualizar_aceleracion(i, j, particles);
                }
            }
        }
    }
}

void Grid::incremento_densidad(int i, int j, ParticleArray& particles) const{
    double modulo_cuadrado = calcular_modulo_cuadrado(i, j, particles);
    if (modulo_cuadrado < calculos.h_cuadrado){

        double const incremento = pow((calculos.h_cuadrado - modulo_cuadrado), 3);
        particles.density[i] += incremento;
        particles.density[j] += incremento;
    }

}
double Grid::calcular_modulo_cuadrado(int i, int j, ParticleArray& particles){
    //No hacemos la raiz cuadrada porque al hacer el cuadrado del modulo, la raiz se va con el cuadrado
    //es una operacion que nos ahorramos
    double modulo;
    modulo = (pow((particles.px[i] - particles.px[j]),2) + pow((particles.py[i] - particles.py[j]),2) + pow((particles.pz[i] - particles.pz[j]),2));
    return modulo;
}

void Grid::transformacion_densidad(int i, ParticleArray& particles) const{

    particles.density[i] = (particles.density[i] + pow(h, 6)) * calculos.trans_densidad;

}

void Grid::actualizar_aceleracion(int i, int j, ParticleArray& particles) const{
    double modulo_cuadrado = calcular_modulo_cuadrado(i, j, particles);
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
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_x(particle_id, particles);
        }
    }
}

void Grid::bucle_bloque_y0(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_y(particle_id, particles);
        }
    }
}

void Grid::bucle_bloque_z0(const std::vector<int> &lista_bloques, ParticleArray& particles){
    for (int block_id: lista_bloques){
        bloque const block = bloques_grid[block_id];
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
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_xmenos1(particle_id, particles);
        }
    }
}
void Grid::bucle_bloque_ymenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_ymenos1(particle_id, particles);
        }
    }
}
void Grid::bucle_bloque_zmenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
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
void Grid::movimiento_particulas(ParticleArray& particles) const{
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
    particles.vx[i] = particles.hvx[i] + (particles.ax[i] * paso_tiempo) * MEDIO;
    particles.vy[i] = particles.hvy[i] + (particles.ay[i] * paso_tiempo) * MEDIO;
    particles.vz[i] = particles.hvz[i] + (particles.az[i] * paso_tiempo) * MEDIO;
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
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_x(particle_id, particles);
        }
    }
}

void Grid::limite_bloque_y0(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_y(particle_id, particles);
        }
    }
}

void Grid::limite_bloque_z0(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
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
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_xmenos1(particle_id, particles);
        }
    }
}
void Grid::limite_bloque_ymenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_ymenos1(particle_id, particles);
        }
    }
}
void Grid::limite_bloque_zmenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
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

    auto ppm_ = static_cast<float>(ppm);
    archivo.write(reinterpret_cast<const char*>(&ppm_), sizeof(float));
    archivo.write(reinterpret_cast<const char*>(&np), sizeof(int));

    for (int i=0; i<np; i++){
        auto px_float = static_cast<float>(particles.px[i]);
        auto py_float = static_cast<float>(particles.py[i]);
        auto pz_float = static_cast<float>(particles.pz[i]);

        auto hvx_float = static_cast<float>(particles.hvx[i]);
        auto hvy_float = static_cast<float>(particles.hvy[i]);
        auto hvz_float = static_cast<float>(particles.hvz[i]);

        auto vx_float = static_cast<float>(particles.vx[i]);
        auto vy_float = static_cast<float>(particles.vy[i]);
        auto vz_float = static_cast<float>(particles.vz[i]);
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