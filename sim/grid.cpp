//
// Created by manu on 8/11/23.
//

#include "grid.hpp"
#include <stdexcept>
#include <sstream>

template <typename T>
T leer_binario(std::ifstream& infile) {
    // Función auxiliar para leer valores binarios del archivo de entrada.
    // Utiliza reinterpret_cast para interpretar los bytes leídos como el tipo T.

    T value;
    infile.read(reinterpret_cast<char*>(&value), sizeof(T));
    return value;
}



Grid::Grid(int time, std::string const& file_input, std::string const& file_output)
        : timeSteps(time), ppm(0), np(0), m(0), h(0), nx(0), ny(0), nz(0), sx(0), sy(0), sz(0){
    // Constructor de la clase Grid.

    // Inicializa las variables con valores predeterminados e igual a 0.
    get_parameters(file_input);  // Obtiene los parámetros iniciales del archivo de entrada.
    grid_properties();  // Calcula las propiedades del grid.
    constantes();  // Calcula constantes utilizadas en la simulación.
    print_grid();  // Imprime información sobre el grid.

    crear_bloques();  // Crea bloques en el grid.

    ParticleArray particles;  // Define un objeto ParticleArray para almacenar las partículas.
    meter_particulas(file_input, particles);  // Lee las partículas del archivo de entrada.
    simulacion(particles);  // Realiza la simulación.
    escribir_informacion(particles, file_output);  // Escribe la información de las partículas en el archivo de salida.
}


void Grid::meter_particulas(const std::string& filename, ParticleArray& particles) {
    // Lee las partículas del archivo binario y las almacena en el objeto ParticleArray.
    std::ifstream file(filename, std::ios::binary);
    // Realiza la lectura de los parámetros.
    if (file){
        file.seekg(8, std::ios::cur); // Colocamos el puntero de lectura del fichero en el byte correcto
        // Resize de todos los parametros en funcion de np
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
            // Loop para obtener todos los valores de los parametros de cada variable y guardarlo en nuestro ParticleArray
            particles.px[count] = static_cast<double>(leer_binario<float>(file));
            particles.i[count] = static_cast<int>((particles.px[count] - bmin[0]) / sx);

            particles.py[count] = static_cast<double>(leer_binario<float>(file));
            particles.j[count] = static_cast<int>((particles.py[count] - bmin[1]) / sy);
            particles.pz[count] = static_cast<double>(leer_binario<float>(file));
            particles.k[count] = static_cast<int>((particles.pz[count] - bmin[2]) / sz);

            particles.hvx[count] = static_cast<double>(leer_binario<float>(file));
            particles.hvy[count] = static_cast<double>(leer_binario<float>(file));
            particles.hvz[count] = static_cast<double>(leer_binario<float>(file));
            particles.vx[count] = static_cast<double>(leer_binario<float>(file));
            particles.vy[count] = static_cast<double>(leer_binario<float>(file));
            particles.vz[count] = static_cast<double>(leer_binario<float>(file));

            count += 1;
        }
        particulas_bloque(particles);
        file.close();
    }

}

void Grid::get_parameters(const std::string& filename){
    // Obtiene parámetros iniciales del archivo binario.
    std::ifstream file(filename, std::ios::binary);
    std::ostringstream errorMessage;
    // Realiza la lectura de los parámetros.
    if (file){
        ppm = static_cast<double>(leer_binario<float>(file)); // Leemos el valor de ppm del fichero y lo pasamos directamente a double
        np = static_cast<int>(leer_binario<int>(file)); // Leemos el valor de np del fichero y lo pasamos directamente a double
        file.close();
    }

}

void Grid::constantes() {
    // Calcula constantes utilizadas en la simulación.
    calculos.h_cuadrado = pow(h, 2); // h^2
    calculos.h_sexta = pow(h, NUMERO6); // h^6
    calculos.trans_densidad = (NUMERO315 / (NUMERO64 * M_PI * pow(h, NUMERO9))) * m; // Cálculo de la densidad

    // Partes de las formulas de actualización de las partículas
    double const aux_calc = 1/(M_PI * calculos.h_sexta);
    calculos.calc1 = NUMERO15 * 3 * m * presion_rigidez * aux_calc * MEDIO;
    calculos.calc2 = 2 * densidad_fluido;
    calculos.calc3 = NUMERO45 * viscosidad * m * aux_calc;
}

void Grid::grid_properties(){
    // Calcula propiedades del grid en función de los parámetros iniciales.
    m = densidad_fluido/pow(ppm, 3); // Masa de una partícula
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
    // Imprime información sobre el grid.
    std::cout << "Numero de particulas: " << np << std::endl;
    std::cout << "Particulas por metro: " << ppm << std::endl;
    std::cout << "Smoothing length: " << h << std::endl;
    std::cout<< "Masa particula: " << m <<std::endl;
    std::cout<< "Tamaño del grid: " << nx << " x " << ny << " x " << nz << std::endl;
    std::cout<< "Numero de bloques: " << nx*ny*nz << std::endl;
    std::cout << "Tamaño de bloque: " << sx << " x " << sy << " x " << sz << std::endl;
}

void Grid::simulacion(ParticleArray& particles){
    // Realiza la simulación en múltiples pasos de tiempo.
    // Hacemos la iteración 0 por separado para ahorrarnos hacer la función de reposicionamiento
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
// Inicializa la aceleración y densidad de las partículas en el grid.
void Grid::inicializacion_aceleracion_densidad(ParticleArray& particles) const{
    for (int i = 0; i < np; i++) {
        // Inicializa la aceleración de las partículas con el vector de gravedad.
        particles.ax[i] = gravedad[0];
        particles.ay[i] = gravedad[1];
        particles.az[i] = gravedad[2];

        // Inicializa la densidad de las partículas a cero.
        particles.density[i] = 0;
    }
}

// Repositiona las partículas en los bloques del grid.
void Grid::reposicionamiento_particulas(ParticleArray& particles){
    // Limpia la lista de índices de partículas en cada bloque.
    for (bloque& bloque : bloques_grid) {
        bloque.index_particle_block.clear();
    }

    // Llena nuevamente la lista de índices de partículas en cada bloque.
    particulas_bloque(particles);
}

// Actualiza la aceleración y densidad de las partículas en función de sus vecinos en bloques adyacentes.
void Grid::actualizar_ac_den(ParticleArray& particles){
    for (int i = 0; i < np; i++) {
        // Calcula el bloque en el que se encuentra la partícula.
        int const bloque = calcular_id_bloque(particles.i[i], particles.j[i], particles.k[i]);

        // Itera sobre los bloques adyacentes al bloque de la partícula.
        for (int const bloque_adyacente : bloques_grid[bloque].adjacent_blocks) {
            // Itera sobre las partículas en el bloque adyacente.
            for (int j : bloques_grid[bloque_adyacente].index_particle_block) {
                // Asegura que no se hagan cálculos redundantes para pares de partículas ya procesados.
                if (i < j) {
                    // Incrementa la densidad de las partículas cercanas.
                    incremento_densidad(i, j, particles);
                }
            }
        }

        // Realiza la transformación de densidad para la partícula actual.
        transformacion_densidad(i, particles);
    }

    // Itera nuevamente para actualizar la aceleración basada en la densidad.
    for (int i = 0; i < np; i++) {
        int const bloque = calcular_id_bloque(particles.i[i], particles.j[i], particles.k[i]);

        for (int const bloque_adyacente : bloques_grid[bloque].adjacent_blocks) {
            for (int j : bloques_grid[bloque_adyacente].index_particle_block) {
                if (i < j) {
                    // Actualiza la aceleración de las partículas basándose en la densidad.
                    actualizar_aceleracion(i, j, particles);
                }
            }
        }
    }
}

// Incrementa la densidad de dos partículas si la distancia al cuadrado entre ellas es menor que el cuadrado del suavizado.
void Grid::incremento_densidad(int i, int j, ParticleArray& particles) const {
    // Calcula el cuadrado de la distancia entre las partículas.
    double modulo_cuadrado = calcular_modulo_cuadrado(i, j, particles);

    // Verifica si la distancia al cuadrado es menor que el suavizado al cuadrado.
    if (modulo_cuadrado < calculos.h_cuadrado) {
        // Calcula el incremento en densidad según el kernel de suavizado.
        double const incremento = pow((calculos.h_cuadrado - modulo_cuadrado), 3);

        // Incrementa la densidad de ambas partículas.
        particles.density[i] += incremento;
        particles.density[j] += incremento;
    }
}

// Calcula el cuadrado del módulo entre dos partículas en el espacio tridimensional.
double Grid::calcular_modulo_cuadrado(int i, int j, ParticleArray& particles) {
    // Calcula la distancia al cuadrado entre las partículas en cada dimensión.
    double modulo = (pow((particles.px[i] - particles.px[j]), 2) +
                     pow((particles.py[i] - particles.py[j]), 2) +
                     pow((particles.pz[i] - particles.pz[j]), 2));

    // No se realiza la raíz cuadrada, ya que se utiliza el módulo al cuadrado en el cálculo de densidad.
    return modulo;
}

// Realiza la transformación de densidad para una partícula específica.
void Grid::transformacion_densidad(int i, ParticleArray& particles) const {
    // Aplica la transformación de densidad a la partícula según el kernel de suavizado.
    particles.density[i] = (particles.density[i] + pow(h, 6)) * calculos.trans_densidad;
}


// Actualiza la aceleración de dos partículas si la distancia al cuadrado entre ellas es menor que el cuadrado del suavizado.
void Grid::actualizar_aceleracion(int i, int j, ParticleArray& particles) const {
    // Calcula el cuadrado de la distancia entre las partículas.
    double modulo_cuadrado = calcular_modulo_cuadrado(i, j, particles);

    // Verifica si la distancia al cuadrado es menor que el suavizado al cuadrado.
    if (modulo_cuadrado < calculos.h_cuadrado) {
        // Calcula la distancia y evita divisiones por cero.
        double const dij = sqrt(std::max(modulo_cuadrado, 1e-12));

        // Variables para las componentes de la aceleración.
        double delta_aij_x;
        double delta_aij_y;
        double delta_aij_z;

        // Cálculos intermedios para la actualización de la aceleración.
        double const parte1 = calculos.calc1 * pow((h - dij), 2) * (particles.density[i] + particles.density[j] - calculos.calc2) * (1 / dij);
        double auxiliar = 1 / (particles.density[i] * particles.density[j]);

        // Actualiza las componentes de la aceleración para las dos partículas.
        delta_aij_x = ((particles.px[i] - particles.px[j]) * parte1 + (particles.vx[j] - particles.vx[i]) * calculos.calc3) * auxiliar;
        delta_aij_y = ((particles.py[i] - particles.py[j]) * parte1 + (particles.vy[j] - particles.vy[i]) * calculos.calc3) * auxiliar;
        delta_aij_z = ((particles.pz[i] - particles.pz[j]) * parte1 + (particles.vz[j] - particles.vz[i]) * calculos.calc3) * auxiliar;

        particles.ax[i] += delta_aij_x;
        particles.ay[i] += delta_aij_y;
        particles.az[i] += delta_aij_z;
        particles.ax[j] -= delta_aij_x;
        particles.ay[j] -= delta_aij_y;
        particles.az[j] -= delta_aij_z;
    }
}

// Realiza las colisiones con los límites del dominio.
void Grid::colisiones(ParticleArray& particles) {
    // Aplica colisiones en los límites de las coordenadas x, y, y z.
    bucle_bloque_x0(nx_0, particles);
    bucle_bloque_y0(ny_0, particles);
    bucle_bloque_z0(nz_0, particles);

    bucle_bloque_xmenos1(nx_menos1, particles);
    bucle_bloque_ymenos1(ny_menos1, particles);
    bucle_bloque_zmenos1(nz_menos1, particles);
}


// Realiza colisiones en el eje x para las partículas en los bloques dados.
void Grid::bucle_bloque_x0(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_x(particle_id, particles);
        }
    }
}

// Realiza colisiones en el eje y para las partículas en los bloques dados.
void Grid::bucle_bloque_y0(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_y(particle_id, particles);
        }
    }
}

// Realiza colisiones en el eje z para las partículas en los bloques dados.
void Grid::bucle_bloque_z0(const std::vector<int> &lista_bloques, ParticleArray& particles){
    for (int block_id: lista_bloques){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_z(particle_id, particles);
        }
    }
}

// Realiza colisiones en el eje x-1 para las partículas en los bloques dados.
void Grid::bucle_bloque_xmenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_xmenos1(particle_id, particles);
        }
    }
}

// Realiza colisiones en el eje y-1 para las partículas en los bloques dados.
void Grid::bucle_bloque_ymenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_ymenos1(particle_id, particles);
        }
    }
}

// Realiza colisiones en el eje z-1 para las partículas en los bloques dados.
void Grid::bucle_bloque_zmenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            colisiones_particulas_eje_zmenos1(particle_id, particles);
        }
    }
}

// Realiza colisiones en el eje x para una partícula específica.
void Grid::colisiones_particulas_eje_x(int id, ParticleArray& particles){
    double const limite_eje_x = particles.px[id] + particles.hvx[id] * paso_tiempo;
    double const incremento_x = t_particula - (limite_eje_x - bmin[0]);

    // Verifica si se debe aplicar la colisión.
    if (incremento_x > comparar_colision){
        double const aux = colision_rigidez * incremento_x - amortiguamiento * particles.vx[id];
        particles.ax[id] += aux;
    }
}

// Realiza colisiones en el eje y para una partícula específica.
void Grid::colisiones_particulas_eje_y(int id, ParticleArray& particles){
    double const limite_eje_y = particles.py[id] + particles.hvy[id] * paso_tiempo;
    double const incremento_y = t_particula - (limite_eje_y - bmin[1]);

    // Verifica si se debe aplicar la colisión.
    if (incremento_y > comparar_colision){
        double const aux = colision_rigidez * incremento_y - amortiguamiento * particles.vy[id];
        particles.ay[id] += aux;
    }
}

// Realiza colisiones en el eje z para una partícula específica.
void Grid::colisiones_particulas_eje_z(int id, ParticleArray& particles){
    double const limite_eje_z = particles.pz[id] + particles.hvz[id] * paso_tiempo;
    double const incremento_z = t_particula - (limite_eje_z - bmin[2]);

    // Verifica si se debe aplicar la colisión.
    if (incremento_z > comparar_colision){
        double const aux = colision_rigidez * incremento_z - amortiguamiento * particles.vz[id];
        particles.az[id] += aux;
    }
}

// Realiza colisiones en el eje x-1 para una partícula específica.
void Grid::colisiones_particulas_eje_xmenos1(int id, ParticleArray& particles){
    double const limite_eje_x = particles.px[id] + particles.hvx[id] * paso_tiempo;
    double const incremento_x = t_particula - (bmax[0] - limite_eje_x);

    // Verifica si se debe aplicar la colisión.
    if (incremento_x > comparar_colision){
        double const aux = colision_rigidez * incremento_x + amortiguamiento * particles.vx[id];
        particles.ax[id] -= aux;
    }
}

// Realiza colisiones en el eje y-1 para una partícula específica.
void Grid::colisiones_particulas_eje_ymenos1(int id, ParticleArray& particles){
    double const limite_eje_y = particles.py[id] + particles.hvy[id] * paso_tiempo;
    double const incremento_y = t_particula - (bmax[1] - limite_eje_y);

    // Verifica si se debe aplicar la colisión.
    if (incremento_y > comparar_colision){
        double const aux = colision_rigidez * incremento_y + amortiguamiento * particles.vy[id];
        particles.ay[id] -= aux;
    }
}

// Realiza colisiones en el eje z-1 para una partícula específica.
void Grid::colisiones_particulas_eje_zmenos1(int id, ParticleArray& particles){
    double const limite_eje_z = particles.pz[id] + particles.hvz[id] * paso_tiempo;
    double const incremento_z = t_particula - (bmax[2] - limite_eje_z);

    // Verifica si se debe aplicar la colisión.
    if (incremento_z > comparar_colision){
        double const aux = colision_rigidez * incremento_z + amortiguamiento * particles.vz[id];
        particles.az[id] -= aux;
    }
}

// Realiza el movimiento de partículas, actualizando posición, velocidad y gradiente.
void Grid::movimiento_particulas(ParticleArray& particles) const{
    for (int i = 0; i<np; i++){
        act_posicion(i, particles);
        act_velocidad(i, particles);
        act_gradiente(i, particles);
    }
}

// Actualiza la posición de una partícula específica.
void Grid::act_posicion(int i, ParticleArray& particles){
    particles.px[i] += particles.hvx[i] * paso_tiempo + particles.ax[i] * pow(paso_tiempo, 2);
    particles.py[i] += particles.hvy[i] * paso_tiempo + particles.ay[i] * pow(paso_tiempo, 2);
    particles.pz[i] += particles.hvz[i] * paso_tiempo + particles.az[i] * pow(paso_tiempo, 2);
}

// Actualiza la velocidad de una partícula específica.
void Grid::act_velocidad(int i, ParticleArray& particles){
    particles.vx[i] = particles.hvx[i] + (particles.ax[i] * paso_tiempo) * MEDIO;
    particles.vy[i] = particles.hvy[i] + (particles.ay[i] * paso_tiempo) * MEDIO;
    particles.vz[i] = particles.hvz[i] + (particles.az[i] * paso_tiempo) * MEDIO;
}

// Actualiza el gradiente de una partícula específica.
void Grid::act_gradiente(int i, ParticleArray& particles){
    particles.hvx[i] += particles.ax[i] * paso_tiempo;
    particles.hvy[i] += particles.ay[i] * paso_tiempo;
    particles.hvz[i] += particles.az[i] * paso_tiempo;
}


// Aplica límites al recinto para evitar que las partículas salgan de él.
void Grid::limites_recinto(ParticleArray& particles){
    limite_bloque_x0(nx_0, particles);
    limite_bloque_y0(ny_0, particles);
    limite_bloque_z0(nz_0, particles);

    limite_bloque_xmenos1(nx_menos1, particles);
    limite_bloque_ymenos1(ny_menos1, particles);
    limite_bloque_zmenos1(nz_menos1, particles);
}

// Aplica límites en el eje X al conjunto de partículas en un bloque específico.
void Grid::limite_bloque_x0(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_x(particle_id, particles);
        }
    }
}

// Aplica límites en el eje Y al conjunto de partículas en un bloque específico.
void Grid::limite_bloque_y0(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_y(particle_id, particles);
        }
    }
}

// Aplica límites en el eje Z al conjunto de partículas en un bloque específico.
void Grid::limite_bloque_z0(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_z(particle_id, particles);
        }
    }
}

// Aplica límites en el eje X a una partícula específica.
void Grid::limite_particulas_eje_x(int id, ParticleArray &particles) {
    double const d_x = particles.px[id] - bmin[0];
    if (d_x < 0){
        particles.px[id] = bmin[0] - d_x;
        particles.vx[id] = -particles.vx[id];
        particles.hvx[id] = -particles.hvx[id];
    }
}

// Aplica límites en el eje Y a una partícula específica.
void Grid::limite_particulas_eje_y(int id, ParticleArray &particles) {
    double const d_y = particles.py[id] - bmin[1];
    if (d_y < 0){
        particles.py[id] = bmin[1] - d_y;
        particles.vy[id] = -particles.vy[id];
        particles.hvy[id] = -particles.hvy[id];
    }
}

// Aplica límites en el eje Z a una partícula específica.
void Grid::limite_particulas_eje_z(int id, ParticleArray &particles) {
    double const d_z = particles.pz[id] - bmin[2];
    if (d_z < 0){
        particles.pz[id] = bmin[2] - d_z;
        particles.vz[id] = -particles.vz[id];
        particles.hvz[id] = -particles.hvz[id];
    }
}

// Aplica límites en el eje X a una partícula específica en el bloque adyacente a x=0.
void Grid::limite_bloque_xmenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_xmenos1(particle_id, particles);
        }
    }
}

// Aplica límites en el eje Y a una partícula específica en el bloque adyacente a y=0.
void Grid::limite_bloque_ymenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_ymenos1(particle_id, particles);
        }
    }
}

// Aplica límites en el eje Z a una partícula específica en el bloque adyacente a z=0.
void Grid::limite_bloque_zmenos1(const std::vector<int> &block_list, ParticleArray& particles){
    for (int block_id: block_list){
        bloque const block = bloques_grid[block_id];
        for (int particle_id: block.index_particle_block){
            limite_particulas_eje_zmenos1(particle_id, particles);
        }
    }
}

// Aplica límites en el eje X a una partícula específica en el bloque adyacente a x=0.
void Grid::limite_particulas_eje_xmenos1(int id, ParticleArray& particles){
    double const d_x = bmax[0] - particles.px[id];
    if (d_x < 0){
        particles.px[id] = bmax[0] + d_x;
        particles.vx[id] = -particles.vx[id];
        particles.hvx[id] = -particles.hvx[id];
    }
}

// Aplica límites en el eje Y a una partícula específica en el bloque adyacente a y=0.
void Grid::limite_particulas_eje_ymenos1(int id, ParticleArray& particles){
    double const d_y = bmax[1] - particles.py[id];
    if (d_y < 0){
        particles.py[id] = bmax[1] + d_y;
        particles.vy[id] = -particles.vy[id];
        particles.hvy[id] = -particles.hvy[id];
    }
}

// Aplica límites en el eje Z a una partícula específica en el bloque adyacente a z=0.
void Grid::limite_particulas_eje_zmenos1(int id, ParticleArray& particles){
    double const d_z= bmax[2] - particles.pz[id];
    if (d_z < 0){
        particles.pz[id] = bmax[2] + d_z;
        particles.vz[id] = -particles.vz[id];
        particles.hvz[id] = -particles.hvz[id];
    }
}


// Escribe la información de las partículas en un archivo binario.
void Grid::escribir_informacion(ParticleArray &particles, const std::string &file_output) {
    // Abre un archivo en modo binario para escritura.
    std::ofstream archivo(file_output, std::ios::binary);

    // Verifica si el archivo se abrió correctamente.
    if (!archivo.is_open()) {
        std::cerr << "Error: No se pudo abrir el archivo para escritura." << std::endl;
        return;
    }

    // Convierte la resolución a formato de punto flotante y escribe en el archivo.
    auto ppm_ = static_cast<float>(ppm);
    archivo.write(reinterpret_cast<const char*>(&ppm_), sizeof(float));

    // Escribe el número total de partículas en el archivo.
    archivo.write(reinterpret_cast<const char*>(&np), sizeof(int));

    // Recorre todas las partículas y escribe sus propiedades en el archivo.
    for (int i=0; i<np; i++){
        // Convierte las coordenadas de posición a formato de punto flotante y las escribe en el archivo.
        auto px_float = static_cast<float>(particles.px[i]);
        auto py_float = static_cast<float>(particles.py[i]);
        auto pz_float = static_cast<float>(particles.pz[i]);

        // Convierte las velocidades medias a formato de punto flotante y las escribe en el archivo.
        auto hvx_float = static_cast<float>(particles.hvx[i]);
        auto hvy_float = static_cast<float>(particles.hvy[i]);
        auto hvz_float = static_cast<float>(particles.hvz[i]);

        // Convierte las velocidades a formato de punto flotante y las escribe en el archivo.
        auto vx_float = static_cast<float>(particles.vx[i]);
        auto vy_float = static_cast<float>(particles.vy[i]);
        auto vz_float = static_cast<float>(particles.vz[i]);

        // Escribe en el archivo las coordenadas y velocidades de la partícula actual.
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

    // Cierra el archivo después de escribir toda la información.
    archivo.close();
}