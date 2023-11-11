//
// Created by manu on 8/11/23.
//

#include "grid.hpp"

void guardar_particulas(ParticleArray& particles, const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);

    if (file) {
        float ppm;
        file.read(reinterpret_cast<char*>(&ppm), sizeof(ppm));

        int32_t intValue;
        file.read(reinterpret_cast<char*>(&intValue), sizeof(intValue));
        int np = *reinterpret_cast<int*>(&intValue);

        std::cout << "Numero de particulas: " << np << std::endl;
        std::cout << "Particulas por metro: " << ppm << std::endl;

        //Calcular h después de leer ppm
        const float h = multiplicador_radio/ppm;
        std::cout << "Smoothing length: " << h << std::endl;

        // Calcular m después de leer ppm
        const float m = densidad_fluido / (ppm * ppm * ppm);
        std::cout<< "Masa particula: " << m <<std::endl;

        const int nx = (bmax[0] - bmin[0])/h;
        const int ny = (bmax[1] - bmin[1])/h;
        const int nz = (bmax[2] - bmin[2])/h;

        std::cout<< "Tamaño del grid: " << nx << " x " << ny << " x " << nz << std::endl;
        std::cout<< "Numero de bloques: " << nx*ny*nz << std::endl;

        const float sx = (bmax[0] - bmin[0])/nx;
        const float sy = (bmax[1] - bmin[1])/ny;
        const float sz = (bmax[2] - bmin[2])/nz;

        std::cout << "Tamaño de bloque: " << sx << " x " << sy << " x " << sz << std::endl;



        // Redimensiona los vectores según el número total de partículas
        particles.px.resize(np);
        particles.py.resize(np);
        particles.pz.resize(np);
        particles.hvx.resize(np);
        particles.hvy.resize(np);
        particles.hvz.resize(np);
        particles.vx.resize(np);
        particles.vy.resize(np);
        particles.vz.resize(np);

        int count = 0;
        while (count < np) {
            file.read(reinterpret_cast<char*>(&particles.px[count]), sizeof(float));
            file.read(reinterpret_cast<char*>(&particles.py[count]), sizeof(float));
            file.read(reinterpret_cast<char*>(&particles.pz[count]), sizeof(float));
            file.read(reinterpret_cast<char*>(&particles.hvx[count]), sizeof(float));
            file.read(reinterpret_cast<char*>(&particles.hvy[count]), sizeof(float));
            file.read(reinterpret_cast<char*>(&particles.hvz[count]), sizeof(float));
            file.read(reinterpret_cast<char*>(&particles.vx[count]), sizeof(float));
            file.read(reinterpret_cast<char*>(&particles.vy[count]), sizeof(float));
            file.read(reinterpret_cast<char*>(&particles.vz[count]), sizeof(float));


            count += 1;
        }

        // Cerrar el archivo después de usarlo
        file.close();

    } else {
        std::cerr << "Error opening the file." << std::endl;
    }
}