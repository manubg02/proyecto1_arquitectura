//
// Created by manu on 8/11/23.
//

#include "grid.hpp"
#include <stdexcept>
#include <sstream>

void guardar_particulas(ParticleArray& particles, const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);

    if (file) {
        float ppm;
        file.read(reinterpret_cast<char*>(&ppm), sizeof(ppm));

        int32_t intValue;
        file.read(reinterpret_cast<char*>(&intValue), sizeof(intValue));
        int np = *reinterpret_cast<int*>(&intValue);

        if (np <= 0){

        }

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
        particles.px.resize(np+1);
        particles.py.resize(np+1);
        particles.pz.resize(np+1);
        particles.hvx.resize(np+1);
        particles.hvy.resize(np+1);
        particles.hvz.resize(np+1);
        particles.vx.resize(np+1);
        particles.vy.resize(np+1);
        particles.vz.resize(np+1);
        particles.i.resize(np+1);
        particles.j.resize(np+1);
        particles.k.resize(np+1);

        int count = 0;
        while (count < np+1) {
            file.read(reinterpret_cast<char*>(&intValue), sizeof(intValue));
            float px = *reinterpret_cast<float *>(&intValue);
            double dpx = static_cast<double>(px);
            particles.px[count] = dpx;
            particles.i[count] = (particles.px[count]-bmin[0])/sx;
            if (particles.i[count] > nx-1){
                particles.i[count] = nx-1;
            }
            file.read(reinterpret_cast<char*>(&intValue), sizeof(intValue));
            float py = *reinterpret_cast<float *>(&intValue);
            double dpy = static_cast<double>(py);
            particles.py[count] = dpy;
            particles.j[count] = (particles.py[count]-bmin[1])/sy;
            file.read(reinterpret_cast<char*>(&intValue), sizeof(intValue));
            float pz = *reinterpret_cast<float *>(&intValue);
            double dpz = static_cast<double>(pz);
            particles.pz[count] = dpz;
            particles.k[count] = (particles.pz[count]-bmin[2])/sz;
            if (particles.k[count] > nz-1){
                particles.k[count] = nz-1;
            }
            file.read(reinterpret_cast<char*>(&intValue), sizeof(intValue));
            float hvx = *reinterpret_cast<float *>(&intValue);
            double dhvx = static_cast<double>(hvx);
            particles.hvx[count] = dhvx;
            file.read(reinterpret_cast<char*>(&intValue), sizeof(intValue));
            float hvy = *reinterpret_cast<float *>(&intValue);
            double dhvy = static_cast<double>(hvy);
            particles.hvy[count] = dhvy;
            file.read(reinterpret_cast<char*>(&intValue), sizeof(intValue));
            float hvz = *reinterpret_cast<float *>(&intValue);
            double dhvz = static_cast<double>(hvz);
            particles.hvz[count] = dhvz;
            file.read(reinterpret_cast<char*>(&intValue), sizeof(intValue));
            float vx = *reinterpret_cast<float *>(&intValue);
            double dvx = static_cast<double>(vx);
            particles.vx[count] = dvx;
            file.read(reinterpret_cast<char*>(&intValue), sizeof(intValue));
            float vy = *reinterpret_cast<float *>(&intValue);
            double dvy = static_cast<double>(vy);
            particles.vy[count] = dvy;
            file.read(reinterpret_cast<char*>(&intValue), sizeof(intValue));
            float vz = *reinterpret_cast<float *>(&intValue);
            double dvz = static_cast<double>(vz);
            particles.vz[count] = dvz;

            if (count > 1){
                std::cout<< particles.px[count] <<" == " << particles.vz[count-1] << "?"  << std::endl;
                if (particles.px[count] == particles.vz[count-1]){
                    std::ostringstream errorMessage;
                    errorMessage << "Number of particles mismatch. Header: " << np << ", Found: " << count;
                    throw std::runtime_error(errorMessage.str());
                }
            }

            count += 1;

        }

        std::cout<<"Prueba"<<std::endl;



        // Cerrar el archivo después de usarlo
        file.close();

    } else {
        std::cerr << "Error opening the file." << std::endl;
    }
}