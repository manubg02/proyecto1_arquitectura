//
// Created by manu on 8/11/23.
//

#ifndef PROYECTO1_ARQUITECTURA_BLOCK_HPP
#define PROYECTO1_ARQUITECTURA_BLOCK_HPP

#include <iostream>
#include <vector>

struct bloque{
    std::vector<int> adjacent_blocks; //Se incluye a el mismo
    std::vector<int> index_particle_block;
};

#endif //PROYECTO1_ARQUITECTURA_BLOCK_HPP
