#!/bin/sh

./etc/profile

module avail
module load gcc/12.1.0

/snap/bin/cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
/snap/bin/cmake --build build -j 4
